#pragma once
#include <cstdio>
#include <exception>
#include <mpi.h>
#include <string>
#include <sstream>
#include <chrono>
#include <mutex>
#include <chrono>
#include <thread>
#include <condition_variable>
#include <array>
#include <cstdlib>
#include <cstring>

#include "log.hpp"
#include "setsm_code.hpp"

static bool requested_custom_async_progress() {
    static const std::array<std::string, 6> enable_strs = {
        "1", "y", "on", "yes", "enable" };

    const char *s = getenv("SETSM_CUSTOM_ASYNC");
    if(!s)
        return false;
    
    // convert to lowercase for comparison
    std::unique_ptr<char []> lower(new char[strlen(s) +1]());
    for(int i = 0; s[i]; i++){
          lower[i] = tolower(s[i]);
    }

    for(auto &s : enable_strs) {
        if(s == lower.get()) {
            return true;
        }
    }
    return false;
}

// adapted from https://stackoverflow.com/a/29775639
class InterruptableTimer {
public:
    /** sleep for time, interruptable
     *
     * Returns false if interrupted, otherwise true
     **/
    template<class R, class P>
    bool sleep( std::chrono::duration<R,P> const& time ) {
        auto is_stopped = [&](){ return _stop; };

        std::unique_lock<std::mutex> lock(m);
        return !cv.wait_for(lock, time, is_stopped);
    }

    void stop() {
        std::unique_lock<std::mutex> lock(m);
        _stop = true;
        cv.notify_all();
    }

private:
    std::mutex m;
    std::condition_variable cv;
    bool _stop = false;

};

template<class R, class P>
void async_progress(InterruptableTimer &timer, std::chrono::duration<R,P> interval) {
    int flag;
    MPI_Status status;
    while(timer.sleep(interval)) {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    }
}

class MPIProgressor {
public:
    template<class R, class P>
    MPIProgressor(std::chrono::duration<R,P> interval) :
        worker(&async_progress<R, P>, std::ref(timer), interval) {
        }

    ~MPIProgressor() {
        timer.stop();
        worker.join();
    }
        

private:
    InterruptableTimer timer;
    std::thread worker;


};


class MPIException: public std::exception
{
public:
    MPIException(const std::string& func, int error_code, const std::string& msg) {
        int len;
        char err_buf[MPI_MAX_ERROR_STRING + 1] = {};
        const char *err_msg = err_buf;

        int ret = MPI_Error_string(error_code, err_buf, &len); 
        if(ret != MPI_SUCCESS)
            err_msg = "Unknown error code";

        std::ostringstream os;
        os << "function " << func << " returned (" << error_code << ") " << err_msg
           << ": " << msg;
        m_msg = os.str();
    }

    virtual const char* what() const throw () { return m_msg.c_str(); }

private:
   std::string m_msg;
};

class MPICounter {
    public:
    MPICounter(int start_value=0) {
        int rank, ret;

        ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Comm_rank", ret, "Could not get comm rank for counter");
        }

        int sz = rank == 0 ? sizeof(int) : 0;

        ret = MPI_Alloc_mem(sz, MPI_INFO_NULL, &local);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Alloc_mem", ret, "Could not allocate mem for counter window");
        }

        // initialize counter to zero
        if(sz)
            local[0] = start_value;

        ret = MPI_Win_create(local, sz, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Win_create", ret, "Could not create counter window");
        }
    }

    ~MPICounter() {
        MPI_Win_free(&win);
        MPI_Free_mem(local);

    }

    MPICounter(const MPICounter&) = delete;
    MPICounter operator=(const MPICounter&) = delete;

    int next() {
        int ret = MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Win_lock", ret, "Failed to lock counter window");
        }

        const int one = 1;
        int val;
        ret = MPI_Fetch_and_op(&one, &val, MPI_INT, 0, 0, MPI_SUM, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Fetch_and_op", ret, "failed to get and increment counter");
        }

        ret = MPI_Win_unlock(0, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_win_unlock", ret, "Failed to unlock window for counter");
        }

        return val;
    }
private:
    int *local;
    MPI_Win win;

};

class MPITileIndexer : public TileIndexer {
private:
    std::shared_ptr<MPICounter> counter;
    int length;
    int rank;
public:
    MPITileIndexer(int length, int rank) : counter(new MPICounter), length(length), rank(rank) {}
    int next() {
        int i = counter->next();
        if(i < length) {
            return i;
        }
        return -1;
    }
};

// check to see if there's async progress with current
// configuration. Do this by creating a counter, and have
// other ranks try to increment it. If they are able to
// without any actions by rank0 (which is where the counter
// memory lives), then we have async progress spport.
// Otherwise, if they're not able to increment the counter
// until rank0 gets to it, then async progress support is
// not there or not good enough.
static bool has_async_support() {
    int rank, ret, size;

    ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(ret != MPI_SUCCESS) {
        throw MPIException("MPI_Comm_rank", ret, "Could not get comm rank");
    }
    ret = MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(ret != MPI_SUCCESS) {
        throw MPIException("MPI_Comm_size", ret, "Could not get comm size");
    }
    if(size <= 1) {
        // async communications don't really happen if
        // there's only 1 rank
        return true;
    }

    MPICounter counter;

    bool supports_async = true;

    // each process other than rank 0 will be incrementing the counter
    // That's (size - 1) processes. If rank zero waits for n seconds
    // before reading from the counter, let's require that the other
    // processes each be able to read from the counter n/2 times
    const int sleep_delay = 10;
    const int min_expected_count = (size - 1) * (sleep_delay / 2);
    // let's set a max here to keep things sane
    const int max_count = (size - 1) * sleep_delay;

    if(rank > 0) {
        for(;;) {
            int val = counter.next();
            if(val > max_count)
                break;
        }
    } else {
        sleep(10);
        int val = counter.next();
        supports_async = (val >= min_expected_count);
    }

    ret = MPI_Bcast(&supports_async, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    if(ret != MPI_SUCCESS) {
        throw MPIException("MPI_Comm_size", ret, "Could not get comm size");
    }

    return supports_async;
}

static void init_mpi(int&argc, char* argv[]) {
    int provided = 1;
    bool wants_custom_async = requested_custom_async_progress();
    int requested = wants_custom_async ? MPI_THREAD_MULTIPLE : MPI_THREAD_FUNNELED;
    MPI_Init_thread(&argc, &argv, requested, &provided);
    if(provided < requested) {
        fprintf(stderr, "ERROR: mpi requested threading support level %d but got %d\n",
                requested, provided);
        MPI_Finalize();
        exit(1);
    }

    if(!wants_custom_async && !has_async_support()) {
        fprintf(stderr, "ERROR: "
                "Async support not detected for current MPI configuration. "
                "Please enable async progress support for your MPI library "
                "or enable setsm's custom async progress feature by setting "
                "SETSM_CUSTOM_ASYNC=1 in the environment. See the REAMDE "
                "for more information\n");
        MPI_Finalize();
        exit(1);
    }
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      printf("MPI: Number of processes: %d\n", size);
    }
}
