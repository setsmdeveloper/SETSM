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

static bool need_custom_async_progress() {
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
        LOG("timer hit, probing MPI\n");
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
    }
    LOG("worker exiting...\n");
}

class MPIProgressor {
public:
    template<class R, class P>
    MPIProgressor(std::chrono::duration<R,P> interval) :
        worker(&async_progress<R, P>, std::ref(timer), interval) {
            LOG("progressor starting\n");
        }

    ~MPIProgressor() {
        LOG("MPIProgressor destructor entrance. stopping timer..\n");
        timer.stop();
        LOG("timer stopped...\n");
        worker.join();
        LOG("worked has joined, MPIProgressor exiting destructor\n");
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
        LOG("freeing window...\n");
        MPI_Win_free(&win);
        LOG("window freed\n");
        LOG("freeing memory...\n");
        MPI_Free_mem(local);
        LOG("memory freed\n");

    }

    MPICounter(const MPICounter&) = delete;
    MPICounter operator=(const MPICounter&) = delete;

    int next() {
        LOG("locking window...\n");
        int ret = MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Win_lock", ret, "Failed to lock counter window");
        }

        LOG("window locked. Fetching and updaating...\n");
        const int one = 1;
        int val;
        ret = MPI_Fetch_and_op(&one, &val, MPI_INT, 0, 0, MPI_SUM, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Fetch_and_op", ret, "failed to get and increment counter");
        }

        LOG("fetch and update done. unlocking...\n");
        ret = MPI_Win_unlock(0, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_win_unlock", ret, "Failed to unlock window for counter");
        }
        LOG("unlocked. returning %d\n", val);

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
        printf("DBG: rank %d calling counter->next()...\n", rank);
        int i = counter->next();
        printf("DBG: rank %d counter->next() returned %d\n", rank, i);
        if(i < length) {
            printf("DBG: rank %d (%d < %d), returning %d from next()\n", rank, i, length, i);
            return i;
        }
        printf("DBG: rank %d (%d >= %d), returning %d from next()\n", rank, i, length, -1);
        return -1;
    }
};

