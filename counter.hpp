#pragma once
#include <cstdio>
#include <exception>
#include <mpi.h>
#include <string>
#include <sstream>
#include <chrono>
#include <cstdarg>

static void log(const char *fmt, ...) {
    static bool started = 0;
    static long start_time;
    static int rank;

    unsigned long now = 
        std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();

    if(!started) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if(rank == 0) {
            start_time = now;
        }
        int ret = MPI_Bcast(&start_time, 1,MPI_LONG, 0, MPI_COMM_WORLD);
        if(ret != MPI_SUCCESS) {
            printf("FAILED to init logging!\n");
            return;
        }
        started = 1;
    }
    double elapsed = now - start_time;
    elapsed /= 1000;

    char buf[1024];

    int n = sprintf(buf, "%02d - %06.2f: ", rank, elapsed);

    va_list ap;
    va_start(ap, fmt);
    vsprintf(buf + n, fmt, ap);
    va_end(ap);

    printf("%s", buf);

}


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
        log("freeing window...\n");
        MPI_Win_free(&win);
        log("window freed\n");
        log("freeing memory...\n");
        MPI_Free_mem(local);
        log("memory freed\n");

    }

    MPICounter(const MPICounter&) = delete;
    MPICounter operator=(const MPICounter&) = delete;

    int next() {
        log("locking window...\n");
        int ret = MPI_Win_lock(MPI_LOCK_SHARED, 0, 0, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Win_lock", ret, "Failed to lock counter window");
        }

        log("window locked. Fetching and updaating...\n");
        const int one = 1;
        int val;
        ret = MPI_Fetch_and_op(&one, &val, MPI_INT, 0, 0, MPI_SUM, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_Fetch_and_op", ret, "failed to get and increment counter");
        }

        log("fetch and update done. unlocking...\n");
        ret = MPI_Win_unlock(0, win);
        if(ret != MPI_SUCCESS) {
            throw MPIException("MPI_win_unlock", ret, "Failed to unlock window for counter");
        }
        log("unlocked. returning %d\n", val);

        return val;
    }
private:
    int *local;
    MPI_Win win;

};
