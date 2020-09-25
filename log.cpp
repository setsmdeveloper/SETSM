#include <cstring>
#include <cstdio>
#include <chrono>

#ifdef BUILDMPI
#include <mpi.h>
#endif

#include <cstdarg>

static bool started = 0;
static long start_time;
static int rank;

static int _init_logging(long now) {
#ifdef BUILDMPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == 0) {
        start_time = now;
    }
    int ret = MPI_Bcast(&start_time, 1,MPI_LONG, 0, MPI_COMM_WORLD);
    if(ret != MPI_SUCCESS) {
        printf("FAILED to init logging!\n");
        return 1;
    }
#else
    rank = 0;
    start_time = now;
#endif
    started = 1;
    return 0;
}

void LOG(const char *fmt, ...) {
    unsigned long now = 
        std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();

    if(!started) {
        if(_init_logging(now))
            return;
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

void init_logging() {
    LOG("logging initialized\n");
}
