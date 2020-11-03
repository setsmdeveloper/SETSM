#include <cstring>
#include <cstdio>
#include <chrono>

#ifdef BUILDMPI
#include <mpi.h>
#endif

#include <cstdarg>

struct SINGLE_FILE {
    FILE *fp;
};

static bool started = 0;
static long start_time;
static int rank;

int init_logging() {
    unsigned long now =
        std::chrono::duration_cast<std::chrono::milliseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();

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
        printf("WARNING: logging not initialized, defaulting to printf\n");
        va_list ap;
        va_start(ap, fmt);
        vprintf(fmt, ap);
        va_end(ap);
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

SINGLE_FILE *single_fopen(const char *path, const char *mode) {
#ifdef BUILDMPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank != 0)
        return nullptr;
#endif
    FILE *fp = fopen(path, mode);
    return new SINGLE_FILE{fp};
}

int fclose(SINGLE_FILE *stream) {
#ifdef BUILDMPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank != 0)
        return 0;
#endif
    int ret = fclose(stream->fp);
    delete stream;
    return ret;
}

int fprintf(SINGLE_FILE *stream, const char *format, ...) {
#ifdef BUILDMPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank != 0)
        return 0;
#endif
    va_list ap;
    va_start(ap, format);
    int ret = vfprintf(stream->fp, format, ap);
    va_end(ap);
    return ret;
}

int single_printf(const char *fmt, ...) {
#ifdef BUILDMPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank != 0)
        return 0;
#endif
    va_list ap;
    va_start(ap, fmt);
    int ret = vprintf(fmt, ap);
    va_end(ap);
    return ret;
}
