#ifndef LOG_H
#define LOG_H

#include<chrono>

class StopWatch
{
public:
    void start() {
        if(is_running)
            stop();

        start_time = std::chrono::high_resolution_clock::now();
        is_running = true;
    }
    void stop() {
        auto stop_time = std::chrono::high_resolution_clock::now();
        duration = stop_time - start_time;
        is_running = false;
    }

    std::chrono::high_resolution_clock::duration get_elapsed_time() {
        return duration;
    }
    long get_elapsed_seconds() {
       return std::chrono::duration_cast<std::chrono::seconds>(duration).count();
    }
    long get_elapsed_milliseconds() {
       return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
    }
    long get_elapsed_microseconds() {
       return std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    }

private:
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::duration duration;
    bool is_running;
};

int init_logging();
void LOG(const char *fmt, ...);
#endif

struct SINGLE_FILE;
SINGLE_FILE *single_fopen(const char *path, const char *mode);
int fclose(SINGLE_FILE *stream);
int fprintf(SINGLE_FILE *stream, const char *format, ...);
int single_printf(const char *fmt, ...);


