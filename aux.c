//
// Created by mm on 11/30/16.
//

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include "aux.h"

void die(const char *format, ...) {
    // DIE DIE DIE
    va_list argptr;
    va_start(argptr, format);
    vfprintf(stderr, format, argptr);
    va_end(argptr);
//    fprintf(stderr, "%s x_x", msg);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void *emalloc(size_t numBytes){ /* Error checked malloc */
    void *buffer = malloc(numBytes);
    if(!buffer){
        fprintf(stderr, "Fatal error: malloc error. Out of memory. \n  x_X\n");
        exit(EXIT_FAILURE);
    }
    return buffer;
}

double **new_2d_double_array(long numRows, long numCols) {
    long i;
    size_t size = sizeof(double *);
    double **array = (double **) emalloc(size * numRows);
    for (i = 0; i < numRows; i++) {
        array[i] = (double *) emalloc(size * numCols);
    }
    return array;
}

FILE *e_fopen(const char *__filename, const char *__modes) {
    /* Error checked file open */
    FILE *infilep = fopen(__filename, __modes);

    if (infilep == NULL) {
        die("Fatal error! Can't open file!\n");
    } else {
        return infilep;
    }
}

#ifdef __linux__

void startTimer(stopwatch_t *stopwatch) {
    clock_gettime(CLOCK_REALTIME, &(stopwatch->startElapsed));
}

void endTimer(stopwatch_t *stopwatch) {
    clock_gettime(CLOCK_REALTIME, &(stopwatch->endElapsed));
}

double getElaspedTime(stopwatch_t *stopwatch) {
    struct timespec now, startElapsed;
    startElapsed = stopwatch->startElapsed;
    clock_gettime(CLOCK_REALTIME, &now);
    __time_t tmp_seconds = (now.tv_sec - startElapsed.tv_sec);
    double seconds = (double) tmp_seconds;
    __time_t tmp_nano = (now.tv_nsec - startElapsed.tv_nsec);
    double nano = tmp_nano / 1e9;

    if (nano > 1e9) {
        tmp_nano = (startElapsed.tv_nsec - now.tv_nsec);
        nano = tmp_nano / 1e9;
    }
    return seconds + nano;
}

double getElaspedns(stopwatch_t *stopwatch) {
    struct timespec now, startElapsed;
    startElapsed = stopwatch->startElapsed;
    clock_gettime(CLOCK_REALTIME, &now);
    __time_t tmp_seconds = (now.tv_sec - startElapsed.tv_sec);
    double seconds = (double) tmp_seconds;
    __time_t tmp_nano = (now.tv_nsec - startElapsed.tv_nsec);
    double nano = tmp_nano / 1e9;
    if (nano > 1e9) {
        tmp_nano = (startElapsed.tv_nsec - now.tv_nsec);
        nano = tmp_nano / 1e9;
    }
    return 1e9 * (seconds + nano);
}

//
//void startTimer() {
//    clock_gettime(CLOCK_REALTIME, &__startElasped);
//}
//
//void endTimer() {
//    clock_gettime(CLOCK_REALTIME, &__endElasped);
//}
//
//
//// Returns total clock time in seconds
//double getElaspedTime() {
//    clock_gettime(CLOCK_REALTIME, &__now);
//    __time_t tmp_seconds = (__now.tv_sec - __startElasped.tv_sec);
//    double seconds = (double) tmp_seconds;
//    __time_t tmp_nano = (__now.tv_nsec - __startElasped.tv_nsec);
//    double nano = tmp_nano / 1e9;
//    if (nano > 1e9) {
//        tmp_nano = (__startElasped.tv_nsec - __now.tv_nsec);
//        nano = tmp_nano / 1e9;
//    }
//    double time = seconds + nano;
//    return time;
//}

#else
void startTimer() {}

void endTimer() {}

double getElaspedTime() {
    return 0;
}
#endif

void wipe(FILE *file, int n) {
    for (; n>0; n--) {
        fprintf(file, "\033[A\033[2K");
    }
}

void progress_bar(long idx, long max, stopwatch_t *stopwatch) {
    double elapsed = getElaspedTime(stopwatch);
    const char SPINNER[5] = "\\|/-";
    int WIDTH = 40;
    int i;
    int subdiv_print = 10;
    double perc = (double) idx / (double) max;
    int fullbars = (int) (perc * WIDTH);
    int halfbar = ((int) (perc * WIDTH*2 )) - 2*fullbars;
    double remain_sec = (elapsed / perc) * (1 - perc);

    if (idx % subdiv_print == 0) {
    wipe(stderr, 1);
    fprintf(stderr, "%c [", SPINNER[((idx/subdiv_print))%4]);
    for (i = 0; i < fullbars; i++) fprintf(stderr, "=");
    if (halfbar) fprintf(stderr, "-");
    for (i = 0; i < (WIDTH-fullbars-halfbar); i++) fprintf(stderr, " ");

    fprintf(stderr, "] % 2.1lf%% % 5ld/%ld", perc*100, idx, max);
    fprintf(stderr, " % 3.0lf:%02.0lf|% 3.0lf:%02.0lf\n",
            remain_sec/60, fmod(remain_sec, 60),
            elapsed/60, fmod(elapsed, 60));
    }
}

