//
// Created by mm on 11/30/16.
//

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "aux.h"


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

