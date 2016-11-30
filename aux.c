//
// Created by mm on 11/30/16.
//

#include <time.h>
#include "aux.h"


#ifdef __linux__
struct timespec __startElasped, __now, __endElasped;

void startTimer() {
    clock_gettime(CLOCK_REALTIME, &__startElasped);
}

void endTimer() {
    clock_gettime(CLOCK_REALTIME, &__endElasped);
}


// Returns total clock time in seconds
double getElaspedTime() {
    clock_gettime(CLOCK_REALTIME, &__now);
    __time_t tmp_seconds = (__now.tv_sec - __startElasped.tv_sec);
    double seconds = (double) tmp_seconds;
    __time_t tmp_nano = (__now.tv_nsec - __startElasped.tv_nsec);
    double nano = tmp_nano / 1e9;
    if (nano > 1e9) {
        tmp_nano = (__startElasped.tv_nsec - __now.tv_nsec);
        nano = tmp_nano / 1e9;
    }
    double time = seconds + nano;
    return time;
}

#else
void startTimer() {}

void endTimer() {}

double getElaspedTime() {
    return 0;
}
#endif

