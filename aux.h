//
// Created by mm on 11/30/16.
//

#ifndef MOL_AUX_H
#define MOL_AUX_H


typedef struct _stopwatch {
    struct timespec startElapsed;
    struct timespec now;
    struct timespec endElapsed;
} stopwatch_t;
//struct timespec __startElasped, __now, __endElasped;


void startTimer(stopwatch_t *stopwatch);
void endTimer(stopwatch_t *stopwatch);
double getElaspedTime(stopwatch_t *stopwatch);
double getElaspedns(stopwatch_t *stopwatch);

//void startTimer();
//void endTimer();
//double getElaspedTime();

void *emalloc(size_t numBytes);
double **new_2d_double_array(long numRows, long numCols);
FILE *e_fopen(const char *__filename, const char *__modes);
void die(const char *format, ...);
void wipe(FILE *file, int n);
void progress_bar(long idx, long max);

#endif //MOL_AUX_H
