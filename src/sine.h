#ifndef SINE_H
#define SINE_H
#include <math.h>

#define PI 3.14159

typedef struct Data {
    double* y;
    int size;
} Data;

void sine(Data *output, double a, double m, double b, double Ts, double t0, double t1) {
    int ARRAY_SIZE = (t1 - t0) / Ts;
    output->size = ARRAY_SIZE;
    output->y = malloc(sizeof(double) * ARRAY_SIZE);
    double t = t0;
    int i = 0;
    while (t0 <= t1 && i < ARRAY_SIZE) {
        double value = a * sin(m *t) + b;
        output->y[i] = value;
        t += Ts;
        i++;
    }
}

#endif
