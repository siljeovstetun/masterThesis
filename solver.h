#include "people.h"

void movePeople(People *people, double *x, double *y, double *condition,
                int numberOfPeople, int L, double v0, double dt, int it, int step,
                double probabiliySerious, double probabilityInfected, double infectionRadius);

void quarantine(People *people, double *x, double *y, double *condition,
                int numberOfPeople, int L, double v0, double dt, double dr, int it, int step,
                double probabilityInfected, double infectionRadius, double probQuarantine, int cases);


void order(People *people, double *x, double *y, double *condition,
                int numberOfPeople, int L, double v0, double dt, double dr, int it, int step,
                double probabilityInfected, double infectionRadius);
