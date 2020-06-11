#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <chrono>

#include "solver.h"
#include "utilities.h"

using namespace std;


int main(){
  srand(time(NULL));
  int time = 250;
  double dt = 0.01;
  int numberOfIterations = (int) time/dt;
  int N = 900;
  int L = 210;
  double v0 = 2;
  double dr = 0.01;
  int step = 10;
  double p = 0.4;
  double infectionRadius = 12.;
  double probQuarantine = 0.7;
  int cases;


  double *x = new double [(int)(N*numberOfIterations/step)];
  double *y = new double [(int)(N*numberOfIterations/step)];
  double *condition = new double [(int)(N*numberOfIterations/step)];

  int numberOfInfected;
  int infectionNumber;
  double *R = new double [(int) numberOfIterations/step];


  People people[N];
  for(int i=0; i<N; i++)
  {
    people[i] = People(L, i, N, v0, dr);
  }


  int person = (double)rand()/RAND_MAX * N;
  people[person].setCondition(1, 0);
  person = (double)rand()/RAND_MAX * N;
  people[person].setCondition(1, 0);

  for(int i=0; i<numberOfIterations; i++)
  {
    cases = 0;
    for(int j=0; j<N; j++)
    {
      if(people[j].getCondition() != 0)
      {
       cases ++;
      }
    }
    quarantine(people, x, y, condition, N, L, v0, dt, dr, i, step, p/2000, infectionRadius, probQuarantine, cases);


    if(i%step==0)
    {
      numberOfInfected = 0;
      infectionNumber = 0;
      for(int j=0; j<N; j++)
      {
        if(people[j].getCondition() == 1 && (i - people[j].getTime()) > 1500)
        {
          numberOfInfected ++;
          infectionNumber += people[j].getNumberOfInfected();
        }
      }

      if(numberOfInfected == 0)
      {
        R[(int) i/step] = 0;
      }
      else
      {
        R[(int) i/step] = (double) infectionNumber/numberOfInfected;
      }
    }
  }

  writeToFileHDF5(x, y, condition, N, numberOfIterations, L, dt, step, R);

  return 0;
}
