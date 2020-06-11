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
  int time = 400;
  double dt = 0.01;
  int numberOfIterations = (int) time/dt;
  int N = 900;
  int L = 90;
  double v0 = 0.3;
  double dr = 0.5;
  int step = 100;
  double infectionRadius = 10.;
  int numberOfProbabilities = 101;
  int numberOfRunnings = 10;

  double *x = new double [(int)(N*numberOfIterations/step)];
  double *y = new double [(int)(N*numberOfIterations/step)];
  double *condition = new double [(int)(N*numberOfIterations/step)];

  int numberOfInfected;
  int infectionNumber;

  double pi[500] = {0};
  double eta[50];

  for(int i=0; i<50; i++){
    eta[i] = (double) i/50;
  }

  for(int middle=0; middle<numberOfRunnings; middle++)
  {
    cout << middle << endl;
    for(int run=0; run<50; run++)
    {
      People people[N];
      for(int i=0; i<N; i++)
      {
        people[i] = People(L, i, N, v0, eta[run]);
      }

      for(int i=0; i<numberOfIterations; i++)
      {
        order(people, x, y, condition, N, L, v0, dt, dr, i, step,
                     0.8, infectionRadius);
      }
      double c = 0;
      double s = 0;
      for(int i=0; i<N; i++)
      {
        s += sin(people[i].getAngle());
        c += cos(people[i].getAngle());
      }
      pi[middle*50 + run] = sqrt(c*c + s*s)/N;
    }
  }

  ofstream outFile1("../data/average_order_10_tenth.dat");
  ofstream outFile2("../data/eta.dat");


  for(int i = 0; i < 500; ++i) {
    outFile1 << pi[i] << '\n';
  }

  for(int i = 0; i < 50; ++i) {
    outFile2 << eta[i] << '\n';
  }

  return 0;
}
