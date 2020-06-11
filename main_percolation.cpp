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
  int time = 500;
  double dt = 0.01;
  int numberOfIterations = (int) time/dt;
  int N = 900;
  int L = 210;
  double v0 = 0.3;
  double dr = 0.5;
  int step = 100;
  double infectionRadius = 10.;
  int numberOfProbabilities = 101;
  int numberOfRunnings = 15;

  double *x = new double [(int)(N*numberOfIterations/step)];
  double *y = new double [(int)(N*numberOfIterations/step)];
  double *condition = new double [(int)(N*numberOfIterations/step)];

  int numberOfInfected;
  int infectionNumber;
  double *R = new double [(int) numberOfIterations/step];

  double infectedAfterRunningTime[numberOfProbabilities*numberOfRunnings] = {0};
  double probabilities[numberOfProbabilities];

  for(int i=0; i<numberOfProbabilities; i++){
    probabilities[i] = (double) i/100;
  }

  for(int middle=0; middle<numberOfRunnings; middle++)
  {
    cout << middle << endl;
    for(int prob=0; prob<numberOfProbabilities; prob++)
    {
      //cout << prob << endl;
      People people[N];
      for(int i=0; i<N; i++)
      {
        people[i] = People(L, i, N, v0, dr);
      }

      /*for(int i=0; i<20; i++)
      {
        int person = (double)rand()/RAND_MAX * N;
        people[person].setv0(5);
        people[person].setDr(0.05);
      }*/
      int person = (double)rand()/RAND_MAX * N;
      people[person].setCondition(1, 0);
      person = (double)rand()/RAND_MAX * N;
      people[person].setCondition(1, 0);


      for(int i=0; i<numberOfIterations; i++)
      {
          quarantine(people, x, y, condition, N, L, v0, dt, dr, i, step,
                     probabilities[prob]/2000, infectionRadius);
      }
      //writeToFileHDF5(x, y, condition, N, numberOfIterations, L, dt, step, R);

      int numberSubceptilbe = 0;
      for(int j=0; j<N; j++)
      {
        if(people[j].getCondition() == 0 )
        {
          numberSubceptilbe++;
        }
      }
      infectedAfterRunningTime[middle * numberOfProbabilities + prob] += 1. - ((double) numberSubceptilbe/N);
    }
  }

  ofstream outFile1("../data/infected_active_15_third.dat");


  for(int i = 0; i < numberOfProbabilities * numberOfRunnings; ++i) {
    outFile1 << infectedAfterRunningTime[i] << '\n';
  }

  return 0;
}








/*
numberOfInfected = 0;
infectionNumber = 0;
for(int j=0; j<numberOfPeople; j++)
{
  if(people[j].getCondition() == 1 && (i - people[j].getTime()) > 1500)
  {
    numberOfInfected ++;
    infectionNumber += people[j].getNumberOfInfected();
  }
}
*/

/*
if(i%step==0)
{
  if(numberOfInfected == 0)
  {
    R[(int) i/step] = 0;
  }
  else
  {
    R[(int) i/step] = (double) infectionNumber/numberOfInfected;
  }
}
*/
