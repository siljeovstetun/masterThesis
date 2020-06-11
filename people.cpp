#include <iostream>
#include <cmath>

#include "people.h"

using namespace std;

People::People(){}

People::People(int L, int i, int N, double v, double d)
{
  condition = 0;
  xPosition = L / ((int) sqrt(N)) * (i % (int) sqrt(N));
  yPosition = L / ((int) sqrt(N)) * (int) (i / (int) sqrt(N));
  angle = -M_PI + (double)rand()/RAND_MAX * 2 * M_PI;
  v0 = v;
  dr = d;
  numberOfInfected = 0;
}


void People::setCondition(int cond, int time)
{
  condition = cond;

  if(cond == 1) // condition = 1 means that the person get infected
  {
    infectedTime = time;
  }
}

void People::setPosition(double x, double y, double a){
  xPosition = x;
  yPosition = y;
  angle =  a;
}

void People::setNumberOfInfected(){
  numberOfInfected++;
}

void People::setv0(double v)
{
  v0 = v;
}

void People::setDr(double d)
{
  dr = d;
}


int People::getCondition()
{
    return condition;
}

double People::getxPosition()
{
    return xPosition;
}

double People::getyPosition()
{
    return yPosition;
}

double People::getAngle()
{
    return angle;
}

int People::getTime()
{
  return infectedTime;
}

double People::getv0()
{
  return v0;
}

double People::getDr()
{
  return dr;
}

int People::getNumberOfInfected(){
  return numberOfInfected;
}
