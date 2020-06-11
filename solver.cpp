#include <iostream>
#include <cmath>
#include <omp.h>

#include "people.h"

using namespace std;

double newCoordinate(double oldCoordinate, double changeCoordinate, int L){
  double newCoordinate = oldCoordinate + changeCoordinate;

  if(newCoordinate > L){
    newCoordinate -= L;
  }
  else if (newCoordinate < 0) {
    newCoordinate += L;
  }
  return newCoordinate;
}


void movePeople(People *people, double *x, double *y, double *condition,
                int numberOfPeople, int L, double v0, double dt, int it, int step,
                double probabiliySerious, double probabilityInfected, double infectionRadius)
{
  double dx;
  double dy;
  double xi;
  double xj;
  double yi;
  double yj;
  double anglei;
  double dist;
  double force = 0;
  double forcex;
  double forcey;
  double newx;
  double newy;
  double newAngle;
  double serious;
  double infected;

  #pragma omp parallel for private(dx, dy, force, forcex, forcey, dist, newx, newy, newAngle, xi, xj, yi, yj, anglei)
  for(int i=0; i<numberOfPeople; i++)
  {
      xi = people[i].getxPosition();
      yi = people[i].getyPosition();
      anglei = people[i].getAngle();
      forcex = 0;
      forcey = 0;

      if(people[i].getCondition() == 1 && (it - people[i].getTime()) > 40000){
        people[i].setCondition(2, 0);
      }
      if(people[i].getCondition() == 3 && (it - people[i].getTime()) > 25000){
        people[i].setCondition(2, 0);
      }
      if(people[i].getCondition() == 1 && (it - people[i].getTime()) == 4000)
      {
        serious = (double)rand()/RAND_MAX;
        if(serious < probabiliySerious)
        {
          people[i].setCondition(3, 0);
        }
      }

      for(int j=0; j<numberOfPeople; j++)
      {
        if(i!=j)
        {
          xj = people[j].getxPosition();
          yj = people[j].getyPosition();
          dx = xj - xi;
          dy = yj - yi;

          if(dx > L/2. || dx < -L/2.)
          {
            dx = -(L - abs(xj - xi)) * (xj - xi) / abs(xj - xi);
          }

          if(dy > L/2. || dy < -L/2)
          {
            dy = -(L - abs(yj - yi)) * (yj - yi) / abs(yj - yi);
          }

          dist = sqrt(dx*dx+dy*dy);

          if(dist < infectionRadius)
          {
            if(people[i].getv0() != 0){
              force = -10 * 0.2227 / (dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist);
              forcex += force * dx / dist;
              forcey += force * dy / dist;
            }

            if(people[j].getCondition() == 1 && dist < infectionRadius && people[i].getCondition() == 0)
            {
              //cout << infectionRadius << "     " << dist;
              infected = (double)rand()/RAND_MAX;
              if(infected < probabilityInfected)
              {
                people[i].setCondition(1, it);
                people[j].setNumberOfInfected();
                //cout << "and here" << endl;
              }
            }
          }
        }
      }
      newx = newCoordinate(xi, people[i].getv0() * cos(anglei) * dt + forcex * dt, L);
      newy = newCoordinate(yi, people[i].getv0() * sin(anglei) * dt + forcey * dt, L);
      newAngle = anglei + (-M_PI + (double)rand()/RAND_MAX * 2 * M_PI);
      people[i].setPosition(newx, newy, newAngle);
    }

  if(it%step==0)
  {
    for(int person=0; person<numberOfPeople; person++)
    {
      x[(int) it/step * numberOfPeople + person] = people[person].getxPosition();
      y[(int) it/step * numberOfPeople + person] = people[person].getyPosition();
      condition[(int) it/step * numberOfPeople + person] = people[person].getCondition();
    }
  }
}


void quarantine(People *people, double *x, double *y, double *condition,
                int numberOfPeople, int L, double v0, double dt, double dr, int it, int step,
              double probabilityInfected, double infectionRadius, double probQuarantine, int cases)
{
  double dx;
  double dy;
  double xi;
  double xj;
  double yi;
  double yj;
  double anglei;
  double dist;
  double force;
  double forcex;
  double forcey;
  double newx;
  double newy;
  double newAngle;
  double infected;
  double diffAngle;
  int numOfNeighbours;
  double prob;


  #pragma omp parallel for private(dx, dy, xi, xj, yi, anglei, dist, force, forcex, forcey, newx, newy, newAngle, infected, diffAngle, numOfNeighbours, prob)
  for(int i=0; i<numberOfPeople; i++)
  {
      xi = people[i].getxPosition();
      yi = people[i].getyPosition();
      anglei = people[i].getAngle();
      forcex = 0;
      forcey = 0;
      diffAngle = 0;
      numOfNeighbours = 0;


      if(people[i].getCondition() == 1 && (it - people[i].getTime()) > 2000){
        people[i].setCondition(2, 0);
      }

      if(people[i].getCondition() == 1 && (it - people[i].getTime()) == 400 && cases > 50){
        prob = (double)rand()/RAND_MAX;
        if(prob < probQuarantine)
        {
          people[i].setCondition(3, 0);
        }
      }

      if(people[i].getCondition() == 3 && (it - people[i].getTime()) > 2000){
        people[i].setCondition(2, 0);
      }

      for(int j=0; j<numberOfPeople; j++)
      {
        if(i!=j)
        {
          xj = people[j].getxPosition();
          dx = xj - xi;

          if(abs(dx) < infectionRadius)
          {
            yj = people[j].getyPosition();
            dy = yj - yi;
            if(abs(dy) < infectionRadius)
            {
              dist = sqrt(dx*dx+dy*dy);
              if(dist < infectionRadius)
              {

                if(dx > L/2. || dx < -L/2.)
                {
                  dx = -(L - abs(xj - xi)) * (xj - xi) / abs(xj - xi);
                }

                if(dy > L/2. || dy < -L/2)
                {
                  dy = -(L - abs(yj - yi)) * (yj - yi) / abs(yj - yi);
                }

                force =  -2.227 / (dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist);
                //force -= 1 * exp(-dist/2);// * (dist<12);
                //force += 4 * exp(-dist);// * (dist<8);
                forcex += force * dx / dist;
                forcey += force * dy / dist;


                if(people[j].getCondition() == 1 && dist < infectionRadius && people[i].getCondition() == 0)
                {
                  infected = (double)rand()/RAND_MAX;
                  if(infected < probabilityInfected)
                  {
                    people[i].setCondition(1, it);
                    people[j].setNumberOfInfected();
                  }
                }
                if(dist<4){
                  //diffAngle += 3/2 * sin(2*people[j].getAngle() - 2*anglei);
                  numOfNeighbours ++;
                }
              }
            }
          }
        }
      }

      newx = newCoordinate(xi, people[i].getv0() * cos(anglei) * dt + forcex * dt, L);
      newy = newCoordinate(yi, people[i].getv0() * sin(anglei) * dt + forcey * dt, L);
      newAngle = anglei + sqrt(2*people[i].getDr()*dt)*(-M_PI + (double)rand()/RAND_MAX * 2 * M_PI);
      if(numOfNeighbours != 0)
      {
        newAngle += diffAngle * dt / numOfNeighbours;
      }
      people[i].setPosition(newx, newy, newAngle);

    }

  if(it%step==0)
  {
    for(int person=0; person<numberOfPeople; person++)
    {
      x[(int) it/step * numberOfPeople + person] = people[person].getxPosition();
      y[(int) it/step * numberOfPeople + person] = people[person].getyPosition();
      condition[(int) it/step * numberOfPeople + person] = people[person].getCondition();
    }
  }

}




void order(People *people, double *x, double *y, double *condition,
                int numberOfPeople, int L, double v0, double dt, double dr, int it, int step,
                double probabilityInfected, double infectionRadius)
{
  double dx;
  double dy;
  double xi;
  double xj;
  double yi;
  double yj;
  double anglei;
  double dist;
  double force;
  double forcex;
  double forcey;
  double newx;
  double newy;
  double newAngle;
  double infected;
  double diffAngle;
  int numOfNeighbours;
  double prob;


  #pragma omp parallel for private(dx, dy, xi, xj, yi, anglei, dist, force, forcex, forcey, newx, newy, newAngle, infected, diffAngle, numOfNeighbours, prob)
  for(int i=0; i<numberOfPeople; i++)
  {
      xi = people[i].getxPosition();
      yi = people[i].getyPosition();
      anglei = people[i].getAngle();
      forcex = 0;
      forcey = 0;
      diffAngle = 0;
      numOfNeighbours = 0;


      for(int j=0; j<numberOfPeople; j++)
      {
        if(i!=j)
        {
          xj = people[j].getxPosition();
          dx = xj - xi;

          if(abs(dx) < infectionRadius)
          {
            yj = people[j].getyPosition();
            dy = yj - yi;
            if(abs(dy) < infectionRadius)
            {
              dist = sqrt(dx*dx+dy*dy);
              if(dist < infectionRadius)
              {

                if(dx > L/2. || dx < -L/2.)
                {
                  dx = -(L - abs(xj - xi)) * (xj - xi) / abs(xj - xi);
                }

                if(dy > L/2. || dy < -L/2)
                {
                  dy = -(L - abs(yj - yi)) * (yj - yi) / abs(yj - yi);
                }

                force =  -2.227 / (dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist*dist);
                //force -= 1 * exp(-dist/2);// * (dist<12);
                //force += 4 * exp(-dist);// * (dist<8);
                forcex += force * dx / dist;
                forcey += force * dy / dist;

                if(dist<4){
                  //diffAngle += 3/2 * sin(2*people[j].getAngle() - 2*anglei);
                  diffAngle += 3/2 * sin(people[j].getAngle() - anglei);
                  numOfNeighbours ++;
                }
              }
            }
          }
        }
      }

      newx = newCoordinate(xi, people[i].getv0() * cos(anglei) * dt + forcex * dt, L);
      newy = newCoordinate(yi, people[i].getv0() * sin(anglei) * dt + forcey * dt, L);
      //newAngle = anglei + sqrt(2*people[i].getDr()*dt)*(-M_PI + (double)rand()/RAND_MAX * 2 * M_PI);
      newAngle = anglei + people[i].getDr()*sqrt(dt)*(-M_PI + (double)rand()/RAND_MAX * 2 * M_PI);
      if(numOfNeighbours != 0)
      {
        newAngle += diffAngle * dt / numOfNeighbours;
      }
      people[i].setPosition(newx, newy, newAngle);

    }

  if(it%step==0)
  {
    for(int person=0; person<numberOfPeople; person++)
    {
      x[(int) it/step * numberOfPeople + person] = people[person].getxPosition();
      y[(int) it/step * numberOfPeople + person] = people[person].getyPosition();
      condition[(int) it/step * numberOfPeople + person] = people[person].getCondition();
    }
  }

}
