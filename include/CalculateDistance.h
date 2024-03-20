#ifndef CALCULATEDISTANCE_H
#define CALCULATEDISTANCE_H
#include "MathFunc.h"

double root3Sigmoid(double x, double a, double b, double c, double d, double m)
{
  return 1. / (1. + exp(-b * (x - m))) * a * pow(x, 1. / 3.) + (1. - 1. / (1 + exp(-b * (x - m)))) * c + d;
}
double Sigmoid(double x, double a, double b, double c, double d, double N)
{
  return a / (b + exp(d * x - c)) + N;
}
double root3Fit(double x, double a, double b)
{
  return a * pow(x, 1. / 3.) + b;
}
double GetRandomDistance(double mt, double Mult)
{
  //double AParam = root3Sigmoid(Mult, 0.62613191, 0.47410809, 1.2921366, -0.49796381, 9.82385661); /*without shifted multiplicities 5 TeV + HM, medium anticorrelation*/
  double AParam = root3Sigmoid(Mult, 0.59780856, 0.58983674, 1.26480697, -0.48262374, 10.88898618); /*without shifted multiplicities 5 TeV + HM, medium anticorrelation*/
  double BParam = Sigmoid(Mult, -0.46688056, 0.39714293, 5.94004927, 0.32703649, 0.81354578);       /*without shifted multiplicities 5 TeV + HM, medium anticorrelation*/

  double sigma = AParam * pow(mt, -BParam); //5.7% global uncertainty on the A parameter (proton+deuteron)
  //This is to prevent the source size becomming 0 and the distance draw being stuck in an infite loop.
  if (sigma < 0.01) {
    sigma = 0.01;
  }
  double maxVal = GausSource(2 * sigma, sigma);
  TRandom3* RandomR = new TRandom3(0);
  while (true) {
    double dstGuess = gRandom->Uniform(0, 20 * sigma);
    double probguess = gRandom->Uniform(0, maxVal);
    if (GausSource(dstGuess, sigma) > probguess) {
      delete RandomR;
      return dstGuess;
    }
  }
  return 1;
}

double GetRandomDistanceSigma(double sigma)
{
  double max = GausSource(2 * sigma, sigma);
  TRandom3* RandomR = new TRandom3(0);
  while (true) {
    double dstGuess = gRandom->Uniform(0, 5. * sigma);
    double probguess = gRandom->Uniform(0, max);
    if (GausSource(dstGuess, sigma) > probguess) {
      delete RandomR;
      return dstGuess;
    }
  }
  return 2 * sigma;
}
double GetRandomDistanceSigma(double sigma, double Mult, int NParticles)
{ /// Use Doenigus parameterization as a proxy
  //double RRef = RDoen(Mult);
  //double RReal = RDoen(NParticles / 1.192);
  //sigma *= RReal / RRef;
  sigma *= pow(NParticles / 1.192, 1. / 3.) / (pow(Mult, 1. / 3.));
  double max = GausSource(2 * sigma, sigma);
  bool found = false;
  TRandom3* RandomR = new TRandom3(0);
  while (true) {
    double dstGuess = gRandom->Uniform(0, 5 * sigma);
    double probguess = gRandom->Uniform(0, max);
    if (GausSource(dstGuess, sigma) > probguess) {
      delete RandomR;
      return dstGuess;
    }
  }
  return 2 * sigma;
}

#endif