#ifndef MULTIPLCIITY_H
#define MULTIPLCIITY_H
#include "MathFunc.h"

/*Parameterization for Erlang function for Multiplicity distributions*/
double LamErl(double Mult)
{
  return PowerLaw(Mult, 0., 0.72653, -0.23343); //0.72653 * pow(Mult, -0.23343);
}
double kErl(double Mult) //Modify Kappa using WildFunction if Multiplicity is below 10
{
  double kpar = PowerLaw(Mult, 0., 0.77544, 0.74669);
  if (Mult < 10.) {
    kpar /= WildFunction(Mult, 6.81004681e-03, -4.72757094e+01, 4.77281562e-01, 6.94118921e+00, 6.87846225e+00, 9.97532044e-01, -8.07181789e-01);
  }
  return kpar; //PowerLaw(Mult,0.,0.77544,0.74669);// 0.77544 * pow(Mult, 0.74669);
}
double Erlang(double Mult, double Kappa, double Lambda)
{
  return pow(Lambda, Kappa) * pow(Mult, Kappa - 1) * exp(-Lambda * Mult) / TMath::Gamma(Kappa);
}
int SampleFromErlang(double Mult)
{
  TRandom3* RandomR = new TRandom3(0);
  while (true) {
    int MultGuess = RandomR->Integer(150);
    double MultD = MultGuess;
    //if (Mult < 2) {
    if (MultGuess == 0) {
      MultD = 0.1;
    }

    double ProbGuess = RandomR->Uniform(0, 0.5);
    double k = kErl(Mult);   //Modify Kappa using WildFunction if Multiplicity is below 10
    double l = LamErl(Mult); //k/l=Mult
    double guess = Erlang(MultD * 0.98, k, k / Mult);
    if (guess > ProbGuess) {
      delete RandomR;
      return MultGuess;
    }
  }
}
double Yield(double Mult)
{
  return 0.06152 * pow(Mult, 0.95621); //this is used in the latest results, Central value
  //return 0.03267894 * 2 * pow(Mult, 0.94494328); //Systematic variation, upper bound fit
  //return 0.02885286 * 2 * pow(Mult, 0.96828117); //Systematic variation, lower bound fit
}

int GetNParticles(int type, double Multiplicity, double MultScaleFactor)
{
  int NParticles = 0;
  // Find the way to determine the multiplicity
  // First option: Fixed value, second option: Poissionian distributions
  if (type == 1) {
    NParticles = gRandom->Poisson(MultScaleFactor * (double)Multiplicity);
  } else if (type == 2) {
    NParticles = (int)floor(MultScaleFactor * (double)Multiplicity);
  } else if (type == 3) /*this one needs a bit more work*/
  {
    //    NParticles = SampleFromErlang(Multiplicity * (double)MultScaleFactor);
    NParticles = SampleFromErlang(Multiplicity * (double)MultScaleFactor);

    /*int ErlangInt = (int)ErlangDouble;
        double Decimals = ErlangDouble - ErlangInt;

        NParticles = (gRandom->Uniform(0, 1) < Decimals) ? (ErlangInt + 1) : (ErlangInt); // factor 1.02 introduced after the fact, most likely a binning effect when calculating the mean.
        */
  } else if (type > 4) {
    cout << "Please specify a Multiplicity type! chosing default type 1 (poissonian)" << endl;
    NParticles = gRandom->Poisson(MultScaleFactor * (double)Multiplicity);
  }
  return NParticles;
}
int GetNumberOfParticlesFromHistogram(TH1* MultDist)
{

  bool found = false;
  int iter = 0;
  double max = MultDist->GetMaximum();
  int nbins = MultDist->GetNbinsX();
  int maxMult = MultDist->GetBinCenter(nbins);
  while (found == false) {
    iter++;
    double PartsGuess = gRandom->Integer(maxMult);
    double probguess = gRandom->Uniform(0, max);
    if (MultDist->GetBinContent(PartsGuess) > probguess) {
      found = true;
      return (int)floor(PartsGuess);
    }
    if (iter > 1000) {
      return 10;
    }
  }

  return 10;
}
#endif