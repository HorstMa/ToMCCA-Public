#ifndef KINEMATICS_H
#define KINEMATICS_H
#include "Multiplicity.h"
#include "MathFunc.h"
#include "TH1.h"
#include "TF1.h"
#define M_P 0.938272

//TRandom3* RandomH = new TRandom3(time(NULL));
double CalcpZ(double px, double py, double m, double rap)
{
  double sign = rap / std::abs(rap);
  return sign * sqrt(px * px + py * py + m * m) * 0.5 * sqrt(-2 + exp(-2 * rap) + exp(2 * rap));
}
double CalcpZFromEta(double px, double py, double eta)
{
  double sign = eta / std::abs(eta);
  return sign * sqrt(px * px + py * py) * 0.5 * sqrt(-2 + exp(-2 * eta) + exp(2 * eta));
}

double CalcEngy(double px, double py, double pz, double m)
{
  return sqrt(px * px + py * py + pz * pz + m * m);
}
double CalcDeltaP(TLorentzVector p1, TLorentzVector p2)
{
  TLorentzVector pPair = p1 + p2;
  TVector3 b = -pPair.BoostVector();
  p1.Boost(b);
  p2.Boost(b);
  return p1.Vect().Mag();
}

double GetRandompTMult(double Mult, int MultType, std::string HadMethod)
{

  double dNdy = Yield(Mult);                          /// Fit was done to p+pbar!
  //dNdy *= 2*(1. / (1 + exp(2.28 * Mult - 5.)) + 0.89); //suppression at low Mult due to edge effects
  //double n = PowerLaw(Mult,7.32306908,2.08554186,-1.); // This was used in the most recent run
  //double C = PowerLaw(Mult,0.,0.11693,0.36305); // This was used in the most recent run
  double n = 7.3;
  double C = PowerLaw(Mult, 0., 0.11263, 0.35086);
  if (MultType == 3) {
    if (HadMethod == "QuarkRecombination") {
      //recalibration Nov08 for Quark recombination
      C *= PowerLaw(Mult, 0.99534332, -0.73804899, -2.00693344);
    } else if (HadMethod == "StringFragmentation") {
      //C *= 0.92817133;
      C *= 1.;
    }
    //dNdy *= PowerLaw(Mult,0.9,0.41517883,-0.55482612);
    //C *= PowerLaw(Mult,1,-0.35195,-1.),//1 - 0.35195 / Mult; // 0.89706 * pow(Mult, 0.0275); // this is needed since the Erlang Multiplicity produces harder spectra since it has a larger tail compared to poisson
    //dNdy *= PowerLaw(Mult+3.07859,0.95740,0.55989,-1.);//0.95740 + 0.55989 / (Mult + 3.07859);
  }
  TF1* LevyFunc = new TF1("Levy", LevyTsallisTF1, 0.0, 8, 4);
  LevyFunc->SetParameters(dNdy, n, C, 0.938271);
  double MaxProb = LevyFunc->GetMaximum();

  while (true) {
    double pTGuess = gRandom->Uniform(0, 8);
    double probguess = gRandom->Uniform(0, MaxProb);

    double guess = LevyFunc->Eval(pTGuess);
    if (guess > probguess) {
      delete LevyFunc;
      return pTGuess;
    }
  }
  return 0;
}
#endif