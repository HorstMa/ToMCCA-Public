#ifndef ANGULARCORRELATIONS_H
#define ANGULARCORRELATIONS_H

double GetRandomdPhi(TH1* distHist)
{
  double DistMax = distHist->GetMaximum();
  double DistMin = distHist->GetMinimum();
  double width = distHist->GetBinWidth(1);
  double FirstBin = distHist->GetBinCenter(1);
  double LastBin = distHist->GetBinCenter(distHist->GetNbinsX() + 1);
  while (true) {
    double PhiGuess = gRandom->Uniform(FirstBin - width / 2, LastBin + width / 2);
    double probguess = gRandom->Uniform(0, DistMax);
    if (distHist->GetBinContent(distHist->FindBin(PhiGuess)) > probguess) {
      return PhiGuess;
    }
  }
  return 0;
}

double MixedFunction(double dphi, double a, double b, double N1, double N2)
{
  return N2 * exp(-pow(dphi, 2) / pow(a, 2)) + (N1 * abs(dphi) + b);
}

double DeltaPhiFunction(double DPhi, double Mult, double pt)
{
  //double N=1-0.57250742*pow(Mult,-0.4075921);
  // Retune Nov10 '23
  //double N=PowerLaw(Mult,-0.04756259,-0.45095972,-0.41476098); some old stuff?
  // Jan 18
  //N*sin(a*x-b)+c
  //double N = PowerLaw(Mult, 0., 0.35661963, -0.32215926) + 0.056614 * pt - 0.08437; //make it pt dependent
  //double N = PowerLaw(Mult, 0., 0.36755048, -0.35345157) * (0.33823262 * pt + 0.49543104); //make it pt dependent

  //N = 0.24099052 * exp(-0.02763867 * Mult) * (0.33823262 * pt + 0.49543104);

  double N = 0.16633569; //fuck Multiplicity and pT dependence just use 20-40%
  double c = PowerLaw(Mult, 0.86744736, 0.08185635, 0.13604509);
  c = 1.00;
  //double c = 0.0010407 * Mult + 0.96645;

  //double N = 0.16738;
  //double c = 1;

  // 0.865 from the shift pp->pL
  return (DeltaPhiSinFit(DPhi, 1., PI / 2., c, 0.865 * N)) * MixedFunction(DPhi, 0.85747601, 5.7197689, -0.89304215, 0.26041062); //(-1.09355734 * abs(DPhi) + 5.99787991); //0.865 for pp->pL
}

double GetRandomdPhiPar(double Multiplicity, double pt)
{

  double FirstBin = -1.41;
  double LastBin = 4.88;
  while (true) {
    double PhiGuess = gRandom->Uniform(FirstBin, LastBin);
    double probguess = gRandom->Uniform(0, 6.3);
    if (DeltaPhiFunction(PhiGuess, Multiplicity, pt) > probguess) {
      return PhiGuess;
    }
  }
  return 0;
}
double GetRandomdEta(double EtaP1)
{
  return gRandom->Uniform(-1, 1);
}

#endif