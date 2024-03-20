#ifndef DEUTERONPARAMETER_H
#define DEUTERONPARAMETER_H
#include "MathFunc.h"
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> GetDeuteronTarget(std::string ROOTName, std::string DataName, std::string ErrorName, std::string Mode)
{
  std::vector<double> DataVec;
  std::vector<double> ErrorVec;
  std::vector<double> edges;
  TFile* fin = new TFile(Form("Ressources/%s", ROOTName.c_str()));
  TH1* DataHist = (TH1D*)fin->Get(DataName.c_str());
  int nbins = DataHist->GetXaxis()->GetNbins();
  for (int i = 0; i < nbins; i++) {
    edges.push_back(DataHist->GetBinCenter(i + 1) - DataHist->GetBinWidth(i + 1) / 2);
    DataVec.push_back(DataHist->GetBinContent(i + 1));
    if (Mode == "Same") {
      ErrorVec.push_back(DataHist->GetBinError(i + 1));
    }
  }
  edges.push_back(DataHist->GetBinCenter(nbins) + DataHist->GetBinWidth(nbins) / 2);
  if (Mode == "Separate") {
    TH1* ErrorHist = (TH1D*)fin->Get(ErrorName.c_str());
    for (int i = 0; i < ErrorHist->GetXaxis()->GetNbins(); i++) {
      ErrorVec.push_back(ErrorHist->GetBinContent(i + 1));
    }
  }
  return std::make_tuple(DataVec, ErrorVec, edges);
}
double DeuteronErr(double pt)
{
  return 0.04276827 * pow(pt, 0.27378294);
}
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> GetDeuteronTargetParam(double mult)
{

  std::vector<double> DataVec;
  std::vector<double> edges;
  std::vector<double> ErrorVec;
  edges.push_back(0);
  /*    double dN = 2.30657910e-07 * pow(mult, 2) + 5.59888568e-05 * (mult)-6.73488159e-05;
    /// double n = 0.0268 * pow(mult, 2) - 0.25525 * (mult) + 5.44044;
    double n = -0.0162396427 * pow(mult - 3.32367935, 2) + 205.146797 * mult + 7.07789373 - 206.957870 * pow(mult, 0.996191333);
    double C = PowerLaw(mult,0,0.04165,0.72037);//0.04165 * pow(mult, 0.72037);
*/

  double dN = 0;
  double n = 0;
  double C = 0;
  if (mult == 18.5) {
    dN = 1.0763e-03;
    n = 17.19;
    C = 0.37;
  } else if (mult == 14.5) {
    dN = 8.1017e-04;
    n = 29.12;
    C = 0.34;
  } else if (mult == 11.9) {
    dN = 6.3748e-04;
    n = 14.33;
    C = 0.29;
  } else if (mult == 9.7) {
    dN = 4.9134e-04;
    n = 11.69;
    C = 0.26;
  } else if (mult == 7.8) {
    dN = 3.602e-04;
    n = 16.6;
    C = 0.25;
  } else if (mult == 6.3) {
    dN = 2.635e-04;
    n = 16.88;
    C = 0.24;
  } else if (mult == 5.2) {
    dN = 1.982e-04;
    n = 25.99;
    C = 0.24;
  } else if (mult == 3.9) {
    dN = 1.280e-04;
    n = 20.63;
    C = 0.21;
  } else if (mult == 2.4) {
    dN = 0.47528e-04;
    n = 29.85;
    C = 0.18;
  } else if (mult == 30.1) {
    dN = 1.84513e-03;
    n = 17.731;
    C = 0.48716;
  } else if (mult == 32.2) {
    dN = 1.9876e-03;
    n = 21.669;
    C = 0.5112;
  } else if (mult == 35.8) {
    dN = 2.2184e-03;
    n = 28.754;
    C = 0.55;
  } else {
    cout << "You messed up dawg! " << mult << endl;
  }

  int nBins = 40;
  double maxpT = 4;
  double stepsize = maxpT / nBins;
  /// now assume a binning: [20,0,4]
  for (int i = 1; i <= nBins; i++) {
    edges.push_back(i * stepsize);
    double binCenter = i * stepsize - stepsize / 2.;
    /// double LevyTsallis(double x, double Norm, double n, double T, double m)
    DataVec.push_back(LevyTsallis(binCenter, dN, n, C, 1.875612));
    ErrorVec.push_back(DeuteronErr(binCenter) * LevyTsallis(binCenter, dN, n, C, 1.875612));
    cout << LevyTsallis(binCenter, dN, n, C, 1.875612) << endl;
  }
  return std::make_tuple(DataVec, ErrorVec, edges);
}

double Chi2(std::vector<double> model, std::vector<double> sim, std::vector<double> Err, int N = 0)
{
  double x2 = 0;
  int ndf = 0;
  int lim = N;
  int si = model.size();
  if (N == 0) {
    lim = model.size();
  }
  for (int i = 0; i < si; i++) {
    if (model.at(i) > 0) {
      if (sim.at(i) > 0) {
        ndf++;
        x2 += pow((model.at(i) - sim.at(i)) / Err.at(i), 2);
      }
    }
  }
  return x2 / ndf;
}
#endif