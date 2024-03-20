#ifndef LOADSETUP_H
#define LOADSETUP_H
#include "TH1.h"
// default 7 TeV MB angular distribution
TH1* LoadAngularDistribution(std::string type = "ALICE") /*options: ALICE, EPOS, Mixed*/
{
  /*Load correlation function measured by ALICE*/
  TFile* fin0 = new TFile("Ressources/AngularCorrelationsALICE.root");
  TH1* dphiCorrelatorHist = (TH1D*)fin0->Get("ALICE_dphi");
  dphiCorrelatorHist->SetDirectory(0);
  fin0->Close();
  /// We dont need the Delta Phi correlation function, but the same event distribution. we "regain" it by multiplying the Correlation function with the mixed event from EPOS 3, hoping its similar enough to the real data
  TFile* fin1 = new TFile("Ressources/AngularCorrelationsEPOS.root");
  TH1* dPhiMixedHist = (TH1D*)fin1->Get("ProjectionsDPhi_Mixed");
  dPhiMixedHist->SetDirectory(0);

  TH1* dPhiSameALICE = (TH1D*)dphiCorrelatorHist->Clone("dPhiSameALICE");
  for (int i = 1; i < 30; i++) {
    dPhiSameALICE->SetBinContent(i, dphiCorrelatorHist->GetBinContent(i) * dPhiMixedHist->GetBinContent(i));
  }
  dPhiSameALICE->SetDirectory(0);

  TH1* dPhiSameEPOS = (TH1D*)fin1->Get("ProjectionsDPhi_Same");
  dPhiSameEPOS->SetDirectory(0);
  fin1->Close();

  if (type == "ALICE") {
    return dPhiSameALICE;
  } else if (type == "EPOS") {
    return dPhiSameEPOS;
  }
  return dPhiMixedHist;
}

// When provided a multiplicity
TH1* LoadAngularDistribution(std::string type = "ALICE", double Multiplicity = 35.8) /*options: ALICE, EPOS, Mixed*/
{
  /*Correlation function parameterization from analysis note https://alice-notes.web.cern.ch/system/files/notes/analysis/1247/2023-03-12-DanielaRuggiano_AnalysisNotes.pdf*/
  double N = PowerLaw(Multiplicity, 0.29732885, 0.41280528, 0.11923053); // 0.41280528 * pow(Multiplicity, 0.11923053) + 0.29732885;
  // DeltaPhiSinFit(x,a,b,c,N)
  /// We dont need the Delta Phi correlation function, but the same event distribution. we "regain" it by multiplying the Correlation function with the mixed event from EPOS 3, hoping its similar enough to the real data
  TFile* fin1 = new TFile("Ressources/AngularCorrelationsEPOS.root");
  TH1* dPhiMixedHist = (TH1D*)fin1->Get("ProjectionsDPhi_Mixed");
  dPhiMixedHist->SetDirectory(0);

  TH1* dPhiSameALICE = (TH1D*)dPhiMixedHist->Clone("dPhiSameALICE");
  for (int i = 1; i < 30; i++) {
    dPhiSameALICE->SetBinContent(i, DeltaPhiSinFit(dPhiMixedHist->GetBinCenter(i), 1, PI / 2., 1., 1 - N) * dPhiMixedHist->GetBinContent(i));
  }
  dPhiSameALICE->SetDirectory(0);

  TH1* dPhiSameEPOS = (TH1D*)fin1->Get("ProjectionsDPhi_Same");
  dPhiSameEPOS->SetDirectory(0);
  fin1->Close();

  if (type == "ALICE") {
    return dPhiSameALICE;
  } else if (type == "EPOS") {
    return dPhiSameEPOS;
  }
  return dPhiMixedHist;
}

void LoadAngularDistributionPar(double Multiplicity, TH1* dPhiSameALICE) /*options: ALICE, EPOS, Mixed*/
{
  // dPhiSameALICE->Reset();
  /*Correlation function parameterization from analysis note https://alice-notes.web.cern.ch/system/files/notes/analysis/1247/2023-03-12-DanielaRuggiano_AnalysisNotes.pdf*/
  //double N = 0.32015428 * pow(Multiplicity, 0.10940088) + 0.43196852;
  // this uses the non-published data from private communication. from these the projections are made. FIS oeaks are excluded
  double N = 1 - 0.57250741 / pow(Multiplicity, 0.4075921);
  // DeltaPhiSinFit(x,a,b,c,N)
  /// We dont need the Delta Phi correlation function, but the same event distribution. we "regain" it by multiplying the Correlation function with the mixed event from EPOS 3, hoping its similar enough to the real data

  double dphi;
  for (int i = 1; i <= dPhiSameALICE->GetNbinsX(); i++) {
    dphi = dPhiSameALICE->GetBinCenter(i);
    // -1 |phi| + 6 is the fit to the mixed event distribution from EPOS
    dPhiSameALICE->SetBinContent(i, DeltaPhiSinFit(dphi, 1, PI / 2., 1., 1 - N) * (-1.09355734 * std::abs(dphi) + 5.99787991));
  }
  // return dPhiSameALICE;
}

TH2* LoadCoalescenceProbability(std::string wavefunction = "Argonne", int type = 2) /*options: wf so far only argoone, type =0,1,2 for S,D,S+D wave, note that S and D wave are not normalized*/
{
  TH2* ProbabilityHistogram;
  if (!strcmp(wavefunction.c_str(), "Argonne")) {
    /// This is the probability histogram for Argonne v18 with added S and D Wave. There are versions for only S wave and only D wave available in the same root file
    TFile* ProbabilityHistogramF = new TFile("Ressources/ArgonneProbability.root");
    switch (type) {
      case 0:
        ProbabilityHistogram = (TH2D*)ProbabilityHistogramF->Get("SWave");
        break;
      case 1:
        ProbabilityHistogram = (TH2D*)ProbabilityHistogramF->Get("DWave");
        break;
      case 2:
        ProbabilityHistogram = (TH2D*)ProbabilityHistogramF->Get("AddedSDWave");
        break;
      default:
        cout << "no wavefunction found!" << endl;
        abort();
        break;
    }
    ProbabilityHistogram->SetDirectory(0);
    ProbabilityHistogramF->Close();
  } else {
    cout << "no wavefunction found!" << endl;
    abort();
  }
  return ProbabilityHistogram;
}
#endif