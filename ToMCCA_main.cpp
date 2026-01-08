#include "TH1.h"
#include "TH2.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <math.h>

#include "include/MathFunc.h"
#include "include/AngularCorrelations.h"
#include "include/CalculateDistance.h"
#include "include/DeuteronParameter.h"
//#include "include/Doenigus.h"
#include "include/HadronizationModels.h"
#include "include/Kinematics.h"
#include "include/LoadSetup.h"
#include "include/Multiplicity.h"

#define MP 0.938272

void ToMCCA_main(int modeI = 1, double Multiplicity = 310.0, int MultiplicityType = 0, int64_t NEvents = 1000000, std ::string HadronizationMode = "QuarkRecombination")
{

  /* initiallize random seed */
  TRandom3* RandomR = new TRandom3(0);

  /*These two parameters are very experimental, use with Caution!!!!!!*/
  double IsospinRatio = 1;    // more Protons --> smaller Value, no Neutrons --> 0 >>>>>>> N/P // this can lead to unexpected behaviour! use with caution!!
  double AntimatterRatio = 1; // more Matter --> smaller Value, no antimatter --> 0  >>>>> AM/M
  /*#################################################*/
  double MultScaleFactor = 1.192;
  Multiplicity /= 10.; // divide by 10 to get correct Multiplicity with one digit precision
  cout << Multiplicity << endl;
  bool Devt = false;
  /*determine which mode we are in*/
  std::string mode;
  switch (modeI) {
    case 1:
      mode = "Default";
      break;
  }
  /*Load physical input*/
  TH2* CoalescenceProbability = LoadCoalescenceProbability("Argonne", 2);
  
  /*Debug: load Multiplicity distribution from file (EPOS HM)*/
  //// EPOS HM

  TFile* fin3 = new TFile("Ressources/EPOS_HM_MultDistdNdy.root");
  TCanvas* MultCanvas = (TCanvas*)fin3->Get("c_5592926ced90_projection_6301");
  TH1* MultHist = (TH1D*)MultCanvas->GetPrimitive("slice_px_of_dNdyV0");
  MultHist->SetDirectory(0);
  fin3->Close();
  cout << "starting" << endl;
  auto start = std::chrono::high_resolution_clock::now();

  TH1* Files = new TH1D("Files", "Files", 2, 0, 2); /*this histogram is interesting when merging files with hadd later on so we know how many files we merged*/
  Files->Fill(1);
  double distance = 0;

  /*define Histograms. Remember to also add them to the output!*/
  TH1* ptDeuteronArgHist = new TH1D("ptDeuteronArg", "ptDeuteronArg", 40, 0, 4);
  TH1* ptDeuteronArgHistFine = new TH1D("ptDeuteronArgFine", "ptDeuteronArgFine", 1000, 0, 10);
  TH1* DeuteronB2Hist = new TH1D("DeuteronB2Hist", "DeuteronB2Hist;pt/A[GeV];dN/dydpt", 50, 0, 8);

  TH1* ptNHist = new TH1D("ptN", "ptN", 100, 0, 8);
  TH1* ptHist = new TH1D("pt", "pt", 100, 0, 8);
  TH1* DoverPMult = new TH1D("DoverPMult", "d/p;<Nch>", 101, -0.5, 100.5);
  TH1* B2pt = new TH1D("B2pt", "B2;pT[GeV];B2", 100, 0, 8);

  /*Debug output*/
  TH1* NParticlesH = new TH1D("NParticles", "NParticles;dN/dy", 201, -0.5, 200.5);
  TH1* DeuteronEvtMult = new TH1D("DeuteronEvtMult", "Multiplicity of events with at least one deuteron;dN/dy", 201, -0.5, 200.5);
  TH2* QvsR = new TH2D("QvsR", "QvsR;q[GeV];r[fm]", 500, 0, 5, 250, 0, 10);
  TH1* Events = new TH1D("Events", "Events", 2, 0, 2);
  TH2* RvspT = new TH2D("RvspT", "RvspT;r[fm];pT[GeV/c]", 500, 0, 5, 500, 0, 5);
  TH2* RvsMt = new TH2D("RvsMt", "RvsMt;r[fm];pT[GeV/c]", 125, 0, 7.5, 16, 0.8, 2.4);
  TH1* EtaHistP = new TH1D("EtaP", "EtaP;#eta;N", 400, -2, 2);
  TH1* EtaHistD = new TH1D("EtaD", "EtaD;#eta;N", 400, -2, 2);
  TH1* DeltaPhi = new TH1D("DeltaPhi", "DeltaPhi;#Delta#phi;P(#Delte#phi)", 29, -1.41, 4.88);

  for (int64_t evt = 0; evt < NEvents; evt++) {
    Devt = false;
    bool EvtChecked = false;
    if (evt % 100000 == 0) {
      cout << "accessing event " << evt << " / " << NEvents << " (" << (double)evt / NEvents * 100 << "%)        "
           << "\r" << std::flush;
    }

    Events->Fill(1);

    int NParticles = GetNParticles(MultiplicityType, Multiplicity, MultScaleFactor);
    if (MultiplicityType == 4) {
      NParticles = GetNumberOfParticlesFromHistogram(MultHist);
    }
    NParticlesH->Fill(NParticles / MultScaleFactor);
    std::tuple<int, int> NumberOfProtNeut = GetNumberNucleons(NParticles, Multiplicity, HadronizationMode, MultScaleFactor, AntimatterRatio, IsospinRatio);
    int NProtons = std::get<0>(NumberOfProtNeut);  // number of protons in the event
    int NNeutrons = std::get<1>(NumberOfProtNeut); // number of neutrons in the event

    for (int pr = 0; pr < NProtons; pr++) {
      double pT = GetRandompTMult((double)NParticles / MultScaleFactor, MultiplicityType, HadronizationMode);
      double rap = RandomR->Uniform(-0.5, 0.5);
      double m = MP;
      double phi = RandomR->Uniform(0, 2 * PI);
      double px = pT * cos(phi);
      double py = pT * sin(phi);
      double pz = CalcpZ(px, py, m, rap);
      double E = CalcEngy(px, py, pz, m);
      double ptot = sqrt(px * px + py * py + pz * pz);
      double theta = acos(pz / ptot);
      double Eta = 0.5 * log((ptot + pz) / (ptot - pz));
      ptHist->Fill(pT);
      EtaHistP->Fill(Eta);

      for (int ne = 0; ne < NNeutrons; ne++) {
        double rapN = RandomR->Uniform(-0.5, 0.5);
        double mN = 0.939565;
        double ptN = GetRandompTMult((double)NParticles / MultScaleFactor, MultiplicityType, HadronizationMode); // scaling because NParticles is dN/dy but Mult classes are dN/deta

        double dPhi = GetRandomdPhiPar((double)NParticles / MultScaleFactor, (ptN + pT) / 2.);
        double phiN = phi + dPhi;
        DeltaPhi->Fill(dPhi);

        double pxN = ptN * cos(phiN);
        double pyN = ptN * sin(phiN);
        double pzN = CalcpZ(pxN, pyN, mN, rapN);
        double EN = CalcEngy(pxN, pyN, pzN, mN);
        double ptotN = sqrt(pxN * pxN + pyN * pyN + pzN * pzN);
        double thetaN = acos(pzN / ptotN);
        double EtaN = 0.5 * log((ptotN + pzN) / (ptotN - pzN));
        ptNHist->Fill(ptN);
        double dEta = Eta - EtaN;
        double mT = sqrt(pow((pT + ptN) / 2, 2) + pow(MP, 2));
        double DeltaP = CalcDeltaP(TLorentzVector(px, py, pz, E), TLorentzVector(pxN, pyN, pzN, EN));
        double distance = 0;
        distance = GetRandomDistance(mT, (double)NParticles / MultScaleFactor);

        double RandomCoalescence = RandomR->Uniform(0, 1);

        int Binx = CoalescenceProbability->GetXaxis()->FindBin(distance / (4. / sqrt(PI)));
        int Biny = CoalescenceProbability->GetYaxis()->FindBin(DeltaP);

        QvsR->Fill(DeltaP, distance);
        RvspT->Fill(distance, pT + ptN);
        RvsMt->Fill(distance, mT);
        if (!EvtChecked) {
          DeuteronEvtMult->Fill(NParticles);
          EvtChecked = true;
        }
        if (CoalescenceProbability->GetBinContent(Binx, Biny) > RandomCoalescence) {
          double RandomSpin = RandomR->Uniform(0, 1);
          if (RandomSpin < 3. / 8.) {
            double ptD = pT + ptN;
            double EtaD = TLorentzVector(px + pxN, py + pyN, pz + pzN, E + EN).Eta();
            EtaHistD->Fill(EtaD);
            NNeutrons--;
            ptDeuteronArgHist->Fill(ptD);
            DeuteronB2Hist->Fill(ptD);
            if (Devt == false) {
              Devt = true;
            }
            break;
          }
        }
      }
    }
    // delete DeltaPhiCorrelation;
  } // NEvts

  DoverPMult->Fill(Multiplicity, ptDeuteronArgHistFine->GetEntries() / ptHist->GetEntries());

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  cout << "\n Finished in " << duration.count() / 1000000 << "s! \n";

  /*Normalize yield histograms to number of events*/
  ptNHist->Scale(1. / (double)(NEvents), "width");
  ptDeuteronArgHist->Scale(1. / (double)(NEvents), "width");
  B2pt->Scale(1. / (double)(NEvents), "width");
  DeuteronB2Hist->Scale(1. / (double)(NEvents), "width");
  ptHist->Scale(1. / (double)NEvents, "width");
  /*Calculate B2*/
  for (int i = 1; i <= 100; i++) {
    B2pt->SetBinContent(i, DeuteronB2Hist->GetBinContent(i) / pow(ptHist->GetBinContent(i), 2) * PI * B2pt->GetBinCenter(i));
  }
  /*Write output*/
  TFile* fout = new TFile("Output.root", "UPDATE");
  fout->cd();
  ptHist->Write();
  ptNHist->Write();
  ptDeuteronArgHist->Write();
  QvsR->Write();
  RvspT->Write();
  RvsMt->Write();
  DoverPMult->Write();
  B2pt->Write();
  EtaHistP->Write();
  EtaHistD->Write();
  DeltaPhi->Write();
  NParticlesH->Write();
  DeuteronEvtMult->Write();
  Files->Write();
  Events->Write();
  fout->Close();

  delete ptHist;
  delete ptNHist;
  delete ptDeuteronArgHist;
  delete QvsR;
  delete RvspT;
  delete RvsMt;
  delete DoverPMult;
  delete B2pt;
  delete EtaHistP;
  delete EtaHistD;
  delete DeltaPhi;
  delete NParticlesH;
  delete DeuteronEvtMult;
  delete Files;
  delete Events;
}
