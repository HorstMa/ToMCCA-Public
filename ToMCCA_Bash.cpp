#include "TH1.h"
#include "ToMCCA_func.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <math.h>
#define PI 3.14159265
#define MP 0.938272

void ToMCCA_Bash(double Multiplicity = 31.0, int MultiplicityType = 0, int64_t NEvents = 1000000)
{
    gRandom = new TRandom3(time(NULL));

    double MultScaleFactor = 1.192; // 1.16 from EPOS // 1.192;
    Multiplicity /= 10.;            // divide by 10 to get correct Multiplicity with one digit precision

    TFile *fin0 = new TFile("Ressources/MCComparison.root");
    TH1 *dphiCorrelatorHist = (TH1D *)fin0->Get("ALICE_dphi");
    dphiCorrelatorHist->SetDirectory(0);
    fin0->Close();

    /// We dont need the Delta Phi correlation function, but the same event distribution. we "regain" it by multiplying the Correlation function with the mixed event from EPOS 3, hoping its similar enough to the real data
    TFile *fin1 = new TFile("Ressources/Correlations_HM.root");
    TH1 *dphiMixedHist = (TH1D *)fin1->Get("ProjectionsDPhi_Mixed");
    dphiMixedHist->SetDirectory(0);
    fin1->Close();

    TH1 *dPhiSameALICE = (TH1D *)dphiCorrelatorHist->Clone("dPhiSameALICE");
    for (int i = 1; i < 30; i++)
    {
        dPhiSameALICE->SetBinContent(i, dphiCorrelatorHist->GetBinContent(i) * dphiMixedHist->GetBinContent(i));
    }
    dPhiSameALICE->SetDirectory(0);

    TFile *fin = new TFile("Ressources/Correlations_HM.root");
    TH1 *dphiSameEPOS = (TH1D *)fin->Get("ProjectionsDPhi_Same");
    dphiSameEPOS->SetDirectory(0);
    fin->Close();
    /// This is the probability histogram for Argonne v18 with added S and D Wave. There are versions for only S wave and only D wave available in the same root file
    TFile *ArgonneProbabilityHistogramF = new TFile("Ressources/ArgonneProbabilityFinal3001.root");
    TH1 *ArgonneProbabilityHistogram = (TH1D *)ArgonneProbabilityHistogramF->Get("AddedSDWave");
    ArgonneProbabilityHistogram->SetDirectory(0);
    ArgonneProbabilityHistogramF->Close();

    TH1 *ptDeuteronArgHist = new TH1D("ptDeuteronArg", "ptDeuteronArg", 40, 0, 4);
    TH1 *ptDeuteronArgHistFine = new TH1D("ptDeuteronArgFine", "ptDeuteronArgFine", 500, 0, 5);
    TH1 *ptNHist = new TH1D("ptN", "ptN", 100, 0, 8);
    TH1 *ptHist = new TH1D("pt", "pt", 100, 0, 8);

    TH1 *NParticlesH = new TH1D("NParticles", "NParticles;dN/dy", 200, 0, 200);
    TH2 *QvsR = new TH2D("QvsR", "QvsR;q[GeV];r[fm]", 500, 0, 5, 500, 0, 5);
    TH2 *QvspT = new TH2D("QvspT", "QvspT pt<1;q[GeV];pT[GeV/c]", 500, 0, 5, 500, 0, 5);
    TH1 *Events = new TH1D("Events", "Events", 2, 0, 2);
    TH2 *RvspT = new TH2D("RvspT", "RvspT;r[fm];pT[GeV/c]", 500, 0, 5, 500, 0, 5);
    TH2 *RvsMt = new TH2D("RvsMt", "RvsMt;r[fm];pT[GeV/c]", 125, 0, 7.5, 16, 0.8, 2.4);
    TH1 *EtaHistP = new TH1D("EtaP", "EtaP;#eta;N", 400, -2, 2);
    TH1 *EtaHistD = new TH1D("EtaD", "EtaD;#eta;N", 400, -2, 2);
    TH1 *AngularDist = new TH1D("AngularDist", "AngularDistribution;#Delta#phi;N", 1000, -6.4, 6.4);
    TH1 *AngularDistCorr = (TH1D *)dphiMixedHist->Clone("AngularCorrelationToMCCA");
    AngularDistCorr->Reset("ICE");
    auto start = std::chrono::high_resolution_clock::now();

    for (int64_t evt = 0; evt < NEvents; evt++)
    {
        if (evt % 500000 == 0)
        {
            cout << "accessing event " << evt << " / " << NEvents << " (" << (double)evt / NEvents * 100 << "%)        "
                 << "\r" << std::flush;
        }

        Events->Fill(1);
        int NProtons = 0;
        int NNeutrons = 0;
        int NParticles = 0;
        // Find the way to determine the multiplicity
        // First option: Fixed value, second option: Poissionian distributions
        if (MultiplicityType == 2)
        {
            NParticles = (int)floor(MultScaleFactor * (double)Multiplicity);
        }
        else if (MultiplicityType == 1)
        {
            NParticles = gRandom->Poisson(MultScaleFactor * (double)Multiplicity);
        }
        else
        {
            cout << "Please specify a Multiplicity type! chosing default type 1 (poissonian)" << endl;
            NParticles = gRandom->Poisson(MultScaleFactor * (double)Multiplicity);
        }
        NParticlesH->Fill(NParticles);

        double NucleonProbability = Yield() / NParticles;

        for (int r = 0; r < NParticles; r++)
        {
            if (gRandom->Rndm() < NucleonProbability)
            {
                NProtons++;
            }
            if (gRandom->Rndm() < NucleonProbability)
            {
                NNeutrons++;
            }
        }

        for (int pr = 0; pr < NProtons; pr++)
        {
            double pT = GetRandompTHM();
            double rap = gRandom->Rndm() - 0.5;
            double m = MP;
            double phi = gRandom->Rndm() * 2 * PI;
            double px = pT * cos(phi);
            double py = pT * sin(phi);
            double pz = CalcpZ(px, py, m, rap);
            double E = CalcEngy(px, py, pz, m);
            double ptot = sqrt(px * px + py * py + pz * pz);
            double theta = acos(pz / ptot);
            double Eta = 0.5 * log((ptot + pz) / (ptot - pz));
            ptHist->Fill(pT);
            EtaHistP->Fill(Eta);

            for (int ne = 0; ne < NNeutrons; ne++)
            {
                double dPhi = GetRandomdPhi(dPhiSameALICE);
                double phiN = phi + dPhi;
                double rapN = gRandom->Rndm() - 0.5;
                double mN = 0.939565;
                double ptN = GetRandompTHM(); // scaling because NParticles is dN/dy but Mult classes are dN/deta --> we scale from 8.222 to 6.9

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
                double distance = GetRandomDistance(mT);
                double DeltaP = CalcDeltaP(TLorentzVector(px, py, pz, E), TLorentzVector(pxN, pyN, pzN, EN));

                double RandomCoalescence = gRandom->Rndm();

                int Binx = ArgonneProbabilityHistogram->GetXaxis()->FindBin(distance / 2.26);
                int Biny = ArgonneProbabilityHistogram->GetYaxis()->FindBin(DeltaP);

                QvspT->Fill(DeltaP, pT + ptN);
                QvsR->Fill(DeltaP, distance);
                RvspT->Fill(distance, pT + ptN);
                RvsMt->Fill(distance, mT);
                AngularDist->Fill(dPhi);
                AngularDistCorr->Fill(dPhi);
                if (ArgonneProbabilityHistogram->GetBinContent(Binx, Biny) > RandomCoalescence)
                {
                    double RandomSpin = gRandom->Rndm();
                    if (RandomSpin < 3. / 8.)
                    {
                        double ptD = pT + ptN;
                        double EtaD = TLorentzVector(px + pxN, py + pyN, pz + pzN, E + EN).Eta();
                        EtaHistD->Fill(EtaD);
                        NNeutrons--;
                        ptDeuteronArgHist->Fill(ptD);
                        ptDeuteronArgHistFine->Fill(ptD);
                        break;
                    }
                }
            }
        }
    } // NEvts

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    cout << "\n Finished in " << duration.count() / 1000000 << "s! \n";
    ptNHist->Scale(1. / (double)(NEvents), "width");
    ptDeuteronArgHist->Scale(1. / (double)(NEvents), "width");
    ptDeuteronArgHistFine->Scale(1. / (double)(NEvents), "width");
    ptHist->Scale(1. / (double)NEvents, "width");

    AngularDistCorr->Scale(dphiMixedHist->Integral() / AngularDistCorr->Integral());
    AngularDistCorr->Divide(dphiMixedHist);

    TFile *fout = new TFile("Output.root", "RECREATE");
    fout->cd();
    ptHist->Write();
    ptNHist->Write();
    ptDeuteronArgHist->Write();
    ptDeuteronArgHistFine->Write();
    QvspT->Write();
    QvsR->Write();
    RvspT->Write();
    RvsMt->Write();
    NParticlesH->Write();
    Events->Write();
    EtaHistP->Write();
    AngularDist->Write();
    AngularDistCorr->Write();
    dPhiSameALICE->Write();
    fout->Close();

    delete ptDeuteronArgHist;
    delete ptDeuteronArgHistFine;
    delete ptNHist;
    delete ptHist;
    delete NParticlesH;
    delete QvsR;
    delete QvspT;
    delete RvsMt;
    delete Events;
    delete RvspT;
}
