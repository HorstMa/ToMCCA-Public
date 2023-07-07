#include "TH1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <random>

double LevyTsallis(double x, double Norm, double n, double T, double m)
{
    return (x * Norm * (n - 1) * (n - 2)) / (n * T * (n * T + m * (n - 2))) * pow(1 + (sqrt(m * m + x * x) - m) / (n * T), -n);
}
// If you want to make a custom pT dependence, change this yield also!
double Yield()
{
    return 0.938;
}
/// Use this function for custom pT distributions
double GetRandompT()
{
    bool found = false;
    while (found == false)
    {
        double pTGuess = gRandom->Uniform(0, 5);
        double probguess = gRandom->Uniform(0, 0.6);
        // INEL 0.13850372 6.64502561 0.23982105
        // INEL>0 0.43670269, 8.83619089, 0.24308461
        if (LevyTsallis(pTGuess, Yield(), 3.66063241, 0.11909471, 1.95716802) > probguess)
        {
            found = true;
            return pTGuess;
        }
    }
}
/// This function is hard-coded for HM collisions
double GetRandompTHM()
{
    bool found = false;
    while (found == false)
    {
        double pTGuess = gRandom->Uniform(0, 5);
        double probguess = gRandom->Uniform(0, 0.9);

        if (LevyTsallis(pTGuess, 0.938, 7.43277365, 0.40303088, 0.938271) > probguess)
        {
            found = true;
            return pTGuess;
        }
    }
}
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

double GetRandomdPhi(TH1 *distHist)
{
    double DistMax = distHist->GetMaximum();
    double DistMin = distHist->GetMinimum();
    double FirstBin = distHist->GetBinCenter(1);
    double LastBin = distHist->GetBinCenter(distHist->GetNbinsX() + 1);
    bool found = false;
    while (found == false)
    {
        double PhiGuess = gRandom->Rndm() * (LastBin - FirstBin) + FirstBin;
        double probguess = gRandom->Rndm() * DistMax;
        if (distHist->GetBinContent(distHist->FindBin(PhiGuess)) > probguess)
        {
            found = true;
            return PhiGuess;
        }
    }
}
double GetRandomdEta(double EtaP1)
{
    return gRandom->Rndm() * 2 - 1;
}
double CalcDeltaP(TLorentzVector p1, TLorentzVector p2)
{
    TLorentzVector pPair = p1 + p2;
    TVector3 b = -pPair.BoostVector();
    p1.Boost(b);
    p2.Boost(b);
    return p1.Vect().Mag();
}

double GausSource(double r, double s)
{
    return 4 * 3.1415 * pow(r, 2) / pow(4 * 3.1415 * pow(s, 2), 3. / 2.) * exp(-pow(r, 2) / (4 * pow(s, 2)));
}
double GetRandomDistance(double mt)
{

    double sigma = 1.419 * pow(mt, -0.536); /// HM Source
    bool found = false;
    while (found == false)
    {
        double dstGuess = gRandom->Rndm() * 15;
        double probguess = gRandom->Rndm() * 2;
        if (GausSource(dstGuess, sigma) > probguess)
        {
            found = true;
            return dstGuess;
        }
    }
}

double Chi2(std::vector<double> model, std::vector<double> sim, std::vector<double> Err, int N = 0)
{
    double x2 = 0;
    int ndf = 0;
    int lim = N;
    int si = model.size();
    if (N == 0)
    {
        lim = model.size();
    }
    for (int i = 0; i < si; i++)
    {
        if (model.at(i) > 0)
        {
            if (sim.at(i) > 0)
            {
                ndf++;
                x2 += pow((model.at(i) - sim.at(i)) / Err.at(i), 2);
            }
        }
    }
    return x2 / ndf;
}