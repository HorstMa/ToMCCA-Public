#ifndef HADRONIZATIONMODELS_H
#define HADRONIZATIONMODELS_H

#include "Multiplicity.h"

std::tuple<int, int> GetNumberNucleons(int NParticles = 20, double Multiplicity = 20, std::string HadronizationMethod = "CorrelatedEmission", double MultScaleFactor = 1.192, double AntimatterRatio = 1., double IsospinRatio = 1.)
{
  double NucleonProbability = Yield((double)NParticles / MultScaleFactor) / NParticles;
  /* This produces positively correlated Proton and Neutron Yields */
  int NumberOfProtons = 0;
  int NumberOfNeutrons = 0;
  if (HadronizationMethod == "CorrelatedEmission") {
    for (int r = 0; r < NParticles; r++) {
      if (gRandom->Rndm() < NucleonProbability / 2 * AntimatterRatio / IsospinRatio) /// Isospin Proton=Antineutron, Neutron=Antiproton --> if Proton>Neutron >>> Antiproton<Antineutron
      {
        NumberOfProtons++;
      }
    }
    for (int r = 0; r < NParticles; r++) {
      if (gRandom->Rndm() < NucleonProbability / 2 * AntimatterRatio / IsospinRatio) /// Isospin Proton=Antineutron, Neutron=Antiproton --> if Proton>Neutron >>> Antiproton<Antineutron
      {
        NumberOfNeutrons++;
      }
    }
  }

  if (HadronizationMethod == "QuarkRecombination") { // this is off by ~ factor 2
    int nUp = NParticles;
    int nDown = NParticles;
    int nPion = 0; // mainly for debugging purpose!
    while (nUp + nDown > 2) {
      if (gRandom->Rndm() <= NucleonProbability * PowerLaw((double)NParticles / MultScaleFactor, 0.87627309, 1.01566945, -0.54537776)) // // (1.2880868 / (Multiplicity + 0.96228454) + 1.03304187) *
      {
        nUp--;
        nDown--;
        if (gRandom->Rndm() <= 1 / (IsospinRatio + 1)) {
          nUp--;
          NumberOfProtons++;
        } else {
          nDown--;
          NumberOfNeutrons++;
        }
      } else {
        if (gRandom->Rndm() <= 0.5) {
          nUp--;
          nPion++;
        } else {
          nDown--;
          nPion++;
        }
      }
    }
  }

  if (HadronizationMethod == "StringFragmentation") {
    int NFractions = NParticles;
    int start = gRandom->Integer(2) ? -1 : 1;
    int order = -start * 1;
    int StartingQuark = gRandom->Integer(2) ? start : start * 2;
    std::vector<int> String{ StartingQuark };
    bool LastWasDiquark = false;
    double prob = NucleonProbability;
    // * (0.23431261 / (Multiplicity + 0.6806183) + 1.04231348)
    double DiquarkProbability = 2 * NucleonProbability * PowerLaw((double)NParticles / MultScaleFactor, 1.06978999, -1.92787585, -2.0928103);
    for (int i = 0; i < NFractions; i++) {
      if (i % 2) {
        order *= -1;
      }
      prob = (LastWasDiquark) ? 0 : DiquarkProbability;

      if (gRandom->Uniform(0, 1) > prob) // dont add a diquark
      {
        LastWasDiquark = false;
        if (gRandom->Uniform(0, 1) < 0.5) // add ddbar pair
        {
          String.emplace_back(-1 * order);
          String.emplace_back(order);
        } else // add uubar pair
        {
          String.emplace_back(-2 * order);
          String.emplace_back(2 * order);
        }
      } else /// add a diquark
      {
        LastWasDiquark = true;
        double diquark = gRandom->Uniform(0, 1);
        if (diquark < 1. / 3.) // add a dd-dbardbar pair (3-3)
        {
          String.emplace_back(-3 * order);
          String.emplace_back(3 * order);
        } else if (diquark > 2. / 3.) // add a uu-ubarubar pair (4-4)
        {
          String.emplace_back(-4 * order);
          String.emplace_back(4 * order);
        } else // add a ud-ubardbar pair (4-4)
        {
          String.emplace_back(-5 * order);
          String.emplace_back(5 * order);
        }
      }
    } /// for (size_t i = 0; i < NFractions; i++)
    String.push_back(-1 * String.at(0));

    for (size_t j = 0; j < String.size() / 2; j++) {
      if (String.at(2 * j) > 2 || String.at(2 * j + 1) > 2) // one of the pair is a diquark
      {
        std::vector<int> partArray;
        //(String.at(2 * j) < String.at(2 * j + 1)) ? partArray = {String.at(2 * j), String.at(2 * j + 1)} : partArray = {String.at(2 * j + 1), String.at(2 * j)};
        int particle = (String.at(2 * j) < String.at(2 * j + 1)) ? stoi(std::to_string(String.at(2 * j)) + std::to_string(String.at(2 * j + 1))) : stoi(std::to_string(String.at(2 * j + 1)) + std::to_string(String.at(2 * j)));
        // int particle = stoi(std::to_string(partArray.at(0)) + std::to_string(partArray.at(1)));

        switch (particle) {
          case 15:
          case 14:
          case 24:
            NumberOfProtons++;
            break;
          case 13:
          case 25:
          case 23:
            NumberOfNeutrons++;
            break;
          default:
            break;
        }
      }
    }
  }
  if (HadronizationMethod == "TunedHadronization") {
    //NucleonProbability *= 1.1;
    double NucleonMean = NucleonProbability * NParticles;
    for (int r = 0; r < NParticles; r++) {
      if (gRandom->Rndm() < NucleonProbability / 2 * AntimatterRatio / IsospinRatio) /// Isospin Proton=Antineutron, Neutron=Antiproton --> if Proton>Neutron >>> Antiproton<Antineutron
      {
        NumberOfProtons++;
      }
    }
    double suppression = 1.508 * pow(NParticles, -1.) * (NucleonMean - NumberOfProtons) + 1.; /*The amount of suppression needed to reproduce data*/
    for (int r = 0; r < NParticles; r++) {
      if (gRandom->Rndm() < NucleonProbability / 2 * AntimatterRatio / IsospinRatio * suppression) /// Isospin Proton=Antineutron, Neutron=Antiproton --> if Proton>Neutron >>> Antiproton<Antineutron
      {
        NumberOfNeutrons++;
      }
    }
  }

  return std::make_tuple(NumberOfProtons, NumberOfNeutrons);
}
#endif