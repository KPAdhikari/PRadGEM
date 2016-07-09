#ifndef _GEMHISTOMANAGER_H__
#define _GEMHISTOMANAGER_H__

/*
 * This class fill all physics histograms
 */

#include <TH2F.h>
#include <TH1F.h>
#include <TString.h>
#include <TCanvas.h>

#include <cassert>

class PRadGEMSystem;
class PRadDataHandler;
class PRDMapping;

class GEMHistoManager
{
public:
  GEMHistoManager();
  ~GEMHistoManager();

  virtual void BookHistos();

protected:
  PRDMapping * fMapping;

public:
  int nbDetector;
  TH2F* hhClusterDist[10];
  TH2F* hhStripHitDist[10];
  TH2F* hhClusterDist_pRad[10];
  TH2F* hhChargeRatio[10];
  TH1F* hChargeRatio[10];
  TH1F* hNbOfClusterXSide[10];
  TH1F* hNbOfClusterYSide[10];
  TH1F* hClusterADCDistXSide[10];
  TH1F* hClusterADCDistYSide[10];
  
  TH1F* hNbClusterPerPlaneX[10];
  TH1F* hNbClusterPerPlaneY[10];
  TH1F* hClusterDistX[10];
  TH1F* hClusterDistY[10];

public:
  TH2F* hhGEMClusterMap;
  TH2F* hhGEMClusterMapSingleCluster;
  TH2F* hhGEMClusterMapMaxTwoCluster;
  TH2F* hhHyCalClusterMap;
  TH1F* hGEMClusterMul;
  TH1F* hHyCalEnergy;
  TH1F* hHyCalEnergyEp;
  TH1F* hHyCalEnergyMoller;
  TH1F* hHyCalClusterMul;

  TH1F* hThetaDistribution;
  TH1F* hThetaDistributionEp;
  TH1F* hThetaDistributionMoller;
  TH1F* hThetaDistributionMollerLarge;
  TH1F* hThetaDistributionMollerSmall;
  TH1F* hThetaDistributionMollerSlice;
  TH1F* hLinearDeviationMoller;
  TH1F* hQSquareEp;
  TH1F* hQSquareEp1;
  TH1F* hQSquareEp2;
  TH1F* hQSquareEp3;
  TH1F* hQSquareEp4;
  TH1F* hQSquareEp5;
  TH2F* hhQSquareScattAngleEp;


  TH2F* hhEnergyVsAngle;
  TH2F* hhEnergyVsAngleEp;
  TH2F* hhEnergyVsAngleMoller;
  TH2F* hhEnergyVsAngleHyCal2CGEM1C;
  //2d cluster map for hycal 2 gem 1
  TH2F * hhGEMClusterMapHyCal2GEM1;
  TH2F * hhHyCalClusterMapHyCal2GEM1;
  TH2F * hhHyCalClusterMap2ClusterBeforeMatch;
  TH2F * hhGEMClusterMap2ClusterAfterMatch;

  // after offset correction with beam line
  TH2F* hhEnergyVsAngleBeamLineCorrection;
  TH2F* hhEnergyVsAngleMollerBeamLineCorrection;
  TH1F* hThetaDistributionMollerBeamLineCorrection;

  // x, y offset between two GEM chambers thrugh ep
  TH1F* hXOffsetGEM;
  TH1F* hYOffsetGEM;
  TH1F* hXOffsetGEMAfterCorrection;
  TH1F* hYOffsetGEMAfterCorrection;
 
  // x, y offset through moller
  TH1F* hXOffsetFromMoller;
  TH1F* hYOffsetFromMoller;
  TH2F* hhMollerCenter;
  TH1F* hXOffsetFromMollerAfterCorrection;
  TH1F* hYOffsetFromMollerAfterCorrection;
  TH2F* hhMollerCenterAfterCorrection;
  TH1F* hLinearDeviationMollerAfterCorrection;

  TH1F* hXOffsetFromSymmetricMoller;
  TH1F* hYOffsetFromSymmetricMoller;
  TH2F* hhSymmetricMollerCenter;
  TH1F* hSymmetricMollerTheta;

  // moller 2d position
  TH2F * hhMoller2DRing;
  TH2F * hhMoller2DRing2;
  TH2F * hhMoller2DRingBeforeMatch;
  TH2F * hhMoller2DRingBeforeMatch2;

  // hycal gem positon match 
  TH1F * hXDiffMatch;
  TH1F * hYDiffMatch;
  TH2F * hhXDiffMatchVsEnergy;
  TH2F * hhYDiffMatchVsEnergy;
  TH1F * hNbPointsMatch;
 
};

#endif
