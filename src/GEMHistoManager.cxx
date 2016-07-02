#include "GEMHistoManager.h"

using namespace std;

GEMHistoManager::GEMHistoManager()
{
  fMapping = PRDMapping::GetInstance();
  nbDetector = fMapping->GetNbOfDetectors(); 
}

GEMHistoManager::~GEMHistoManager()
{
}

void GEMHistoManager::BookHistos()
{
  for(int i=0;i<nbDetector;i++)
  {
    TString detectorName = fMapping->GetDetectorFromID(i);
    list<TString> plane = fMapping->GetPlaneListFromDetector(detectorName);
    assert(plane.size() == 2);
    Float_t ySide = fMapping->GetPlaneSize(plane.front() );
    Float_t xSide = fMapping->GetPlaneSize(plane.back() );
    
    //2D cluster distribution
    TString detectorTitle = detectorName + " 2D Cluster Distribution";
    hhClusterDist[i] = new TH2F(detectorName.Data(), detectorTitle.Data(), 2000, -xSide/2-50, xSide/2+50, 1200, -ySide/2-50, ySide/2+50);
    TString stripHit = detectorName+"_stripHit";
    TString stripHitTitle = stripHit + " Distribution";
    hhStripHitDist[i] = new TH2F(stripHit.Data(), stripHitTitle.Data(), 601, -300.5, 300.5, 1241, -620.5, 620.5);
    
    //2D cluster distribution HyCal Coordiante
    TString detectorTitlepRad = detectorName + " 2D Cluster Distribution pRad";
    TString detectorNamepRad = detectorName + "_pRad";
    hhClusterDist_pRad[i] = new TH2F(detectorNamepRad.Data(), detectorTitlepRad.Data(), 600, -xSide/2-100, xSide/2+100, 1200, -ySide/2-50, ySide/2+50);
    

    // 2d Charge Ratio Distribution
    TString chargeRatio = detectorName + " Charge Ratio";
    TString hchargeratio = detectorName+"charge_ratio";
    hhChargeRatio[i] = new TH2F(hchargeratio, chargeRatio, 2000, 0, 4000, 2000, 0, 4000);
    // 1d Charge Ratio Distribution
    TString chargeRatio1d = detectorName + " Charge Ratio 1D";
    TString hchargeratio1d = detectorName+"charge_ratio_1d";
    hChargeRatio[i] = new TH1F(hchargeratio1d, chargeRatio1d, 1000, -1, 10);

    // 1d cluster adc distribution
    TString clusteradcNamex = detectorName+"XClusterADC";
    TString clusteradcNamey = detectorName+"YClusterADC";
    TString clusteradctitlex = detectorName+"X side Cluster ADC";
    hClusterADCDistXSide[i] = new TH1F(clusteradcNamex, clusteradctitlex, 1000,0,8000);
    TString clusteradctitley = detectorName+"Y Side Cluster ADC";
    hClusterADCDistYSide[i] = new TH1F(clusteradcNamey, clusteradctitley, 1000, 0, 8000);

    // 1d cluster position distribution
    TString clusterdistNameX = detectorName+"_XClusterPos";
    TString clusterdistNameY = detectorName+"_YClusterPos";
    TString clusterdistTitleX = detectorName+"_X Cluster Dist";
    TString clusterdistTitleY = detectorName+"_Y Cluster Dist";
    TString nbclusterNameX = detectorName+"_XClusterSize";
    TString nbclusterNameY = detectorName+"_YClusterSize";
    TString nbclusterTitleX = detectorName+"_X Cluster Size";
    TString nbclusterTitleY = detectorName+"_Y Cluster Size";

    hNbClusterPerPlaneX[i] = new TH1F(nbclusterNameX, nbclusterTitleX, 100, 0, 100);
    hNbClusterPerPlaneY[i] = new TH1F(nbclusterNameY, nbclusterTitleY, 100, 0, 100);
    hClusterDistX[i] = new TH1F(clusterdistNameX, clusterdistTitleX,  600, -xSide/2-100, xSide/2+100);
    hClusterDistY[i] = new TH1F(clusterdistNameY, clusterdistTitleY, 1200, -ySide/2-50, ySide/2+50);
    // other histograms
  }

  // single histograms physics histogram
  hGEMClusterMul = new TH1F("gemnbClusterPerEvent", "# of Clusters Per Event GEM", 100, 0, 100);
  hHyCalClusterMul = new TH1F("hycalnbClusterPerEvent", "# of Clusters Per Event HyCal", 100, 0, 100);

  hhGEMClusterMap = new TH2F("hhGEMClusterMap","Cluster Map", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMapSingleCluster = new TH2F("hhGEMClusterMapSingleCluster","Cluster Map single cluster", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMapMaxTwoCluster = new TH2F("hhGEMClusterMapMaxTwoCluster","Cluster Map max two cluster", 1000, -650, 650, 1000, -650, 650);
  hhHyCalClusterMap = new TH2F("hhHyCalClusterMap","Cluster Map", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMapHyCal2GEM1 = new TH2F("hhGEMClusterMapHyCal2GEM1","gem Cluster Map hycal 2 gem 1", 1000, -650, 650, 1000, -650, 650);
  hhHyCalClusterMapHyCal2GEM1 = new TH2F("hhHyCalClusterMapHyCal2GEM1","hycal Cluster Map hycal 2 gem 1", 1000, -650, 650, 1000, -650, 650);
  hhHyCalClusterMap2ClusterBeforeMatch = new TH2F("hhHyCalClusterMap2ClusterBeforeMatch","Cluster Map hycal 2 cluster before match", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMap2ClusterAfterMatch = new TH2F("hhGEMClusterMap2ClusterBeforeMatch","Cluster Map gem 2 cluster after match", 1000, -650, 650, 1000, -650, 650);

  hHyCalEnergy = new TH1F("hHyCalEnergy", "HyCal Energy [MeV]", 5100, -100, 5000);
  hHyCalEnergyEp = new TH1F("hHyCalEnergyEp", "HyCal Energy Ep [MeV]", 5100, -100, 5000);
  hHyCalEnergyMoller = new TH1F("hHyCalEnergyMoller", "HyCal Energy Moller [MeV]", 5100, -100, 5000);

  hThetaDistribution = new TH1F("thetadistribution", "Theta Distribution", 100, 0, 10);
  hThetaDistributionEp = new TH1F("thetadistribution_singleCluster", "theta distribution Ep", 1000, 0, 10 );
  hThetaDistributionMoller = new TH1F("thetadistribution_moller", "theta distribution moller", 1000, 0, 10 );
  hThetaDistributionMollerBeamLineCorrection = new TH1F("thetadistribution_mollerbeamline", "theta distribution moller beam line", 1000, 0, 10 );

  hThetaDistributionMollerLarge = new TH1F("thetadistribution_mollerLarge", "theta distribution moller large angle", 1000, 0, 10 );
  hThetaDistributionMollerSmall = new TH1F("thetadistribution_mollerSmall", "theta distribution moller small angle", 1000, 0, 10 );
  hThetaDistributionMollerSlice = new TH1F("thetadistribution_mollerSlice", "theta distribution moller slice", 1000, 0, 10 );

  hLinearDeviationMoller = new TH1F("LinearDeviationMoller", "Moller Clusters Linear Deviation", 10000, -180, 180);
  hQSquareEp = new TH1F("hQSquareEp", "hQSquareEp / GeV^2 all",     1000, 0,0.1);
  hQSquareEp1 = new TH1F("hQSquareEp1", "hQSquareEp / GeV^2 0.7-0.8", 400, 0.00008,0.01);
  hQSquareEp2 = new TH1F("hQSquareEp2", "hQSquareEp / GeV^2 1.0-1.1", 400, 0.0002,0.01);
  hQSquareEp3 = new TH1F("hQSquareEp3", "hQSquareEp / GeV^2 1.5-1.6", 400, 0.0008,0.01);
  hQSquareEp4 = new TH1F("hQSquareEp4", "hQSquareEp / GeV^2 2.0-2.1", 400, 0.002,0.01);
  hQSquareEp5 = new TH1F("hQSquareEp5", "hQSquareEp / GeV^2 2.2-", 400, 0.008,0.01);
  hhQSquareScattAngleEp = new TH2F("hhQSquareScattAngleEp", "Q^2 vs Scatt Angle", 1000, 0,8, 1000, 0, 0.1);

  hhEnergyVsAngle = new TH2F("hhEnergyVsAngle", "Energy Angle", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleBeamLineCorrection = new TH2F("hhEnergyVsAngleBeamLineCorrection", "Energy Angle beam line", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleEp = new TH2F("hhEnergyVsAngleEp", "Energy Angle Ep", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleMoller = new TH2F("hhEnergyVsAngleMoller", "Energy Angle Moller", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleMollerBeamLineCorrection = new TH2F("hhEnergyVsAngleMollerBeamLineCorrection", "Energy Angle Moller beam line", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleHyCal2CGEM1C = new TH2F("hhEnergyVsAngleHyCal2CGEM1C", "Energy Angle > 2 cluster hycal, 1 cluster on gem ", 1000, 0, 10, 4000, 0, 2500);

 
  hXOffsetFromMoller = new TH1F("hXOffsetFromMoller", "x offset using moller events", 1000, -50, 50);
  hYOffsetFromMoller = new TH1F("hYOffsetFromMoller", "y offset using moller events", 1000, -50, 50);
  hhMollerCenter = new TH2F("hhMollerCenter", "Moller Center", 1000, -50, 50, 1000, -50, 50);
  hXOffsetFromSymmetricMoller = new TH1F("hXOffsetFromSymmetricMoller", "x offset using moller events", 1000, -50, 50);
  hYOffsetFromSymmetricMoller = new TH1F("hYOffsetFromSymmetricMoller", "y offset using moller events", 1000, -50, 50);
  hhSymmetricMollerCenter = new TH2F("hhSymmetricMollerCenter", "Moller Center", 1000, -50, 50, 1000, -50, 50);
  hSymmetricMollerTheta = new TH1F("hSymmetricMollerTheta", "Symmetric Moller Theta", 1000, -10, 10);

  hXOffsetFromMollerAfterCorrection = new TH1F("hXOffsetFromMollerAfterCorrection", "x offset using moller events after correction", 1000, -50, 50);
  hYOffsetFromMollerAfterCorrection = new TH1F("hYOffsetFromMollerAfterCorrection", "y offset using moller events after correction", 1000, -50, 50);
  hhMollerCenterAfterCorrection = new TH2F("hhMollerCenterAfterCorrection", "Moller Center after correction", 1000, -50, 50, 1000, -50, 50);
  hLinearDeviationMollerAfterCorrection = new TH1F("LinearDeviationMollerAfterCorrection", "Moller Clusters Linear Deviation after correction", 10000, -180, 180);
 
  hXOffsetGEM = new TH1F("hXOffsetGEM", "x offset between two GEM chambers", 1000, -20, 20);
  hYOffsetGEM = new TH1F("hYOffsetGEM", "y offset between two GEM chambers", 1000, -20, 20);
  hXOffsetGEMAfterCorrection = new TH1F("hXOffsetGEMAfterCorrection", "x offset between two GEM chambers after correction", 1000, -20, 20);
  hYOffsetGEMAfterCorrection = new TH1F("hYOffsetGEMAfterCorrection", "y offset between two GEM chambers after correction", 1000, -20, 20);

  // moller 2d position
  hhMoller2DRing = new TH2F("hhMoller2DRing"," Moller ", 1000, -650, 650, 1000, -650, 650);
  hhMoller2DRing2 = new TH2F("hhMoller2DRing2"," Moller ", 1000, -650, 650, 1000, -650, 650);
  hhMoller2DRingBeforeMatch = new TH2F("hhMoller2DRingBeforeMatch"," 2 clusters ", 1000, -650, 650, 1000, -650, 650);
  hhMoller2DRingBeforeMatch2 = new TH2F("hhMoller2DRingBeforeMatch2"," 2 clusters  ", 1000, -650, 650, 1000, -650, 650);

  // hycal gem position match
  hXDiffMatch = new TH1F("hXDiffMatch", "x_gem - x_hycal", 1000, -60, 60);
  hYDiffMatch = new TH1F("hYDiffMatch", "y_gem - y_hycal", 1000, -60, 60);
  hhXDiffMatchVsEnergy = new TH2F("hhXDiffMatch", "x diff. vs energy", 2500, 0, 2500, 1000, -60, 60);
  hhYDiffMatchVsEnergy = new TH2F("hhYDiffMatch", "y diff. vs energy", 2500, 0, 2500, 1000, -60, 60);
  hNbPointsMatch = new TH1F("hNbPointsMatch", "# cluster match", 100, -1.5, 98.5);
}
