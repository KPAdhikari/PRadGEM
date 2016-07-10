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
    Float_t xSide = fMapping->GetPlaneSize(plane.front() ); 
    Float_t ySide = fMapping->GetPlaneSize(plane.back() );
    
    //2D cluster distribution
    TString detectorTitle = detectorName + " 2D Cluster Distribution";
    hhClusterDist[i] = new TH2F(detectorName.Data(), detectorTitle.Data(), 2000, -xSide/2-50, xSide/2+50, 1200, -ySide/2-50, ySide/2+50);
    hhClusterDist[i] -> GetXaxis()->SetTitle(" X [mm] ");
    hhClusterDist[i] -> GetYaxis()->SetTitle(" Y [mm] ");

    TString stripHit = detectorName+"_stripHit";
    TString stripHitTitle = stripHit + " Distribution";
    hhStripHitDist[i] = new TH2F(stripHit.Data(), stripHitTitle.Data(), 601, -300.5, 300.5, 1241, -620.5, 620.5);
    hhStripHitDist[i] -> GetXaxis()->SetTitle(" X [mm] "); 
    hhStripHitDist[i] -> GetXaxis()->SetTitle(" Y [mm] "); 

    //2D cluster distribution HyCal Coordiante
    TString detectorTitlepRad = detectorName + " 2D Cluster Distribution pRad";
    TString detectorNamepRad = detectorName + "_pRad";
    hhClusterDist_pRad[i] = new TH2F(detectorNamepRad.Data(), detectorTitlepRad.Data(), 600, -xSide/2-100, xSide/2+100, 1200, -ySide/2-50, ySide/2+50);
    hhClusterDist_pRad[i] -> GetXaxis()->SetTitle(" X [mm] ");
    hhClusterDist_pRad[i] -> GetYaxis()->SetTitle(" Y [mm] ");

    // 2d Charge Ratio Distribution
    TString chargeRatio = detectorName + " Charge Ratio";
    TString hchargeratio = detectorName+"charge_ratio";
    hhChargeRatio[i] = new TH2F(hchargeratio, chargeRatio, 2000, 0, 4000, 2000, 0, 4000);
    hhChargeRatio[i] -> GetXaxis()->SetTitle("X ADC");
    hhChargeRatio[i] -> GetYaxis()->SetTitle("Y ADC");
 
    // 1d Charge Ratio Distribution
    TString chargeRatio1d = detectorName + " Charge Ratio 1D";
    TString hchargeratio1d = detectorName+"charge_ratio_1d";
    hChargeRatio[i] = new TH1F(hchargeratio1d, chargeRatio1d, 1000, -1, 10);
    hChargeRatio[i] -> GetXaxis() -> SetTitle("Charge Ratio x/y"); 
    hChargeRatio[i] -> GetYaxis() -> SetTitle("entries"); 
 
    // 1d cluster adc distribution
    TString clusteradcNamex = detectorName+"XClusterADC";
    TString clusteradcNamey = detectorName+"YClusterADC";
    TString clusteradctitlex = detectorName+"X side Cluster ADC";
    hClusterADCDistXSide[i] = new TH1F(clusteradcNamex, clusteradctitlex, 1000,0,8000);
    hClusterADCDistXSide[i] -> GetXaxis()->SetTitle(" X side cluster adc");
    hClusterADCDistXSide[i] -> GetYaxis()->SetTitle(" entries ");
    TString clusteradctitley = detectorName+"Y Side Cluster ADC";
    hClusterADCDistYSide[i] = new TH1F(clusteradcNamey, clusteradctitley, 1000, 0, 8000);
    hClusterADCDistYSide[i] -> GetXaxis()->SetTitle(" Y side cluster adc");
    hClusterADCDistYSide[i] -> GetYaxis()->SetTitle(" entries ");
 
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
    hNbClusterPerPlaneX[i] -> GetXaxis() -> SetTitle("cluster quantity per event on x side");
    hNbClusterPerPlaneX[i] -> GetYaxis() -> SetTitle("entries");
    hNbClusterPerPlaneY[i] = new TH1F(nbclusterNameY, nbclusterTitleY, 100, 0, 100);
    hNbClusterPerPlaneY[i] -> GetXaxis() -> SetTitle("cluster quantity per event on y side");
    hNbClusterPerPlaneY[i] -> GetYaxis() -> SetTitle("entries");
 
    hClusterDistX[i] = new TH1F(clusterdistNameX, clusterdistTitleX,  600, -xSide/2-100, xSide/2+100);
    hClusterDistX[i] -> GetXaxis()->SetTitle(" x [mm] ");
    hClusterDistX[i] -> GetYaxis()->SetTitle(" entries ");
 
    hClusterDistY[i] = new TH1F(clusterdistNameY, clusterdistTitleY, 1200, -ySide/2-50, ySide/2+50);
    hClusterDistY[i] -> GetXaxis()->SetTitle(" y [mm] ");
    hClusterDistY[i] -> GetYaxis()->SetTitle(" entries ");
 
    // other histograms
  }

  // single histograms physics histogram
  hGEMClusterMul = new TH1F("gemnbClusterPerEvent", "# of Clusters Per Event GEM", 100, 0, 100);
  hGEMClusterMul -> GetXaxis() -> SetTitle("cluster quantity per event on GEM");
  hGEMClusterMul -> GetYaxis() -> SetTitle("entries");
 
  hHyCalClusterMul = new TH1F("hycalnbClusterPerEvent", "# of Clusters Per Event HyCal", 100, 0, 100);
  hHyCalClusterMul -> GetXaxis() -> SetTitle("cluster quantity per event on HyCal");
  hHyCalClusterMul -> GetYaxis() -> SetTitle("entries");
 
  hhGEMClusterMap = new TH2F("hhGEMClusterMap","Cluster Map", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMap -> GetXaxis() -> SetTitle(" x [mm] ");
  hhGEMClusterMap -> GetYaxis() -> SetTitle(" y [mm] ");
 
  hhGEMClusterMapSingleCluster = new TH2F("hhGEMClusterMapSingleCluster","Cluster Map single cluster", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMapSingleCluster -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhGEMClusterMapSingleCluster -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  hhGEMClusterMapMaxTwoCluster = new TH2F("hhGEMClusterMapMaxTwoCluster","Cluster Map max two cluster", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMapMaxTwoCluster -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhGEMClusterMapMaxTwoCluster -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  hhHyCalClusterMap = new TH2F("hhHyCalClusterMap","Cluster Map", 1000, -650, 650, 1000, -650, 650);
  hhHyCalClusterMap -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhHyCalClusterMap -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  hhGEMClusterMapHyCal2GEM1 = new TH2F("hhGEMClusterMapHyCal2GEM1","gem Cluster Map hycal 2 gem 1", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMapHyCal2GEM1 -> GetXaxis()->SetTitle(" x [mm] "); 
  hhGEMClusterMapHyCal2GEM1 -> GetYaxis()->SetTitle(" y [mm] "); 
 
  hhHyCalClusterMapHyCal2GEM1 = new TH2F("hhHyCalClusterMapHyCal2GEM1","hycal Cluster Map hycal 2 gem 1", 1000, -650, 650, 1000, -650, 650);
  hhHyCalClusterMapHyCal2GEM1 -> GetXaxis() -> SetTitle(" x [mm] ");
  hhHyCalClusterMapHyCal2GEM1 -> GetYaxis() -> SetTitle(" y [mm] ");
 
  hhHyCalClusterMap2ClusterBeforeMatch = new TH2F("hhHyCalClusterMap2ClusterBeforeMatch","Cluster Map hycal 2 cluster before match", 1000, -650, 650, 1000, -650, 650);
  hhHyCalClusterMap2ClusterBeforeMatch -> GetXaxis()->SetTitle(" x [mm] "); 
  hhHyCalClusterMap2ClusterBeforeMatch -> GetYaxis()->SetTitle(" y [mm] "); 
 
  hhGEMClusterMap2ClusterAfterMatch = new TH2F("hhGEMClusterMap2ClusterBeforeMatch","Cluster Map gem 2 cluster after match", 1000, -650, 650, 1000, -650, 650);
  hhGEMClusterMap2ClusterAfterMatch -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhGEMClusterMap2ClusterAfterMatch -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  hHyCalEnergy = new TH1F("hHyCalEnergy", "HyCal Energy [MeV]", 5100, -100, 5000);
  hHyCalEnergy -> GetXaxis() -> SetTitle(" Energy [MeV] ");
  hHyCalEnergy -> GetYaxis() -> SetTitle("entries ");

  hHyCalEnergyEp = new TH1F("hHyCalEnergyEp", "HyCal Energy ep [MeV]", 5100, -100, 5000);
  hHyCalEnergyEp -> GetXaxis() -> SetTitle(" Energy [MeV] ");
  hHyCalEnergyEp -> GetYaxis() -> SetTitle("entries ");

  hHyCalEnergyMoller = new TH1F("hHyCalEnergyMoller", "HyCal Energy Moller [MeV]", 5100, -100, 5000);
  hHyCalEnergyMoller -> GetXaxis() -> SetTitle(" Energy [MeV] ");
  hHyCalEnergyMoller -> GetYaxis() -> SetTitle("entries ");

  hThetaDistribution = new TH1F("thetadistribution", "#Theta Distribution", 100, 0, 10);
  hThetaDistribution -> GetXaxis()->SetTitle(" #theta [degree]"); 
  hThetaDistribution -> GetYaxis()->SetTitle(" entries "); 
 
  hThetaDistributionEp = new TH1F("thetadistribution_singleCluster", "#theta distribution Ep", 1000, 0, 10 );
  hThetaDistributionEp -> GetXaxis()->SetTitle(" #theta [degree]"); 
  hThetaDistributionEp -> GetYaxis()->SetTitle(" entries "); 
 
  hThetaDistributionMoller = new TH1F("thetadistribution_moller", "#theta distribution moller", 1000, 0, 10 );
  hThetaDistributionMoller -> GetXaxis()->SetTitle(" #theta [degree]"); 
  hThetaDistributionMoller -> GetYaxis()->SetTitle(" entries "); 
 
  hThetaDistributionMollerBeamLineCorrection = new TH1F("thetadistribution_mollerbeamline", "#theta distribution moller beam line", 1000, 0, 10 );
  hThetaDistributionMollerBeamLineCorrection -> GetXaxis()->SetTitle(" #theta [degree]"); 
  hThetaDistributionMollerBeamLineCorrection -> GetYaxis()->SetTitle(" entries "); 

  hThetaDistributionMollerLarge = new TH1F("thetadistribution_mollerLarge", "#theta distribution moller large angle", 1000, 0, 10 );
  hThetaDistributionMollerLarge -> GetXaxis()->SetTitle(" #theta [degree]"); 
  hThetaDistributionMollerLarge -> GetYaxis()->SetTitle(" entries "); 
 
  hThetaDistributionMollerSmall = new TH1F("thetadistribution_mollerSmall", "#theta distribution moller small angle", 1000, 0, 10 );
  hThetaDistributionMollerSmall -> GetXaxis()->SetTitle(" #theta [degree]"); 
  hThetaDistributionMollerSmall -> GetYaxis()->SetTitle(" entries "); 
 
  hThetaDistributionMollerSlice = new TH1F("thetadistribution_mollerSlice", "#theta distribution moller slice", 1000, 0, 10 );
  hThetaDistributionMollerSlice -> GetXaxis()->SetTitle(" #theta [degree]"); 
  hThetaDistributionMollerSlice -> GetYaxis()->SetTitle(" entries "); 
 
  hLinearDeviationMoller = new TH1F("LinearDeviationMoller", "Moller Clusters coplanarity", 10000, -180, 180);
  hLinearDeviationMoller -> GetXaxis() -> SetTitle("coplanarity [degree] "); 
  hLinearDeviationMoller -> GetYaxis() -> SetTitle("entries "); 
 
  hQSquareEp = new TH1F("hQSquareEp",   "Q^{2} ep / GeV^{2} all",     1000, 0,0.1);
  hQSquareEp -> GetXaxis() -> SetTitle("Q^{2} [ GeV^{2} / c ]"); 
  hQSquareEp -> GetYaxis() -> SetTitle(" entries "); 
  hQSquareEp1 = new TH1F("hQSquareEp1", "Q^{2} ep / GeV^{2} 0.7-0.8", 400, 0.00008,0.01);
  hQSquareEp1 -> GetXaxis() -> SetTitle("Q^{2} [ GeV^{2} / c ]"); 
  hQSquareEp1 -> GetYaxis() -> SetTitle(" entries "); 
 
  hQSquareEp2 = new TH1F("hQSquareEp2", "Q^{2} ep / GeV^{2} 1.0-1.1", 400, 0.0002,0.01);
  hQSquareEp2 -> GetXaxis() -> SetTitle("Q^{2} [ GeV^{2} / c ]"); 
  hQSquareEp2 -> GetYaxis() -> SetTitle(" entries "); 
 
  hQSquareEp3 = new TH1F("hQSquareEp3", "Q^{2} ep / GeV^{2} 1.5-1.6", 400, 0.0008,0.01);
  hQSquareEp3 -> GetXaxis() -> SetTitle("Q^{2} [ GeV^{2} / c ]"); 
  hQSquareEp3 -> GetYaxis() -> SetTitle(" entries "); 
 
  hQSquareEp4 = new TH1F("hQSquareEp4", "Q^{2} ep / GeV^{2} 2.0-2.1", 400, 0.002,0.01);
  hQSquareEp4 -> GetXaxis() -> SetTitle("Q^{2} [ GeV^{2} / c ]"); 
  hQSquareEp4 -> GetYaxis() -> SetTitle(" entries "); 
 
  hQSquareEp5 = new TH1F("hQSquareEp5", "Q^{2} ep / GeV^{2} 2.2-", 400, 0.008,0.01);
  hQSquareEp5 -> GetXaxis() -> SetTitle("Q^{2} [ GeV^{2} / c ]"); 
  hQSquareEp5 -> GetYaxis() -> SetTitle(" entries "); 
 
  hhQSquareScattAngleEp = new TH2F("hhQSquareScattAngleEp", "Q^{2} vs Scatt Angle", 1000, 0,8, 1000, 0, 0.1);
  hhQSquareScattAngleEp -> GetXaxis()->SetTitle(" scattering angle [degree] "); 
  hhQSquareScattAngleEp -> GetYaxis()->SetTitle(" Q^{2} [ GeV^{2} / c ]"); 
 
  hhEnergyVsAngle = new TH2F("hhEnergyVsAngle", "Energy Angle", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngle -> GetXaxis() -> SetTitle("  scattering angle [degree] "); 
  hhEnergyVsAngle -> GetYaxis() -> SetTitle("  energy [MeV] "); 

  hhEnergyVsAngleBeamLineCorrection = new TH2F("hhEnergyVsAngleBeamLineCorrection", "Energy Angle beam line correction", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleBeamLineCorrection -> GetXaxis() -> SetTitle("  scattering angle [degree] "); 
  hhEnergyVsAngleBeamLineCorrection -> GetYaxis() -> SetTitle("  energy [MeV] "); 


  hhEnergyVsAngleEp = new TH2F("hhEnergyVsAngleEp", "Energy Angle ep", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleEp -> GetXaxis() -> SetTitle("  scattering angle [degree] "); 
  hhEnergyVsAngleEp -> GetYaxis() -> SetTitle("  energy [MeV] "); 

  hhEnergyVsAngleMoller = new TH2F("hhEnergyVsAngleMoller", "Energy Angle Moller", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleMoller -> GetXaxis() -> SetTitle("  scattering angle [degree] "); 
  hhEnergyVsAngleMoller -> GetYaxis() -> SetTitle("  energy [MeV] "); 

  hhEnergyVsAngleMollerBeamLineCorrection = new TH2F("hhEnergyVsAngleMollerBeamLineCorrection", "Energy Angle Moller beam line correction", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleMollerBeamLineCorrection -> GetXaxis() -> SetTitle("  scattering angle [degree] "); 
  hhEnergyVsAngleMollerBeamLineCorrection  -> GetYaxis() -> SetTitle("  energy [MeV] "); 

  hhEnergyVsAngleHyCal2CGEM1C = new TH2F("hhEnergyVsAngleHyCal2CGEM1C", "Energy Angle > 2 cluster on hycal, 1 cluster on gem ", 1000, 0, 10, 4000, 0, 2500);
  hhEnergyVsAngleHyCal2CGEM1C -> GetXaxis() -> SetTitle("  scattering angle [degree] "); 
  hhEnergyVsAngleHyCal2CGEM1C -> GetYaxis() -> SetTitle("  energy [MeV] "); 
 
  hXOffsetFromMoller = new TH1F("hXOffsetFromMoller", "x offset using moller events", 1000, -50, 50);
  hXOffsetFromMoller -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hXOffsetFromMoller -> GetYaxis()->SetTitle(" entries "); 
 
  hYOffsetFromMoller = new TH1F("hYOffsetFromMoller", "y offset using moller events", 1000, -50, 50);
  hYOffsetFromMoller -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hYOffsetFromMoller -> GetYaxis()->SetTitle(" entries "); 
 
  hhMollerCenter = new TH2F("hhMollerCenter", "Moller Center", 1000, -50, 50, 1000, -50, 50);
  hhMollerCenter -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhMollerCenter -> GetYaxis() -> SetTitle(" x [mm] "); 
 
  hXOffsetFromSymmetricMoller = new TH1F("hXOffsetFromSymmetricMoller", "x offset using moller events", 1000, -50, 50);
  hXOffsetFromSymmetricMoller -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hXOffsetFromSymmetricMoller -> GetYaxis()->SetTitle(" entries "); 
 
  hYOffsetFromSymmetricMoller = new TH1F("hYOffsetFromSymmetricMoller", "y offset using moller events", 1000, -50, 50);
  hYOffsetFromSymmetricMoller -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hYOffsetFromSymmetricMoller -> GetYaxis()->SetTitle(" entries "); 
 
  hhSymmetricMollerCenter = new TH2F("hhSymmetricMollerCenter", "Moller Center", 1000, -50, 50, 1000, -50, 50);
  hhSymmetricMollerCenter -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhSymmetricMollerCenter -> GetYaxis() -> SetTitle(" x [mm] "); 
 
  hSymmetricMollerTheta = new TH1F("hSymmetricMollerTheta", "Symmetric Moller #Theta", 1000, -10, 10);
  hSymmetricMollerTheta -> GetXaxis() -> SetTitle(" #theta [ degree ] "); 
  hSymmetricMollerTheta -> GetYaxis() -> SetTitle(" entries "); 
 
  hXOffsetFromMollerAfterCorrection = new TH1F("hXOffsetFromMollerAfterCorrection", "x offset using moller events after correction", 1000, -50, 50);
  hXOffsetFromMollerAfterCorrection -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hXOffsetFromMollerAfterCorrection -> GetYaxis()->SetTitle(" entries "); 
 
  hYOffsetFromMollerAfterCorrection = new TH1F("hYOffsetFromMollerAfterCorrection", "y offset using moller events after correction", 1000, -50, 50);
  hYOffsetFromMollerAfterCorrection -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hYOffsetFromMollerAfterCorrection -> GetYaxis()->SetTitle(" entries "); 
 
  hhMollerCenterAfterCorrection = new TH2F("hhMollerCenterAfterCorrection", "Moller Center after correction", 1000, -50, 50, 1000, -50, 50);
  hhMollerCenterAfterCorrection -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhMollerCenterAfterCorrection -> GetYaxis() -> SetTitle(" x [mm] "); 
 
  hLinearDeviationMollerAfterCorrection = new TH1F("LinearDeviationMollerAfterCorrection", "Moller Clusters coplanarity after correction", 10000, -180, 180);
  hLinearDeviationMollerAfterCorrection -> GetXaxis() -> SetTitle("coplanarity [degree]");
  hLinearDeviationMollerAfterCorrection -> GetYaxis() -> SetTitle("entries");
 
  hXOffsetGEM = new TH1F("hXOffsetGEM", "x offset between two GEM chambers", 1000, -20, 20);
  hXOffsetGEM -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hXOffsetGEM -> GetYaxis()->SetTitle(" entries "); 
 
  hYOffsetGEM = new TH1F("hYOffsetGEM", "y offset between two GEM chambers", 1000, -20, 20);
  hYOffsetGEM -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hYOffsetGEM -> GetYaxis()->SetTitle(" entries "); 
 
  hXOffsetGEMAfterCorrection = new TH1F("hXOffsetGEMAfterCorrection", "x offset between two GEM chambers after correction", 1000, -20, 20);
  hXOffsetGEMAfterCorrection -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hXOffsetGEMAfterCorrection -> GetYaxis()->SetTitle(" entries "); 
 
  hYOffsetGEMAfterCorrection = new TH1F("hYOffsetGEMAfterCorrection", "y offset between two GEM chambers after correction", 1000, -20, 20);
  hYOffsetGEMAfterCorrection -> GetXaxis()->SetTitle(" x offset [mm] "); 
  hYOffsetGEMAfterCorrection -> GetYaxis()->SetTitle(" entries "); 
 
  // moller 2d position
  hhMoller2DRing = new TH2F("hhMoller2DRing"," Moller ", 1000, -650, 650, 1000, -650, 650);
  hhMoller2DRing -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhMoller2DRing -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  hhMoller2DRing2 = new TH2F("hhMoller2DRing2"," Moller ", 1000, -650, 650, 1000, -650, 650);
  hhMoller2DRing2 -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhMoller2DRing2 -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  hhMoller2DRingBeforeMatch = new TH2F("hhMoller2DRingBeforeMatch"," 2 clusters ", 1000, -650, 650, 1000, -650, 650);
  hhMoller2DRingBeforeMatch -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhMoller2DRingBeforeMatch -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  hhMoller2DRingBeforeMatch2 = new TH2F("hhMoller2DRingBeforeMatch2"," 2 clusters  ", 1000, -650, 650, 1000, -650, 650);
  hhMoller2DRingBeforeMatch2 -> GetXaxis() -> SetTitle(" x [mm] "); 
  hhMoller2DRingBeforeMatch2 -> GetYaxis() -> SetTitle(" y [mm] "); 
 
  // hycal gem position match
  hXDiffMatch = new TH1F("hXDiffMatch", "x_gem - x_hycal", 1000, -60, 60);
  hXDiffMatch -> GetXaxis() -> SetTitle(" #Delta X [mm] "); 
  hXDiffMatch -> GetYaxis() -> SetTitle(" entries "); 

  hYDiffMatch = new TH1F("hYDiffMatch", "y_gem - y_hycal", 1000, -60, 60);
  hYDiffMatch -> GetXaxis() -> SetTitle(" #Delta Y [mm] "); 
  hYDiffMatch -> GetYaxis() -> SetTitle(" entries "); 

  hhXDiffMatchVsEnergy = new TH2F("hhXDiffMatch", "x diff. vs energy", 2500, 0, 2500, 1000, -60, 60);
  hhXDiffMatchVsEnergy -> GetXaxis() -> SetTitle(" #Delta X [mm] "); 
  hhXDiffMatchVsEnergy -> GetYaxis() -> SetTitle(" energy [MeV]  "); 

  hhYDiffMatchVsEnergy = new TH2F("hhYDiffMatch", "y diff. vs energy", 2500, 0, 2500, 1000, -60, 60);
  hhYDiffMatchVsEnergy -> GetXaxis() -> SetTitle(" #Delta Y [mm] "); 
  hhYDiffMatchVsEnergy -> GetYaxis() -> SetTitle(" energy [MeV]  "); 

  hNbPointsMatch = new TH1F("hNbPointsMatch", "# cluster match", 100, -1.5, 98.5);
  hNbPointsMatch -> GetXaxis() -> SetTitle("clusuter quantity after match per event"); 
  hNbPointsMatch -> GetYaxis() -> SetTitle("entries"); 

}
