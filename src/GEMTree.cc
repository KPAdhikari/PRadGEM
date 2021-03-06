//========================================================================//
//                                                                        //
// Xinzhan Bai                                                            //
// 08/05/2016                                                             //
//========================================================================//
#include "GEMTree.h"
#include "TTree.h"
#include "TFile.h"
#include <cassert>
#include <iostream>
#include "GEMDataStruct.h"
#include <vector>
#include "PRadMoller.h"
#include "PRadEP.h"
#include "TH1F.h"
#include "MollerGEMSpatialRes.h"
#include "GEMConfigure.h"

using namespace std;

// ------------------------------- common setting --------------------------------------
GEMTree::GEMTree()
{
    TH1::AddDirectory(kFALSE);
    file = new TFile("root_file/res.root", "recreate");
    event_id = 0;
    empty_moller_event = true;
    empty_ep_event = true;
    empty_epics_event = true;
    empty_gem_event = true;
}

void GEMTree::WriteToDisk()
{
    file->Write();
    //file->Save();
}

GEMTree::~GEMTree()
{
}

void GEMTree::SetGEMConfigure(GEMConfigure * c)
{
    configure = c;

    // init res tree
    this -> InitGEMTree(2); // 2 gem detectors
    this -> InitPhysicsTree();
    this -> InitEpicsTree();
    this -> InitCaliOffsetTree();
    this -> InitProdOffsetTree();
    this -> InitOverlapTree();
}


// ------------------------------- gem tree -------------------------------------
void GEMTree::SetEventID(unsigned int &id)
{
    event_id = (int) id;
}

void GEMTree::InitGEMTree(int ndet)
{
    gem_tree = new TTree("T", "res");
    gem_tree->SetDirectory(file);

    assert( ndet <= NDETECTOR);

    // event id
    gem_tree->Branch("event_id", &evt_id, "event_id/I");

    // gem information
    for(int i=0;i<ndet;i++)
    {
	gem_tree->Branch(Form("nhits%d",i),&nhits[i], Form("nhits%d/I",i));
	gem_tree->Branch(Form("x%d",i), &x[i], Form("x%d[nhits%d]/F",i,i));
	gem_tree->Branch(Form("y%d",i), &y[i], Form("y%d[nhits%d]/F",i,i));
	gem_tree->Branch(Form("x%d_charge",i), &x_charge[i], Form("x%d_charge[nhits%d]",i,i));
	gem_tree->Branch(Form("y%d_charge",i), &y_charge[i], Form("y%d_charge[nhits%d]",i,i));
	gem_tree->Branch(Form("energy%d",i), &energy[i], Form("energy%d[nhits%d]",i,i));
	gem_tree->Branch(Form("z%d",i), &z[i], Form("z%d[nhits%d]",i,i));
	gem_tree->Branch(Form("x%d_size",i), &x_size[i], Form("x%d_size[nhits%d]",i,i));
	gem_tree->Branch(Form("y%d_size",i), &y_size[i], Form("y%d_size[nhits%d]",i,i));
    }

    // tdc information
    gem_tree->Branch("n_S2", &n_S2, "n_S2/I");
    gem_tree->Branch("n_S1", &n_S1, "n_S1/I");
    gem_tree->Branch("TDCS2", &TDCS2, "TDCS2[n_S2]/F");
    gem_tree->Branch("TDCS1", &TDCS1, "TDCS1[n_S1]/F");
    gem_tree->Branch("n_hycal", &n_hycal, "n_hycal/I");
    gem_tree->Branch("TDCHyCal", &TDCHyCal, "TDCHyCal[n_hycal]/F");

    // hycal information
    gem_tree -> Branch("hycal_nhits", &hycal_nhits, "hycal_nhits/I");
    gem_tree -> Branch("hycal_x", &hycal_x, "hycal_x[n_hycal]/F");
    gem_tree -> Branch("hycal_y", &hycal_y, "hycal_y[n_hycal]/F");
    gem_tree -> Branch("hycal_e", &hycal_e, "hycal_e[n_hycal]/F");

    InitTDCGroup();
}

void GEMTree::InitTDCGroup()
{
    hycal_group_q = configure->TDC_Quan;

    for(int i=0;i<hycal_group_q;i++)
    {
	hycal_group[i] = configure->TDC[i];
    }

    scin_group[0] = "S1";
    scin_group[1] = "S2";
}

void GEMTree::PushCalibrationData(vector<GEMClusterStruct> &gem1,
	vector<GEMClusterStruct> &gem2,
	unordered_multimap<string, double> &tdc_map,
	vector<HyCalHit> *hycalhit)
{
    PushDetector(0, gem1);
    PushDetector(1, gem2);
    PushTDCValue(tdc_map);
    PushHyCalData(hycalhit);
}


void GEMTree::PushDetector(int nthDet, std::vector<GEMClusterStruct> Gem)
{
    int index = 0;
    int n = nthDet;
    int nh = Gem.size();
    if(nh != 0)
	empty_gem_event = false;

    nhits[n] = nh;
    for(auto &i: Gem)
    {
	x[n][index]=i.x;
	y[n][index] = i.y;
	x_charge[n][index] = i.x_charge;
	y_charge[n][index] = i.y_charge;
	energy[n][index] = i.energy;
	z[n][index] = i.z;
	x_size[n][index] = i.x_size;
	y_size[n][index] = i.y_size;
	index++;
    }
}

void GEMTree::PushTDCValue(unordered_multimap<string, double> &tdc_map)
{
    vector<float> scin_val[2];

    for(int i=0;i<2;i++)
    {
	auto range = tdc_map.equal_range( (string)scin_group[i] );
	for(auto it=range.first;it!=range.second;++it)
	    scin_val[i].push_back( (*it).second );
    }

    vector<float> hycal_val[4];

    for(int i=0;i<hycal_group_q;i++)
    {
	auto range = tdc_map.equal_range( (string)hycal_group[i] );
	for(auto it=range.first;it!=range.second;++it)
	    hycal_val[i].push_back( (*it).second );
    }

    int index = 0;
    // fill scintillator tdc
    n_S1 = scin_val[0].size();
    for(auto &ii: scin_val[0])
    {
	TDCS1[index] = ii;
	index++;
    }
    index = 0;
    n_S2 = scin_val[1].size();
    for(auto &ii: scin_val[1])
    {
	TDCS2[index] = ii;
	index++;
    }
    // fill hycal tdc
    index = 0; n_hycal = 0;
    for(int i=0;i<hycal_group_q;i++)
    {
	n_hycal += hycal_val[i].size();
	for(auto &ii: hycal_val[i])
	{
	    TDCHyCal[index] = ii;
	    index++;
	}
    }
}

void GEMTree::PushHyCalData(vector<HyCalHit> *hycal_hit)
{
    int index = 0;
    hycal_nhits = hycal_hit->size();
    for(auto &i: *hycal_hit)
    {
        hycal_x[index] = i.x;
	hycal_y[index] = i.y;
	hycal_e[index] = i.E;
	index++;
    }
}

void GEMTree::FillGEMTree()
{
    if( ! empty_gem_event){
	gem_tree->Fill();
    }
    empty_gem_event = true;
}


// ------------------------------- physics tree -------------------------------------
void GEMTree::InitPhysicsTree()
{
    moller_tree = new TTree("moller_tree", "moller tree");
    moller_tree -> SetDirectory(file);
    ep_tree = new TTree("ep_tree", "ep tree");
    ep_tree -> SetDirectory(file);
    ep_moller_tree = new TTree("ep_moller_tree", "ep moller tree");
    ep_moller_tree->SetDirectory(file);
    // moller 
    moller_tree->Branch("evt_id", &moller_data.event_id, "evt_id/I");
    moller_tree->Branch("chamber_id1", &moller_data.chamber_id1, "chamber_id1/I");
    moller_tree->Branch("x1", &moller_data.x1, "x1/F");
    moller_tree->Branch("y1", &moller_data.y1, "y1/F");
    moller_tree->Branch("x1_hycal", &moller_data.x1_hycal, "x1_hycal/F");
    moller_tree->Branch("y1_hycal", &moller_data.y1_hycal, "y1_hycal/F");
    moller_tree->Branch("e1", &moller_data.e1, "e1/F");
    moller_tree->Branch("angle1", &moller_data.angle1, "angle1/F");
    moller_tree->Branch("chamber_id2", &moller_data.chamber_id2, "chamber_id2/I");
    moller_tree->Branch("x2", &moller_data.x2, "x2/F");
    moller_tree->Branch("y2", &moller_data.y2, "y2/F");
    moller_tree->Branch("x2_hycal", &moller_data.x2_hycal, "x2_hycal/F");
    moller_tree->Branch("y2_hycal", &moller_data.y2_hycal, "y2_hycal/F");
    moller_tree->Branch("e2", &moller_data.e2, "e2/F");
    moller_tree->Branch("angle2", &moller_data.angle2, "angle2/F");
    moller_tree->Branch("coplanarity", &moller_data.coplanarity, "coplanarity/F");
    moller_tree->Branch("moller_center_x", &moller_center_x, "moller_center_x/F");
    moller_tree->Branch("moller_center_y", &moller_center_y, "moller_center_y/F");
    moller_tree->Branch("moller_pos_res_dx", &moller_pos_res_dx, "moller_pos_res_dx/F");
    moller_tree->Branch("moller_pos_res_dy", &moller_pos_res_dy, "moller_pos_res_dy/F");
    // ep
    ep_tree -> Branch("evt_id", &ep_data.event_id, "evt_id/I");
    ep_tree -> Branch("chamber_id", &ep_data.chamber_id, "chamber_id/I");
    ep_tree -> Branch("x", &ep_data.x, "x/F");
    ep_tree -> Branch("y", &ep_data.y, "y/F");
    ep_tree -> Branch("x_hycal", &ep_data.x_hycal, "x_hycal/F");
    ep_tree -> Branch("y_hycal", &ep_data.y_hycal, "y_hycal/F");
    ep_tree -> Branch("e", &ep_data.e, "e/F");
    ep_tree -> Branch("angle", &ep_data.angle, "angle/F");
    ep_tree -> Branch("q_square", &q_square, "q_square/F");
    // moller + ep
    ep_moller_tree -> Branch("evt_id", &evt_id, "evt_id/I");
    ep_moller_tree -> Branch("nCluster", &nCluster, "nCluster/I");
    ep_moller_tree -> Branch("det_id", detector_id, "det_id[nCluster]/I");
    ep_moller_tree -> Branch("scatt_energy", scatt_energy, "scatt_energy[nCluster]/F");
    ep_moller_tree -> Branch("scatt_angle", scatt_angle, "scatt_angle[nCluster]/F");
    ep_moller_tree -> Branch("scatt_x", scatt_x, "scatt_x[nCluster]/F");
    ep_moller_tree -> Branch("scatt_y", scatt_y, "scatt_y[nCluster]/F");
}

void GEMTree::PushMoller(PRadMoller * moller)
{
    vector<pair<double, double> > ea = moller->EnergyAngle();
    vector<HyCalHit> hycal_hit = moller->GetHyCalMatch();

    if(ea.size() == 2) {
	empty_moller_event = false;
    }
    else {
	return;
    }

    // moller
    moller_data.event_id = moller->GetEvtID();
    moller_data.coplanarity = moller->Coplanarity();

    vector<pair<int, pair<double, double> > > pos = moller->Positions();
    moller_data.chamber_id1 = pos[0].first;
    moller_data.x1 = pos[0].second.first;
    moller_data.y1 = pos[0].second.second;
    moller_data.x1_hycal = hycal_hit[0].x;
    moller_data.y1_hycal = hycal_hit[0].y;
    moller_data.e1 = ea[0].first;
    moller_data.angle1 = ea[0].second;
    moller_data.chamber_id2 = pos[1].first;
    moller_data.x2 = pos[1].second.first;
    moller_data.y2 = pos[1].second.second;
    moller_data.x2_hycal = hycal_hit[1].x;
    moller_data.y2_hycal = hycal_hit[1].y;
    moller_data.e2 = ea[1].first;
    moller_data.angle2 = ea[1].second;

    moller_center_x = moller->MollerCenter().first;
    moller_center_y = moller->MollerCenter().second;
    pair<double, double> reso = moller->GetSpatialResHandler()->GetDiff();
    moller_pos_res_dx = reso.first;
    moller_pos_res_dy = reso.second;

    // moller + ep
    nCluster = 2;
    for(int i=0;i<nCluster;i++)
    {
	detector_id[i] = pos[i].first;
	scatt_x[i] = pos[i].second.first;
	scatt_y[i] = pos[i].second.second;
	scatt_energy[i] = ea[i].first;
	scatt_angle[i] = ea[i].second;
    }
}

void GEMTree::PushEP(PRadEP * ep)
{
    vector<pair<double, double> > ea = ep->EnergyAngle();
    vector<HyCalHit> hycal_hit = ep->GetHyCalMatch();

    if(ea.size() == 1 ) {
	empty_ep_event = false;
    }
    else {
	return;
    }

    ep_data.event_id = ep->GetEvtID();
    ep_data.chamber_id = ep->GetChamberID();
    vector<pair<double, double> > pos = ep->Positions();
    ep_data.x = pos[0].first;
    ep_data.y = pos[0].second;
    ep_data.x_hycal = hycal_hit[0].x;
    ep_data.y_hycal = hycal_hit[0].y;
    ep_data.e = ea[0].first;
    ep_data.angle = ea[0].second;

    q_square = ep -> QSquare();

    // moller + ep
    nCluster = 1;
    detector_id[0] = ep->GetChamberID();
    scatt_x[0] = pos[0].first;
    scatt_y[0] = pos[0].second;
    scatt_energy[0] = ea[0].first;
    scatt_angle[0] = ea[0].second;
}

void GEMTree::FillMollerTree()
{
    if( ! empty_moller_event) {
	moller_tree->Fill();
	ep_moller_tree->Fill();
    }
    empty_moller_event = true;
}

void GEMTree::FillEPTree()
{
    if( ! empty_ep_event) {
	ep_tree->Fill();
	ep_moller_tree->Fill();
    }
    empty_ep_event = true;
}


// ------------------------------- epics tree -------------------------------------
void GEMTree::InitEpicsTree()
{
    //file = new TFile("root_file/epics_res.root", "recreate");
    epics_tree = new TTree("epics_tree", "epics tree");
    epics_tree -> SetDirectory(file);

    string expression("");

    expression+="MFC_Flow/D:Cell_Gas_T:Tank_P_P:Chamber_P:";
    expression+="Cell_P:Cell_Body_T:hallb_ptrans_y2_encoder:";
    expression+="hallb_ptrans_y1_encoder: hallb_ptrans_x_encoder:";
    expression+="ptrans_y:ptrans_x:hallb_IPM2H01_CUR:";
    expression+="hallb_IPM2H01_YPOS:hallb_IPM2H01_XPOS:";
    expression+="hallb_IPM2C24A_CUR:hallb_IPM2C24A_XPOS:";
    expression+="hallb_IPM2C21A_CUR:hallb_IPM2C21A_YPOS:";
    expression+="hallb_IPM2C21A_XPOS:hallb_IPM2C24A_YPOS:";
    expression+="scaler_calc1:VCG2H02A:VCG2H01A:VCG2H00A:";
    expression+="VCG2C24A:VCG2C21A:VCG2C21:MBSY2C_energy";

    epics_tree->Branch("epics", &_epics_tree, expression.c_str());
    epics_tree->Branch("evtID", &event_id, "evtID/I");
}

void GEMTree::PushEpics(unordered_map<string, double> & epics_map)
{
    //memset(&_epics_tree, 0, sizeof _epics_tree);
    _epics_tree.Clear();
    if(epics_map.size() == 0)
	return;

    if(epics_map.find("TGT:PRad:MFC_Flow") != epics_map.end())
	_epics_tree.MFC_Flow = epics_map["TGT:PRad:MFC_Flow"];
    if(epics_map.find("TGT:PRad:Cell_Gas_T") != epics_map.end() )
	_epics_tree.Cell_Gas_T = epics_map["TGT:PRad:Cell_Gas_T"];
    if(epics_map.find("TGT:PRad:Tank_P_P") != epics_map.end() )
	_epics_tree.Tank_P_P = epics_map["TGT:PRad:Tank_P_P"];
    if(epics_map.find("TGT:PRad:Chamber_P") != epics_map.end() )
	_epics_tree.Chamber_P = epics_map["TGT:PRad:Chamber_P"];
    if(epics_map.find("TGT:PRad:Cell_P") != epics_map.end() )
	_epics_tree.Cell_P = epics_map["TGT:PRad:Cell_P"];
    if(epics_map.find("TGT:PRad:Cell_Body_T") != epics_map.end() )
	_epics_tree.Cell_Body_T = epics_map["TGT:PRad:Cell_Body_T"];
    if(epics_map.find("hallb_ptrans_y2_encoder") != epics_map.end() )
	_epics_tree.hallb_ptrans_y2_encoder = epics_map["hallb_ptrans_y2_encoder"];
    if(epics_map.find("hallb_ptrans_y1_encoder") != epics_map.end() )
	_epics_tree.hallb_ptrans_y1_encoder = epics_map["hallb_ptrans_y1_encoder"];
    if(epics_map.find("hallb_ptrans_x_encoder") != epics_map.end() )
	_epics_tree.hallb_ptrans_x_encoder = epics_map["hallb_ptrans_x_encoder"];
    if(epics_map.find("ptrans_y") != epics_map.end() )
	_epics_tree.ptrans_y = epics_map["ptrans_y"];
    if(epics_map.find("ptrans_x") != epics_map.end() )
	_epics_tree.ptrans_x = epics_map["ptrans_x"];
    if(epics_map.find("hallb_IPM2H01_CUR") != epics_map.end() )
	_epics_tree.hallb_IPM2H01_CUR = epics_map["hallb_IPM2H01_CUR"];
    if(epics_map.find("hallb_IPM2H01_YPOS") != epics_map.end() )
	_epics_tree.hallb_IPM2H01_YPOS = epics_map["hallb_IPM2H01_YPOS"];
    if(epics_map.find("hallb_IPM2H01_XPOS") != epics_map.end() )
	_epics_tree.hallb_IPM2H01_XPOS = epics_map["hallb_IPM2H01_XPOS"];
    if(epics_map.find("hallb_IPM2C24A_CUR") != epics_map.end() )
	_epics_tree.hallb_IPM2C24A_CUR = epics_map["hallb_IPM2C24A_CUR"];
    if(epics_map.find("hallb_IPM2C24A_XPOS") != epics_map.end() )
	_epics_tree.hallb_IPM2C24A_XPOS = epics_map["hallb_IPM2C24A_XPOS"];
    if(epics_map.find("hallb_IPM2C21A_CUR") != epics_map.end() )
	_epics_tree.hallb_IPM2C21A_CUR = epics_map["hallb_IPM2C21A_CUR"];
    if(epics_map.find("hallb_IPM2C21A_YPOS") != epics_map.end() )
	_epics_tree.hallb_IPM2C21A_YPOS = epics_map["hallb_IPM2C21A_YPOS"];
    if(epics_map.find("hallb_IPM2C21A_XPOS") != epics_map.end() )
	_epics_tree.hallb_IPM2C21A_XPOS = epics_map["hallb_IPM2C21A_XPOS"];
    if(epics_map.find("hallb_IPM2C24A_YPOS") != epics_map.end() )
	_epics_tree.hallb_IPM2C24A_YPOS = epics_map["hallb_IPM2C24A_YPOS"];
    if(epics_map.find("scaler_calc1") != epics_map.end() )
	_epics_tree.scaler_calc1 = epics_map["scaler_calc1"];
    if(epics_map.find("VCG2H02A") != epics_map.end() )
	_epics_tree.VCG2H02A = epics_map["VCG2H02A"];
    if(epics_map.find("VCG2H01A") != epics_map.end() )
	_epics_tree.VCG2H01A = epics_map["VCG2H01A"];
    if(epics_map.find("VCG2H00A") != epics_map.end() )
	_epics_tree.VCG2H00A = epics_map["VCG2H00A"];
    if(epics_map.find("VCG2C24A") != epics_map.end() )
	_epics_tree.VCG2C24A = epics_map["VCG2C24A"];
    if(epics_map.find("VCG2C21A") != epics_map.end() )
	_epics_tree.VCG2C21A = epics_map["VCG2C21A"];
    if(epics_map.find("VCG2C21") != epics_map.end() )
	_epics_tree.VCG2C21 = epics_map["VCG2C21"];
    if(epics_map.find("MBSY2C_energy") != epics_map.end() )
	_epics_tree.MBSY2C_energy = epics_map["MBSY2C_energy"];
}

void GEMTree::FillEpicsTree()
{
    if( _epics_tree.MBSY2C_energy == -1000 )
	return;

    epics_tree->Fill();
}

// ------------------------------- calibration offset tree -------------------------------------
void GEMTree::InitCaliOffsetTree()
{
    cali_offset_tree = new TTree("cali_offset_tree", "calibration offset tree");
    cali_offset_tree->SetDirectory(file);
    cali_offset_tree->Branch("cali_x_offset", &cali_x_offset, "cali_x_offset/D");
    cali_offset_tree->Branch("cali_y_offset", &cali_y_offset, "cali_y_offset/D");
}

void GEMTree::PushCaliOffset(double x, double y)
{
    cali_x_offset = x;
    cali_y_offset = y;
}

void GEMTree::FillCaliOffsetTree()
{
    cali_offset_tree->Fill();
}

// ------------------------------- producton offset tree -------------------------------------
// overlapping area method
void GEMTree::InitProdOffsetTree()
{
    prod_offset_tree1 = new TTree("prod_offset_tree1", "production offset tree 1");
    prod_offset_tree2 = new TTree("prod_offset_tree2", "production offset tree 2");
    prod_offset_tree_res = new TTree("prod_offset_tree_res", "production offset tree res");
    prod_offset_tree1 -> SetDirectory(file);
    prod_offset_tree2 -> SetDirectory(file);
    prod_offset_tree_res -> SetDirectory(file);

    empty_prod_gem1 = true;
    empty_prod_gem2 = true;

    prod_offset_tree1 -> Branch("evt_id", &evt_id, "evt_id/I");
    prod_offset_tree1 -> Branch("prod_gem1_dx", &prod_gem1_dx, "prod_gem1_dx/D");
    prod_offset_tree1 -> Branch("prod_gem1_dy", &prod_gem1_dy, "prod_gem1_dy/D");
    prod_offset_tree1 -> Branch("prod_gem1_ncluster", &prod_gem1_ncluster, "prod_gem1_ncluster/I");
    prod_offset_tree1 -> Branch("prod_gem1_scattx", prod_gem1_scattx, "prod_gem1_scattx[prod_gem1_ncluster]/D");
    prod_offset_tree1 -> Branch("prod_gem1_scatty", prod_gem1_scatty, "prod_gem1_scatty[prod_gem1_ncluster]/D");
    prod_offset_tree1 -> Branch("prod_gem1_energy1", &prod_gem1_energy1, "prod_gem1_energy1/D");
    prod_offset_tree1 -> Branch("prod_gem1_energy2", &prod_gem1_energy2, "prod_gem1_energy2/D");
    prod_offset_tree1 -> Branch("prod_gem1_coplanarity", &prod_gem1_coplanarity, "prod_gem1_coplanarity/D");
    prod_offset_tree1 -> Branch("prod_gem1_open_angle", &prod_gem1_open_angle, "prod_gem1_open_angle/D");
    prod_offset_tree1 -> Branch("prod_gem1_angle1", &prod_gem1_angle1, "prod_gem1_angle1/D");
    prod_offset_tree1 -> Branch("prod_gem1_angle2", &prod_gem1_angle2, "prod_gem1_angle2/D");

    prod_offset_tree2 -> Branch("evt_id", &evt_id, "evt_id/I");
    prod_offset_tree2 -> Branch("prod_gem2_dx", &prod_gem2_dx, "prod_gem2_dx/D");
    prod_offset_tree2 -> Branch("prod_gem2_dy", &prod_gem2_dy, "prod_gem2_dy/D");
    prod_offset_tree2 -> Branch("prod_gem2_ncluster", &prod_gem2_ncluster, "prod_gem2_ncluster/I");
    prod_offset_tree2 -> Branch("prod_gem2_scattx", prod_gem2_scattx, "prod_gem2_scattx[prod_gem2_ncluster]/D");
    prod_offset_tree2 -> Branch("prod_gem2_scatty", prod_gem2_scatty, "prod_gem2_scatty[prod_gem2_ncluster]/D");
    prod_offset_tree2 -> Branch("prod_gem2_energy1", &prod_gem2_energy1, "prod_gem2_energy1/D");
    prod_offset_tree2 -> Branch("prod_gem2_energy2", &prod_gem2_energy2, "prod_gem2_energy2/D");
    prod_offset_tree2 -> Branch("prod_gem2_coplanarity", &prod_gem2_coplanarity, "prod_gem2_coplanarity/D");
    prod_offset_tree2 -> Branch("prod_gem2_open_angle", &prod_gem2_open_angle, "prod_gem2_open_angle/D");
    prod_offset_tree2 -> Branch("prod_gem2_angle1", &prod_gem2_angle1, "prod_gem2_angle1/D");
    prod_offset_tree2 -> Branch("prod_gem2_angle2", &prod_gem2_angle2, "prod_gem2_angle2/D");

    prod_offset_tree_res -> Branch("evt_id", &evt_id, "evt_id/I");
    prod_offset_tree_res -> Branch("prod_offset_x", &prod_offset_x, "prod_offset_x/D");
    prod_offset_tree_res -> Branch("prod_offset_y", &prod_offset_y, "prod_offset_y/D");
}

void GEMTree::PushProdOffset(int i, PRadMoller * moller)
{
    vector<pair<double, double> > ea = moller->EnergyAngle();
    if( ea.size() == 2 ){
	evt_id = moller -> GetEvtID();
	if( i==0){
	    empty_prod_gem1 = false;
	    prod_gem1_ncluster=2;
	    prod_gem1_dx = moller->MollerCenter().first;
	    prod_gem1_dy = moller->MollerCenter().second;
	    prod_gem1_energy1 = ea[0].first;
	    prod_gem1_energy2 = ea[1].first;
	    prod_gem1_angle1 = ea[0].second;
	    prod_gem1_angle2 = ea[1].second;
	    prod_gem1_coplanarity = moller->Coplanarity();
	    prod_gem1_open_angle = moller->OpenAngle();
	    vector<pair<int, pair<double, double> > > pos = moller->Positions();
	    prod_gem1_scattx[0] = pos[0].second.first;
	    prod_gem1_scatty[0] = pos[0].second.second;
	    prod_gem1_scattx[1] = pos[1].second.first;
	    prod_gem1_scatty[1] = pos[1].second.second;
	}
	else if( i==1){
	    empty_prod_gem2 = false;
	    prod_gem2_ncluster=2;
	    prod_gem2_dx = moller->MollerCenter().first;
	    prod_gem2_dy = moller->MollerCenter().second;
	    prod_gem2_energy1 = ea[0].first;
	    prod_gem2_energy2 = ea[1].first;
	    prod_gem2_angle1 = ea[0].second;
	    prod_gem2_angle2 = ea[1].second;
	    prod_gem2_coplanarity = moller->Coplanarity();
	    prod_gem2_open_angle = moller->OpenAngle();
	    vector<pair<int, pair<double, double> > > pos = moller->Positions();
	    prod_gem2_scattx[0] = pos[0].second.first;
	    prod_gem2_scatty[0] = pos[0].second.second;
	    prod_gem2_scattx[1] = pos[1].second.first;
	    prod_gem2_scatty[1] = pos[1].second.second;
	}
    }
}

void GEMTree::FillProdOffsetTree()
{
    if( (!empty_prod_gem1) && (!empty_prod_gem2)){
	prod_offset_x = -(prod_gem1_dx - prod_gem2_dx);
	prod_offset_y = -(prod_gem1_dy - prod_gem2_dy);
	prod_offset_tree_res->Fill();
    }
    if( ! empty_prod_gem1)
	prod_offset_tree1 -> Fill();
    if( ! empty_prod_gem2)
	prod_offset_tree2 -> Fill();

    empty_prod_gem1 = true;
    empty_prod_gem2 = true;
}


// ------------------------------- overlap area tree -------------------------------------
void GEMTree::InitOverlapTree()
{
    olp_empty_moller_event = true;
    olp_empty_ep_event = true;

    overlap_moller_tree = new TTree("overlap_moller_tree", "overlap moller tree");
    overlap_moller_tree -> SetDirectory(file);
    overlap_ep_tree = new TTree("overlap_ep_tree", "overlap ep tree");
    overlap_ep_tree -> SetDirectory(file);

    //gem1
    overlap_moller_tree -> Branch("gem1.event_id", &olp_moller_data1.event_id, "gem1.event_id/I");
    overlap_moller_tree -> Branch("gem1.chamber_id1", &olp_moller_data1.chamber_id1, "gem1.chamber_id1/I");
    overlap_moller_tree -> Branch("gem1.x1", &olp_moller_data1.x1, "gem1.x1/F");
    overlap_moller_tree -> Branch("gem1.y1", &olp_moller_data1.y1, "gem1.y1/F");
    overlap_moller_tree -> Branch("gem1.x1_hycal", &olp_moller_data1.x1_hycal, "gem1.x1_hycal/F");
    overlap_moller_tree -> Branch("gem1.y1_hycal", &olp_moller_data1.y1_hycal, "gem1.y1_hycal/F");
    overlap_moller_tree -> Branch("gem1.e1", &olp_moller_data1.e1, "gem1.e1/F");
    overlap_moller_tree -> Branch("gem1.angle1", &olp_moller_data1.angle1, "gem1.angle1/F");
    overlap_moller_tree -> Branch("gem1.chamber_id2", &olp_moller_data1.chamber_id2, "gem1.chamber_id2/I");
    overlap_moller_tree -> Branch("gem1.x2", &olp_moller_data1.x2, "gem1.x2/F");
    overlap_moller_tree -> Branch("gem1.y2", &olp_moller_data1.y2, "gem1.y2/F");
    overlap_moller_tree -> Branch("gem1.x2_hycal", &olp_moller_data1.x2_hycal, "gem1.x2_hycal/F");
    overlap_moller_tree -> Branch("gem1.y2_hycal", &olp_moller_data1.y2_hycal, "gem1.y2_hycal/F");
    overlap_moller_tree -> Branch("gem1.e2", &olp_moller_data1.e2, "gem1.e2/F");
    overlap_moller_tree -> Branch("gem1.angle2", &olp_moller_data1.angle2, "gem1.angle2/F");
    overlap_moller_tree -> Branch("gem1.coplanarity", &olp_moller_data1.coplanarity, "gem1.coplanarity/F");
    //gem2
    overlap_moller_tree -> Branch("gem2.event_id", &olp_moller_data2.event_id, "gem2.event_id/I");
    overlap_moller_tree -> Branch("gem2.chamber_id1", &olp_moller_data2.chamber_id1, "gem2.chamber_id1/I");
    overlap_moller_tree -> Branch("gem2.x1", &olp_moller_data2.x1, "gem2.x1/F");
    overlap_moller_tree -> Branch("gem2.y1", &olp_moller_data2.y1, "gem2.y1/F");
    overlap_moller_tree -> Branch("gem2.x1_hycal", &olp_moller_data2.x1_hycal, "gem2.x1_hycal/F");
    overlap_moller_tree -> Branch("gem2.y1_hycal", &olp_moller_data2.y1_hycal, "gem2.y1_hycal/F");
    overlap_moller_tree -> Branch("gem2.e1", &olp_moller_data2.e1, "gem2.e1/F");
    overlap_moller_tree -> Branch("gem2.angle1", &olp_moller_data2.angle1, "gem2.angle1/F");
    overlap_moller_tree -> Branch("gem2.chamber_id2", &olp_moller_data2.chamber_id2, "gem2.chamber_id2/I");
    overlap_moller_tree -> Branch("gem2.x2", &olp_moller_data2.x2, "gem2.x2/F");
    overlap_moller_tree -> Branch("gem2.y2", &olp_moller_data2.y2, "gem2.y2/F");
    overlap_moller_tree -> Branch("gem2.x2_hycal", &olp_moller_data2.x2_hycal, "gem2.x2_hycal/F");
    overlap_moller_tree -> Branch("gem2.y2_hycal", &olp_moller_data2.y2_hycal, "gem2.y2_hycal/F");
    overlap_moller_tree -> Branch("gem2.e2", &olp_moller_data2.e2, "gem2.e2/F");
    overlap_moller_tree -> Branch("gem2.angle2", &olp_moller_data2.angle2, "gem2.angle2/F");
    overlap_moller_tree -> Branch("gem2.coplanarity", &olp_moller_data2.coplanarity, "gem2.coplanarity/F");
    //gem 1
    overlap_ep_tree -> Branch("gem1.event_id", &olp_ep_data1.event_id, "gem1.event_id/I");
    overlap_ep_tree -> Branch("gem1.chamber_id", &olp_ep_data1.chamber_id, "gem1.chamber_id/I");
    overlap_ep_tree -> Branch("gem1.x", &olp_ep_data1.x, "gem1.x/F");
    overlap_ep_tree -> Branch("gem1.y", &olp_ep_data1.y, "gem1.y/F");
    overlap_ep_tree -> Branch("gem1.x_hycal", &olp_ep_data1.x_hycal, "gem1.x_hycal/F");
    overlap_ep_tree -> Branch("gem1.y_hycal", &olp_ep_data1.y_hycal, "gem1.y_hycal/F");
    overlap_ep_tree -> Branch("gem1.e", &olp_ep_data1.e, "gem1.e/F");
    overlap_ep_tree -> Branch("gem1.angle", &olp_ep_data1.angle, "gem1.angle/F");
    //gem2
    overlap_ep_tree -> Branch("gem2.event_id", &olp_ep_data2.event_id, "gem2.event_id/I");
    overlap_ep_tree -> Branch("gem2.chamber_id", &olp_ep_data2.chamber_id, "gem2.chamber_id/I");
    overlap_ep_tree -> Branch("gem2.x", &olp_ep_data2.x, "gem2.x/F");
    overlap_ep_tree -> Branch("gem2.y", &olp_ep_data2.y, "gem2.y/F");
    overlap_ep_tree -> Branch("gem2.x_hycal", &olp_ep_data2.x_hycal, "gem2.x_hycal/F");
    overlap_ep_tree -> Branch("gem2.y_hycal", &olp_ep_data2.y_hycal, "gem2.y_hycal/F");
    overlap_ep_tree -> Branch("gem2.e", &olp_ep_data2.e, "gem2.e/F");
    overlap_ep_tree -> Branch("gem2.angle", &olp_ep_data2.angle, "gem2.angle/F");
}

void GEMTree::PushOverlapMollerTree(PRadMoller *moller1, PRadMoller * moller2)
{
    vector<pair<double, double> > ea1 = moller1->EnergyAngle();
    vector<pair<double, double> > ea2 = moller2->EnergyAngle();
    vector<HyCalHit> hycal_gem1 = moller1->GetHyCalMatch();
    vector<HyCalHit> hycal_gem2 = moller2->GetHyCalMatch();

    if(ea1.size() == 2 && ea2.size()==2) {
	olp_empty_moller_event = false;
    }
    else {
	return;
    }
    // gem1
    olp_moller_data1.event_id = moller1->GetEvtID();
    olp_moller_data1.coplanarity = moller1->Coplanarity();

    vector<pair<int, pair<double, double> > > pos1 = moller1->Positions();
    olp_moller_data1.chamber_id1 = pos1[0].first;
    olp_moller_data1.x1 = pos1[0].second.first;
    olp_moller_data1.y1 = pos1[0].second.second;
    olp_moller_data1.x1_hycal = hycal_gem1[0].x;
    olp_moller_data1.y1_hycal = hycal_gem1[0].y;
    olp_moller_data1.e1 = ea1[0].first;
    olp_moller_data1.angle1 = ea1[0].second;

    olp_moller_data1.chamber_id2 = pos1[1].first;
    olp_moller_data1.x2 = pos1[1].second.first;
    olp_moller_data1.y2 = pos1[1].second.second;
    olp_moller_data1.x2_hycal = hycal_gem1[1].x;
    olp_moller_data1.y2_hycal = hycal_gem1[1].y;
    olp_moller_data1.e2 = ea1[1].first;
    olp_moller_data1.angle2 = ea1[1].second;
    // gem2
    olp_moller_data2.event_id = moller2->GetEvtID();
    olp_moller_data2.coplanarity = moller2->Coplanarity();

    vector<pair<int, pair<double, double> > > pos2 = moller2->Positions();
    olp_moller_data2.chamber_id1 = pos2[0].first;
    olp_moller_data2.x1 = pos2[0].second.first;
    olp_moller_data2.y1 = pos2[0].second.second;
    olp_moller_data2.x1_hycal = hycal_gem2[0].x;
    olp_moller_data2.y1_hycal = hycal_gem2[0].y;
    olp_moller_data2.e1 = ea2[0].first;
    olp_moller_data2.angle1 = ea2[0].second;

    olp_moller_data2.chamber_id2 = pos2[1].first;
    olp_moller_data2.x2 = pos2[1].second.first;
    olp_moller_data2.y2 = pos2[1].second.second;
    olp_moller_data2.x2_hycal = hycal_gem2[1].x;
    olp_moller_data2.y2_hycal = hycal_gem2[1].y;
    olp_moller_data2.e2 = ea2[1].first;
    olp_moller_data2.angle2 = ea2[1].second;
}

void GEMTree::PushOverlapEpTree(PRadEP *ep1, PRadEP *ep2)
{
    vector<pair<double, double> > ea1 = ep1->EnergyAngle();
    vector<pair<double, double> > ea2 = ep2->EnergyAngle();
    vector<HyCalHit> hycal_gem1 = ep1 -> GetHyCalMatch();
    vector<HyCalHit> hycal_gem2 = ep2 -> GetHyCalMatch();

    if(ea1.size() == 1  && ea2.size() ==1) {
	olp_empty_ep_event = false;
    }
    else {
	return;
    }
    //gem1
    olp_ep_data1.event_id = ep1->GetEvtID();
    olp_ep_data1.chamber_id = ep1->GetChamberID();
    vector<pair<double, double> > pos1 = ep1->Positions();
    olp_ep_data1.x = pos1[0].first;
    olp_ep_data1.y = pos1[0].second;
    olp_ep_data1.x_hycal = hycal_gem1[0].x;
    olp_ep_data1.y_hycal = hycal_gem1[0].y;
    olp_ep_data1.e = ea1[0].first;
    olp_ep_data1.angle = ea1[0].second;
    //gem2
    olp_ep_data2.event_id = ep2->GetEvtID();
    olp_ep_data2.chamber_id = ep2->GetChamberID();
    vector<pair<double, double> > pos2 = ep2->Positions();
    olp_ep_data2.x = pos2[0].first;
    olp_ep_data2.y = pos2[0].second;
    olp_ep_data2.x_hycal = hycal_gem2[0].x;
    olp_ep_data2.y_hycal = hycal_gem2[0].y;
    olp_ep_data2.e = ea2[0].first;
    olp_ep_data2.angle = ea2[0].second;
}

void GEMTree::FillOverlapTree()
{
    if( !olp_empty_moller_event)
	overlap_moller_tree->Fill();
    if( !olp_empty_ep_event)
	overlap_ep_tree->Fill();
    olp_empty_moller_event = true;
    olp_empty_ep_event = true;
}
