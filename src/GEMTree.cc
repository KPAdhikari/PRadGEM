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

using namespace std;

GEMTree::GEMTree()
{
    TH1::AddDirectory(kFALSE);
    file = new TFile("root_file/res.root", "recreate");
    event_id = 0;
    empty_moller_event = true;
    empty_ep_event = true;
    empty_epics_event = true;
    empty_gem_event = true;

    // init res tree
    this -> InitGEMTree(2); // 2 gem detectors
    this -> InitPhysicsTree();
    this -> InitEpicsTree();
    this -> InitCaliOffsetTree();
    this -> InitProdOffsetTree();
}

GEMTree::~GEMTree()
{
}

void GEMTree::SetEventID(unsigned int &id)
{
    event_id = (int) id;
}

void GEMTree::InitGEMTree(int ndet)
{
    gem_tree = new TTree("T", "res");
    gem_tree->SetDirectory(file);

    assert( ndet < NDETECTOR);
    nDet = ndet;
    for(int i=0;i<2;i++)
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
}

void GEMTree::InitPhysicsTree()
{
    //physics_file = new TFile("root_file/physics.root", "recreate");
    moller_tree = new TTree("moller_tree", "moller tree");
    moller_tree -> SetDirectory(file);
    ep_tree = new TTree("ep_tree", "ep tree");
    ep_tree -> SetDirectory(file);
    ep_moller_tree = new TTree("ep_moller_tree", "ep moller tree");
    ep_moller_tree->SetDirectory(file);

    moller_tree->Branch("nCluster", &nCluster, "nCluster/I");
    moller_tree->Branch("moller_scatt_angle1", &moller_scatt_angle1, "moller_scatt_angle1/F");
    moller_tree->Branch("moller_scatt_angle2", &moller_scatt_angle2, "moller_scatt_angle2/F");
    moller_tree->Branch("moller_scatt_energy1", &moller_scatt_energy1, "moller_scatt_energy1/F");
    moller_tree->Branch("moller_scatt_energy2", &moller_scatt_energy2, "moller_scatt_energy2/F");
    moller_tree->Branch("scatt_x", scatt_x, "scatt_x[nCluster]/F");
    moller_tree->Branch("scatt_y", scatt_y, "scatt_y[nCluster]/F");
    moller_tree->Branch("coplanarity", &coplanarity, "coplanarity/F");
    moller_tree->Branch("open_angle", &open_angle, "open_angle/F");
    moller_tree->Branch("scatt_energy", scatt_energy, "scatt_energy[nCluster]/F");
    moller_tree->Branch("scatt_angle", scatt_angle, "scatt_angle[nCluster]/F");
    moller_tree->Branch("moller_center_x", &moller_center_x, "moller_center_x/F");
    moller_tree->Branch("moller_center_y", &moller_center_y, "moller_center_y/F");
    moller_tree->Branch("moller_pos_res_dx", &moller_pos_res_dx, "moller_pos_res_dx/F");
    moller_tree->Branch("moller_pos_res_dy", &moller_pos_res_dy, "moller_pos_res_dy/F");

    ep_tree -> Branch("cCluster", &nCluster, "nCluster/I");
    ep_tree -> Branch("scatt_energy", scatt_energy, "scatt_energy[nCluster]/F");
    ep_tree -> Branch("scatt_angle", scatt_angle, "scatt_angle[nCluster]/F");
    ep_tree->Branch("scatt_x", scatt_x, "scatt_x[nCluster]/F");
    ep_tree->Branch("scatt_y", scatt_y, "scatt_y[nCluster]/F");
    ep_tree -> Branch("q_square", &q_square, "q_square/F");

    ep_moller_tree -> Branch("cCluster", &nCluster, "nCluster/I");
    ep_moller_tree -> Branch("scatt_energy", scatt_energy, "scatt_energy[nCluster]/F");
    ep_moller_tree -> Branch("scatt_angle", scatt_angle, "scatt_angle[nCluster]/F");
}

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

void GEMTree::PushMoller(PRadMoller * moller)
{
    nCluster = 2;
    vector<pair<double, double> > ea = moller->EnergyAngle();
    if(ea.size() == 2) {
	empty_moller_event = false;
    }
    else {
	return;
    }
    vector<pair<double, double> > pos = moller->Positions();
    open_angle = moller->OpenAngle();
    coplanarity = moller->Coplanarity();
    int ii = 0;
    for(auto &i: ea)
    {
	scatt_angle[ii] = i.second;
	scatt_energy[ii] = i.first;
	//temp
	scatt_x[ii] = pos[ii].first;
	scatt_y[ii] = pos[ii].second;
	ii++;
    }
    moller_scatt_angle1 = scatt_angle[0];
    moller_scatt_angle2 = scatt_angle[1];
    moller_scatt_energy1 = scatt_energy[0];
    moller_scatt_energy2 = scatt_energy[1];
    moller_center_x = moller->MollerCenter().first;
    moller_center_y = moller->MollerCenter().second;
    pair<double, double> reso = moller->GetSpatialResHandler()->GetDiff();
    moller_pos_res_dx = reso.first;
    moller_pos_res_dy = reso.second;
}

void GEMTree::PushEP(PRadEP * ep)
{
    nCluster = 1;
    vector<pair<double, double> > ea = ep->EnergyAngle();
    if(ea.size() == 1 ) {
	empty_ep_event = false;
    }
    else {
	return;
    }
    vector<pair<double, double> > pos = ep->Positions();
    int ii = 0;
    for(auto &i: ea)
    {
	scatt_angle[ii] = i.second;
	scatt_energy[ii] = i.first;
	//temp
	scatt_x[ii] = pos[ii].first;
	scatt_y[ii] = pos[ii].second;
	ii++;
    }
    q_square = ep -> QSquare();
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

void GEMTree::FillGEMTree()
{
    if( ! empty_gem_event){
	gem_tree->Fill();
    }
    empty_gem_event = true;
}

void GEMTree::WriteToDisk()
{
    file->Write();
    //file->Save();
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

    prod_offset_tree_res -> Branch("prod_offset_x", &prod_offset_x, "prod_offset_x/D");
    prod_offset_tree_res -> Branch("prod_offset_y", &prod_offset_y, "prod_offset_y/D");
}

void GEMTree::PushProdOffset(int i, PRadMoller * moller)
{
    vector<pair<double, double> > ea = moller->EnergyAngle();
    if( ea.size() == 2 ){
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
	    vector<pair<double, double> > pos = moller->Positions();
	    prod_gem1_scattx[0] = pos[0].first;
	    prod_gem1_scatty[0] = pos[0].second;
	    prod_gem1_scattx[1] = pos[1].first;
	    prod_gem1_scatty[1] = pos[1].second;
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
	    vector<pair<double, double> > pos = moller->Positions();
	    prod_gem2_scattx[0] = pos[0].first;
	    prod_gem2_scatty[0] = pos[0].second;
	    prod_gem2_scattx[1] = pos[1].first;
	    prod_gem2_scatty[1] = pos[1].second;
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
