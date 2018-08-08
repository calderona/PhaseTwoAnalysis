/*
 * ntupler.cpp
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#include "interface/ntupler.h"
#include "interface/scaleFactors.h"
//dirty hack
#include "../../../NTupler/src/MiniEvent.cc"
#include "TDirectory.h"

#include "TH1F.h"

void ntupler::analyze(size_t childid /* this info can be used for printouts */){

        d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
	d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
	d_ana::dBranchHandler<GenParticle> genpart(tree(),"Particle");
	d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>         jet(tree(),"JetPUPPI");
        d_ana::dBranchHandler<Jet>         taujet(tree(),"Jet");
	d_ana::dBranchHandler<Muon>        muonloose(tree(),"MuonLooseCHS");
	d_ana::dBranchHandler<Photon>      photon(tree(),"PhotonLoose");
	d_ana::dBranchHandler<MissingET>   met(tree(),"PuppiMissingET");
	size_t nevents=tree()->entries();
	if(isTestMode())
		nevents/=100;


	//create output
	TString chilidstr="";
	chilidstr+=childid;
	TFile * outfile= new TFile(getOutDir()+"/p2ntuple_"+(TString)getLegendName()+"_"+chilidstr+".root","RECREATE");
	TDirectory *counterdir = outfile->mkdir("weightCounter");
	counterdir->cd();
	TH1F * h_event_weight = new TH1F("Event_weight","Event_weight",1,0,1);

	outfile->cd();
	TDirectory *ntupledir = outfile->mkdir("ntuple");
	ntupledir->cd();

	MiniEvent_t ev_;
	TTree * t_event_        = new TTree("Event","Event");
	TTree * t_genParts_     = new TTree("Particle","Particle");
	TTree * t_genPhotons_   = new TTree("GenPhoton","GenPhoton");
	TTree * t_vertices_     = new TTree("Vertex","Vertex");
	TTree * t_genJets_      = new TTree("GenJet","GenJet");
	TTree * t_looseElecs_   = new TTree("ElectronLoose","ElectronLoose");
	TTree * t_mediumElecs_  = new TTree("ElectronMedium","ElectronMedium");
	TTree * t_tightElecs_   = new TTree("ElectronTight","ElectronTight");
	TTree * t_looseMuons_   = new TTree("MuonLoose","MuonLoose");
	TTree * t_tightMuons_   = new TTree("MuonTight","MuonTight");
	TTree * t_allTaus_   = new TTree("TauAll","TauAll");
	TTree * t_puppiJets_    = new TTree("JetPUPPI","JetPUPPI");
	TTree * t_puppiMET_     = new TTree("PuppiMissingET","PuppiMissingET");
	TTree * t_loosePhotons_ = new TTree("PhotonLoose","PhotonLoose");
	TTree * t_tightPhotons_ = new TTree("PhotonTight","PhotonTight");

	//createMiniEventTree(t_event_, t_genParts_, t_vertices_, t_genJets_, t_looseMuons_, t_puppiJets_, t_puppiMET_, ev_);

	createMiniEventTree(t_event_, t_genParts_, t_vertices_, t_genJets_, t_genPhotons_, t_looseElecs_, t_mediumElecs_,t_tightElecs_, t_looseMuons_, t_tightMuons_, t_allTaus_, t_puppiJets_, t_puppiMET_, t_loosePhotons_,t_tightPhotons_, ev_);



	for(size_t eventno=0;eventno<nevents;eventno++){
		/*
		 * The following two lines report the status and set the event link
		 * Do not remove!
		 */
		reportStatus(eventno,nevents);
		tree()->setEntry(eventno);

		if(event.size()<1)continue;

		h_event_weight->Fill(0.,(double)event.at(0)->Weight);


		std::vector<Muon*>selectedMuons;
		for(size_t i=0;i<muonloose.size();i++){
                  if(muonloose.at(i)->PT<5)continue;
                  selectedMuons.push_back(muonloose.at(i));
		}


		std::vector<Jet*>selectedjets;
		for(size_t i=0;i<jet.size();i++){
			if(jet.at(i)->PT<10)continue;
			selectedjets.push_back(jet.at(i));
		}


            
		ev_.event = event.at(0)->Number;
		ev_.g_nw = 1;
		ev_.g_w[0] = event.at(0)->Weight;


            ev_.ngl=0;

            for(size_t i=0;i<genpart.size();i++){
                  if(ev_.ngl>=MiniEvent_t::maxpart)break;

                  int pid= fabs(genpart.at(i)->PID);

                  if( pid != 13  || pid != -13  ) continue;

                  ev_.gl_pid[ev_.ngl]=genpart.at(i)->PID;
                  ev_.gl_ch[ev_.ngl]=genpart.at(i)->Charge;
                  ev_.gl_st[ev_.ngl]=genpart.at(i)->Status;
                  ev_.gl_p[ev_.ngl]=genpart.at(i)->P;
                  ev_.gl_pz[ev_.ngl]=genpart.at(i)->Pz;
                  ev_.gl_pt[ev_.ngl]=genpart.at(i)->PT;
                  ev_.gl_eta[ev_.ngl]=genpart.at(i)->Eta;
                  ev_.gl_phi[ev_.ngl]=genpart.at(i)->Phi;
                  ev_.gl_mass[ev_.ngl]=genpart.at(i)->Mass;
                  ev_.ngl++;

                  //std::cout<<i<<"   "<<pid<<"    "<<genpart.at(i)->Status<<"  ->"<<genpart.at(i)->M1<<"  "<<genpart.at(i)->M2<<" ;  "<<genpart.at(i)->D1<<"  "<<genpart.at(i)->D2<<std::endl;

            }

            std::cout<<std::endl;

		ev_.ntm=0;
		for(size_t i=0;i<selectedMuons.size();i++){
			if(ev_.ntm>=MiniEvent_t::maxpart)break;
			ev_.tm_pt    [ev_.ntm] =selectedMuons.at(i)->PT;
			ev_.tm_eta   [ev_.ntm]=selectedMuons.at(i)->Eta;
			ev_.tm_phi   [ev_.ntm]=selectedMuons.at(i)->Phi;
			ev_.tm_mass  [ev_.ntm]=0.105;
			ev_.tm_relIso[ev_.ntm]=selectedMuons.at(i)->IsolationVarRhoCorr; // /selectedMuons.at(i)->PT;
			//ev_.tm_g     [ev_.ntm] =selectedMuons.at(i)->Particle.PID;
			ev_.ntm++;
		}

		ev_.nj=0;
		for(size_t i=0;i<selectedjets.size();i++){
			if(ev_.nj>=MiniEvent_t::maxjets)break;
			ev_.j_pt  [ev_.nj] =selectedjets.at(i)->PT;
			ev_.j_eta [ev_.nj]=selectedjets.at(i)->Eta;
			ev_.j_phi [ev_.nj]=selectedjets.at(i)->Phi;
			ev_.j_mass[ev_.nj]=selectedjets.at(i)->Mass;

			ev_.j_hadflav[ev_.nj]=selectedjets.at(i)->Flavor;

			ev_.j_deepcsv[ev_.nj]=0;
			ev_.j_mvav2[ev_.nj]=0;
			if(selectedjets.at(i)->BTag){
				ev_.j_deepcsv[ev_.nj]=0b00000111;
				ev_.j_mvav2[ev_.nj]=0b00000111;
			}
			ev_.nj++;
		}

		ev_.nmet=0;
		for(size_t i=0;i<met.size();i++){
			if(ev_.nmet>=MiniEvent_t::maxpart) break;
			ev_.met_eta[ev_.nmet]=met.at(i)->Eta ;
			ev_.met_pt [ev_.nmet]=met.at(i)->MET ;
			ev_.met_phi[ev_.nmet]=met.at(i)->Phi ;
			ev_.nmet++;
		}



		t_event_->Fill();
		t_genParts_->Fill();
		//t_genPhotons_->Fill();
		t_vertices_->Fill();
		t_genJets_->Fill();
		//t_looseElecs_->Fill();
		//t_mediumElecs_->Fill();
		//t_tightElecs_->Fill();
		t_looseMuons_->Fill();
		//t_tightMuons_->Fill();
		//t_allTaus_->Fill();
		t_puppiJets_->Fill();
		t_puppiMET_->Fill();
		//t_loosePhotons_->Fill();
		//t_tightPhotons_->Fill();

	}

	counterdir->cd();
	h_event_weight->Write();

	ntupledir->cd();
	t_event_        ->Write();
	t_genParts_     ->Write();
	//t_genPhotons_   ->Write();
	t_vertices_     ->Write();
	t_genJets_      ->Write();
	//t_looseElecs_   ->Write();
	//t_mediumElecs_   ->Write();
	//t_tightElecs_   ->Write();
	t_looseMuons_   ->Write();
	//t_tightMuons_   ->Write();
	//t_allTaus_      ->Write();
	t_puppiJets_    ->Write();
	t_puppiMET_     ->Write();
	//t_loosePhotons_ ->Write();
	//t_tightPhotons_ ->Write();

	outfile->Close();
	/*
	 * Must be called in the end, takes care of thread-safe writeout and
	 * call-back to the parent process
	 */
	processEndFunction();
}



void ntupler::postProcess(){

	/* empty */

}



