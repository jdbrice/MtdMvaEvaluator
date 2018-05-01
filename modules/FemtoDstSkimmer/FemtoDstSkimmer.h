#ifndef FEMTO_DST_SKIMMER_H
#define FEMTO_DST_SKIMMER_H

#include "TreeAnalyzer.h"
#include "HistoBins.h"

// FemtoDstFormat
#include "FemtoDstFormat/BranchReader.h"
#include "FemtoDstFormat/TClonesArrayReader.h"
#include "FemtoDstFormat/FemtoEvent.h"
#include "FemtoDstFormat/FemtoTrack.h"
#include "FemtoDstFormat/FemtoMcTrack.h"
#include "FemtoDstFormat/FemtoTrackHelix.h"
#include "FemtoDstFormat/FemtoMtdPidTraits.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

// #include "Filters/MuonBDTFilter.h"
// #include "Filters/MuonMLPFilter.h"
#include "Filters/MuonMVAFilter.h"




#define LOGURU_WITH_STREAMS 1
#include "vendor/loguru.h"

#include "XmlHistogram.h"
#include "XmlFunction.h"

#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"

class FemtoDstSkimmer : public TreeAnalyzer
{
protected:
	FemtoEvent *_event;

	BranchReader<FemtoEvent> _rEvent;
	TClonesArrayReader<FemtoTrack> _rTracks;
	TClonesArrayReader<FemtoMcTrack> _rMcTracks;
	TClonesArrayReader<FemtoTrackHelix> _rHelices;
	TClonesArrayReader<FemtoMtdPidTraits> _rMtdPid;

	int ptIndex = 0;
	vector <MuonMVAFilter> mlps;
	vector <MuonMVAFilter> bdts;

	XmlHistogram bg_deltaTOF;
	XmlHistogram sig_deltaTOF_2d;
	vector<TH1*> sig_deltaTOF_pt;
	TH1 * sig_deltaTOF_pt_all;

	XmlHistogram kaon_deltaTOF_2d;
	vector<TH1*> kaon_deltaTOF_pt;
	TH1 * kaon_deltaTOF_pt_all;

	XmlHistogram proton_deltaTOF_2d;
	vector<TH1*> proton_deltaTOF_pt;
	TH1 * proton_deltaTOF_pt_all;

	XmlHistogram bg_dca;
	XmlHistogram sig_dca;

	HistoBins bins_pt_mva;
	HistoBins bins_pt_dtof;
	HistoBins bins_kaon_pt_dtof;

	TNtuple *tuple;
	TRandom3 r3;

	bool use_bdt = false;
public:
	virtual const char* classname() const {return "FemtoDstSkimmer";}
	FemtoDstSkimmer() {}
	~FemtoDstSkimmer() {}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		_rEvent.setup( chain, "Event" );
		_rTracks.setup( chain, "Tracks" );
		_rMcTracks.setup( chain, "McTracks" );
		_rMtdPid.setup( chain, "MtdPidTraits" );

		
		// bdt.load( config, nodePath + ".MuonBDTFilter" );

		string path_bg = config.q( nodePath+".DeltaTOF.XmlHistogram{name==bg}" );
		string path_kaon = config.q( nodePath+".DeltaTOF.XmlHistogram{name==kaon_dtof_vs_pt}" );
		string path_proton = config.q( nodePath+".DeltaTOF.XmlHistogram{name==proton_dtof_vs_pt}" );
		string path_sig = config.q( nodePath+".DeltaTOF.XmlHistogram{name==hsignalPdf}" );


		
		LOG_F( INFO, "Loading bg from : %s", path_bg.c_str()  );
		LOG_F( INFO, "Loading kaon from : %s", path_kaon.c_str()  );
		LOG_F( INFO, "Loading proton from : %s", path_proton.c_str()  );
		LOG_F( INFO, "Loading sig from : %s", path_sig.c_str()  );

		bg_deltaTOF.load( config, path_bg );
		sig_deltaTOF_2d.load( config, path_sig );
		kaon_deltaTOF_2d.load( config, path_kaon );
		proton_deltaTOF_2d.load( config, path_proton );

		bg_dca.load( config, nodePath + ".DCA.XmlHistogram[0]" );
		sig_dca.load( config, nodePath + ".DCA.XmlHistogram[1]" );

		book->cd();
		// bg_dca.getTH1()->Write();
		// sig_dca.getTH1()->Write();

		tuple= new TNtuple( "MvaTree", "mva variables", "classId:dY:dZ:dTof:nsp:nhf:dca:cell:mod:bl:pt:c:bdt:mlp:gid" );

		int seed= config.getInt( "seed", 0 );
		r3.SetSeed( seed );
		gRandom->SetSeed(r3.Integer( 32000 ));
		LOG_F( INFO, "RANDOM SEED : %d", seed );

		bins_pt_mva.load( config, "bins.mva_pt" );
		string mlp_template_str = config.getString( nodePath + ".MuonMVAFilter.WeightsMLP" );
		string bdt_template_str = config.getString( nodePath + ".MuonMVAFilter.WeightsBDT", "" );

		/* LOAD THE MLPs */
		// use this first one to setup the sinlgeton reader and load the variables
		MuonMVAFilter m_setup;
		m_setup.loadVars( config, nodePath + ".MuonMVAFilter" );

		if ( mlp_template_str.find( "%" ) != string::npos ){
			for ( size_t i = 0; i < bins_pt_mva.nBins(); i++ ){
				MuonMVAFilter m;
				m.load( string( TString::Format( mlp_template_str.c_str(), i ) ), string( TString::Format( "mlp_%zu", i ) ) );
				mlps.push_back( m );
			}
		} else {
			MuonMVAFilter m;
			m.load( mlp_template_str, "DNN" );
			mlps.push_back( m );
		}

		/* LOAD THE BDTs */
		if ( bdt_template_str.find( "%" ) != string::npos ){
			if ( "" != bdt_template_str ){
				use_bdt = true;
				for ( size_t i = 0; i < bins_pt_mva.nBins(); i++ ){
					MuonMVAFilter m;
					m.load( string( TString::Format( bdt_template_str.c_str(), i ) ), string( TString::Format( "bdt_%zu", i ) ) );
					bdts.push_back( m );
				}
			}
		} else {
			if ( "" != bdt_template_str ){
				use_bdt = true;
				MuonMVAFilter m;
				m.load( bdt_template_str, "BDT" );
				bdts.push_back( m );
			}
		}


		// build the signal dtof 1d histos
		TH2 * h2 = (TH2*)sig_deltaTOF_2d.getTH1().get();
		bins_pt_dtof.load( config, "bins.dtof_pt" );
		size_t nPtBinsDtof = bins_pt_dtof.nBins();
		for ( size_t i = 0; i < nPtBinsDtof; i++ ){
			int ib1 = h2->GetXaxis()->FindBin( bins_pt_dtof.bins[i] );
			int ib2 = h2->GetXaxis()->FindBin( bins_pt_dtof.bins[i+1] );
			LOG_F( 3, "Projecting ( %f -> %f )", bins_pt_dtof.bins[i], bins_pt_dtof.bins[i+1] );
			TH1 * h1 = h2->ProjectionY( TString::Format( "sig_dtof_pt_%zu", i ), ib1, ib2 );
			sig_deltaTOF_pt.push_back( h1 );
			h1->SetTitle( TString::Format( "%0.2f < pT < %0.2f (GeV/c)", bins_pt_dtof.bins[i], bins_pt_dtof.bins[i+1] ) );
			h1->Write();
		}
		sig_deltaTOF_pt_all = h2->ProjectionY( "sig_dtof_pt_all", 1, -1 );


		// build the kaon dtof 1d histos
		TH2 * kh2 = (TH2*)kaon_deltaTOF_2d.getTH1().get();
		LOG_F( INFO, "kh2=%p", kh2 );
		bins_kaon_pt_dtof.load( config, "bins.kaon_dtof_pt" );
		size_t nKaonPtBinsDtof = bins_kaon_pt_dtof.nBins();
		for ( size_t i = 0; i < nKaonPtBinsDtof; i++ ){
			int ib1 = kh2->GetXaxis()->FindBin( bins_kaon_pt_dtof.bins[i] );
			int ib2 = kh2->GetXaxis()->FindBin( bins_kaon_pt_dtof.bins[i+1] );
			LOG_F( 3, "Projecting ( %f -> %f )", bins_kaon_pt_dtof.bins[i], bins_kaon_pt_dtof.bins[i+1] );
			TH1 * h1 = kh2->ProjectionY( TString::Format( "kaon_dtof_pt_%zu", i ), ib1, ib2 );
			kaon_deltaTOF_pt.push_back( h1 );
			h1->SetTitle( TString::Format( "%0.2f < pT < %0.2f (GeV/c)", bins_kaon_pt_dtof.bins[i], bins_kaon_pt_dtof.bins[i+1] ) );
			h1->Write();
		}
		kaon_deltaTOF_pt_all = kh2->ProjectionY( "kaon_dtof_pt_all", 1, -1 );

		// build the proton dtof 1d histos
		TH2 * ph2 = (TH2*)proton_deltaTOF_2d.getTH1().get();
		LOG_F( INFO, "ph2=%p", ph2 );
		for ( size_t i = 0; i < nKaonPtBinsDtof; i++ ){
			int ib1 = ph2->GetXaxis()->FindBin( bins_kaon_pt_dtof.bins[i] );
			int ib2 = ph2->GetXaxis()->FindBin( bins_kaon_pt_dtof.bins[i+1] );
			LOG_F( 3, "Projecting ( %f -> %f )", bins_kaon_pt_dtof.bins[i], bins_kaon_pt_dtof.bins[i+1] );
			TH1 * h1 = ph2->ProjectionY( TString::Format( "proton_dtof_pt_%zu", i ), ib1, ib2 );
			proton_deltaTOF_pt.push_back( h1 );
			h1->SetTitle( TString::Format( "%0.2f < pT < %0.2f (GeV/c)", bins_kaon_pt_dtof.bins[i], bins_kaon_pt_dtof.bins[i+1] ) );
			h1->Write();
		}
		proton_deltaTOF_pt_all = kh2->ProjectionY( "kaon_dtof_pt_all", 1, -1 );

		// LOG_F( INFO, "[%d] => pT Range : %f <= pT < %f", ptIndex, bins_pt_mva.bins[ptIndex], bins_pt_mva.bins[ptIndex+1] );

	}


protected:


	float sample_bg_dca( FemtoTrackProxy &proxy ){
		float pt = proxy._track->mPt;
		if ( pt > 14 ) pt = 14;
		TH2 * hbgdca = ((TH2*)bg_dca.getTH1().get());
		int ptbin = hbgdca->GetXaxis()->FindBin( pt );
		TH1 * hbgdcaslice = hbgdca->ProjectionY( "tmp", ptbin, ptbin );
		return hbgdcaslice->GetRandom();
	}

	float sample_sig_dca( FemtoTrackProxy &proxy ){
		float pt = proxy._track->mPt;
		if ( pt > 14 ) pt = 14;
		TH2 * hsigdca = ((TH2*)sig_dca.getTH1().get());
		int ptbin = hsigdca->GetXaxis()->FindBin( pt );
		TH1 * hsigdcaslice = hsigdca->ProjectionY( "tmp", ptbin, ptbin );
		return hsigdcaslice->GetRandom();
	}

	float sample_bg_deltaTOF( ){
		return bg_deltaTOF.getTH1().get()->GetRandom();
	} 
	float sample_sig_deltaTOF( FemtoTrackProxy &proxy ){
		LOG_F( 3, "pT=%f ", proxy._track->mPt );

		int ib = bins_pt_dtof.findBin( proxy._track->mPt );
		LOG_F( 3, "ib=%d", ib );
		if ( ib < sig_deltaTOF_pt.size() && ib >= 0){
			LOG_F( 3, "h1: %s", sig_deltaTOF_pt[ib]->GetName() );
			return sig_deltaTOF_pt[ib]->GetRandom();
		}

		
		return sig_deltaTOF_pt_all->GetRandom();
	} 

	float sample_kaon_deltaTOF( FemtoTrackProxy &proxy ){
		LOG_F( 3, "pT=%f ", proxy._track->mPt );

		int ib = bins_kaon_pt_dtof.findBin( proxy._track->mPt );
		LOG_F( 3, "ib=%d", ib );
		if ( ib < kaon_deltaTOF_pt.size() && ib >= 0){
			LOG_F( 3, "h1: %s", kaon_deltaTOF_pt[ib]->GetName() );
			return kaon_deltaTOF_pt[ib]->GetRandom();
		}
		return kaon_deltaTOF_pt_all->GetRandom();
	} 

	float sample_proton_deltaTOF( FemtoTrackProxy &proxy ){
		LOG_F( 3, "pT=%f ", proxy._track->mPt );

		int ib = bins_kaon_pt_dtof.findBin( proxy._track->mPt );
		LOG_F( 3, "ib=%d", ib );
		if ( ib < proton_deltaTOF_pt.size() && ib >= 0){
			LOG_F( 3, "h1: %s", proton_deltaTOF_pt[ib]->GetName() );
			return proton_deltaTOF_pt[ib]->GetRandom();
		}
		return proton_deltaTOF_pt_all->GetRandom();
	} 

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();
		book->cd();
	}

	virtual void analyzeTrack( FemtoTrackProxy &_proxy ){
		
		int ipt = bins_pt_mva.findBin( _proxy._track->mPt );
		// LOG_F( INFO, "pt=%f, ipt = %d", _proxy._track->mPt, ipt );
		if ( ipt < 0 ) return;
		// if ( ipt != ptIndex ) return;

		MuonMVAFilter &_mlp = mlps[ipt];
		MuonMVAFilter *_bdt = &mlps[ipt];

		if ( use_bdt )
			_bdt = &bdts[ipt];

		float classId = 0; // bg
		
		if ( isDecayMuonOutsideTPC( _proxy._mtdPid, _proxy._mcTrack ) ){
			classId = 2;
		}
		if ( isDecayMuonInsideTPC( _proxy._mcTrack ) ){
			return;
			classId = 1;
		}
		
		// if ( cleanPunchThrough( _proxy._mtdPid, track, mcTrack ) ){
		// 	classId = 3;
		// }
		
		if ( isSignal( _proxy._mcTrack ) ){

			classId = 5;
		}


		/************************************
		* Sample dTOF distribution for sig/bg
		************************************/
		double deltaTof = -999;
		if ( isSignal( _proxy._mcTrack  )){
			float sdtof = sample_sig_deltaTOF( _proxy );
			// book->fill( "signal_dtof", sdtof );
			deltaTof = sdtof;
		} else if ( isKaon( _proxy._mcTrack ) ){
			float kdtof = sample_kaon_deltaTOF( _proxy );
			deltaTof = kdtof;
			classId = 10;
		} else if ( isProton( _proxy._mcTrack ) ){
			float pdtof = sample_proton_deltaTOF( _proxy );
			deltaTof = pdtof;
			classId = 11;
		} else {
			deltaTof = sample_bg_deltaTOF( );
		}


		/************************************
		* Sample dca distribution for sig/bg
		************************************/
		if ( isSignal( _proxy._mcTrack ) ){
			_proxy._track->gDCA( sample_sig_dca( _proxy ) );
		}
		 else {
			_proxy._track->gDCA( sample_bg_dca( _proxy ) );
		}


		_proxy._mtdPid->mDeltaTimeOfFlight = deltaTof;

		MuonMVAFilter::fillVars( _proxy );
		float mlpr = _mlp.evaluate( _proxy );
		float bdtr = 0;
		if ( use_bdt  )
			bdtr = _bdt->evaluate( _proxy );



		float data[] = {
			classId,
			_proxy._mtdPid->mDeltaY,
			_proxy._mtdPid->mDeltaZ,
			_proxy._mtdPid->mDeltaTimeOfFlight,
			_proxy._track->nSigmaPion(),
			(Float_t)fabs(_proxy._track->mNHitsFit),
			_proxy._track->gDCA(),
			(Float_t)_proxy._mtdPid->cell(),
			(Float_t)_proxy._mtdPid->module(),
			(Float_t)_proxy._mtdPid->backleg(),
			_proxy._track->mPt,
			(Float_t)_proxy._track->charge(),
			bdtr,
			mlpr,
			(float)_proxy._mcTrack->mGeantPID
		};
		tuple->Fill( data );
	}

	virtual void analyzeEvent(){
		// return;
		_event = _rEvent.get();

		size_t nTracks = _rTracks.N();
		FemtoTrackProxy _proxy;

		for (size_t i = 0; i < nTracks; i++ ){
			
			FemtoTrack        *track   = _rTracks.get(i);
			FemtoMcTrack      *mcTrack = nullptr;
			FemtoMtdPidTraits *mtdPid  = nullptr;

			// get MTD PidTraits
			if ( track->mMtdPidTraitsIndex >= 0) 
				mtdPid = _rMtdPid.get( track->mMtdPidTraitsIndex );

			// get MC Track 
			if ( track->mMcIndex >= 0 )
				mcTrack = _rMcTracks.get( track->mMcIndex );

			// Reject if one is missing
			if ( nullptr == mtdPid || nullptr == mcTrack )
				continue;

			FemtoTrackProxy _proxy;
			_proxy._track = track;
			_proxy._mtdPid = mtdPid;
			_proxy._mcTrack = mcTrack;

			analyzeTrack( _proxy );

		}

		

	}


	/* Post Event Loop
	 * Adds config and writes tuple
	 */
	virtual void postEventLoop(){
		tuple->Write();
		TreeAnalyzer::postEventLoop();

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}
	} // postEventLoop

	bool isMuon( FemtoMcTrack * mcTrack ){
		if ( nullptr == mcTrack )
			return false;
		if ( 5 == mcTrack->mGeantPID || 6 == mcTrack->mGeantPID )
			return true;
		return false;
	}

	bool isKaon( FemtoMcTrack *mcTrack ){
		if ( nullptr == mcTrack )
			return false;
		if ( 11 == mcTrack->mGeantPID || 12 == mcTrack->mGeantPID )
			return true;
		return false;
	}

	bool isProton( FemtoMcTrack *mcTrack ){
		if ( nullptr == mcTrack )
			return false;
		if ( 14 == mcTrack->mGeantPID || 15 == mcTrack->mGeantPID )
			return true;
		return false;
	}

	bool isSignal( FemtoMcTrack * mcTrack ){
		if ( isMuon( mcTrack ) && mcTrack->mParentIndex < 0 )
			return true;

		return false;
	}	

	bool isDecayMuon( FemtoMcTrack * mcTrack ){
		if (nullptr == mcTrack) 
			return false;
		if ( 5 == mcTrack->mGeantPID || 6 == mcTrack->mGeantPID ){
			if ( mcTrack->mParentIndex < 0 ) // signal muon (primary)
				return false;
			
			auto parent = _rMcTracks.get( mcTrack->mParentIndex );
			
			if (8 == parent->mGeantPID || 9 == parent->mGeantPID || 11 == parent->mGeantPID || 12 == parent->mGeantPID )
				return true;
			
			return false;
		}

		return false;
	}

	bool isDecayMuonInsideTPC( FemtoMcTrack *mcTrack ){
		return isDecayMuon( mcTrack );
	}

	bool isDecayMuonOutsideTPC( FemtoMtdPidTraits *mtdPid, FemtoMcTrack *mcTrack ){
		if (nullptr == mtdPid) 
			return false;
		if ( mtdPid->mIdTruth < 0 )
			return false;
		auto mtdMcTrack = _rMcTracks.get( mtdPid->mIdTruth );

		if ( mcTrack->mParentIndex >= 0 )
			return false;

		return isDecayMuon( mtdMcTrack );
	}
	bool cleanPunchThrough( FemtoMtdPidTraits *mtdPid, FemtoTrack *track, FemtoMcTrack *mcTrack ){
		
		if ( nullptr == track || nullptr == mtdPid )
			return false;
		
		if ( mcTrack->mGeantPID != 8 && mcTrack->mGeantPID != 9 )
			return false;

		if ( mcTrack->mParentIndex >=0 )
			return false;

		if ( track->mMcIndex == mtdPid->mIdTruth )
			return true;

		return false;
	}
	
};

#endif
