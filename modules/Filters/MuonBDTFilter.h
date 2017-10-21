#ifndef FILTER_MUON_BDT_H
#define FILTER_MUON_BDT_H

#include "XmlConfig.h"
#include "XmlRange.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "vendor/loguru.h"

#include <memory>

class MuonBDTFilter {
public:


	shared_ptr<TMVA::Reader> reader = nullptr;
	// BDT variables
	Float_t MVA_dY;
	Float_t MVA_dZ;
	Float_t MVA_dTof;
	Float_t MVA_nSigmaPion;
	Float_t MVA_nHitsFit;
	Float_t MVA_DCA;
	Float_t MVA_Cell;
	Float_t MVA_Module;
	Float_t MVA_BL;
	Float_t MVA_Pt;
	Float_t MVA_Charge;

	XmlRange signal_range;

	MuonBDTFilter() {}
	MuonBDTFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~MuonBDTFilter() {}

	float evaluate( FemtoTrackProxy &_proxy ){
		
		if ( nullptr == _proxy._mtdPid  ) return -1;

		MVA_dY         = _proxy._mtdPid->mDeltaY;
		MVA_dZ         = _proxy._mtdPid->mDeltaZ;
		MVA_dTof       = _proxy._mtdPid->mDeltaTimeOfFlight;
		MVA_nSigmaPion = _proxy._track->nSigmaPion();
		MVA_nHitsFit   = (Float_t)fabs(_proxy._track->mNHitsFit);
		MVA_DCA        = _proxy._track->gDCA();
		MVA_Cell       = (Float_t)_proxy._mtdPid->cell();
		MVA_Module     = (Float_t)_proxy._mtdPid->module();
		MVA_BL         = (Float_t)_proxy._mtdPid->backleg();
		MVA_Pt         = _proxy._track->mPt;
		MVA_Charge     = (Float_t)_proxy._track->charge();

		float bdt = reader->EvaluateMVA( "BDT" );

		// LOG_F( INFO, "MVA_dY = %f ", MVA_dY );
		// LOG_F( INFO, "MVA_dZ = %f ", MVA_dZ );
		// LOG_F( INFO, "MVA_dTof = %f ", MVA_dTof );
		// LOG_F( INFO, "MVA_nSigmaPion = %f ", MVA_nSigmaPion );
		// LOG_F( INFO, "MVA_nHitsFit = %f ", MVA_nHitsFit );
		// LOG_F( INFO, "MVA_DCA = %f ", MVA_DCA );
		// LOG_F( INFO, "MVA_Cell = %f ", MVA_Cell );
		// LOG_F( INFO, "MVA_Module = %f ", MVA_Module );
		// LOG_F( INFO, "MVA_BL = %f ", MVA_BL );
		// LOG_F( INFO, "MVA_Pt = %f ", MVA_Pt );
		// LOG_F( INFO, "MVA_Charge = %f ", MVA_Charge );
		// LOG_F( INFO, "BDT ======= %f", bdt);
		return bdt;

		
	}

	void load( XmlConfig &_cfg, string _nodePath ){
		
		string weights_xml = _cfg.getString( _nodePath + ".weights" );
		LOG_F( INFO, "Loading BDT weights from : %s", weights_xml.c_str() );

		reader = shared_ptr<TMVA::Reader>(new TMVA::Reader( "!Color:!Silent" ) ); 
		
		reader->AddVariable( "dY := (MtdPidTraits_mDeltaY)", &MVA_dY );
		reader->AddVariable( "dZ := MtdPidTraits_mDeltaZ", &MVA_dZ );
		// reader->AddVariable( "nSigmaPi := Tracks_mNSigmaPion", &MVA_nSigmaPion );
		// reader->AddVariable( "nh := Tracks_mNHitsFit", &MVA_nHitsFit );
		// reader->AddVariable( "dca := Tracks_mDCA", &MVA_DCA );
		reader->AddVariable( "Cell := MtdPidTraits_mCell", &MVA_Cell );
		reader->AddVariable( "Module := MtdPidTraits_mModule", &MVA_Module );
		reader->AddVariable( "BL := MtdPidTraits_mBL", &MVA_BL );
		reader->AddVariable( "pT := Tracks_mPt", &MVA_Pt );
		reader->AddVariable( "charge := Tracks_mCharge", &MVA_Charge );
		reader->AddVariable( "dTof := MtdPidTraits_mDeltaTOF", &MVA_dTof );


		reader->BookMVA( "BDT", weights_xml.c_str() ); 

		signal_range.loadConfig( _cfg, _nodePath + ".Range" );

	}

	bool pass( FemtoTrackProxy &_proxy ){

		if ( nullptr == _proxy._mtdPid  ) return false;

		float lh = evaluate( _proxy );

		if ( !signal_range.inInclusiveRange( lh ) )
			return false;
	
		return true;
	}

	bool fail( FemtoTrackProxy &_proxy ){
		return !pass( _proxy );
	}

};

#endif