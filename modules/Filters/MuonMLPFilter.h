#ifndef FILTER_MUON_MLP_H
#define FILTER_MUON_MLP_H

#include "XmlConfig.h"
#include "XmlRange.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "vendor/loguru.h"

#include <memory>

class MuonMLPFilter {
public:

	string name = "MLP";
	static shared_ptr<TMVA::Reader> reader;
	// MLP variables
	Float_t MVA_dY;
	Float_t MVA_dZ;
	Float_t MVA_dTof;
	
	Float_t MVA_nSigmadY;
	Float_t MVA_nSigmadZ;
	Float_t MVA_nSigmadTof;

	Float_t MVA_nSigmaPion;
	Float_t MVA_nHitsFit;
	Float_t MVA_DCA;
	Float_t MVA_Cell;
	Float_t MVA_Module;
	Float_t MVA_BL;
	Float_t MVA_Pt;
	Float_t MVA_Charge;

	XmlRange signal_range;

	vector<string> vars;

	MuonMLPFilter() {}
	MuonMLPFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~MuonMLPFilter() {}

	double nSigmaDeltaY( float pt, float dy ){
		double sigy = -17.6867 + 18.4528*exp(0.637142/pt);
		return dy / sigy;
	}

	double nSigmaDeltaZ( float pt, float dz ){
		double sigz = -32.6793 + 32.6034 * exp( 0.444217 / pt );
		return dz / sigz;
	}

	double nSigmaDeltaTOF( float pt, float dtof ){
		double sigdtof = 0.100915 + 0.00911555 * exp( 3.47368 / pt );
		return dtof / sigdtof;
	}

	float evaluate( FemtoTrackProxy &_proxy ){
		
		if ( nullptr == _proxy._mtdPid  ) return -1;

		MVA_dY         = _proxy._mtdPid->mDeltaY;
		MVA_dZ         = _proxy._mtdPid->mDeltaZ;
		MVA_dTof       = _proxy._mtdPid->mDeltaTimeOfFlight;

		MVA_nSigmadY   = nSigmaDeltaY( _proxy._track->mPt, _proxy._mtdPid->mDeltaY );
		MVA_nSigmadZ   = nSigmaDeltaZ( _proxy._track->mPt, _proxy._mtdPid->mDeltaZ );
		MVA_nSigmadTof = nSigmaDeltaTOF( _proxy._track->mPt, _proxy._mtdPid->mDeltaTimeOfFlight );

		MVA_nSigmaPion = _proxy._track->nSigmaPion();
		MVA_nHitsFit   = (Float_t)fabs(_proxy._track->mNHitsFit);
		MVA_DCA        = _proxy._track->gDCA();
		MVA_Cell       = (Float_t)_proxy._mtdPid->cell();
		MVA_Module     = (Float_t)_proxy._mtdPid->module();
		MVA_BL         = (Float_t)_proxy._mtdPid->backleg();
		MVA_Pt         = _proxy._track->mPt;
		MVA_Charge     = (Float_t)_proxy._track->charge();

		MVA_dY        *= MVA_Charge;

		return reader->EvaluateMVA( this->name.c_str() );
	}

	bool hasVar( string var ){
		if ( std::find( this->vars.begin(), this->vars.end(), var ) != this->vars.end() )
			return true;
		return false;
	}

	void load( string _weight_file, string _name = "MLP" ){

		
		this->name=_name;
		LOG_F( INFO, "Loading MLP weights from : %s", _weight_file.c_str() );

		if ( nullptr == reader ){
			reader = shared_ptr<TMVA::Reader>(new TMVA::Reader( "!Color:!Silent" ) ); 
		
			if ( hasVar( "qdY" ) )
				reader->AddVariable( "qdY := (MtdPidTraits_mDeltaY * Tracks_mCharge)", &MVA_dY );
			if ( hasVar( "dZ" ) )
				reader->AddVariable( "dZ := MtdPidTraits_mDeltaZ", &MVA_dZ );
			if ( hasVar( "nSigmaqdY" ) )
				reader->AddVariable( "nSigmaqdY := (MtdPidTraits_mNSigDeltaY * Tracks_mCharge)", &MVA_nSigmadY );
			if ( hasVar( "nSigmadZ" ) )
				reader->AddVariable( "nSigmadZ := MtdPidTraits_mNSigDeltaZ", &MVA_nSigmadZ );
			if ( hasVar( "nSigmadTof" ) )
				reader->AddVariable( "dTof := MtdPidTraits_mNSigDeltaTOF", &MVA_nSigmadTof );
			if ( hasVar( "nSigmaPi" ) )
				reader->AddVariable( "nSigmaPi := Tracks_mNSigmaPion", &MVA_nSigmaPion );
			if ( hasVar( "nh" ) )
				reader->AddVariable( "nh := Tracks_mNHitsFit", &MVA_nHitsFit );
			if ( hasVar( "dca" ) )
				reader->AddVariable( "dca := Tracks_mDCA", &MVA_DCA );
			if ( hasVar( "Cell" ) )
				reader->AddVariable( "Cell := MtdPidTraits_mCell", &MVA_Cell );
			if ( hasVar( "Module" ) )
				reader->AddVariable( "Module := MtdPidTraits_mModule", &MVA_Module );
			if ( hasVar( "BL" ) )
				reader->AddVariable( "BL := MtdPidTraits_mBL", &MVA_BL );
			if ( hasVar( "pT" ) )
				reader->AddVariable( "pT := Tracks_mPt", &MVA_Pt );
			if ( hasVar( "charge" ) )
				reader->AddVariable( "charge := Tracks_mCharge", &MVA_Charge );
			if ( hasVar( "dTof" ) )
				reader->AddVariable( "dTof := MtdPidTraits_mDeltaTOF", &MVA_dTof );
			
		}
		reader->BookMVA( this->name.c_str(), _weight_file.c_str() ); 
	}

	void load( XmlConfig &_cfg, string _nodePath ){
		
		this->vars = config.getStringVector( _nodePath + ".vars" );
		LOG_F( "vars: %s", vts(vars).c_str() );
		
		string weights_xml = _cfg.getString( _nodePath + ".weights" );
		load( weights_xml );
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