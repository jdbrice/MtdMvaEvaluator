#ifndef FILTER_MUON_MVA_H
#define FILTER_MUON_MVA_H

#include "XmlConfig.h"
#include "XmlRange.h"
#include "FemtoDstFormat/FemtoTrackProxy.h"

#include "TMVA/Reader.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "vendor/loguru.h"

#include <memory>

class MuonMVAFilter {
public:

	string name = "MLP";

	static shared_ptr<TMVA::Reader> reader;
	static vector<string> vars;
	// MLP variables
	static Float_t MVA_dY;
	static Float_t MVA_dZ;
	static Float_t MVA_dTof;
	
	static Float_t MVA_nSigmadY;
	static Float_t MVA_nSigmadZ;
	static Float_t MVA_nSigmadTof;

	static Float_t MVA_nSigmaPion;
	static Float_t MVA_nHitsFit;
	static Float_t MVA_DCA;
	static Float_t MVA_Cell;
	static Float_t MVA_Module;
	static Float_t MVA_BL;
	static Float_t MVA_Pt;
	static Float_t MVA_Charge;

	
	XmlRange signal_range;

	MuonMVAFilter() {}
	MuonMVAFilter(XmlConfig &_cfg, string _nodePath) {
		load( _cfg, _nodePath );
	}
	~MuonMVAFilter() {}

	static double nSigmaDeltaY( float pt, float dy ){
		double sigy = -17.6867 + 18.4528*exp(0.637142/pt);
		return dy / sigy;
	}

	static double nSigmaDeltaZ( float pt, float dz ){
		double sigz = -32.6793 + 32.6034 * exp( 0.444217 / pt );
		return dz / sigz;
	}

	static double nSigmaDeltaTOF( float pt, float dtof ){
		double sigdtof = 0.100915 + 0.00911555 * exp( 3.47368 / pt );
		return dtof / sigdtof;
	}

	static float fillVars(FemtoTrackProxy &_proxy){
		if (nullptr == _proxy._mtdPid)
			return -1;

		MVA_dY = _proxy._mtdPid->mDeltaY;
		MVA_dZ = _proxy._mtdPid->mDeltaZ;
		MVA_dTof = _proxy._mtdPid->mDeltaTimeOfFlight;

		MVA_nSigmadY = nSigmaDeltaY(_proxy._track->mPt, _proxy._mtdPid->mDeltaY);
		MVA_nSigmadZ = nSigmaDeltaZ(_proxy._track->mPt, _proxy._mtdPid->mDeltaZ);
		MVA_nSigmadTof = nSigmaDeltaTOF(_proxy._track->mPt, _proxy._mtdPid->mDeltaTimeOfFlight);

		MVA_nSigmaPion = _proxy._track->nSigmaPion();
		MVA_nHitsFit = (Float_t)fabs(_proxy._track->mNHitsFit);
		MVA_DCA = _proxy._track->gDCA();
		MVA_Cell = (Float_t)_proxy._mtdPid->cell();
		MVA_Module = (Float_t)_proxy._mtdPid->module();
		MVA_BL = (Float_t)_proxy._mtdPid->backleg();
		MVA_Pt = _proxy._track->mPt;
		MVA_Charge = (Float_t)_proxy._track->charge();

		MVA_dY *= MVA_Charge;

		// LOG_F(INFO, "MVA_dY=%f", MVA_dY);
		// LOG_F(INFO, "MVA_dZ=%f", MVA_dZ);
		// LOG_F(INFO, "MVA_dTof=%f", MVA_dTof);
		// LOG_F(INFO, "MVA_nSigmadY=%f", MVA_nSigmadY);
		// LOG_F(INFO, "MVA_nSigmadZ=%f", MVA_nSigmadZ);
		// LOG_F(INFO, "MVA_nSigmadTof=%f", MVA_nSigmadTof);
		// LOG_F(INFO, "MVA_nSigmaPion=%f", MVA_nSigmaPion);
		// LOG_F(INFO, "MVA_nHitsFit=%f", MVA_nHitsFit);
		// LOG_F(INFO, "MVA_DCA=%f", MVA_DCA);
		// LOG_F(INFO, "MVA_Cell=%f", MVA_Cell);
		// LOG_F(INFO, "MVA_Module=%f", MVA_Module);
		// LOG_F(INFO, "MVA_BL=%f", MVA_BL);
		// LOG_F(INFO, "MVA_Pt=%f", MVA_Pt);
		// LOG_F(INFO, "MVA_Charge=%f", MVA_Charge);
	}

	float evaluate( FemtoTrackProxy &_proxy ){
		return reader->EvaluateMVA( this->name.c_str() );
	}

	bool hasVar( string var ){
		if ( std::find( vars.begin(), vars.end(), var ) != vars.end() )
			return true;
		return false;
	}

	void load( string _weight_file, string _name = "MLP" ){

		this->name=_name;
		LOG_F( INFO, "Loading MVA weights from : %s", _weight_file.c_str() );

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

	void loadVars( XmlConfig &_cfg, string _nodePath ){
		LOG_SCOPE_FUNCTION(INFO);
		vars = _cfg.getStringVector( _nodePath + ".vars" );
		LOG_F( INFO, "vars: %s", vts(vars).c_str() );
	}
	void load( XmlConfig &_cfg, string _nodePath ){
		
		this->vars = _cfg.getStringVector( _nodePath + ".vars" );
		LOG_F( INFO, "vars: %s", vts(vars).c_str() );
		
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