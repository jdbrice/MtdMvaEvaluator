

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>


#include "FemtoDstSkimmer/FemtoDstSkimmer.h"
shared_ptr<TMVA::Reader> MuonMVAFilter::reader = nullptr;
vector<string> MuonMVAFilter::vars;
Float_t MuonMVAFilter::MVA_dY = -999.9;
Float_t MuonMVAFilter::MVA_dZ = -999.9;
Float_t MuonMVAFilter::MVA_dTof = -999.9;
Float_t MuonMVAFilter::MVA_nSigmadY = -999.9;
Float_t MuonMVAFilter::MVA_nSigmadZ = -999.9;
Float_t MuonMVAFilter::MVA_nSigmadTof = -999.9;
Float_t MuonMVAFilter::MVA_nSigmaPion = -999.9;
Float_t MuonMVAFilter::MVA_nHitsFit = -999.9;
Float_t MuonMVAFilter::MVA_DCA = -999.9;
Float_t MuonMVAFilter::MVA_Cell = -999.9;
Float_t MuonMVAFilter::MVA_Module = -999.9;
Float_t MuonMVAFilter::MVA_BL = -999.9;
Float_t MuonMVAFilter::MVA_Pt = -999.9;
Float_t MuonMVAFilter::MVA_Charge = -999.9;

#define LOGURU_IMPLEMENTATION 1
#include "vendor/loguru.h"

int main( int argc, char* argv[] ) {
	// loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);
	
	TaskFactory::registerTaskRunner<FemtoDstSkimmer>( "FemtoDstSkimmer" );

	TaskEngine engine( argc, argv );

	return 0;
}
