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

#include "Filters/MuonBDTFilter.h"
#include "Filters/MuonMLPFilter.h"




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

	MuonMLPFilter mlp;
	MuonBDTFilter bdt;

	XmlHistogram bg_deltaTOF;
	XmlFunction sig_deltaTOF;
	vector<TH1*> bg_deltaTOF_dY;
	TH1 * bg_deltaTOF_dY_all;

	TNtuple *tuple;
	TRandom3 r3;
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

		mlp.load( config, nodePath + ".MuonMLPFilter" );
		bdt.load( config, nodePath + ".MuonBDTFilter" );

		bg_deltaTOF.load( config, nodePath + ".DeltaTOF.XmlHistogram" );
		sig_deltaTOF.set( config, nodePath + ".DeltaTOF.XmlFunction" );

		tuple= new TNtuple( "MvaTree", "mva variables", "classId:dY:dZ:dTof:nsp:nhf:dca:cell:mod:bl:pt:c:bdt:mlp:gid" );

		book->cd();

		int seed= config.getInt( "seed", 0 );
		r3.SetSeed( seed );
		gRandom->SetSeed(r3.Integer( 32000 ));
		LOG_F( INFO, "RANDOM SEED : %d", seed );

		size_t nBgdTofBins = bg_deltaTOF.getTH1()->GetXaxis()->GetNbins() + 1;
		TH2 * h2d = (TH2*) bg_deltaTOF.getTH1().get();
		for ( size_t i = 1; i < nBgdTofBins; i++ ){
			TH1 * h1d = h2d->ProjectionY( TString::Format( "bg_dTOF_%lu", i ), i, i );
			h1d->Write();
			bg_deltaTOF_dY.push_back( h1d );
		}
		bg_deltaTOF_dY_all = h2d->ProjectionY( "bg_dTOF_all" );


	}


protected:


	float sample_bg_deltaTOF( float deltaY ){
		int b = bg_deltaTOF.getTH1()->GetXaxis()->FindBin( deltaY );
		if ( b >= 1 && b < bg_deltaTOF_dY.size() ){
			float sdtof = bg_deltaTOF_dY[ b - 1]->GetRandom() + r3.Gaus( 0.0, 0.2 );
			book->fill( "bg_dtof_vs_dY", deltaY, sdtof );
			return sdtof;
		}

		float sdtof =bg_deltaTOF_dY_all->GetRandom();
		book->fill( "bg_dtof_vs_dY", deltaY, sdtof );
		return sdtof;
	} 

	virtual void preEventLoop(){
		TreeAnalyzer::preEventLoop();
		book->cd();
	}

	bool isMuon( FemtoMcTrack * mcTrack ){
		if ( nullptr == mcTrack )
			return false;
		if ( 5 == mcTrack->mGeantPID || 6 == mcTrack->mGeantPID )
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

	virtual void analyzeEvent(){
	
		_event = _rEvent.get();

		size_t nTracks = _rTracks.N();
		FemtoTrackProxy _proxy;

		for (size_t i = 0; i < nTracks; i++ ){
			
			FemtoTrack * track = _rTracks.get(i);
			FemtoMcTrack * mcTrack = nullptr;
			FemtoMtdPidTraits *mtdPid = nullptr;

			if ( track->mMtdPidTraitsIndex >= 0) 
				mtdPid = _rMtdPid.get( track->mMtdPidTraitsIndex );
			if ( track->mMcIndex >= 0 )
				mcTrack = _rMcTracks.get( track->mMcIndex );

			if ( nullptr == mtdPid || nullptr == mcTrack )
				continue;

			FemtoTrackProxy _proxy;
			_proxy._track = track;
			_proxy._mtdPid = mtdPid;
			_proxy._mcTrack = mcTrack;

			bool inTPC = isDecayMuonInsideTPC(mcTrack);
			bool outTPC = isDecayMuonOutsideTPC( mtdPid, mcTrack );
			// if ( inTPC != decayInsideTPC ){
			// 	continue;
			// }
			// if ( outTPC != true ){
			// 	continue;
			// }

			// _proxy.assemble( i, _rTracks, _rMtdPid );
			// if ( nullptr == _proxy._mtdPid  ) continue;
			// // LOG_F( INFO, "Track Mc Index=%d", _proxy._track->mMcIndex );

			// FemtoMcTrack * mcTrack = nullptr;
			// if ( _proxy._track->mMcIndex >= 0 )
			// 	mcTrack = _rMcTracks.get( _proxy._track->mMcIndex );

			// if ( nullptr == mcTrack  ) continue;

			double deltaTof = -999;
			if ( isMuon( mcTrack  )){
				// LOG_F( INFO, "Signal" );
				float sdtof = sig_deltaTOF.getTF1()->GetRandom();
				book->fill( "signal_dtof", sdtof );
				deltaTof = sdtof;
			} else {
				// LOG_F( INFO, "Background" );
				deltaTof = sample_bg_deltaTOF( _proxy._mtdPid->mDeltaY );// bg_deltaTOF.getTH1()->GetRandom();
			}

			_proxy._mtdPid->mDeltaTimeOfFlight = deltaTof;

			// if ( _proxy._track->mPt < 2.0  ) continue;
			float mlpr = mlp.evaluate( _proxy );
			float bdtr = bdt.evaluate( _proxy );
			
			// book->fill( "mlp", mlpr );
			// book->fill( "bdt", bdtr );

			// book->fill( "pt_vs_mlp", mlpr, _proxy._track->mPt );
			// book->fill( "pt_vs_bdt", bdtr, _proxy._track->mPt );

			// LOG_F( INFO, "mGeantPID=%d", mcTrack->mGeantPID  );

			float classId = 0; // bg
			
			if ( isDecayMuonOutsideTPC( _proxy._mtdPid, mcTrack ) ){
				classId = 2;
			}
			if ( isDecayMuonInsideTPC( mcTrack ) ){
				classId = 1;
			}
			
			// if ( cleanPunchThrough( _proxy._mtdPid, track, mcTrack ) ){
			// 	classId = 3;
			// }
			
			if ( isSignal( mcTrack ) ){

				classId = 5;
			}

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

		

	}


	virtual void postEventLoop(){
		tuple->Write();
		TreeAnalyzer::postEventLoop();

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}
		
	}
	
};

#endif
