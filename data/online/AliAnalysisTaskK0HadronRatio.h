/* Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full copyright notice */
/*
                AliAnalysisTaskK0HadronRatio class
                This task is for determining the k0/hadron (~pion) ratio
                in different kinematic regions using a two-particle azimuthal
                correlation method
                Origin: Ryan Hannigan, July 2021, ryan.hannigan@cern.ch
*/

#ifndef AliAnalysisTaskK0HadronRatio_H
#define AliAnalysisTaskK0HadronRatio_H

// Includes order: standard, ROOT, AliRoot (for objects in this file only)
#include "TString.h"
#include "TLorentzVector.h"

#include "AliAnalysisTaskSE.h"

// Forward declarations order: standard, ROOT, AliRoot (for pointers in this file only)
class TList;
class TH1D;
class TH1F;
class THnSparse;

class AliAODEvent;
class AliEventPoolManager;
class AliPIDResponse;
class AliMultSelection;
class AliEventPool;
class AliAODTrack;


class AliAnalysisTaskK0HadronRatio : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskK0HadronRatio();
  AliAnalysisTaskK0HadronRatio(const char *name);
  virtual ~AliAnalysisTaskK0HadronRatio();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);
  void LoadEfficiencies(TString filePath);

  void SetMultBounds(float multLow, float multHigh);
  void SetTriggerBit(float trigBit);
  void SetAssociatedBit(float assocBit);
  void SetCentEstimator(TString estimator);

  struct AliMotherContainer {
    TLorentzVector particle;
    int daughter1ID;
    int daughter2ID;
  };

private:
  float fMultLow; // lower bound for multiplicity
  float fMultHigh; // upper bound for multiplicity
  float fDaughterBit; // filter bit for daughter particle
  float fAssociatedBit; // filter bit for associated particle
  float fTriggerBit; // filter bit for trigger particle

  TString fCentEstimator;

  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager
  
  TH1D* fTriggerEff; ///> trigger efficiency
  TH1D* fAssociatedEff; ///> associated efficiency
  TH1D* fK0Eff; ///> K0 efficiency

  THnSparse* fTriggerDistEff;  //!>! single particle trigger dist (corrected for efficiency)

  THnSparse* fTriggeredK0Dist;  //!>! single particle K0 dist within a triggered event
  THnSparse* fTriggeredK0DistFilterbit;  //!>! single particle K0 dist where daughters have filter bit 16 within a triggered event

  THnSparse* fDphiHK0Filterbit;  //!>! hadron-K0 correlation hist where daughter has filter bit 16
  THnSparse* fDphiHK0Eff;  //!>! hadron-K0 correlation hist (efficiency corrected)
  THnSparse* fDphiHK0V0;  //!>! hadron-K0 correlation hist (using v0 finder for K0)
  THnSparse* fDphiHHEff;   //!>! hadron-hadron correlation hist (efficiency corrected)
  THnSparse* fDphiHK0Mixed; //!>! hadron-K0 mixed correlation hist
  THnSparse* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparse* fDphiHK0V0Mixed;  //!>! hadron-K0 mixed correlation hist (using v0 finder for K0)

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection


  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  AliMotherContainer RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle);
  AliMotherContainer FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff=false);
  void FillMotherDist(std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist);
  void MakeSameHK0Correlations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHK0Correlations(AliEventPool *fPool, std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list, THnSparse* fDphi, double zVtx, bool eff=true);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx, bool eff=true);
  bool PassDaughterCuts(AliAODTrack *track);
  bool PassTriggerCuts(AliAODTrack *track);
  bool PassAssociatedCuts(AliAODTrack *track);

  ClassDef(AliAnalysisTaskK0HadronRatio, 3);

};
#endif