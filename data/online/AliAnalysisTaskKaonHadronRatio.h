#ifndef AliAnalysisTaskKaonHadronRatio_H
#define AliAnalysisTaskKaonHadronRatio_H

//All relevant AliPhysics includes (this list will continue to grow):
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH3D.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliPID.h"
#include "TChain.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliCFParticle.h"

//These includes probably aren't necessary but I feel like it's proper etiquette
#include <vector>
#include <iostream>

class AliAnalysisTaskKaonHadronRatio : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskKaonHadronRatio();
  AliAnalysisTaskKaonHadronRatio(const char *name);
  virtual ~AliAnalysisTaskKaonHadronRatio();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);

  virtual void Terminate(Option_t* option);

  struct AliMotherContainer {
    TLorentzVector particle;
    int daughter1ID;
    int daughter2ID;
  };

 private:
  AliAODEvent* fAOD; //!>! input event
  TList* fOutputList; //!>! output list

  AliEventPoolManager *fCorPoolMgr; //!>! correlation pool manager

  TH1D* fMultDist; //!>! event multiplicity dist

  TH2D* fTriggersAndKaonsPerEvent_All; //!>! triggers and all Kaons per event
  TH2D* fTriggersAndKaonsPerEvent_2_4; //!>! triggers and 2-4 GeV Kaons per event

  THnSparseF* fLooseDist;  //!>! single particle all hadron dist (no cuts at all)
  THnSparseF* fTriggerDist;  //!>! single particle trigger dist
  THnSparseF* fAssociatedHDist;  //!>! single particle associated hadron dist
  THnSparseF* fKaonDist;  //!>! single particle Kaon dist

  THnSparseF* fDphiHKaon;  //!>! hadron-Kaon correlation hist
  THnSparseF* fDphiHKaonRotated;  //!>! hadron-Kaon correlation hist with rotated pion
  THnSparseF* fDphiHKaonRotatedPi;  //!>! hadron-Kaon correlation hist with daughter rotated by pi+
  THnSparseF* fDphiHKaonRotatedPiMinus;  //!>! hadron-Kaon correlation hist with rotated pi-
  THnSparseF* fDphiHKaonFlipped;  //!>! hadron-Kaon correlation hist with flipped pion
  THnSparseF* fDphiHH;   //!>! hadron-hadron correlation hist
  THnSparseF* fDphiTriggerTrigger;   //!>! trigger-trigger correlation hist
  THnSparseF* fDphiHKaonLS; //!>! hadron-proton+pion like sign correlation hist
  THnSparseF* fDphiHKaonMixed; //!>! hadron-Kaon mixed correlation hist
  THnSparseF* fDphiHHMixed; //!>! hadron-hadron mixed correlation hist
  THnSparseF* fDphiHKaonLSMixed; //!>! hadron-proton+pion like sign mixed correlation hist
  THnSparseF* fDphiTriggerTriggerMixed;   //!>! mixed trigger-trigger correlation hist

  TH3D* fKaonDaughterDCA;

  // THnSparseF* fPid; //!>! histogram to visualize pid cuts
  // THnSparseF* fSignalAnalysis; //!>! histogram to analyze signal with nsigma cuts

  AliPIDResponse *fpidResponse; //!>!pid response
  AliMultSelection *fMultSelection; //!>!mult selection

  //hand written functions:

  AliMotherContainer DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  AliMotherContainer RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle);
  AliMotherContainer FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2);
  void FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist);
  void FillSingleParticleDist(std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist);
  void MakeSameHKaonCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list, THnSparse* fDphi, double zVtx);
  void MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx);
  void MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx);
  void MakeMixedHKaonCorrelations(AliEventPool *fPool, std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list, THnSparse* fDphi, double zVtx);
  void MakeMixedHHCorrelations(AliEventPool *fPool, std::vector<AliAODTrack*> associated_h_list , THnSparse* fDphi, double zVtx);

  ClassDef(AliAnalysisTaskKaonHadronRatio, 0);

};
#endif
