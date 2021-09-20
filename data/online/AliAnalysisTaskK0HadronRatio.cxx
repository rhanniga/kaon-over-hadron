/************************************************************************* 
 * Copyright(c) 1998-2021, ALICE Experiment at CERN, All rights reserved. * 
 *                                                                        * 
 * Author: Ryan Hannigan
 * Contributors are mentioned in the code where appropriate.              * 
 *                                                                        * 
 * Permission to use, copy, modify and distribute this software and its   * 
 * documentation strictly for non-commercial purposes is hereby granted   * 
 * without fee, provided that the above copyright notice appears in all   * 
 * copies and that both the copyright notice and this permission notice   * 
 * appear in the supporting documentation. The authors make no claims     * 
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  * 
 **************************************************************************/

#include <iostream>

#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TParticle.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliEventPoolManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODMCHeader.h"
#include "AliCFParticle.h"
#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskK0HadronRatio.h"

ClassImp(AliAnalysisTaskK0HadronRatio);

AliAnalysisTaskK0HadronRatio::AliAnalysisTaskK0HadronRatio() :
    AliAnalysisTaskSE(),
    fAOD(0x0),
    fOutputList(0x0),
    fCorPoolMgr(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fK0Eff(0x0),
    fTriggerDistEff(0x0),
    fTriggeredK0Dist(0x0),
    fTriggeredK0DistFilterbit(0x0),
    fDphiHK0Filterbit(0x0),
    fDphiHK0Eff(0x0),
    fDphiHK0V0(0x0),
    fDphiHHEff(0x0),
    fDphiHK0Mixed(0x0),
    fDphiHHMixed(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0)
{
}

AliAnalysisTaskK0HadronRatio::AliAnalysisTaskK0HadronRatio(const char *name) :
    AliAnalysisTaskSE(name),
    fAOD(0x0),
    fOutputList(0x0),
    fCorPoolMgr(0x0),
    fTriggerEff(0x0),
    fAssociatedEff(0x0),
    fK0Eff(0x0),
    fTriggerDistEff(0x0),
    fTriggeredK0Dist(0x0),
    fTriggeredK0DistFilterbit(0x0),
    fDphiHK0Filterbit(0x0),
    fDphiHK0Eff(0x0),
    fDphiHK0V0(0x0),
    fDphiHHEff(0x0),
    fDphiHK0Mixed(0x0),
    fDphiHHMixed(0x0),
    fpidResponse(0x0),
    fMultSelection(0x0),
    fMultLow(0.0),
    fMultHigh(0.0),
    fDaughterBit(0.0),
    fAssociatedBit(0.0),
    fTriggerBit(0.0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

AliAnalysisTaskK0HadronRatio::~AliAnalysisTaskK0HadronRatio()
{
    if(fOutputList) delete fOutputList;
    if(fTriggerEff) delete fTriggerEff;
    if(fAssociatedEff) delete fAssociatedEff;
    if(fK0Eff) delete fK0Eff;
}

void AliAnalysisTaskK0HadronRatio::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //Generating the mixed event pools:
    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {fMultLow, fMultHigh};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);


    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {200, 16, 20, 10};
    double dist_mins[4] = {0.0, 0, -1, -10};
    double dist_maxes[4] = {20.0, 6.28, 1, 10};

    fTriggerDistEff = new THnSparseF("fTriggerDistEff", "Efficiency Corrected Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDistEff->Sumw2();
    fOutputList->Add(fTriggerDistEff);

    //Mother distribution axes are: Pt, Phi, Eta, Mass, Event multiplicity
    int mother_dist_bins[5] = {100, 16, 20, 100, 10};
    double mother_dist_mins[5] = {0, -3.14, -1, 0.4, 0};
    double mother_dist_maxes[5] = {15, 3.14, 1, 0.6, 100};

    fTriggeredK0Dist = new THnSparseF("fTriggeredK0Dist", "K0 Distribution (with triggered event)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredK0Dist->Sumw2();
    fOutputList->Add(fTriggeredK0Dist);

    fTriggeredK0DistFilterbit = new THnSparseF("fTriggeredK0DistFilter", "K0 Distribution (with triggered event, filterbit 16 on daughters)", 5, mother_dist_bins, mother_dist_mins, mother_dist_maxes);
    fTriggeredK0DistFilterbit->Sumw2();
    fOutputList->Add(fTriggeredK0DistFilterbit);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Inv Mass, Zvtx
    int hl_cor_bins[6] = {8, 10, 16, 20, 100, 10};
    double hl_cor_mins[6] = {4.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 0.4, -10};
    double hl_cor_maxes[6] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 0.6, 10};

    fDphiHK0Filterbit = new THnSparseF("fDphiHK0Filterbit", "Hadron-K0 Correlation Histogram (daughter has filter bit kTrkGlobalNoDCA) ", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHK0Filterbit->Sumw2();
    fOutputList->Add(fDphiHK0Filterbit);

    fDphiHK0Eff = new THnSparseF("fDphiHK0Eff", "Efficiency-corrected Hadron-K0 Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHK0Eff->Sumw2();
    fOutputList->Add(fDphiHK0Eff);

    fDphiHK0V0 = new THnSparseF("fDphiHK0V0", "Hadron-K0 (using V0) Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHK0V0->Sumw2();
    fOutputList->Add(fDphiHK0V0);

    fDphiHK0Mixed = new THnSparseF("fDphiHK0Mixed", "Mixed Hadron-K0 Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHK0Mixed->Sumw2();
    fOutputList->Add(fDphiHK0Mixed);

    fDphiHK0V0Mixed = new THnSparseF("fDphiHK0V0Mixed", "Mixed Hadron-K0 (using V0) Correlation Histogram", 6, hl_cor_bins, hl_cor_mins, hl_cor_maxes);
    fDphiHK0V0Mixed->Sumw2();
    fOutputList->Add(fDphiHK0V0Mixed);

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Zvtx
    int hh_cor_bins[5] = {20, 20, 16, 20, 10};
    double hh_cor_mins[5] = {2, 2, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {12, 12, 3.0*TMath::Pi()/2.0, 2.0, 10};

    fDphiHHEff = new THnSparseF("fDphiHHEff", "Efficiency corrected Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHEff->Sumw2();
    fOutputList->Add(fDphiHHEff);

    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHHMixed->Sumw2();
    fOutputList->Add(fDphiHHMixed);

    PostData(1, fOutputList);
}


void AliAnalysisTaskK0HadronRatio::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist, bool trig_eff)
{
    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i];
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        dist_points[3] = zVtx;
        bool in_pt_range = (particle->Pt() < 10 && particle->Pt() > 0.5);
        if(trig_eff && in_pt_range) {
            int trigBin = fTriggerEff->FindBin(particle->Pt());
            double trigEff = fTriggerEff->GetBinContent(trigBin);
            double triggerScale = 1.0/trigEff;
            fDist->Fill(dist_points, triggerScale);
        }
        else{
            fDist->Fill(dist_points);
        }

    }
}

void AliAnalysisTaskK0HadronRatio::FillMotherDist(std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> particle_list, float multPercentile, THnSparse* fDist)
{
    double dist_points[5]; //Pt, Phi, Eta, M, event multiplicity
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = particle.M();
        dist_points[4] = multPercentile;
        fDist->Fill(dist_points);
    }
}

AliAnalysisTaskK0HadronRatio::AliMotherContainer AliAnalysisTaskK0HadronRatio::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{
    AliAnalysisTaskK0HadronRatio::AliMotherContainer mom;

    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));

    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;
}

AliAnalysisTaskK0HadronRatio::AliMotherContainer AliAnalysisTaskK0HadronRatio::RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle)
{
    AliAnalysisTaskK0HadronRatio::AliMotherContainer mom;
    // Rotating track1
    TVector3 track1Vector(track1->Px(), track1->Py(), track1->Pz());
    track1Vector.RotateZ(angle);
    mom.particle.SetPx(track1Vector(0) + track2->Px());
    mom.particle.SetPy(track1Vector(1) + track2->Py());
    mom.particle.SetPz(track1Vector(2) + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;
}

void AliAnalysisTaskK0HadronRatio::SetMultBounds(float multLow, float multHigh) {
    fMultLow = multLow;
    fMultHigh = multHigh;
}

void AliAnalysisTaskK0HadronRatio::SetTriggerBit(float trigBit) {
    fTriggerBit = trigBit;
}

void AliAnalysisTaskK0HadronRatio::SetAssociatedBit(float associatedBit) {
    fAssociatedBit = associatedBit;
}

void AliAnalysisTaskK0HadronRatio::SetCentEstimator(TString centEstimator) {
    fCentEstimator = centEstimator;
}

void AliAnalysisTaskK0HadronRatio::LoadEfficiencies(TString filePath) {
    TFile* effFile = TFile::Open(filePath);

    if(!effFile) {
        AliFatal("NULL INPUT FILE WHEN LOADING EFFICIENCIES, EXITING");
    }

    // TODO: CHANGE THIS TO CORRECT EFFICIENCY ONCE CACLULATED
    fK0Eff = (TH1D*) effFile->Get("fLambdaEff")->Clone("fK0EffClone");
    if(!fK0Eff) {
        AliFatal("UNABLE TO FIND K0 EFF, EXITING");
    }
    
    fAssociatedEff = (TH1D*) effFile->Get("fAssociatedEff")->Clone("fAssociatedEffClone");
    if(!fAssociatedEff) {
        AliFatal("UNABLE TO FIND ASSOCIATED EFF, EXITING");
    }

    fTriggerEff = (TH1D*) effFile->Get("fTriggerEff")->Clone("fTriggerEffClone");
    if(!fTriggerEff) {
        AliFatal("UNABLE TO FIND TRIGGER EFF, EXITING");
    }
}

AliAnalysisTaskK0HadronRatio::AliMotherContainer AliAnalysisTaskK0HadronRatio::FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{
    AliAnalysisTaskK0HadronRatio::AliMotherContainer mom;
    // Flipping track1
    mom.particle.SetPx(-track1->Px() + track2->Px());
    mom.particle.SetPy(-track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;
}

void AliAnalysisTaskK0HadronRatio::MakeSameHK0Correlations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)K0_list.size(); i++) {
            auto K0 = K0_list[i];

            //Make sure trigger isn't one of the daughters of K0
            if((trigger->GetID() == K0.daughter1ID) || (trigger->GetID() == K0.daughter2ID)) continue;

            dphi_point[1] = K0.particle.Pt();
            dphi_point[2] = trigger->Phi() - K0.particle.Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - K0.particle.Eta();
            dphi_point[4] = K0.particle.M();
            dphi_point[5] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (K0.particle.Pt() < 10 && K0.particle.Pt() > 0.5));

            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int K0Bin = fK0Eff->FindBin(K0.particle.Pt());
                double K0Eff = fK0Eff->GetBinContent(K0Bin);
                double K0Scale = 1.0/K0Eff;
                double totalScale = triggerScale*K0Scale;
                fDphi->Fill(dphi_point, totalScale);

            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskK0HadronRatio::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)associated_h_list.size(); i++) {
            auto associate = associated_h_list[i];

            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (associate->Pt() < 10 && associate->Pt() > 0.5));

            if(eff && in_pt_range) {

                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;

                int associatedBin = fAssociatedEff->FindBin(associate->Pt());
                double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                double associatedScale = 1.0/associatedEff;

                double totalScale = triggerScale*associatedScale;

                fDphi->Fill(dphi_point, totalScale);

            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskK0HadronRatio::MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[5];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];

        dphi_point[0] = trigger->Pt();

        for(int i = j+1; i < (int)trigger_list.size(); i++) {
            auto associate = trigger_list[i];

            dphi_point[1] = associate->Pt();
            dphi_point[2] = trigger->Phi() - associate->Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - associate->Eta();
            dphi_point[4] = zVtx;

            bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                               && (associate->Pt() < 10 && associate->Pt() > 0.5));

            if(eff && in_pt_range) {
                int trigBin = fTriggerEff->FindBin(trigger->Pt());
                double trigEff = fTriggerEff->GetBinContent(trigBin);
                double triggerScale = 1.0/trigEff;
                int associatedBin = fTriggerEff->FindBin(associate->Pt());
                double associatedEff = fTriggerEff->GetBinContent(associatedBin);
                double associatedScale = 1.0/associatedEff;
                double totalScale = triggerScale*associatedScale;
                fDphi->Fill(dphi_point, totalScale);
            }
            else{
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskK0HadronRatio::MakeMixedHK0Correlations(AliEventPool* fPool, std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list , THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[6];
    int numEvents = fPool->GetCurrentNEvents();
    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        tracks->SetName(Form("%d_Zvtx", (int)zVtx));
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            if(!trigger) continue;
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)K0_list.size(); j++) {
                auto K0 = K0_list[j];

                dphi_point[1] = K0.particle.Pt();
                dphi_point[2] = trigger->Phi() - K0.particle.Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - K0.particle.Eta();
                dphi_point[4] = K0.particle.M();
                dphi_point[5] = zVtx;
                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                                && (K0.particle.Pt() < 10 && K0.particle.Pt() > 0.5));
                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    int K0Bin = fK0Eff->FindBin(K0.particle.Pt());
                    double K0Eff = fK0Eff->GetBinContent(K0Bin);
                    double K0Scale = 1.0/K0Eff;
                    double totalScale = triggerScale*K0Scale;
                    fDphi->Fill(dphi_point, totalScale);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }
}

void AliAnalysisTaskK0HadronRatio::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx, bool eff)
{
    double dphi_point[5];

    int numEvents = fPool->GetCurrentNEvents();

    for(int iEvent = 0; iEvent < numEvents; iEvent++) {
        TObjArray *tracks = fPool->GetEvent(iEvent);
        int numTracks = tracks->GetEntriesFast();

        for(int i = 0; i < numTracks; i++) {
            AliCFParticle *trigger = (AliCFParticle*) tracks->At(i);
            dphi_point[0] = trigger->Pt();

            for(int j = 0; j < (int)associated_h_list.size(); j++) {
                auto associate = associated_h_list[j];

                dphi_point[1] = associate->Pt();
                dphi_point[2] = trigger->Phi() - associate->Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - associate->Eta();
                dphi_point[4] = zVtx;

                bool in_pt_range = ((trigger->Pt() < 10 && trigger->Pt() > 0.5) 
                                && (associate->Pt() < 10 && associate->Pt() > 0.5));

                if(eff && in_pt_range) {
                    int trigBin = fTriggerEff->FindBin(trigger->Pt());
                    double trigEff = fTriggerEff->GetBinContent(trigBin);
                    double triggerScale = 1.0/trigEff;
                    int associatedBin = fAssociatedEff->FindBin(associate->Pt());
                    double associatedEff = fAssociatedEff->GetBinContent(associatedBin);
                    double associatedScale = 1.0/associatedEff;
                    double totalScale = triggerScale*associatedScale;
                    fDphi->Fill(dphi_point, totalScale);
                }
                else{
                    fDphi->Fill(dphi_point);
                }
            }
        }
    }
}

bool AliAnalysisTaskK0HadronRatio::PassDaughterCuts(AliAODTrack *track){
    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && (track->IsOn(AliAODTrack::kTPCrefit));

    pass = pass && (track->GetTPCCrossedRows() > 70);

    float ratio = (track->GetTPCNclsF() > 0)  ? track->GetTPCCrossedRows()/track->GetTPCNclsF() : 0;
    pass = pass && (ratio > 0.8);

    return pass;
}

bool AliAnalysisTaskK0HadronRatio::PassAssociatedCuts(AliAODTrack *track){
    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestFilterMask(fAssociatedBit);

    return pass;
}

Bool_t AliAnalysisTaskK0HadronRatio::PassTriggerCuts(AliAODTrack *track){
    Bool_t pass = kTRUE;

    pass = pass && (TMath::Abs(track->Eta()) <= 0.8);
    pass = pass && (track->Pt() >= 0.15);

    pass = pass && track->TestBit(fTriggerBit);

    return pass;
}

void AliAnalysisTaskK0HadronRatio::UserExec(Option_t*)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        AliFatal("THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!");
    }


    fpidResponse = fInputHandler->GetPIDResponse();

    //Event cuts
    TString cent_estimator = fCentEstimator;
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;

    if(multPercentile < fMultLow || multPercentile > fMultHigh) return;

    AliVVertex *prim = fAOD->GetPrimaryVertex();
    int NcontV = prim->GetNContributors();
    if(NcontV < 3) return;

    double primZ = prim->GetZ();
    if(primZ < -10 || primZ > 10) return;


    int numTracks = fAOD->GetNumberOfTracks();

    std::vector<AliAODTrack*> piPlus_list;
    std::vector<AliAODTrack*> piMinus_list;
    std::vector<AliAODTrack*> filterbit_piPlus_list;
    std::vector<AliAODTrack*> filterbit_piMinus_list;
    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    // Bool to keep track if the event has a high-pt (> 4 GeV) trigger
    bool is_triggered_event = false;


    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //Filter for trigger particles
        if(PassTriggerCuts(track)) {
            trigger_list.push_back(track);
            AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
            fMixedTrackObjArray->Add(triggerPart);
            if(triggerPart->Pt() > 4) is_triggered_event = true;
        }

        if(PassAssociatedCuts(track)) {
            associated_h_list.push_back(track);
        }

        if(track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA)) {

            double filterbit_TPCNSigmaPion = 1000;
            double filterbit_TOFNSigmaPion = 1000;

            filterbit_TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            filterbit_TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

            if(TMath::Abs(filterbit_TPCNSigmaPion) <= 3 && (TMath::Abs(filterbit_TOFNSigmaPion) <= 3 || filterbit_TOFNSigmaPion == 1000)) {

                if(track->Charge() == 1){
                    filterbit_piPlus_list.push_back(track);
                }
                else {
                    filterbit_piMinus_list.push_back(track);
                }
            }
        } 
            
        if(PassDaughterCuts(track)) {
            double TPCNSigmaPion = 1000;
            double TOFNSigmaPion = 1000;

            TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
            TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

            double time = track->GetTOFsignal() - fpidResponse->GetTOFResponse().GetStartTime(track->P());
            double length = track->GetIntegratedLength();
            double v = length/time;
            double c = 0.0288782;
            double beta = v/c;


            if(TMath::Abs(TPCNSigmaPion) <= 3 && (TMath::Abs(TOFNSigmaPion) <= 3 || TOFNSigmaPion == 1000)) {

                if(track->Charge() == 1){
                    piPlus_list.push_back(track);
                }
                else {
                    piMinus_list.push_back(track);
                }
            }
        }
    }

    //Making list of possible K0s (have to do +/- for proton or pi):

    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_filterbit_daughters;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_v0;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_signal_region;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_signal_region_2_4;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_RotatedPion;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_RotatedPiMinus;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_RotatedPi;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_Flipped;
    std::vector<AliAnalysisTaskK0HadronRatio::AliMotherContainer> K0_list_LS;

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) piMinus_list.size(); j++) {
            AliMotherContainer K0 = DaughtersToMother(piPlus_list[i], piMinus_list[j], 0.1396, 0.1396);
            K0_list.push_back(K0);
        }
    }

    for(int i = 0; i < (int)filterbit_piPlus_list.size(); i++) {
        for(int j = 0; j < (int) filterbit_piMinus_list.size(); j++) {
            AliMotherContainer filterbit_K0 = DaughtersToMother(filterbit_piPlus_list[i], filterbit_piMinus_list[j], 0.1396, 0.1396);
            K0_list_filterbit_daughters.push_back(filterbit_K0);
        }
    }


    // V0 SECTION

    int numV0s = fAOD->GetNumberOfV0s();
    for(int i = 0; i < numV0s; i++) {
        AliAODv0 *v0 = fAOD->GetV0(i);

        AliAODTrack* posTrack = (AliAODTrack*) v0->GetDaughter(0);
        AliAODTrack* negTrack = (AliAODTrack*) v0->GetDaughter(1);

        // Occasionally returns null, not quite sure why...
        if(!posTrack || !negTrack) continue;
        if(!(PassDaughterCuts(posTrack) && PassDaughterCuts(negTrack))) continue;

        double TPCNSigmaPiPlus = 1000;
        double TOFNSigmaPiPlus = 1000;

        double TPCNSigmaPiMinus = 1000;
        double TOFNSigmaPiMinus = 1000;

        TPCNSigmaPiPlus = fpidResponse->NumberOfSigmasTPC(posTrack, AliPID::kPion);
        TOFNSigmaPiPlus = fpidResponse->NumberOfSigmasTOF(posTrack, AliPID::kPion);

        TPCNSigmaPiMinus = fpidResponse->NumberOfSigmasTPC(negTrack, AliPID::kPion);
        TOFNSigmaPiMinus = fpidResponse->NumberOfSigmasTOF(negTrack, AliPID::kPion);

        bool isPosTrackPion = TMath::Abs(TPCNSigmaPiPlus) <= 3 && (TMath::Abs(TOFNSigmaPiPlus) <= 3 || TOFNSigmaPiPlus == 1000);
        bool isNegTrackPion = TMath::Abs(TPCNSigmaPiMinus) <= 3 && (TMath::Abs(TOFNSigmaPiMinus) <= 3 || TOFNSigmaPiMinus == 1000);

        if(isNegTrackPion && isPosTrackPion) {
            auto K0 = DaughtersToMother(negTrack, posTrack, 0.1396, 0.1396);
            K0_list_v0.push_back(K0);
        }
    }


    // Filling all of our single particle distribution histograms:
    FillSingleParticleDist(trigger_list, primZ, fTriggerDistEff, true);

    // Filling our single particle K0 distribution histogram:
    if(is_triggered_event) FillMotherDist(K0_list, multPercentile, fTriggeredK0Dist);
    if(is_triggered_event) FillMotherDist(K0_list_filterbit_daughters, multPercentile, fTriggeredK0DistFilterbit);

    // Filling all of our correlation histograms
    MakeSameHK0Correlations(trigger_list, K0_list_filterbit_daughters, fDphiHK0Filterbit, primZ, false);
    MakeSameHK0Correlations(trigger_list, K0_list, fDphiHK0Eff, primZ, true);
    MakeSameHK0Correlations(trigger_list, K0_list_v0, fDphiHK0V0, primZ);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHHEff, primZ, true);

    if(K0_list.size() > 0 && associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }


        else {
            if(fCorPool->IsReady()) {
                MakeMixedHK0Correlations(fCorPool, K0_list, fDphiHK0Mixed, primZ);
                MakeMixedHK0Correlations(fCorPool, K0_list_v0, fDphiHK0V0Mixed, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);
}

void AliAnalysisTaskK0HadronRatio::Terminate(Option_t *option)
{
}