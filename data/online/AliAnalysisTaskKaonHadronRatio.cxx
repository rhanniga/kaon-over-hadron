#include "AliAnalysisTaskKaonHadronRatio.h"

//BOTH OF THESE ARE WEIRD, BUT APPARENTLY NECESSARRY
class AliAnalysisTaskKaonHadronRatio;
ClassImp(AliAnalysisTaskKaonHadronRatio);

static const int centLow = 20;
static const int centHigh = 50;

AliAnalysisTaskKaonHadronRatio::AliAnalysisTaskKaonHadronRatio() :
    AliAnalysisTaskSE(),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHKaon{0},
    fDphiHKaonRotated{0},
    fDphiHKaonRotatedPiMinus{0},
    fDphiHKaonRotatedPi{0},
    fDphiHKaonFlipped{0},
    fDphiHH{0},
    fDphiHKaonLS{0},
    fDphiHKaonMixed{0},
    fDphiHHMixed{0},
    fDphiHKaonLSMixed{0},
    fCorPoolMgr{0}
    // fPid{0},
    // fSignalAnalysis{0}
{

}

AliAnalysisTaskKaonHadronRatio::AliAnalysisTaskKaonHadronRatio(const char *name) :
    AliAnalysisTaskSE(name),
    fpidResponse{0},
    fAOD{0},
    fOutputList{0},
    fTriggerDist{0},
    fAssociatedHDist{0},
    fDphiHKaon{0},
    fDphiHKaonRotated{0},
    fDphiHKaonRotatedPiMinus{0},
    fDphiHKaonRotatedPi{0},
    fDphiHKaonFlipped{0},
    fDphiHH{0},
    fDphiHKaonLS{0},
    fDphiHKaonMixed{0},
    fDphiHHMixed{0},
    fDphiHKaonLSMixed{0},
    fCorPoolMgr{0}
    // fPid{0},
    // fSignalAnalysis{0}
{

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());

}

AliAnalysisTaskKaonHadronRatio::~AliAnalysisTaskKaonHadronRatio()
{

    if(fOutputList) delete fOutputList;

}

void AliAnalysisTaskKaonHadronRatio::UserCreateOutputObjects()
{

    fOutputList = new TList();
    fOutputList->SetOwner(true);

    //Generating the mixed event pools:

    int poolSize = 500;
    int trackDepth = 1000;

    int numMultBins = 1;
    double multBins[2] = {centLow, centHigh};

    int numzVtxBins = 10;
    double zVtxBins[11] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

    fCorPoolMgr = new AliEventPoolManager(poolSize, trackDepth, numMultBins, multBins, numzVtxBins, zVtxBins);
    fCorPoolMgr->SetTargetValues(trackDepth, 0.1, 5);

    fTriggersAndKaonsPerEvent_All = new TH2D("fTriggersAndKaonsPerEvent_All", "Triggers and Kaons per event (all p_{T})", 10, 0, 10, 10, 0, 10);
    fTriggersAndKaonsPerEvent_2_4 = new TH2D("fTriggersAndKaonsPerEvent_2_4", "Triggers and Kaons per event (2-4 p_{T})", 10, 0, 10, 10, 0, 10);
    fOutputList->Add(fTriggersAndKaonsPerEvent_All);
    fOutputList->Add(fTriggersAndKaonsPerEvent_2_4);


    //Distribution axes are: Pt, Phi, Eta, zVtx
    int dist_bins[4] = {100, 16, 20, 10};
    double dist_mins[4] = {0, 0, -1, -10};
    double dist_maxes[4] = {15, 6.28, 1, 10};

    fLooseDist = new THnSparseF("fLooseDist", "All Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fTriggerDist = new THnSparseF("fTriggerDist", "Trigger Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fAssociatedHDist = new THnSparseF("fAssociatedHDist", "Associated Hadron Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fKaonDist = new THnSparseF("fKaonDist", "Kaon Distribution", 4, dist_bins, dist_mins, dist_maxes);
    fOutputList->Add(fLooseDist);
    fOutputList->Add(fTriggerDist);
    fOutputList->Add(fAssociatedHDist);
    fOutputList->Add(fKaonDist);
    

    //Correlation axes are: Trigger Pt, Associated Pt, dPhi, dEta, Associated Inv Mass (Kaon only), Zvtx
    int hk_cor_bins[6] = {8, 10, 16, 20, 100, 10};
    double hk_cor_mins[6] = {4.0, 1, -1.0*TMath::Pi()/2.0, -2.0, 0.4, -10};
    double hk_cor_maxes[6] = {12.0, 6, 3.0*TMath::Pi()/2.0, 2.0, 0.6, 10};

    int hh_cor_bins[5] = {20, 20, 16, 20, 10};
    double hh_cor_mins[5] = {2, 2, -1.0*TMath::Pi()/2.0, -2.0, -10};
    double hh_cor_maxes[5] = {12, 12, 3.0*TMath::Pi()/2.0, 2.0, 10};



    fDphiHKaon = new THnSparseF("fDphiHKaon", "Hadron-Kaon Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHKaonRotated = new THnSparseF("fDphiHKaonRotated", "Hadron-Kaon (rotated) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHKaonRotatedPi = new THnSparseF("fDphiHKaonRotatedPi", "Hadron-Kaon (rotated pi+) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHKaonRotatedPiMinus = new THnSparseF("fDphiHKaonRotatedPiMinus", "Hadron-Kaon (rotated pi-) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHKaonFlipped = new THnSparseF("fDphiHKaonFlipped", "Hadron-Kaon (flipped) Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHH = new THnSparseF("fDphiHH", "Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiTriggerTrigger = new THnSparseF("fDphiTriggerTrigger", "Trigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHKaonLS = new THnSparseF("fDphiHKaonLS", "Hadron-Kaon LS Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHKaonMixed = new THnSparseF("fDphiHKaonMixed", "Mixed Hadron-Kaon Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fDphiHHMixed = new THnSparseF("fDphiHHMixed", "Mixed Hadron-Hadron Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiTriggerTriggerMixed = new THnSparseF("fDphiTriggerTriggerMixed", "MixedTrigger-Trigger Correlation Histogram", 5, hh_cor_bins, hh_cor_mins, hh_cor_maxes);
    fDphiHKaonLSMixed = new THnSparseF("fDphiHKaonLSMixed", "Mixed Hadron-Kaon LS Correlation Histogram", 6, hk_cor_bins, hk_cor_mins, hk_cor_maxes);
    fOutputList->Add(fDphiHKaon);
    fOutputList->Add(fDphiHKaonRotated);
    fOutputList->Add(fDphiHKaonRotatedPi);
    fOutputList->Add(fDphiHKaonRotatedPiMinus);
    fOutputList->Add(fDphiHKaonFlipped);
    fOutputList->Add(fDphiHH);
    fOutputList->Add(fDphiTriggerTrigger);
    fOutputList->Add(fDphiHKaonLS);
    fOutputList->Add(fDphiHKaonMixed);
    fOutputList->Add(fDphiHHMixed);
    fOutputList->Add(fDphiTriggerTriggerMixed);
    fOutputList->Add(fDphiHKaonLSMixed);

    //axes are pion DCA proton DCA Kaon mass
    fKaonDaughterDCA = new TH3D("fKaonDaughterDCA", "K^{0}_{S} daughter DCA dist", 100, -5, 5, 100, -5, 5, 100, 1.08, 1.16);
    fOutputList->Add(fKaonDaughterDCA);

    // int pid_bins[7] = {500, 500, 500, 50, 50, 50, 50};
    // double pid_mins[7] = {0, 0, 10000, -10, -10, -10, -10};
    // double pid_maxes[7] = {10, 160, 30000, 10, 10, 10, 10};

    // fPid = new THnSparseF("fPid", "PID Histogram", 7, pid_bins, pid_mins, pid_maxes);
    // fOutputList->Add(fPid);

    // int signal_bins[6] = {50, 50, 50, 50, 100, 50};
    // double signal_mins[6] = {-10, -10, -10, -10, 1.06, 0};
    // double signal_maxes[6] = {10, 10, 10, 10, 1.16, 10};

    // fSignalAnalysis = new THnSparseF("fSignalAnalysis", "Signal Analysis Histogram", 6, signal_bins, signal_mins, signal_maxes);
    // fOutputList->Add(fSignalAnalysis);

    fMultDist = new TH1D("fMultDist", "Event Multiplicty Distribution", 100, 0, 100);
    fOutputList->Add(fMultDist);

    PostData(1, fOutputList);

}

void AliAnalysisTaskKaonHadronRatio::FillSingleParticleDist(std::vector<AliAODTrack*> particle_list, double zVtx, THnSparse* fDist)
{

    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i];
        dist_points[0] = particle->Pt();
        dist_points[1] = particle->Phi();
        dist_points[2] = particle->Eta();
        dist_points[3] = zVtx;
        fDist->Fill(dist_points);
    }

}

void AliAnalysisTaskKaonHadronRatio::FillSingleParticleDist(std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> particle_list, double zVtx, THnSparse* fDist)
{

    double dist_points[4]; //Pt, Phi, Eta, zVtx
    for(int i = 0; i < (int)particle_list.size(); i++) {
        auto particle = particle_list[i].particle;
        dist_points[0] = particle.Pt();
        dist_points[1] = particle.Phi();
        dist_points[2] = particle.Eta();
        dist_points[3] = zVtx;
        fDist->Fill(dist_points);
    }

}

AliAnalysisTaskKaonHadronRatio::AliMotherContainer AliAnalysisTaskKaonHadronRatio::DaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskKaonHadronRatio::AliMotherContainer mom;
    mom.particle.SetPx(track1->Px() + track2->Px());
    mom.particle.SetPy(track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

AliAnalysisTaskKaonHadronRatio::AliMotherContainer AliAnalysisTaskKaonHadronRatio::RotatedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2, double angle)
{

    AliAnalysisTaskKaonHadronRatio::AliMotherContainer mom;
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

AliAnalysisTaskKaonHadronRatio::AliMotherContainer AliAnalysisTaskKaonHadronRatio::FlippedDaughtersToMother(AliAODTrack* track1, AliAODTrack* track2, double mass1, double mass2)
{

    AliAnalysisTaskKaonHadronRatio::AliMotherContainer mom;
    // Flipping track1
    mom.particle.SetPx(-track1->Px() + track2->Px());
    mom.particle.SetPy(-track1->Py() + track2->Py());
    mom.particle.SetPz(track1->Pz() + track2->Pz());
    mom.particle.SetE(track1->E(mass1) + track2->E(mass2));
    mom.daughter1ID = track1->GetID();
    mom.daughter2ID = track2->GetID();
    return mom;

}

void AliAnalysisTaskKaonHadronRatio::MakeSameHKaonCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list, THnSparse* fDphi, double zVtx)
{

    double dphi_point[6];

    for(int j = 0; j < (int)trigger_list.size(); j++) {
        auto trigger = trigger_list[j];
        dphi_point[0] = trigger->Pt();

        for(int i = 0; i < (int)kaon_list.size(); i++) {
            auto kaon = kaon_list[i];

            //Make sure trigger isn't one of the daughters of kaon
            if((trigger->GetID() == kaon.daughter1ID) || (trigger->GetID() == kaon.daughter2ID)) continue;

            dphi_point[1] = kaon.particle.Pt();
            dphi_point[2] = trigger->Phi() - kaon.particle.Phi();

            if(dphi_point[2] < -TMath::Pi()/2.0) {
                dphi_point[2] += 2.0*TMath::Pi();
            }
            else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                dphi_point[2] -= 2.0*TMath::Pi();
            }

            dphi_point[3] = trigger->Eta() - kaon.particle.Eta();
            dphi_point[4] = kaon.particle.M();
            dphi_point[5] = zVtx;
            fDphi->Fill(dphi_point);
        }
    }

}

void AliAnalysisTaskKaonHadronRatio::MakeSameHHCorrelations(std::vector<AliAODTrack*> trigger_list, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx)
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
            fDphi->Fill(dphi_point);
        }
    }

}

void AliAnalysisTaskKaonHadronRatio::MakeSameTriggerTriggerCorrelations(std::vector<AliAODTrack*> trigger_list, THnSparse* fDphi, double zVtx)
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
            fDphi->Fill(dphi_point);
        }
    }

}

void AliAnalysisTaskKaonHadronRatio::MakeMixedHKaonCorrelations(AliEventPool* fPool, std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list , THnSparse* fDphi, double zVtx)
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

            for(int j = 0; j < (int)kaon_list.size(); j++) {
                auto kaon = kaon_list[j];

                dphi_point[1] = kaon.particle.Pt();
                dphi_point[2] = trigger->Phi() - kaon.particle.Phi();

                if(dphi_point[2] < -TMath::Pi()/2.0) {
                    dphi_point[2] += 2.0*TMath::Pi();
                }
                else if(dphi_point[2] > 3.0*TMath::Pi()/2.0) {
                    dphi_point[2] -= 2.0*TMath::Pi();
                }

                dphi_point[3] = trigger->Eta() - kaon.particle.Eta();
                dphi_point[4] = kaon.particle.M();
                dphi_point[5] = zVtx;
                fDphi->Fill(dphi_point);
            }
        }
    }
}

void AliAnalysisTaskKaonHadronRatio::MakeMixedHHCorrelations(AliEventPool* fPool, std::vector<AliAODTrack*> associated_h_list, THnSparse* fDphi, double zVtx)
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
                fDphi->Fill(dphi_point);
            }
        }
    }

}

void AliAnalysisTaskKaonHadronRatio::UserExec(Option_t*)
{

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
        std::cout << "THERE IS NO AOD EVENT, CHECK EVENT HANDLER... ALSO WHERE DOES STANDARD OUT GO WHEN I RUN ON THE GRID??? also is it a good idea to use abort??? Probably not!!" << std::endl;
        std::abort();
    }

    fpidResponse = fInputHandler->GetPIDResponse();


    //Event cuts
    TString cent_estimator = "V0A";
    double multPercentile = 0;

    fMultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(fMultSelection) multPercentile = fMultSelection->GetMultiplicityPercentile(cent_estimator.Data());
    else return;

    fMultDist->Fill(multPercentile);

    if(multPercentile < centLow || multPercentile > centHigh) return;

    AliVVertex *prim = fAOD->GetPrimaryVertex();
    int NcontV = prim->GetNContributors();
    if(NcontV < 3) return;

    double primZ = prim->GetZ();
    if(primZ < -10 || primZ > 10) return;


    int numTracks = fAOD->GetNumberOfTracks();

    std::vector<AliAODTrack*> unlikelyPiPlus_list;
    std::vector<AliAODTrack*> unlikelyPiMinus_list;
    std::vector<AliAODTrack*> piPlus_list;
    std::vector<AliAODTrack*> piMinus_list;
    std::vector<AliAODTrack*> trigger_list;
    std::vector<AliAODTrack*> associated_h_list;
    std::vector<AliAODTrack*> all_hadron_list;

    //Trigger list used for event mixing
    TObjArray* fMixedTrackObjArray = new TObjArray;
    fMixedTrackObjArray->SetOwner(kTRUE);

    for(int trackNum = 0; trackNum < numTracks; trackNum++) {
    

        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(trackNum));
        if(!track) continue;

        //List for comparison with cuts/filter bits
        all_hadron_list.push_back(track);

        //Filter for trigger particles
        bool trigFilter = track->TestBit(AliAODTrack::kIsHybridGCG);
        if(trigFilter) {

            if(track->Pt() > 4 && track->Pt() < 8 && TMath::Abs(track->Eta()) < 1) {
                trigger_list.push_back(track);
                AliCFParticle *triggerPart = new AliCFParticle(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0);
                fMixedTrackObjArray->Add(triggerPart);
           }
        }

        //Filter for associated particles

        bool assocFilter = track->TestFilterBit(AliAODTrack::kTrkGlobalNoDCA);
        if(assocFilter) {

            if(track->Pt() > 0.15 && TMath::Abs(track->Eta()) < 1) {

                associated_h_list.push_back(track);

                double TPCNSigmaPion = 1000;
                double TOFNSigmaPion = 1000;

                TPCNSigmaPion = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
                TOFNSigmaPion = fpidResponse->NumberOfSigmasTOF(track, AliPID::kPion);

                if(TOFNSigmaPion != 1000 && track->Charge() == 1) unlikelyPiPlus_list.push_back(track);
                if(TOFNSigmaPion != 1000 && track->Charge() != 1) unlikelyPiMinus_list.push_back(track);

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
    }

    //Making list of possible kaons (have to do +/- for proton or pi):

    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list;
    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list_signal_region;
    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list_signal_region_2_4;
    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list_RotatedPiPlus;
    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list_RotatedPiMinus;
    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list_RotatedPi;
    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list_Flipped;
    std::vector<AliAnalysisTaskKaonHadronRatio::AliMotherContainer> kaon_list_LS;


    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = 0; j < (int) piMinus_list.size(); j++) {
            AliMotherContainer kaon = DaughtersToMother(piPlus_list[i], piMinus_list[j], 0.1396, 0.1396);
            
            double piplus_dz[2];
            double piplus_covar[3];

            double piminus_dz[2];
            double piminus_covar[3];

            bool is_piplusDCA = piPlus_list[i]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., piplus_dz, piplus_covar);
            bool is_piminusDCA = piMinus_list[j]->PropagateToDCA(prim, fAOD->GetMagneticField(), 20., piminus_dz, piminus_covar);

            if(is_piplusDCA && is_piminusDCA) {
                fKaonDaughterDCA->Fill(piplus_dz[0], piminus_dz[0], kaon.particle.M());
            }

            AliMotherContainer kaon_RotatedPi = RotatedDaughtersToMother(piPlus_list[i], piMinus_list[j], 0.1396, 0.1396, TMath::Pi());
            AliMotherContainer kaon_Flipped = FlippedDaughtersToMother(piPlus_list[i], piMinus_list[j], 0.1396, 0.1396);
            kaon_list.push_back(kaon);
            kaon_list_RotatedPi.push_back(kaon_RotatedPi);
            kaon_list_Flipped.push_back(kaon_Flipped);

            AliMotherContainer kaon_Rotated;
            AliMotherContainer kaon_RotatedPiPlus;
            for(int k = 1; k < 12; k++) {
                kaon_Rotated = RotatedDaughtersToMother(piPlus_list[i], piMinus_list[j], 0.1396, 0.1396, (2*TMath::Pi()*k)/12);
                kaon_RotatedPiPlus = RotatedDaughtersToMother(piPlus_list[i], piMinus_list[j], 0.1396, 0.1396, (2*TMath::Pi()*k)/12);
                kaon_list_RotatedPiPlus.push_back(kaon_Rotated);
                kaon_list_RotatedPiMinus.push_back(kaon_RotatedPiPlus);

            }
        }
    }

    for(int i = 0; i < (int)piPlus_list.size(); i++) {
        for(int j = i + 1; j < (int) piPlus_list.size(); j++) {
            AliMotherContainer kaon = DaughtersToMother(piPlus_list[i], piPlus_list[j], 0.1396, 0.1396);
            kaon_list_LS.push_back(kaon);
        }
    }


    for(int i = 0; i < (int)kaon_list.size(); i++) {
        if(kaon_list[i].particle.M() < 0.520 && kaon_list[i].particle.M() > 0.480) {
            kaon_list_signal_region.push_back(kaon_list[i]);
            if(kaon_list[i].particle.Pt() < 4 && kaon_list[i].particle.Pt() > 2) {
                kaon_list_signal_region_2_4.push_back(kaon_list[i]);
            }
        }
    }



    // Filling all of our single particle distribution histograms:
    FillSingleParticleDist(trigger_list, primZ, fTriggerDist);
    FillSingleParticleDist(associated_h_list, primZ, fAssociatedHDist);
    FillSingleParticleDist(all_hadron_list, primZ, fLooseDist);
    FillSingleParticleDist(kaon_list_signal_region, primZ, fKaonDist);

    // Filling all of our correlation histograms
    MakeSameHKaonCorrelations(trigger_list, kaon_list, fDphiHKaon, primZ);
    MakeSameHKaonCorrelations(trigger_list, kaon_list_RotatedPi, fDphiHKaonRotatedPi, primZ);
    MakeSameHKaonCorrelations(trigger_list, kaon_list_Flipped, fDphiHKaonFlipped, primZ);
    MakeSameHKaonCorrelations(trigger_list, kaon_list_RotatedPiPlus, fDphiHKaonRotated, primZ);
    MakeSameHKaonCorrelations(trigger_list, kaon_list_RotatedPiMinus, fDphiHKaonRotatedPiMinus, primZ);
    MakeSameHHCorrelations(trigger_list, associated_h_list, fDphiHH, primZ);
    MakeSameTriggerTriggerCorrelations(trigger_list, fDphiTriggerTrigger, primZ);
    MakeSameHKaonCorrelations(trigger_list, kaon_list_LS, fDphiHKaonLS, primZ);


    fTriggersAndKaonsPerEvent_All->Fill(trigger_list.size(), kaon_list_signal_region.size());
    fTriggersAndKaonsPerEvent_2_4->Fill(trigger_list.size(), kaon_list_signal_region_2_4.size());


    if(kaon_list.size() > 0 && associated_h_list.size() > 0) {
        AliEventPool *fCorPool = fCorPoolMgr->GetEventPool(multPercentile, primZ);
        if(!fCorPool) {
            AliFatal(Form("No pool found for multiplicity = %f, zVtx = %f", multPercentile, primZ));
        }


        else {
            if(fCorPool->IsReady()) {
                MakeMixedHKaonCorrelations(fCorPool, kaon_list, fDphiHKaonMixed, primZ);
                MakeMixedHKaonCorrelations(fCorPool, kaon_list_LS, fDphiHKaonLSMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, associated_h_list, fDphiHHMixed, primZ);
                MakeMixedHHCorrelations(fCorPool, trigger_list, fDphiTriggerTriggerMixed, primZ);
            }
            if(fMixedTrackObjArray->GetEntries() > 0) {
                fCorPool->UpdatePool(fMixedTrackObjArray);
            }
        }
    }

    PostData(1, fOutputList);

}

void AliAnalysisTaskKaonHadronRatio::Terminate(Option_t *option)
{
}
