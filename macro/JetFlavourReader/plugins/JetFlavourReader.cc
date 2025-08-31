// -*- C++ -*-
//
// Package:    JetFlavourReader/JetFlavourReader
// Class:      JetFlavourReader
//
/**\class JetFlavourReader JetFlavourReader.cc JetFlavourReader/JetFlavourReader/plugins/JetFlavourReader.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:
//         Created:  Mon, 28 Jul 2025 01:41:04 GMT
//
//

// system include files
#include <memory>
#include <stack>
#include <fstream>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

// my information need
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

// class JetFlavourReader : public edm::one::EDAnalyzer<edm::one::SharedResources>
class JetFlavourReader : public edm::one::EDAnalyzer<>
{
public:
    explicit JetFlavourReader(const edm::ParameterSet &);
    ~JetFlavourReader() override;

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
    void beginJob() override;
    void analyze(const edm::Event &, const edm::EventSetup &) override;
    void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;
    edm::EDGetTokenT<std::vector<reco::GenJet>> genJetToken_;
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfoToken_;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;
    // vector<pat::Jet>                      "slimmedJets"               ""                "PAT"
    // vector<reco::GenJet>                  "slimmedGenJets"            ""                "PAT"
    // reco::JetFlavourInfoMatchingCollection   "slimmedGenJetsFlavourInfos"   ""                "PAT"

    // functions to assist information reading
    int findBSource(const reco::GenParticle *genBMeson);

    // for output
    void outputCSV(const std::vector<float> &event);
    std::ofstream outFile_;
    // edm::FileInPath outFileName_;
    const std::string defaultFileName_ = "events_output.csv";
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JetFlavourReader::JetFlavourReader(const edm::ParameterSet &iConfig)
    : jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("slimmedJets"))),
      genJetToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("slimmedGenJets"))),
      jetFlavourInfoToken_(consumes<reco::JetFlavourInfoMatchingCollection>(iConfig.getParameter<edm::InputTag>("slimmedGenJetsFlavourInfos"))),
      muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("slimmedMuons"))),
      vertexToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlineSlimmedPrimaryVertices")))
{
}

JetFlavourReader::~JetFlavourReader()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // please remove this method altogether if it would be left empty
    if (outFile_.is_open())
    {
        outFile_.close();
    }
}

//
// member functions
//
int JetFlavourReader::findBSource(const reco::GenParticle *genBMeson)
{
    // Define a stack to branch deep search
    std::stack<const reco::GenParticle *> particleStack;
    std::set<const reco::GenParticle *> visitedParticles;
    particleStack.push(genBMeson);

    while (!particleStack.empty())
    {
        // each time start from the stack top particle
        const reco::GenParticle *currentParticle = particleStack.top();
        particleStack.pop();
        if (currentParticle == nullptr)
            continue;
        if (visitedParticles.find(currentParticle) != visitedParticles.end())
            continue;
        visitedParticles.insert(currentParticle);

        // target: find Higgs/top/W/Z
        auto currentPID = currentParticle->pdgId();
        if (currentPID == 2212 || currentPID == -2212)
            continue;
        if (currentPID == 25 || currentPID == 24 || currentPID == -24 || currentPID == 23 || currentPID == 6 || currentPID == -6)
            return currentPID;

        // save each b-meson/b-quark in mothers on the stack
        auto numMothers = currentParticle->numberOfMothers();
        for (size_t index = 0; index < numMothers; index++)
        {
            const reco::GenParticle *motherPar = dynamic_cast<const reco::GenParticle *>(currentParticle->mother(index));
            if (motherPar == nullptr)
                continue;

            auto motherPID = motherPar->pdgId();
            if ((motherPID >= 500 && motherPID < 600) || (motherPID > -600 && motherPID <= -500) ||
                (motherPID >= 5000 && motherPID < 6000) || (motherPID > -6000 && motherPID <= -5000) ||
                (motherPID == 5) || (motherPID == -5) ||
                (motherPID == 25 || motherPID == 24 || motherPID == -24 || motherPID == 23 || motherPID == 6 || motherPID == -6))
            {
                particleStack.push(motherPar);
            }
        }
    }

    std::cout << "No B quark source found." << std::endl;
    return 0;
}

void JetFlavourReader::outputCSV(const std::vector<float> &saveInfo)
{
    for (size_t i = 0; i < saveInfo.size(); i++)
    {
        outFile_ << saveInfo[i];
        if (i != saveInfo.size() - 1)
        {
            outFile_ << ",";
        }
    }
    outFile_ << std::endl;
}

// ------------ method called for each event  ------------
void JetFlavourReader::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup)
{
    using namespace edm;

    int verbose = 0;

    // loading
    edm::Handle<std::vector<pat::Jet>> jets;
    iEvent.getByToken(jetToken_, jets);

    edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlaMap;
    iEvent.getByToken(jetFlavourInfoToken_, jetFlaMap);

    edm::Handle<std::vector<reco::GenJet>> genjets;
    iEvent.getByToken(genJetToken_, genjets);

    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonToken_, muons);

    edm::Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(vertexToken_, vertices);

    // if jet container not pass, not use the events
    if (!jets.isValid())
        return;

    // add muon selection condition
    int passDimuon(0);
    int passSameSign(0);
    if (vertices.isValid() && !vertices->empty())
    {
        auto primaryVtx = vertices->front();
        std::vector<float> goodMuonChargeList;
        if (muons.isValid())
        {
            for (const auto &muon : *muons)
            {
                if (muon.pt() < 15)
                    break;
                if (muon.eta() > 2.4 || muon.eta() < -2.4)
                    continue;
                if (!muon.isTightMuon(primaryVtx))
                    continue;
                if (muon.PFIsoTight == 0)
                    continue;
                goodMuonChargeList.push_back(muon.charge());
            }
            passDimuon = goodMuonChargeList.size();
            if (goodMuonChargeList.size() == 2)
            {
                passSameSign = (goodMuonChargeList[0] * goodMuonChargeList[1] > 0.0);
            }
        }
    }
    if (verbose)
        std::cout << "pass Dimuon condition: " << passDimuon << ", pass SS condition: " << passSameSign << std::endl;

    // extract genjet to hadron mapping
    std::map<unsigned, reco::JetFlavourInfo> flavMap;
    if (jetFlaMap.isValid())
    {
        for (auto const &pair : *jetFlaMap)
        {
            edm::RefToBase<reco::Jet> ref = pair.first;
            flavMap[ref.key()] = pair.second;
        }
        if (flavMap.size() < jetFlaMap->size())
        {
            std::cout << "[warning] inside reco::JetFlavourInfo, the iterator entry size " << flavMap.size() << " is smaller than container entry size " << jetFlaMap->size();
        }
    }

    std::vector<float> bjetPt;
    std::vector<float> bjetEta;
    std::vector<float> bjetPhi;
    std::vector<int> bjetSourceLabel;
    for (const auto &jet : *jets)
    {
        // select b-jet condition
        if (jet.pt() < 20)
            break;
        if (jet.eta() > 2.5 || jet.eta() < -2.5)
            continue;
        if (!(jet.neutralHadronEnergyFraction() < 0.99 &&
              jet.neutralEmEnergyFraction() < 0.90 &&
              jet.nConstituents() > 1 &&
              jet.chargedHadronEnergyFraction() > 0.01 &&
              jet.chargedMultiplicity() > 0))
            continue;
        // 2023 threshold 0.1917, 2023BPix threshold 0.1919
        auto pnetBtagValue = jet.bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralDiscriminatorsJetTags:BvsAll");
        if (pnetBtagValue < 0.1917)
            continue;
        if (verbose)
            std::cout << "pass bjet selection" << std::endl;

        // record reco information
        bjetPt.push_back(jet.pt());
        bjetEta.push_back(jet.eta());
        bjetPhi.push_back(jet.phi());

        // source the jet
        auto genjet = jet.genJet();
        int sourceType(0);
        if (genjet != nullptr)
        {
            // record genjet orientation for multi matching
            float genjetEta = genjet->eta();
            float genjetPhi = genjet->phi();

            // read gen jet flavour info if possible
            edm::Ref<std::vector<reco::GenJet>> genjetRef(genjets, genjet - &genjets->front());
            auto genjetIdx = genjetRef.key();
            if (flavMap.find(genjetIdx) != flavMap.end())
            {
                auto genjetflav = flavMap.at(genjetIdx);

                // if genjet is also bjet, then try to find the source
                if (genjetflav.getHadronFlavour() == 5)
                {
                    auto bRefs = genjetflav.getbHadrons();
                    // if more than 1 bHadron, use delta R selection
                    std::vector<int> sourceTypes;
                    std::vector<float> sourceEtas;
                    std::vector<float> sourcePhis;
                    float smallestDR(999);
                    for (const auto &bref : bRefs)
                    {
                        if (!bref.isAvailable())
                            continue;
                        auto bPart = bref.get();
                        sourceTypes.push_back(findBSource(bPart));
                        sourceEtas.push_back(bPart->eta());
                        sourcePhis.push_back(bPart->phi());
                    }

                    // determine the best matching
                    if (sourceTypes.size() == 1)
                        sourceType = sourceTypes[0];
                    else
                    {
                        int tempType(0);
                        for (size_t srcIdx = 0; srcIdx < sourceTypes.size(); srcIdx++)
                        {
                            if (sourceTypes[srcIdx] == 0)
                                continue;
                            float tempDR = std::pow(std::pow((sourceEtas[srcIdx] - genjetEta), 2.0) + std::pow((sourcePhis[srcIdx] - genjetPhi), 2.0), 0.5);
                            if (tempDR < smallestDR)
                            {
                                smallestDR = tempDR;
                                tempType = sourceTypes[srcIdx];
                            }
                        }
                        sourceType = tempType;
                    }
                }
            }
        }
        bjetSourceLabel.push_back(sourceType);
    }
    int bjetMultiplicity = bjetPt.size();

    if (bjetMultiplicity < 4)
    {
        if (verbose)
            std::cout << "Fail 4 bjet requirement" << std::endl;
        return;
    }

    if (verbose)
    {
        std::cout << "pass " << bjetPt.size() << " b-jets" << std::endl;
        for (size_t bjetIdx = 0; bjetIdx < bjetPt.size(); bjetIdx++)
            std::cout << "\t bjet id: " << bjetIdx << ", pt: " << bjetPt[bjetIdx] << ", eta: " << bjetEta[bjetIdx] << ", phi: " << bjetPhi[bjetIdx] << ", truth source label: " << bjetSourceLabel[bjetIdx] << std::endl;
    }

    // compute dR among b-jets
    std::vector<float> leadingBjetsDR;
    std::vector<std::pair<int, int>> pairingOrder = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (const auto &bPair : pairingOrder)
    {
        leadingBjetsDR.push_back(std::pow(std::pow((bjetEta[bPair.first] - bjetEta[bPair.second]), 2.0) + std::pow((bjetPhi[bPair.first] - bjetPhi[bPair.second]), 2.0), 0.5));
    }
    if (verbose)
    {
        std::cout << "b-jet dR:" << std::endl;
        for (size_t i = 0; i < leadingBjetsDR.size(); i++)
        {
            std::cout << "\t dR (" << pairingOrder[i].first << ", " << pairingOrder[i].second << "): " << leadingBjetsDR[i] << std::endl;
        }
    }

    // converting the b matching infor into a label
    int bjetAssignmentLabel = -1;
    if (bjetSourceLabel[0] == 25 && bjetSourceLabel[1] == 25)
        bjetAssignmentLabel = 0;
    else if (bjetSourceLabel[0] == 25 && bjetSourceLabel[2] == 25)
        bjetAssignmentLabel = 1;
    else if (bjetSourceLabel[0] == 25 && bjetSourceLabel[3] == 25)
        bjetAssignmentLabel = 2;
    else if (bjetSourceLabel[1] == 25 && bjetSourceLabel[2] == 25)
        bjetAssignmentLabel = 3;
    else if (bjetSourceLabel[1] == 25 && bjetSourceLabel[3] == 25)
        bjetAssignmentLabel = 4;
    else if (bjetSourceLabel[2] == 25 && bjetSourceLabel[3] == 25)
        bjetAssignmentLabel = 5;

    // final output
    std::vector<float> storedInfo;
    storedInfo.push_back(bjetMultiplicity);
    storedInfo.push_back(bjetPt[0]);
    storedInfo.push_back(bjetEta[0]);
    storedInfo.push_back(bjetPhi[0]);
    storedInfo.push_back(bjetPt[1]);
    storedInfo.push_back(bjetEta[1]);
    storedInfo.push_back(bjetPhi[1]);
    storedInfo.push_back(bjetPt[2]);
    storedInfo.push_back(bjetEta[2]);
    storedInfo.push_back(bjetPhi[2]);
    storedInfo.push_back(bjetPt[3]);
    storedInfo.push_back(bjetEta[3]);
    storedInfo.push_back(bjetPhi[3]);
    storedInfo.insert(storedInfo.end(), leadingBjetsDR.begin(), leadingBjetsDR.end());
    storedInfo.push_back(bjetAssignmentLabel);
    storedInfo.push_back(passDimuon);
    storedInfo.push_back(passSameSign);

    outputCSV(storedInfo);
}

// // ------------ method called once each job just before starting event loop  ------------
void JetFlavourReader::beginJob()
{
    outFile_.open(defaultFileName_, std::ios::out);
    if (!outFile_.is_open())
    {
        std::cerr << "Failed to open file for writing: events_output.csv" << std::endl;
        return;
    }

    outFile_ << "Nbjet, b1pt, b1eta, b1phi, b2pt, b2eta, b2phi, b3pt, b3eta, b3phi, b4pt, b4eta, b4phi, b12dr, b13dr, b14dr, b23dr, b24dr, b34dr, assignmentlabel, passdimuon, passSS" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void JetFlavourReader::endJob()
{
    if (outFile_.is_open())
    {
        outFile_.close();
    }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void JetFlavourReader::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
    // The following says we do not know what parameters are allowed so do no validation
    //  Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    // Specify that only 'tracks' is allowed
    // To use, remove the default given above and uncomment below
    // edm::ParameterSetDescription desc;
    // desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
    // descriptions.addWithDefaultLabel(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(JetFlavourReader);
