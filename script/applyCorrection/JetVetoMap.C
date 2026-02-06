#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

// combining passJetVetoMap
std::shared_ptr<const correction::CorrectionSet> JVM_cset;
std::shared_ptr<const correction::Correction> JVM_corr;

void JVM_init(const std::string& era)
{
    std::map<std::string, std::string> JVMFiles = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},
        {"Run2024C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},
        {"Run2024D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},
        {"Run2024E", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},
        {"Run2024F", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},
        {"Run2024G", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},
        {"Run2024H", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},
        {"Run2024I", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2024_Summer24/jetvetomaps.json"},

        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23/jetvetomaps.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23BPix/jetvetomaps.json"},
        {"Run2023C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23/jetvetomaps.json"},
        {"Run2023D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2023_Summer23BPix/jetvetomaps.json"},

        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2022_Summer22/jetvetomaps.json"},
        {"Run2022C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2022_Summer22/jetvetomaps.json"},
        {"Run2022D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2022_Summer22/jetvetomaps.json"},

        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2022_Summer22EE/jetvetomaps.json"},
        {"Run2022E", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2022_Summer22EE/jetvetomaps.json"},
        {"Run2022F", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2022_Summer22EE/jetvetomaps.json"},
        {"Run2022G", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/2022_Summer22EE/jetvetomaps.json"},
    };
    JVM_cset = correction::CorrectionSet::from_file(JVMFiles[era]);

    std::map<std::string, std::string> JVMTab = {
        {"RunIII2024Summer24NanoAODv15", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024C", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024D", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024E", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024F", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024G", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024H", "Summer24Prompt24_RunBCDEFGHI_V1"},
        {"Run2024I", "Summer24Prompt24_RunBCDEFGHI_V1"},
    
        {"Run3Summer23NanoAODv12", "Summer23Prompt23_RunC_V1"},
        {"Run3Summer23BPixNanoAODv12", "Summer23BPixPrompt23_RunD_V1"},
        {"Run2023C", "Summer23Prompt23_RunC_V1"},
        {"Run2023D", "Summer23BPixPrompt23_RunD_V1"},


        {"Run3Summer22NanoAODv12", "Summer22_23Sep2023_RunCD_V1"},
        {"Run2022C", "Summer22_23Sep2023_RunCD_V1"},
        {"Run2022D", "Summer22_23Sep2023_RunCD_V1"},

        {"Run3Summer22EENanoAODv12", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run2022E", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run2022F", "Summer22EE_23Sep2023_RunEFG_V1"},
        {"Run2022G", "Summer22EE_23Sep2023_RunEFG_V1"},
    };
    auto it = JVMTab.find(era);
    if (it == JVMTab.end()) {
        throw std::runtime_error("Unknown era: " + era);
    }
    JVM_corr = JVM_cset->at(it->second);
}

ROOT::VecOps::RVec<int> passJetVetoFunc(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi) 
{
    ROOT::VecOps::RVec<int> pass_flags(eta.size());
    for (size_t i = 0; i < eta.size(); ++i) {
        // exclude non coverage etas
        if (abs(eta[i])>2.5) {
            pass_flags[i]=0;
            continue;
        }
        float safe_phi = std::clamp(phi[i], -3.1415f, 3.1415f);
        pass_flags[i] = (JVM_corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", eta[i], safe_phi}) == 0);
    }
    return pass_flags;
}

ROOT::VecOps::RVec<float> minDistanceFromMuon(const ROOT::VecOps::RVec<int>& isGoodJet, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const int& leadingMuonIdx, const int& subleadingMuonIdx, const ROOT::VecOps::RVec<float>& Muon_eta, const ROOT::VecOps::RVec<float>& Muon_phi) 
{
    auto jetSize = isGoodJet.size();
    ROOT::VecOps::RVec<float> minDR(jetSize);
    auto leadingMuonEta = Muon_eta[leadingMuonIdx];
    auto leadingMuonPhi = Muon_phi[leadingMuonIdx];
    auto subleadingMuonEta = Muon_eta[subleadingMuonIdx];
    auto subleadingMuonPhi = Muon_phi[subleadingMuonIdx];
    for (auto i=0; i < jetSize; i++)
    {
        // only need to compute for the goodjet
        if (!isGoodJet[i]) 
        {
            minDR[i]=0.0;
            continue;
        }
        float dRleading = ROOT::VecOps::DeltaR(leadingMuonEta, eta[i], leadingMuonPhi, phi[i]);
        float dRsubleading = ROOT::VecOps::DeltaR(subleadingMuonEta, eta[i], subleadingMuonPhi, phi[i]);
        minDR[i] = std::min(dRleading, dRsubleading);
    }
    return minDR;
}