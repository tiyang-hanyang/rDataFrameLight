#include <functional>
#include <stdexcept>
#include <string>
#include <map>
#include <memory>
#include <algorithm>

// combining passJetVetoMap
std::shared_ptr<const correction::CorrectionSet> JVM_oneMuon_cset;
std::shared_ptr<const correction::Correction> JVM_oneMuon_corr;

void JVM_oneMuon_init(const std::string& era)
{
    std::map<std::string, std::string> JVMFiles = {
        {"RunIII2024Summer24NanoAODv15", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},
        {"Run2024C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},
        {"Run2024D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},
        {"Run2024E", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},
        {"Run2024F", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},
        {"Run2024G", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},
        {"Run2024H", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},
        {"Run2024I", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/jetvetomaps.json"},

        {"Run3Summer23NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23CSep23-Summer23-NanoAODv12/jetvetomaps.json"},
        {"Run3Summer23BPixNanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/jetvetomaps.json"},
        {"Run2023C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23CSep23-Summer23-NanoAODv12/jetvetomaps.json"},
        {"Run2023D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/jetvetomaps.json"},

        {"Run3Summer22NanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22CDSep23-Summer22-NanoAODv12/jetvetomaps.json"},
        {"Run2022C", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22CDSep23-Summer22-NanoAODv12/jetvetomaps.json"},
        {"Run2022D", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22CDSep23-Summer22-NanoAODv12/jetvetomaps.json"},

        {"Run3Summer22EENanoAODv12", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jetvetomaps.json"},
        {"Run2022E", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jetvetomaps.json"},
        {"Run2022F", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jetvetomaps.json"},
        {"Run2022G", "/home/tiyang/public/rDataFrameLight_git/correction/POGCorr/POG/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/jetvetomaps.json"},
    };
    JVM_oneMuon_cset = correction::CorrectionSet::from_file(JVMFiles[era]);

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
    JVM_oneMuon_corr = JVM_oneMuon_cset->at(it->second);
}

ROOT::VecOps::RVec<int> passJetVetoFunc_oneMuon(const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi)
{
    ROOT::VecOps::RVec<int> pass_flags(eta.size());
    for (size_t i = 0; i < eta.size(); ++i) {
        float safe_eta = std::clamp(eta[i], -5.191f, 5.191f);
        float safe_phi = std::clamp(phi[i], -3.1415f, 3.1415f);
        pass_flags[i] = (JVM_oneMuon_corr->evaluate(std::vector<correction::Variable::Type>{"jetvetomap", safe_eta, safe_phi}) == 0);
    }
    return pass_flags;
}

ROOT::VecOps::RVec<float> minDistanceFromSingleMuon(const ROOT::VecOps::RVec<int>& isGoodJet, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const int& leadingMuonIdx, const ROOT::VecOps::RVec<float>& Muon_eta, const ROOT::VecOps::RVec<float>& Muon_phi)
{
    auto jetSize = isGoodJet.size();
    ROOT::VecOps::RVec<float> minDR(jetSize);
    auto leadingMuonEta = Muon_eta[leadingMuonIdx];
    auto leadingMuonPhi = Muon_phi[leadingMuonIdx];
    for (auto i = 0; i < jetSize; i++)
    {
        if (!isGoodJet[i])
        {
            minDR[i] = 0.0;
            continue;
        }
        minDR[i] = ROOT::VecOps::DeltaR(leadingMuonEta, eta[i], leadingMuonPhi, phi[i]);
    }
    return minDR;
}
