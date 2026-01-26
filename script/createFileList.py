import os
import argparse
import json

def name_simplification(folderName):
    name_short = {
        "DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": "DY2Mu",
        "DYto2Mu_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8": "DY2Mu_low",
        "DYto2Tau-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": "DY2Tau",
        "DYto2Tau_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8": "DY2Tau_low",

        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": "DY2L",
        "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": "DY2L_low",
    
        "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8": "ttbarDL",
        "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8": "ttbarSL",
        "TTto4Q_TuneCP5_13p6TeV_powheg-pythia8": "ttbar4Q",

        "TTBBto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8": "TTBB_DL",
        "TTBBtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8": "TTBB_SL",
        "TTBBto4Q_TuneCP5_13p6TeV_powheg-pythia8": "TTBB_4Q",
    
        "TTHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8": "TTHBB",
        "TTHtoNon2B_M-125_TuneCP5_13p6TeV_powheg-pythia8": "TTHnonBB",
        "TTH-Hto2B_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8": "TTHBB",
        "TTH-HtoNon2B_Par-M-125_TuneCP5_13p6TeV_powheg-pythia8": "TTHnonBB",
        "TTZH_TuneCP5_13p6TeV_madgraph-pythia8": "TTZH",
        "TTWH_TuneCP5_13p6TeV_madgraph-pythia8": "TTWH",
        "TTZH-ZHto4B_TuneCP5_13p6TeV_madgraph-pythia8": "TTZH",
        "THQ_ctcvcp_HIncl_M-125_4FS_TuneCP5_13p6TeV_madgraph-pythia8": "THQ",
        "THW_ctcvcp_HIncl_M-125_5FS_TuneCP5_13p6TeV_madgraph-pythia8": "THW",

        "WWW_4F_TuneCP5_13p6TeV_amcatnlo-madspin-pythia8": "WWW",
        "WWW-4F_TuneCP5_13p6TeV_amcatnlo-pythia8": "WWW",
        "WWZ-4F_TuneCP5_13p6TeV_amcatnlo-pythia8": "WWZ",
        "WZZ-5F_TuneCP5_13p6TeV_amcatnlo-pythia8": "WZZ",
        "ZZZ-5F_TuneCP5_13p6TeV_amcatnlo-pythia8": "ZZZ",
        "WZZ_TuneCP5_13p6TeV_amcatnlo-pythia8": "WZZ",
        "ZZZ_TuneCP5_13p6TeV_amcatnlo-pythia8": "ZZZ",

        "QCD_PT-15to20_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_15_20_mu",
        "QCD_PT-20to30_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_20_30_mu",
        "QCD_PT-30to50_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_30_50_mu",
        "QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_50_80_mu",
        "QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_80_120_mu",
        "QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_120_170_mu",
        "QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_170_300_mu",
        "QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_300_470_mu",
        "QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_470_600_mu",
        "QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_600_800_mu",
        "QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_800_1000_mu",
        "QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": "QCD_1000_mu",
        "QCD_Bin-PT-15to20_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_15_20_mu",
        "QCD_Bin-PT-20to30_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_20_30_mu",
        "QCD_Bin-PT-30to50_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_30_50_mu",
        "QCD_Bin-PT-50to80_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_50_80_mu",
        "QCD_Bin-PT-80to120_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_80_120_mu",
        "QCD_Bin-PT-120to170_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_120_170_mu",
        "QCD_Bin-PT-170to300_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_170_300_mu",
        "QCD_Bin-PT-300to470_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_300_470_mu",
        "QCD_Bin-PT-470to600_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_470_600_mu",
        "QCD_Bin-PT-600to800_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_600_800_mu",
        "QCD_Bin-PT-800to1000_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_800_1000_mu",
        "QCD_Bin-PT-1000_Fil-MuEnriched_TuneCP5_13p6TeV_pythia8": "QCD_1000_mu",

        "TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8": "TTZ_low",
        "TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8": "TTZ_high",
        "TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": "TTW",
        "TTTT_TuneCP5_13p6TeV_amcatnlo-pythia8": "TTTT",
        "TTLL_Bin-MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8": "TTZ_low",
        "TTLL_Bin-MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8": "TTZ_high",

        "WW_TuneCP5_13p6TeV_pythia8": "WW",
        "WZ_TuneCP5_13p6TeV_pythia8": "WZ",
        "ZZ_TuneCP5_13p6TeV_pythia8": "ZZ",
        "WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8": "WW2L2Nu",
        "WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8": "WZ2L2Q",
        "WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8": "WZ3l",
        "ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8": "ZZ2L2Nu",
        "ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8": "ZZ2L2Q",
        "ZZto4L_TuneCP5_13p6TeV_powheg-pythia8": "ZZ4l",

        "TTWW_TuneCP5_13p6TeV_madgraph-madspin-pythia8": "TTWW",
        "TTZZ_TuneCP5_13p6TeV_madgraph-madspin-pythia8": "TTZZ",
        "TTWW_TuneCP5_13p6TeV_madgraph-pythia8": "TTWW",
        "TTWZ_TuneCP5_13p6TeV_madgraph-pythia8": "TTWZ",

        "WtoLNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8": "WJets",
        "WtoLNu-4Jets_Bin-1J_TuneCP5_13p6TeV_madgraphMLM-pythia8": "WJet_1J",
        "WtoLNu-4Jets_Bin-2J_TuneCP5_13p6TeV_madgraphMLM-pythia8": "WJet_2J",
        "WtoLNu-4Jets_Bin-3J_TuneCP5_13p6TeV_madgraphMLM-pythia8": "WJet_3J",
        "WtoLNu-4Jets_Bin-4J_TuneCP5_13p6TeV_madgraphMLM-pythia8": "WJet_4J",

        "TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8": "TbarWp2L",
        "TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8": "TbarWp1L",
        "TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8": "TbarWp4Q",
        "TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8": "TWm2L",
        "TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8": "TWm1L",
        "TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8": "TWm4Q",
        "TZQB-Zto2L-4FS_Bin-MLL-30_TuneCP5_13p6TeV_amcatnlo-pythia8": "TZQB",

        "TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8": "TBbarQ",
        "TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8": "TbarBQ",
        "TBbarQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8": "TBbarQ2Q",
        "TBbarQtoLNu-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8": "TBbarQ1L",
        "TbarBQto2Q-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8": "TbarBQ2Q",
        "TbarBQtoLNu-t-channel-4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8": "TbarBQ1L",

        "TBbartoLplusNuBbar-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8": "TBbar1L",
        "TbarBtoLminusNuB-s-channel-4FS_TuneCP5_13p6TeV_amcatnlo-pythia8": "TbarB1L",
        "TBbartoLNu-s-channel_TuneCP5_13p6TeV_powheg-pythia8": "TBbar1L",
    }

    if folderName in name_short.keys():
        return name_short[folderName]
    else:
        return ""

# extracting to the last common folder
def find_single_chain_dir(start_dir):
    current = start_dir
    while True:
        if not os.path.isdir(current):
            return os.path.abspath(current)
        entries = list(os.scandir(current))
        subdirs = [e for e in entries if e.is_dir()]
        has_root_files = any(e.is_file() and e.name.endswith(".root") for e in entries)
        # in case inside only one folder and does not have root files inside, going on
        if len(subdirs) == 1 and not has_root_files:
            current = subdirs[0].path
            continue
        # when inside begins to have multiple folders or has root file, return and remain further level to files
        return os.path.abspath(current)

# getting everything under the common folder
def collect_root_files(base_dir):
    root_files = []
    for root, _, files in os.walk(base_dir):
        for fname in files:
            # whenever walk to root file, record the path after the base
            if fname.endswith(".root"):
                full_path = os.path.join(root, fname)
                rel_path = os.path.relpath(full_path, base_dir).replace(os.sep, "/")
                root_files.append(rel_path)
    return sorted(root_files)


def saveFiles(outName, fileDir):
    channels = [
        "Muon0",
        "Muon1",
        "Muon",

        "DY2Mu",
        "DY2Mu_low",

        "ttbarDL",
        "ttbarSL",
        "ttbar4Q",

        "TTBB_DL",
        "TTBB_SL",
        "TTBB_4Q",

        "QCD_120_170_mu",
        "QCD_20_30_mu",
        "QCD_470_600_mu",
        "QCD_800_1000_mu",
        "QCD_15_20_mu",
        "QCD_300_470_mu",
        "QCD_50_80_mu",
        "QCD_80_120_mu",
        "QCD_1000_mu",
        "QCD_170_300_mu",
        "QCD_30_50_mu",
        "QCD_600_800_mu",

        "TTHBB",
        "TTHnonBB",
        "TTWH",
        "TTZH",
        "THQ",
        "THW",

        "WJet_1J",
        "WJet_2J",
        "WJet_3J",
        "WJet_4J",

        "WW2L2Nu",
        "WZ2L2Q",
        "ZZ2L2Q",
        "ZZ2L2Nu",
        "WZ3l",
        "ZZ4l",

        "WWW",
        "WWZ",
        "WZZ",
        "ZZZ",

        "TbarWp2L",
        "TWm2L",
        "TbarWp1L",
        "TWm1L",
        "TZQB",
        
        "TTW",
        "TTZ_low",
        "TTZ_high",

        "TTTT"
    ]

    dirBlk=["{", "    \"dir\": {"]
    fileBlk=["    },", "    \"file\": {"]

    for folderName in os.listdir(fileDir):
        if folderName in channels:
            commonDir = find_single_chain_dir(os.path.join(fileDir, folderName))
            dirBlk.append("        \""+folderName+"\": \""+commonDir+"/\",")
            allFiles = collect_root_files(commonDir)
            if len(allFiles)>0:
                fileBlk.append("        \""+folderName+"\": [")
                for fName in allFiles:
                    fileBlk.append("            \""+fName+"\",")
                fileBlk[-1] = fileBlk[-1][:-1]
                fileBlk.append("        ],")
            else:
                fileBlk.append("        \""+folderName+"\": [],")
            continue
        symbol = name_simplification(folderName)
        if symbol in channels:
            commonDir = find_single_chain_dir(os.path.join(fileDir, folderName))
            dirBlk.append("        \""+symbol+"\": \""+commonDir+"/\",")
            allFiles = collect_root_files(commonDir)
            if len(allFiles)>0:
                fileBlk.append("        \""+symbol+"\": [")
                for fName in allFiles:
                    fileBlk.append("            \""+fName+"\",")
                fileBlk[-1] = fileBlk[-1][:-1]
                fileBlk.append("        ],")
            else:
                fileBlk.append("        \""+symbol+"\": [],")
    dirBlk[-1] = dirBlk[-1][:-1]
    fileBlk[-1] = fileBlk[-1][:-1]
    fileBlk.append("    }")
    fileBlk.append("}")

    with open(outName, "w") as f:
        for line in dirBlk:
            f.write(line+"\n")
        for line in fileBlk:
            f.write(line+"\n")

def build_args():
    parser = argparse.ArgumentParser(description="Create file list JSON.")
    parser.add_argument("output_json", help="Output JSON filename, e.g. RunIII2024Summer24NanoAODv15_raw.json")
    parser.add_argument("base_dir",help="Base directory that contains dataset folders, e.g. /data2/common/NanoAOD/mc/v15/RunIII2024Summer24NanoAODv15/")
    return parser.parse_args()

if __name__ == "__main__":
    args = build_args()
    saveFiles(args.output_json, args.base_dir)