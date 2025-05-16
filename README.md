# rDataFrameLight

This package is based on the RDataFrame, providing an easy-to-operate analysis code containing sample skimming, variable defining, histogram plotting, graph making functions with json configs.
The code only need the c++ 17 and ROOT installed for operating, CMSSW would be optional. The json parsing uses the header from nlohmann: https://github.com/nlohmann/json.

## setup

The minimum practise of getting things ready

```
mkdir source run build
cd source/
git clone \<your forked code\> .
mkdir external
cd external/
wget https://raw.githubusercontent.com/nlohmann/json/develop/single_include/nlohmann/json.hpp
cd ../../build/
cmake ../source
make
source setup.sh
```

In case of dealing with 2024 data, the CMSSW would be needed to avoid large memory leakage due to the NANOAOD structure.
The release setup can be done as follows on a machine with the cvmfs installed. Note that if the compile has been done previously, the cache need to be cleaned first.
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_14_2_2
cd CMSSW_14_2_2
cmsenv
scram b
cd \<Your build directory\>
rm CMakeCache.txt
cmake ../source/
make
source setup.sh
```

## preparing the correctionlib
```
cd ../
mkdir correction/
cd correction
# downloading Rochester correction
git clone https://gitlab.cern.ch/cms-muonPOG/muonscarekit.git RochesterCorr
# downloading PU, jetvetomap, Muon SF corrections
git clone https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration.git POGCorr
cd POGCorr/POG/JME/2023_Summer23/
gunzip jetvetomaps.json.gz
cd ../2023_Summer23BPix/
gunzip jetvetomaps.json.gz
cd ../../LUM/2023_Summer23/
gunzip puWeights.json.gz
cd ../2023_Summer23BPix
gunzip puWeights.json.gz
cd ../../MUO/2023_Summer23/
gunzip muon_Z.json.gz
cd ../2023_Summer23BPix
gunzip muon_Z.json.gz
cd ../../../../
# downloading Drell-Yan Z pt corrections
git clone https://gitlab.cern.ch/cms-higgs-leprare/hleprare.git DYZptCorr
cd DYZptCorr/DYweightCorrlib/
gunzip DY_pTll_weights_2023preBPix_v2.json.gz
gunzip DY_pTll_weights_2023postBPix_v2.json.gz
cd ../../../run/
```
The names of the CMS correction packages are very important here for the correction scripts to find the files. Do not change them here.
For the POG json correctionlib integration, the json files are all archived, here shows the method to unzip the 2023 files. If want to do more period, one can unzip corresponding files as well.

## Running Analysis

### skimSamples

The executable to do the skimming jobs. A job options should be created in json format, with reference to json/DiMuonSkim.json as an example. The needed branches, samples, cuts should be also created in the json format for your job option to read. There are also examples in correspondng json/ folder. 
The output would be a folder containing skimmed samples.
Currently, only backup processing the skimming of local samples. The grid run extension would be added in the future.

### collectHists

The executable to draw histogram from samples using RDataFrame API. The job option should be created in json format, with reference to json/MuonHistCreate.json.
The output would be a folder containning root files storing the histograms.
Currently, only backup processing the histogram extraction from local samples. The grid run extension would be added in the future. If your variables does not exist in the NANOAOD form, you can add them by adding a `cut` config just defining it.

### plotHists

The executable to plot the histogram received by the colletHists step. The job option should be created in json format, with reference to json/MuonHistPlot.json.
The output would be a folder of figures of `.png` `.pdf` and `.eps` copies. 
*To fix the folder storage soon.*
The plotting code currently has a segfault issue with releasing the stack hist, please for now only put one variable in the job option each run.
*To fix the segfault soon.*