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

The executable to do the skimming jobs. 
A job options should be created in json format and given to the skimSamples utility. Either absolute path could be provided or the directory before json could be omitted. An example is provided as follows:
```
skimSamples json/skim/Dimuon_medium_RunIII2024Summer24NanoAODv15_noRochester.json
```
The output will be one or several root files, including a Events TTree and a TH1D histogram named genWeightSum. Any non-raw sample file should have a genWeightSum histogram indicating the total genWeight in the corresponding file without any cuts.

In the skim config, the user needs to provide:
- `preliminary`: boolean, default value 1 for the initial skim from raw sample. For the preliminary skim, the code will collect genWeight from Events tree. Else, the code will sum up from the input file's genWeightSum TH1D. If not the first skimming, must indicate!
- `merge`: specifying the running mode
- `year`: necessary for data skimming as this specifies the path of the GoldenJson.
- `era`: for building up the IO filesystem control, just write down the correct one would be fine.
- `isData`: boolean, 1 means data skimming. Data need to take the Goldenjson at preliminary run.
- `datasets`: all the datasets you want to process as a list.
- `outDir`: specifying the output dir, the internal era and channel structure will be automatically handled.
- `sampleConfig`: specifying the json files of the input samples.
- `cutConfig`: specifying the json files of the selections in a list.
- `branchConfig`: specifying the branches you want to keep. The branches inside takes intersection with the branches exists in the input data file.

optional configs:
- `name`: For naming the output files for the in the merging running mode. Default value to be "{era}_skimmed". The non-merging running mode does not need this.
- `maxFiles`: For setting a limit of the numbers of files in each
- `jobType`: For safety check. Better to always turn on avoid running on wrong configs.
- `comments`: or any other names not specified here, could write anything for the comments as json does not backup comments grammar.

This running calls several additional json configs that the user need to write in advance. Here is the explanation of each json config:

#### sampleConfig
The json file providing the files to be used. The `json/samples/rawNanoAOD.json` could be referred to as an example. 
The `dir` block is a map of dataset name and the absolute path. The dataset name must be exactly the same as the one used in the skim config to be understood.
The `file` block is a map of dataset name and the list of sample file paths. The concatenation of the two part should indicate the absolute path of the file. Now only able to read local path, cannot using xRootD.

The `SampleConfig.h` is designed to parse it. With this auxiliary class, the user could get the list of files for any usage. The auto-filter will identify the channel and files not available to avoid crash.

#### cutConfig
The json file providing the definition of variables and doing cuts. The `json/cuts/MetricMuon/MuonMetricSkim_strict_2024.json` could be referred to as an example. 
The files is a list of operations that will be applied on rDataFrame in order. The default functions include:
- define: ["define", branchName, formula]. This is the config used for defining variables to save in addition or for the further usage. 
- cut: ["cut", comments, formula]. This is the config used for doing Filter.

Notice that due to the limitation of string format config in json to avoid re-compiling, the functionality of this type of config is not very efficient.
If the selection and definition need to use complicated definition as lambda function, one need to append to the gInterpreter lines in `Root/CutControl.cxx` and re-compile the code. 
I am considering how to make this usage more convenient.

#### branchConfig
The json file providing the list of branches to save. The `json/branches/MuonMetricBranchNew.json` is used as an example. There is the automatic intersection with available branches inside the root files, so one could write any branches ever need in such configs.

----------

### Di-muon skimming

Currently, the default skimming contains the following requirements:

<details>
<summary> Trigger menu </summary>
<pre>
HLT_IsoMu24
HLT_IsoMu24_eta2p1
HLT_IsoMu27
HLT_Mu50
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8
HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8
</pre>
</details>

<details>
<summary> Noise filter </summary>
<pre>
Flag_goodVertices
Flag_globalSuperTightHalo2016Filter
Flag_EcalDeadCellTriggerPrimitiveFilter
Flag_BadPFMuonFilter
Flag_eeBadScFilter
Flag_BadPFMuonDzFilter
Flag_hfNoisyHitsFilter
</pre>
</details>

Di-muon selections:
- Has primary verices (PV_npvsGood>0)
- Exactly two good muons
    <details>
    <summary> Good Muon selections </summary>
    <pre>
    Muon_pt > 15, abs(Muon_eta) < 2.4
    abs(Muon_dxy) < 0.05,  abs(Muon_dz) < 0.1, abs(Muon_sip3d) < 8
    Muon_mediumId
    Muon_miniPFRelIso_all < 0.4
    Muon_jetDF fail medium b-tag working points (depending on years)
    Muon_promptMVA > 0.64
    </pre>
    </details>
- leading muon with pT larger than 20 GeV to ensure the trigger fully effective.

----------

### sample corrections

The `script/applyCorrection/` folder contains all the corrections. 
Each correction contains a python file and a C macro. The python instead of cpp UI is due to the interface of the corrections need the `correctionlib` reading and the API are highly variational. The C macro is since the correction definition is usually highly complicated.
The corrections applied are explain one by one here:

#### JEC and JVM correction

This correction changes the Jet_pt value, thus should be apply first. The corrected Jet_pt is named as `Jet_pt_JEC`.
The correction refer to the following files:
```
script/applyCorrection/
|--JEC_Data.py
│--JEC_Data.C
--------------------
|--JEC_MC.py
|--JEC_MC.C
|--jet_smear.json
--------------------
|--JetVeoMap.py
|--JetVeoMap.C
```

In data, we need to do the JES based on the L1FastJet, L2Relative, L2L3Residual corrections. In Run 3, the L1 correction is dummy so we only do the L2 and L2L3 corrections for data.
The correction library currently used are specified inside `JEC_Data.C` as follows:

<details>
<summary> JEC data </summary>

<pre>
common:
    pt -> pt * (1-rawFactor)

2024:
    pt -> pt * Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi(eta, phi, pt)
    pt -> pt * Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi(run, eta, pt)

2023BPix:
    pt -> pt * Summer23BPixPrompt23_V3_DATA_L2Relative_AK4PFPuppi(eta, phi, pt)
    pt -> pt * Summer23BPixPrompt23_V3_DATA_L2L3Residual_AK4PFPuppi(run, eta, pt)

2023:
    pt -> pt * Summer23Prompt23_V2_DATA_L2Relative_AK4PFPuppi(eta, pt)
    pt -> pt * Summer23Prompt23_V2_DATA_L2L3Residual_AK4PFPuppi(run, eta, pt)

2022:
    pt -> pt * Summer22_22Sep2023_RunCD_V3_DATA_L2Relative_AK4PFPuppi(eta, pt)
    pt -> pt * Summer22_22Sep2023_RunCD_V3_DATA_L2L3Residual_AK4PFPuppi(eta, pt)

2022EE:
    pt -> pt * Summer22EE_22Sep2023_Run{E|F|G}_V3_DATA_L2Relative_AK4PFPuppi(eta, pt)
    pt -> pt * Summer22EE_22Sep2023_Run{E|F|G}_V3_DATA_L2L3Residual_AK4PFPuppi(eta, pt)
</pre>
</details>

In mc, we does not do the residual correction thus only L2 for JES. The post-JES step has the JER which need the correctionlib as well as the json_smear.json config taken from https://gitlab.cern.ch/cms-analysis/jme/jerc-application-tutorial/-/blob/master/jer_smear.json.gz.

<details>
<summary> JEC MC </summary>

<pre>
common:
    pt -> pt * (1-rawFactor)

2024:
    pt -> pt * Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi(eta, phi, pt)
    resolution = Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi(eta, pt, rho)
    scale = Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi(eta, pt)


2023BPix:
    pt -> pt * Summer23BPixPrompt23_V3_MC_L2Relative_AK4PFPuppi(eta, phi, pt)
    resolution = Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi(eta, pt, rho)
    scale = Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi(eta, pt)


2023:
    pt -> pt * Summer22_22Sep2023_V3_MC_L2Relative_AK4PFPuppi(eta, pt)
    resolution = Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi(eta, pt, rho)
    scale = Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi(eta, pt)


2022:
    pt -> pt * Summer23Prompt23_V2_MC_L2Relative_AK4PFPuppi(eta, pt)
    resolution = Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi(eta, pt, rho)
    scale = Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi(eta, pt)


2022EE:
    pt -> pt * Summer22EE_22Sep2023_V3_MC_L2Relative_AK4PFPuppi(eta, pt)
    resolution = Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi(eta, pt, rho)
    scale = Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi(eta, pt)

</pre>
</details>

After the jet pt correction. the good Jet is defined with the requirements
- Jet_pt_JEC > 30, abs(Jet_eta) < 2.5
- Jet_jetId == 6 (tight (lep-veto) id)
    <details>
    <summary> Jet tight (lep-veto) Id in 2024</summary>
    <pre>
    Jet_neHEF<0.99
    Jet_neEmEF<0.90
    Jet_nConstituents>1
    Jet_muEF<0.80
    Jet_chHEF>0.01
    Jet_chMultiplicity>0
    Jet_chEmEF<0.80
    </pre>
    </details>
- Jet_neEmEF + Jet_chEmEF < 0.9
- Jet_rawFactor < 0.9 (avoiding problematic jets)
- DR > 0.4 from selected muons

Jet-veto Map

The loose jet is defined with relaxation in the threshold `Jet_pt_JEC > 15 GeV` and relaxation of the DR from muons.
If any of the loose jets in the event falls in the exclusion region, drop the events.

<details>
<summary> Jet-veto map exclusion </summary>

<pre>
2024:
    Summer24Prompt24_RunBCDEFGHI_V1(eta,phi)

2023:
    Summer23Prompt23_RunC_V1(eta,phi)

2023BPix
    Summer23BPixPrompt23_RunD_V1(eta,phi)
    
2022:
    Summer22_23Sep2023_RunCD_V1(eta,phi)

2022EE:
    Summer22EE_23Sep2023_RunEFG_V1(eta,phi)

</pre>
</details>

#### Pile-up reweighting

The pile-up correction is from the LUM POG, demonstrated using the files:
```
script/applyCorrection/
|--PUWeight.py
│--PUWeight.C
```

The correction of the pile-up generates one additional weights `PUWeight` for the MC corrections.

<details>
<summary> PU correction correctionlib </summary>

<pre>
RunIII2024Summer24NanoAODv15
    Collisions24_BCDEFGHI_goldenJSON

Run3Summer23NanoAODv12:
    Collisions2023_366403_369802_eraBC_GoldenJson(Pileup_nTrueInt)

Run3Summer23BPixNanoAODv12
    Collisions2023_369803_370790_eraD_GoldenJson(Pileup_nTrueInt)

Run3Summer22NanoAODv12
    Collisions2022_355100_357900_eraBCD_GoldenJson(Pileup_nTrueInt)

Run3Summer22EENanoAODv12
    Collisions2022_359022_362760_eraEFG_GoldenJson(Pileup_nTrueInt)

</pre>
</details>

#### ttbar overlapping removal

The TT samples has the contribution with additional b-jets other than the two from top pair decay during the parton shower. This has the phase space overlapping with TTBB sample contribution. The selection is needed to identify the TT events with 2 additional b-jets and exclude them to avoid double counting.
Since the generator level decay chain tracing is complicated, we deal it with the files:
```
script/applyCorrection/
|--ttbar_overlapping_removal.py
│--ttbar_overlapping_removal.C
```

In these scripts, the last b-hadrons are sourced to top quark, heavy bosons, and not from these (from light flavour QCD shower) using the GenPart branches in the NANOAOD level.
The GenJet with ghost matching to b-hadrons (GenJet_hadronFlavour == 5) are DR matched to these b-hadrons and assigned with the source from top or additional. No selections are applied to these GenJet as this overlapping removal only for remove the fiducial phase space double count.
Events with 2 additional GenJet ghost matching to b-hadron and not from top quark decay are removed.

#### b-tag correction

The b-tagged jets need the correction as well. The first thing needed is the b-tag efficiency collection for each of the channel using `utils/computeBEff.cxx`. To do the correction, the user need to create a new folder and run with specified era. For example:
```
computeBEff RunIII2024Summer24NanoAODv15
```

The output is in a correctionlib json format, which can be read in the same way as the b-tag scale factors from the BTV POG. The correction using files:
```
script/applyCorrection/
|--BTagCorr.py
│--BTagCorr.C
```

Notice that currently the b-tag correction for the 2024 only has the fixed-working point scale factors provided for UParTAK4 tagging algorithm.
The other era using v12 samples have the recommendation for PNetAK4 tagging and the corresponding shape scale factors.
The correction needs as well the tagged efficiency on the actual channel, thus the user need to specify on the `L44` of `BTagCorr.C` if the b-tagged efficiency computed from `computeBEff` is stored in their own path:
```
44    std::string btagEffDir = "/home/tiyang/public/rDataFrameLight_git/SR_medium_muon/btag_eff_"+era+"/";
45    BTagCorr_cseteff = correction::CorrectionSet::from_file(btagEffDir+"/"+channel+"_btag_eff.json");
46    BTagCorr_ceff =  BTagCorr_cseteff->at("UParTAK4_eff_values");
``` 

----------

### Signal selections

The Signal selection in this study is also based on the skimming methods. 
The step of selections are listed as follows:
- S1. Requiring exactly two good muons, as the skimming introduced in the above selection explains. 
- S2. Requiring at least four good jets.
- S3. Requiring at least three b-tagged jets.
- S4. Requiring two muons of the same-sign charges.

The corresponding cuts are stored in 
`json/cuts/MetricMuon/`

----------

### collectHists

After getting the skimmed, corrected and selected events, the histogram serving codes using RDataFrame API are provided in this framework. Using json configs, the histograms could be created easily.

The first step is to collect the histograms as desired using files
```
Root/
|--HistControl.cxx
rDataFrameLight/
|--HistControl.h
utils/
│--collectHists.cxx
```
The running command example is like:
```
collectHists json/MetricMuonHists/medium_dimuon/S1_dimuon_create_muonVar_RunIII2024Summer24NanoAODv15.json
```

In the json config for the histogram creation, the user needs to provide:
- `varName`: Necessary to provide, the name of the variables to plots. Could use variables not defined yet, as long as providing definition in additional cut.
- `era`: Used for extracting the luminosity and forming the output nomenclature .
- `outDir`: The output directory.
- `datasets`: The list of all channels you want to plot, related to the following two configs as well.
- `isData`: To specify if a channel is data or mc. In case of mc there will be the need of the weight applied.
- `needMerge`: A dictionary from the channel string to a list of sub-channels. The term `datasets` could include the composite channels. like merging the different pt ranges' "QCD" together for simplicity. In such case, the user need to specify here or the code would fail to find such channel.
- `MCWeight`: The list to contain all the corrections weights for the mc samples. Only needed by the MC plotting.

- `sampleConfig`: The necessary configs of the skimmed samples where the histograms are drawn from.
- `cutConfig`: A list of all the configs of RDF level Define and Filter needed for the histogram plotting. For the variables need computation from the original branches, those variables must be defined in one of the cut config provided. Could be empty if only need to plot an existing branch without any further cuts.
- `varConfig`: The variable configs. Different from cutConfig, here for all the variables to plot, the binning and axis name properties need to be pre-defined in such configs.
- `lumiConfig`: The configs of the integrated luminosity of different eras and campaigns. The default value is `json/Lumi/Run3.json`. If need to run on Run2, this entry will be necessary.
- `XSConfig`: The XS configs of the channels. The default value is `json/XS/Run3.json`. If need to run on Run2, this entry will be necessary.

optional configs:
- `jobType`: For safety checking, also recommended to always writing as other collectHists configs do.
- `comments`: Any comments to write.

For the configs needed in this for the histogram creation, the introduction is as follows:

#### sampleConfig and cutConfig the same as introduced in the skimming chapter.

#### varConfig
An example is in `json/variables/variables.json`. This config is essentially a map from the variable name that could be drawn to the histogram properties.
- `key`: Name of the variable. After registered here, such variables could be plotted as histograms.
- `label`: For the histogram x-axis title, can use ROOT-style latex. 
- `nBins`: self-documented
- `min`: self-documented
- `max`: self-documented
- `binLabels`: NOT implemented yet.
Future to add more configurable options including variable binning and bin labels.

#### lumiConfig and XSConfig
Examples in `json/Lumi/Run3.json` and `json/XS/Run3.json`. Mapping from the era to the Golden json integrated luminosity in the fb-1 unit, and mapping from the channel (not allow composite channels) to the expected cross-section in pb unit. 
The user could edit the values according to the need, but make sure the keep the unit correct.

----------

### plotHists

The step for visualize the histograms. Back up the plot of different types of comparison and stacked plotting. The styles currently are hard-coded. Later may develop the flexible styles.

##TODO 
- color skim synchronized with the LFV study
- change the data to no line, suitable size of dots with error bar.

The first step is to collect the histograms as desired using files
```
Root/
|--HistControl.cxx
|--PlotControl.cxx
rDataFrameLight/
|--HistControl.h
|--PlotControl.h
utils/
│--plotHists.cxx
```
The running command example is like:
```
plotHists json/MetricMuonHists/medium_dimuon/S1_dimuon_plot_muonVar_Run2024.json
```

In the json config for the histogram plotting, the user needs to provide:
- `era`: A list for all the corresponding errors. Also make sure that all of the era's histograms have been already created in the previous steps.
- `mc_era`: Considering the possibility that not always having the matched campaign MC for the histogram plotting. In such cases, we need this to re-weighting the MC.
- `varName`: Necessary. A list specifying what variables to be plotting. Multiple entries leading multiple plots output.
- `datasets`: The list of all datasets involved for plotting. Either data or mc, either stacked or not, everything wanted to plot must input here.
- `needMerge`: Similar to the one in collectingHists, a second merging could be applied here.
- `isData`: Necessary specification of whether is data or mc, influencing scaling, plotting styles.

- `inDir`: Necessary. The common name of the input histograms's directory (not contain the era).
- `dataWeight`: Not affect histograms, but used to find the input data histograms.
- `MCWeight`: Not affect histograms, but used to find the input MC histograms.
- `outDir`: Necessary. Defines the output directory.
- `name`: Used for the output nomenclature.

- `yLabel`: Visualization config. By default should use "Events", but if turning on the normalization plots, we should use "a.u." then. If there is the need to use other name, need to specify.
- `yRatioLabel`: Similarly visualization config. Used in case of need the ratio plot. The default name would be "Data / MC". If the other ratios are plotted, need to specify.
- `histXSize`: Visualization config, for the canvas width.
- `histYSize`: Visualization config, for the canvas height.
- `doLogPlot`: Visualization config, when turning this on, we are plotting the histograms (not affecting the ratio plots) in log scale.
- `doRatioPlot`: Visualization config, when turning this on, the plots have upper panel for the original histograms and the lower panel for ratio histograms. Use together with `numerator` and `stackOrder` options.
- `doNormPlot`: Visualization config, when turning this one, the plots do normalization or each histograms (THStack considered as one). Do NOT use together with ratio plot or log plots!

- `needCrop`: Visualization config, specifying which variable we want to have the zoom in plot.
- `cropedRange`: pairing with the `needCrop` option, gives the range for that.
- `colorMapping`: Visualization config. Determines the color for each. Considering take out to be one additional config to synchronize the color scheme. Should be for the final channels (if with needMerge, then only post-merge).
- `datasetLabel`: Visualization config. Determines the TLegend texts. Considering take out to be one additional config to synchronize the legend scheme. Should be for the final channels (if with needMerge, then only post-merge).
- `isSignal`: Visualization config. signal processes and other processes have different plotting styles. signal has amplification to make visible.
- `header`: For the text to be plot above the upper panel histograms. Only allow two sentences. Usually fixed texts.
- `texts`: For the TLatex to write in the empty area of the above panel, specifying the information of the histograms.

- `numerator`: Serving as the numerator when plotting the ratio histograms. Usually is data.
- `stackOrder`: Serving as the denominator when plotting the ratio histograms. Order could be fixed as input order if turning off `reOrder` option.
- `reOrder`: Allowing the auto-ordering withing the THStack to make each component visible. Especially used for the log scale plots.

- `varConfig`: As introduced at collectHists.
- `systConfig`: For the syst up/down ratio joining here. Need additional computations to prepare.
- `lumiConfig`: luminosity input again for the mc rescaling when using un-matching data and MC.

Optional configs:
- `jobType`: For safety checking, also recommended to always writing as other plotHists configs do.
- `comments`: Any comments to write.

The config in the plot is a bit complicated now. The further simplification is under consideration.
In case of having the data-driven contribution (themselves from data, but used similar to MC, different weight from MC), we need additional options:
- `dd_era`
- `dd_datasets`

These are used in the same way as using the mc. 

----------

### Additional scripts and utilities for achieving some simple tasks

Most of these scripts need additional modifications to be usable.

#### scirpt/checkTriggerEff.py
A quick script for testing the trigger efficiency on the signal MC samples.

#### scirpt/createFileList.py
A quick script for testing the trigger efficiency on the signal MC samples.
```
createFileLists.py
```

#### script/splitPromptStats.py

#### script/reportPromptStats.py
A quick check for the 

