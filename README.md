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
wget https://raw.githubusercontent.com/nlohmann/json/develop/include/nlohmann/json.hpp
cd ../../build/
source ../source
make
```

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