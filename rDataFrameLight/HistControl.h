#ifndef HISTOGRAM_CONTROLLER_H
#define HISTOGRAM_CONTROLLER_H

#include <string>
#include <map>
#include <tuple>
#include <optional>

#include "ROOT/RDataFrame.hxx"
#include "TH1D.h"

// default structure for controling binning information.
// used in case of no template hist provided.
struct HistBinning
{
    // if nBins==-1 && bins.size()==0, will go to auto bins
    int nBins = -1;
    int defaultNBins = 100;
    double min = -1;
    double max = -1;
    // to add variable binning later
    // std::vector<double> varBins={};
    // if auto bin, by default stripe out extreme values
    bool stripeHigh = false;
    bool stripeLow = false;
};

class HistControl
{
private:
    // can introduce information from json as well
    // need to determine what histograms to make
    // how to store and use them efficiently

    std::string _varName; // designed for histograms from multiple dataset on the same variable
    std::map<std::string, TH1D*> _histograms; // internal stored histograms
    TH1D* _templateHist = nullptr;

    // internal needed functions, very slow, better not rely on this, but using json
    std::tuple<double, double> getMinMax(ROOT::RDF::RNode &rnode, const std::string &varName, bool stripeLow, bool stripeHigh);

    // need always get the template first
    void generateTemplateFromBinning(int bins, double min, double max);

    TH1D *calculateRatio(TH1D *numerator, TH1D *denominator, const std::string &ratioName);

public:
    explicit HistControl() = default;
    explicit HistControl(const std::string &varName): _varName(varName){};
    explicit HistControl(const std::string &varName, TH1D *templateHist);
    ~HistControl();

    // Create a histogram from RDataFrame
    TH1D *createHistogram(ROOT::RDF::RNode &rnode, const std::string &datasetName, const HistBinning *binning=nullptr, const std::string &varName="", const std::string& weightName="one", const std::string& outDir="");
    // Load a histogram from a file
    TH1D *loadHistogram(const std::string &fileName, const std::string &histName, const std::string & histKey, const std::string &varName="", const std::string& additionalName="");
    // Save a histogram to a file
    void saveHistogram(const TH1D *hist, const std::string &fileName, const std::string& outDir="");

    // get required hists
    std::map<std::string, TH1D *> getHists(std::vector<std::string> histKeys);
    // get corped all internal histograms
    //std::map<std::string, TH1D *> cropHistograms(double newMin, double newMax);
    HistControl cropHistograms(double newMin, double newMax);
    // get merged internal histograms
    TH1D *mergeHistograms(const std::vector<std::string> &histNames, const std::string &mergedName);
    // get ratio hists
    std::map<std::string, TH1D *> getRatios(const std::vector<std::string>& numerator, const std::vector<std::string> &referenceNames);

    // to add histogram together
    HistControl addHistograms(const HistControl& toAdd);

    HistControl(const HistControl &other);
    HistControl &operator=(const HistControl &other);
};

#endif