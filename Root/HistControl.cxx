#include "HistControl.h"
#include "Utility.h"

#include "TFile.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <memory>

namespace
{
std::vector<int> getBinsFullyInsideRange(const TH1D *hist, double newMin, double newMax)
{
    std::vector<int> bins;
    if (hist == nullptr)
        return bins;

    const int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i)
    {
        const double lowEdge = hist->GetXaxis()->GetBinLowEdge(i);
        const double upEdge = hist->GetXaxis()->GetBinUpEdge(i);
        if (lowEdge >= newMin && upEdge <= newMax)
            bins.push_back(i);
    }
    return bins;
}

std::vector<double> getTemplateBinEdgesFromBins(const TH1D *hist, const std::vector<int> &bins)
{
    std::vector<double> edges;
    if (hist == nullptr || bins.empty())
        return edges;

    edges.reserve(bins.size() + 1);
    for (const int bin : bins)
        edges.push_back(hist->GetXaxis()->GetBinLowEdge(bin));
    edges.push_back(hist->GetXaxis()->GetBinUpEdge(bins.back()));
    return edges;
}
}

// providing a template to initialize, then the bin setup for filling would not be needed later
HistControl::HistControl(const std::string &varName, TH1D *templateHist)
    : _varName(varName), _templateHist(templateHist)
{
    if (templateHist)
        rdfWS_utility::messageINFO("HistControl", std::string("Initializing with a template provided: ") + templateHist->GetName());
}

// using TH1D* pointer, thus need to release memory when destruct
HistControl::~HistControl()
{
    if (_templateHist)
    {
        delete _templateHist;
    }
    for (auto &[name, hist] : this->_histograms)
    {
        if (hist)
            delete hist;
    }
    _histograms.clear();
}

std::string HistControl::buildHistogramKey(const std::string &datasetName, const std::string &weightName) const
{
    return datasetName + "_" + this->_varName + "_" + weightName;
}

/// @brief The function for giving the minimum and maximum from the DataFrame itself
/// @param rnode input node of RDataFrame
/// @param varName variable name
/// @param stripeLow whether to negelect 1% low extreme values
/// @param stripeHigh whether to negelect 1% high extreme values
/// @return
std::tuple<double, double> HistControl::getMinMax(ROOT::RDF::RNode &rnode, const std::string &varName, bool stripeLow, bool stripeHigh)
{
    // Initialize variables
    double globalMin = std::numeric_limits<double>::max();
    double globalMax = std::numeric_limits<double>::lowest();
    size_t totalEntries = 0;

    // Case 1: No stripping (use only global min and max)
    if (!stripeLow && !stripeHigh)
    {
        rnode.Foreach([&](double value)
                      {
            ++totalEntries;
            if (value < globalMin) globalMin = value;
            if (value > globalMax) globalMax = value; }, {varName});

        if (totalEntries == 0)
        {
            throw std::runtime_error("No data available in variable: " + varName);
        }

        return {globalMin, globalMax};
    }

    // Case 2: Only strip low extremes
    if (stripeLow && !stripeHigh)
    {
        std::priority_queue<double> highHeap; // Max-heap to track smallest values
        rnode.Foreach([&](double value)
                      {
            ++totalEntries;
            if (value > globalMax) globalMax = value;

            highHeap.push(value);
            if (highHeap.size() > totalEntries * 0.01) {
                highHeap.pop(); // Keep only the smallest 1%
            } }, {varName});

        if (totalEntries == 0)
        {
            throw std::runtime_error("No data available in variable: " + varName);
        }

        double minValue = highHeap.top();
        return {minValue, globalMax};
    }

    // Case 3: Only strip high extremes
    if (!stripeLow && stripeHigh)
    {
        std::priority_queue<double, std::vector<double>, std::greater<>> lowHeap; // Min-heap to track largest values
        rnode.Foreach([&](double value)
                      {
            ++totalEntries;
            if (value < globalMin) globalMin = value;

            lowHeap.push(value);
            if (lowHeap.size() > totalEntries * 0.01) {
                lowHeap.pop(); // Keep only the largest 1%
            } }, {varName});

        if (totalEntries == 0)
        {
            throw std::runtime_error("No data available in variable: " + varName);
        }

        double maxValue = lowHeap.top();
        return {globalMin, maxValue};
    }

    // Case 4: Strip both low and high extremes
    std::priority_queue<double> highHeap;                                     // Max-heap to track smallest values
    std::priority_queue<double, std::vector<double>, std::greater<>> lowHeap; // Min-heap to track largest values

    rnode.Foreach([&](double value)
                  {
        ++totalEntries;
        highHeap.push(value);
        if (highHeap.size() > totalEntries * 0.01) {
            highHeap.pop(); // Keep only the smallest 1%
        }

        lowHeap.push(value);
        if (lowHeap.size() > totalEntries * 0.01) {
            lowHeap.pop(); // Keep only the largest 1%
        } }, {varName});

    if (totalEntries == 0)
    {
        throw std::runtime_error("No data available in variable: " + varName);
    }

    double minValue = highHeap.top();
    double maxValue = lowHeap.top();
    return {minValue, maxValue};
}

/// @brief Internal function for directly setting up a template from binning
/// @param nBins
/// @param min
/// @param max
void HistControl::generateTemplateFromBinning(int nBins, double min, double max)
{
    if (_templateHist)
    {
        throw std::runtime_error("Template histogram already exists.");
    }
    if (nBins <= 0)
    {
        throw std::runtime_error("Invalid number of bins specified.");
    }

    this->_templateHist = new TH1D("templateHist", "Template Histogram", nBins, min, max);
    rdfWS_utility::messageINFO("HistControl", "Template histogram dynamically created from binning: bins = " + std::to_string(nBins) + ", range = [" + std::to_string(min) + ", " + std::to_string(max) + "].");
}

void HistControl::generateTemplateFromBinning(const std::vector<double> &binEdges)
{
    if (_templateHist)
    {
        throw std::runtime_error("Template histogram already exists.");
    }
    if (binEdges.size() < 2)
    {
        throw std::runtime_error("Invalid variable binning specified.");
    }

    this->_templateHist = new TH1D("templateHist", "Template Histogram", static_cast<int>(binEdges.size()) - 1, binEdges.data());
    rdfWS_utility::messageINFO("HistControl", "Template histogram dynamically created from variable binning.");
}

/// @brief Internal function to compute the ratio of two histograms. In principle, all histograms should be processed through the hist controller and the TH1D* easily causing leakage, so that this function is not public. May change upon to the need
/// @param numerator
/// @param denominator
/// @param ratioName
/// @return
TH1D *HistControl::calculateRatio(TH1D *numerator, TH1D *denominator, const std::string &ratioName)
{
    if (!numerator || !denominator)
        rdfWS_utility::messageERROR("HistControl", "Invalid numerator or denominator for ratio calculation.");

    TH1D *ratioHist = (TH1D *)numerator->Clone(ratioName.c_str());
    ratioHist->SetDirectory(0);
    std::string ratioHistTitle = numerator->GetName() + std::string(" / ") + denominator->GetName();
    ratioHist->SetTitle(ratioHistTitle.c_str());

    // settting bin contents
    double maxBinContent(-1e9);
    double minBinContent(1e9);
    for (int i = 1; i <= numerator->GetNbinsX(); ++i)
    {
        double N_num = numerator->GetBinContent(i);
        double N_denom = denominator->GetBinContent(i);
        double sigma_num = numerator->GetBinError(i);
        double sigma_denom = denominator->GetBinError(i);

        if (N_denom == 0)
        {
            rdfWS_utility::messageWARN("HistControl", "The " + std::to_string(i) + "-th bin of denominator is empty. Set the ratio bin to 0");
            ratioHist->SetBinContent(i, 0);
            ratioHist->SetBinError(i, 0);
        }
        else
        {
            double ratio = N_num / N_denom;
            double ratioError = std::sqrt((sigma_num / N_denom) * (sigma_num / N_denom) +
                                          (sigma_denom * N_num / (N_denom * N_denom)) * (sigma_denom * N_num / (N_denom * N_denom)));

            ratioHist->SetBinContent(i, ratio);
            ratioHist->SetBinError(i, ratioError);
            if (ratio > maxBinContent)
                maxBinContent = ratio;
            if (ratio < minBinContent)
                minBinContent = ratio;
        }
    }
    // in case of all empty
    if (maxBinContent < -1e8)
        maxBinContent = 1.2;
    if (minBinContent > 1e8)
        minBinContent = 0.8;
    ratioHist->SetMaximum(maxBinContent);
    ratioHist->SetMinimum(minBinContent);

    return ratioHist;
}

//////////////////////////////////////////////////
/// User interface of the HistControl class
//////////////////////////////////////////////////

/// @brief Creating histograms by filling from a RDataFrame. This function will also save the corresponding histogram generated in the internal map to ENABLE REUSING and AVOID OVERWRITING
/// @param rnode
/// @param datasetName internal identity
/// @param binning binning information of the histogram, can neglect or input nullptr, then the code will find one from dataset
/// @param varName: Can omit the varName, in case there is already varName provided. Note if there is no varName stored initially, you must provide it!
/// @param weightName: By default will take "one".
/// @param ifSave
/// @return histogram extracted
TH1D *HistControl::createHistogram(ROOT::RDF::RNode &rnode, const std::string &datasetName, const HistBinning *binning, const std::string &varName, const std::string &weightName, const std::string &outDir, bool ifSave)
{
    if (this->_varName.empty())
    {
        if (varName.empty())
            rdfWS_utility::messageERROR("HistControl", "Variable name for histograms must be explicitly provided as no internal variable name is set.");
        this->_varName = varName;
    }
    else if (!varName.empty() && this->_varName != varName)
        rdfWS_utility::messageERROR("HistControl", "Variable name mismatch. Expected: " + this->_varName + ", but got: " + varName);

    bool tempBinning = false;
    if (binning == nullptr)
    {
        binning = new HistBinning();
        tempBinning = true;
    }

    if (!this->_templateHist)
    {
        if (!binning->varBins.empty())
        {
            if (binning->nBins != -1 && static_cast<int>(binning->varBins.size()) != binning->nBins + 1)
            {
                throw std::runtime_error("Variable binning size does not match nBins + 1.");
            }
            generateTemplateFromBinning(binning->varBins);
        }
        else
        {
            double min = binning->min;
            double max = binning->max;
            int bins = binning->nBins;
            if (bins == -1)
            {
                std::tie(min, max) = this->getMinMax(rnode, this->_varName, binning->stripeLow, binning->stripeHigh);
                bins = binning->defaultNBins;
            }
            generateTemplateFromBinning(bins, min, max);
        }
    }

    if (tempBinning)
    {
        delete binning;
        binning = nullptr;
    }

    if (!rnode.HasColumn("one"))
    {
        rnode = rnode.Define("one", "1.0");
    }

    ROOT::RDF::TH1DModel histModel("templateHistModel", datasetName.c_str(), this->_templateHist->GetNbinsX(), this->_templateHist->GetXaxis()->GetXmin(), this->_templateHist->GetXaxis()->GetXmax());
    if (this->_templateHist->GetXaxis()->GetXbins()->GetSize() > 0)
    {
        histModel = ROOT::RDF::TH1DModel("templateHistModel", datasetName.c_str(), this->_templateHist->GetNbinsX(), this->_templateHist->GetXaxis()->GetXbins()->GetArray());
    }

    auto hist = rnode.Histo1D(histModel, this->_varName, weightName);
    std::string histName = datasetName + "_" + this->_varName + "_" + weightName;
    TH1D *histClone = (TH1D *)hist.GetPtr()->Clone(histName.c_str());
    histClone->SetDirectory(0);
    const std::string histKey = buildHistogramKey(datasetName, weightName);
    auto [it, inserted] = this->_histograms.emplace(histKey, histClone);
    if (!inserted)
    {
        delete histClone;
        throw std::runtime_error("Histogram with key " + histKey + " already exists.");
    }
    if (ifSave)
        saveHistogram(histClone, histName, outDir);
    return histClone;
}

/// @brief Directly load histograms from files, to save re-filling time from RDataFrame
/// @param fileName
/// @param histName
/// @param histKey internal identity
/// @param varName: can omit if already initialized
/// @param additionalName: can omit, just to avoid hist name collisiton
/// @return
void HistControl::loadHistogram(const std::string &fileName, const std::string &histName, const std::string &histKey, float scaling, const std::string &varName, const std::string &additionalName)
{
    if (_varName.empty())
    {
        if (varName.empty())
        {
            throw std::runtime_error("Variable name for histograms must be explicitly provided as no internal variable name is set.");
        }
        _varName = varName;
    }
    else if (!varName.empty() && _varName != varName)
    {
        throw std::runtime_error("Variable name mismatch. Expected: " + _varName + ", but got: " + varName);
    }

    if (!std::filesystem::exists(std::filesystem::path(fileName)))
    {
        throw std::runtime_error("Failed to open hist TFile: " + fileName);
    }

    TFile *f1 = TFile::Open(fileName.c_str(), "READ");
    if (!f1 || f1->IsZombie())
    {
        delete f1;
        throw std::runtime_error("Failed to open hist TFile: " + fileName);
    }

    auto hist = (TH1D *)f1->Get(histName.c_str());
    if (!hist)
    {
        f1->Close();
        delete f1;
        throw std::runtime_error("Histogram " + histName + " not found in file " + fileName);
    }
    hist->SetDirectory(0);
    f1->Close();
    delete f1;
    hist->Scale(scaling);

    const int maxBin = hist->GetNbinsX();
    const double underflowContent = hist->GetBinContent(0);
    const double firstContent = hist->GetBinContent(1);
    const double underflowError = hist->GetBinError(0);
    const double firstError = hist->GetBinError(1);
    hist->SetBinContent(1, firstContent + underflowContent);
    hist->SetBinError(1, std::hypot(firstError, underflowError));

    const double overflowContent = hist->GetBinContent(maxBin + 1);
    const double lastContent = hist->GetBinContent(maxBin);
    const double overflowError = hist->GetBinError(maxBin + 1);
    const double lastError = hist->GetBinError(maxBin);
    hist->SetBinContent(maxBin, lastContent + overflowContent);
    hist->SetBinError(maxBin, std::hypot(lastError, overflowError));

    if (this->_templateHist == nullptr)
    {
        this->_templateHist = (TH1D *)hist->Clone(("template_" + this->_varName + additionalName).c_str());
        this->_templateHist->SetDirectory(0);
        this->_templateHist->Reset();
    }
    else
    {
        bool sameBinning = this->_templateHist->GetNbinsX() == hist->GetNbinsX();
        if (sameBinning)
        {
            for (int i = 1; i <= hist->GetNbinsX() + 1; ++i)
            {
                if (std::abs(this->_templateHist->GetXaxis()->GetBinLowEdge(i) - hist->GetXaxis()->GetBinLowEdge(i)) > 1e-9)
                {
                    sameBinning = false;
                    break;
                }
            }
        }
        if (!sameBinning)
        {
            delete hist;
            throw std::runtime_error("Binning of histogram " + histName + " does not match the template histogram.");
        }
    }

    std::string histFullName = hist->GetName();
    histFullName += additionalName;
    hist->SetName(histFullName.c_str());
    if (_histograms.find(histKey) != _histograms.end())
    {
        delete hist;
        throw std::runtime_error("Histogram with key " + histKey + " already exists.");
    }
    this->_histograms.emplace(histKey, hist);
    return;
}

void HistControl::addHistogram(const TH1D *hist, const std::string &histKey, const std::string &varName, const std::string &additionalName)
{
    if (!hist)
    {
        throw std::runtime_error("Cannot add a null histogram.");
    }

    if (_varName.empty())
    {
        if (varName.empty())
        {
            throw std::runtime_error("Variable name for histograms must be explicitly provided as no internal variable name is set.");
        }
        _varName = varName;
    }
    else if (!varName.empty() && _varName != varName)
    {
        throw std::runtime_error("Variable name mismatch. Expected: " + _varName + ", but got: " + varName);
    }

    TH1D *histClone = static_cast<TH1D *>(hist->Clone());
    histClone->SetDirectory(0);

    if (this->_templateHist == nullptr)
    {
        this->_templateHist = static_cast<TH1D *>(histClone->Clone(("template_" + this->_varName + additionalName).c_str()));
        this->_templateHist->SetDirectory(0);
        this->_templateHist->Reset();
    }
    else
    {
        bool sameBinning = this->_templateHist->GetNbinsX() == histClone->GetNbinsX();
        if (sameBinning)
        {
            for (int i = 1; i <= histClone->GetNbinsX() + 1; ++i)
            {
                if (std::abs(this->_templateHist->GetXaxis()->GetBinLowEdge(i) - histClone->GetXaxis()->GetBinLowEdge(i)) > 1e-9)
                {
                    sameBinning = false;
                    break;
                }
            }
        }
        if (!sameBinning)
        {
            delete histClone;
            throw std::runtime_error("Binning of histogram " + std::string(hist->GetName()) + " does not match the template histogram.");
        }
    }

    std::string histFullName = histClone->GetName();
    histFullName += additionalName;
    histClone->SetName(histFullName.c_str());
    if (_histograms.find(histKey) != _histograms.end())
    {
        delete histClone;
        throw std::runtime_error("Histogram with key " + histKey + " already exists.");
    }
    this->_histograms.emplace(histKey, histClone);
}

/// @brief Save histograms when created from histograms, maybe should move to internal
/// @param hist
/// @param fileName
void HistControl::saveHistogram(const TH1D *hist, const std::string &fileName, const std::string &outDir)
{
    std::string outFilePath = fileName;
    if (outDir != "")
    {
        rdfWS_utility::creatingFolder("HistControl", outDir);
        outFilePath = outDir + "/" + fileName + ".root";
    }
    TFile *saveFile = TFile::Open(outFilePath.c_str(), "RECREATE");
    if (!saveFile || saveFile->IsZombie())
    {
        delete saveFile;
        throw std::runtime_error("Failed to create output ROOT file: " + outFilePath);
    }
    hist->Write();
    saveFile->Close();
    delete saveFile;
}

void HistControl::removeHistogram(const std::string &datasetName, const std::string &weightName)
{
    const std::string histKey = buildHistogramKey(datasetName, weightName);
    auto it = this->_histograms.find(histKey);
    if (it == this->_histograms.end())
        return;
    delete it->second;
    this->_histograms.erase(it);
}

void HistControl::clearHistograms()
{
    for (auto &[key, hist] : this->_histograms)
    {
        delete hist;
    }
    this->_histograms.clear();
}

// take some modification to restrict extraction range
// so that the merged hist storage would not cause trouble
std::map<std::string, TH1D *> HistControl::getHists(std::vector<std::string> histKeys)
{
    if (this->_histograms.size() == 0)
    {
        throw std::runtime_error("No histograms prepared");
    }
    std::map<std::string, TH1D *> extractedHist;
    // no given keys means extract everything
    if (histKeys.size() == 0)
    {
        for (const auto &[key, hist] : this->_histograms)
        {
            extractedHist.emplace(key, (TH1D *)hist->Clone((key + "_out").c_str()));
            extractedHist[key]->SetDirectory(0);
        }
        return extractedHist;
    }
    // with given keys
    for (auto key : histKeys)
    {
        extractedHist.emplace(key, (TH1D *)this->_histograms[key]->Clone((key + "_out").c_str()));
        extractedHist[key]->SetDirectory(0);
    }
    return extractedHist;
}

/// @brief update cropHistograms function into create a new histcontrol, as the croped hist might also be used for ratio computation and merging operations
/// @param newMin
/// @param newMax
/// @return
HistControl HistControl::cropHistograms(double newMin, double newMax)
{
    // checking crop range
    if (!this->_templateHist)
    {
        throw std::runtime_error("No template histogram available for binning check.");
    }

    double origMin = this->_templateHist->GetXaxis()->GetXmin();
    double origMax = this->_templateHist->GetXaxis()->GetXmax();

    if (newMin < origMin || newMax > origMax)
    {
        throw std::runtime_error("Requested range [" + std::to_string(newMin) + ", " + std::to_string(newMax) +
                                 "] exceeds original range [" + std::to_string(origMin) + ", " + std::to_string(origMax) + "].");
    }

    const std::vector<int> selectedBins = getBinsFullyInsideRange(this->_templateHist, newMin, newMax);
    if (selectedBins.empty())
    {
        throw std::runtime_error("Requested crop range [" + std::to_string(newMin) + ", " + std::to_string(newMax) +
                                 "] does not fully contain any histogram bin.");
    }

    const std::vector<double> croppedEdges = getTemplateBinEdgesFromBins(this->_templateHist, selectedBins);
    const int nNewBins = static_cast<int>(selectedBins.size());

    TH1D *newTemplate = new TH1D("croppedTemplate", "Cropped Template Histogram",
                                 nNewBins, croppedEdges.data());
    newTemplate->SetDirectory(0);
    HistControl croppedControl(this->_varName, newTemplate);

    for (const auto &[key, hist] : this->_histograms)
    {
        std::string newName = std::string(hist->GetName()) + "_cropped_" + std::to_string(newMin) + "_" + std::to_string(newMax);

        TH1D *croppedHist = new TH1D(newName.c_str(), hist->GetTitle(),
                                     nNewBins, croppedEdges.data());
        croppedHist->SetDirectory(0);

        for (int outputBin = 1; outputBin <= nNewBins; ++outputBin)
        {
            const int inputBin = selectedBins[outputBin - 1];
            const double content = hist->GetBinContent(inputBin);
            const double error = hist->GetBinError(inputBin);
            croppedHist->SetBinContent(outputBin, content);
            croppedHist->SetBinError(outputBin, error);
        }

        croppedControl._histograms[key] = croppedHist;
    }

    return croppedControl;
}

/// @brief Function to merge internal stored histograms. Indeed this might be needed when plotting histograms anywhere, thus I set this public
/// @param histNames
/// @param mergedName
/// @return
TH1D *HistControl::mergeHistograms(const std::vector<std::string> &histKeys, const std::string &mergedKey)
{
    if (histKeys.empty())
    {
        rdfWS_utility::messageERROR("HistControl", "No histograms specified for merging.");
    }
    rdfWS_utility::messageINFO("HistControl", "Merging histograms.");

    TH1D *mergedHist = nullptr;
    std::string mergedHistName = mergedKey + "_mergedFrom";
    for (auto key : histKeys)
    {
        mergedHistName += ("_" + key);
    }

    // summing from internal stored hists
    for (const auto &name : histKeys)
    {
        auto it = this->_histograms.find(name);
        if (it == _histograms.end())
        {
            throw std::runtime_error("Histogram " + name + " not found for merging.");
        }
        rdfWS_utility::messageINFO("HistControl", std::string("Hist name: ") + it->second->GetName());

        // at the first loop
        if (!mergedHist)
        {
            mergedHist = (TH1D *)it->second->Clone(mergedHistName.c_str());
            mergedHist->SetDirectory(0);
            mergedHist->Reset();
        }

        mergedHist->Add(it->second);
    }

    // consider add it the map as well
    // as I would need some manually merged things like Diboson
    if (mergedHist)
        this->_histograms.emplace(mergedKey, mergedHist);

    return mergedHist;
}

/// @brief Function to get ratio histograms. The computation takes one internal stored reference histograms. The reference can be summed up of multiple histograms
/// @param numerator The neumerators need to take ratio with up to the summed up denominator
/// @param referenceName all reference histograms, to sum up
/// @return
std::map<std::string, TH1D *> HistControl::getRatios(const std::vector<std::string> &numerator, const std::vector<std::string> &referenceNames, bool doNormalize)
{
    std::map<std::string, TH1D *> ratioHists;

    // preparing reference hists, say MC, or SM, ... depending on your need
    TH1D *mergedRef = nullptr;
    if (!referenceNames.empty())
    {
        mergedRef = mergeHistograms(referenceNames, "merged_reference");
    }
    std::string mergedRefName;
    for (size_t i = 0; i < referenceNames.size(); ++i)
    {
        mergedRefName += referenceNames[i];
        if (i < referenceNames.size() - 1)
        {
            mergedRefName += "+";
        }
    }

    double mergedInt = mergedRef->Integral();

    // only get ratio plots for the channels needed
    for (auto key : numerator)
    {
        if (this->_histograms.find(key) == this->_histograms.end())
        {
            throw std::runtime_error("ratio hist numerator channel " + key + " does not exist in internal histograms");
        }
        double numerInt = this->_histograms[key]->Integral();
        std::string ratioName = key + "_ratio_to_" + mergedRefName;
        TH1D *ratioHist = calculateRatio(this->_histograms[key], mergedRef, ratioName);
        if (doNormalize && numerInt > 0)
        {
            ratioHist->Scale(mergedInt / numerInt);
        }
        ratioHists[ratioName] = ratioHist;
    }

    return ratioHists;
}

/// @brief add function to merge (addup) two HistControl set
/// @param toAdd
/// @return
HistControl HistControl::addHistograms(const HistControl &toAdd)
{
    if (this->_histograms.size() == 0)
    {
        HistControl result = toAdd;
        return result;
    }

    if (this->_varName != toAdd._varName)
    {
        throw std::runtime_error("Variable name mismatch between the two HistControl objects.");
    }

    auto sameBinning = [](const TH1D *lhs, const TH1D *rhs)
    {
        if (lhs == nullptr || rhs == nullptr)
            return lhs == rhs;
        if (lhs->GetNbinsX() != rhs->GetNbinsX())
            return false;
        for (int i = 1; i <= lhs->GetNbinsX() + 1; ++i)
        {
            if (std::abs(lhs->GetXaxis()->GetBinLowEdge(i) - rhs->GetXaxis()->GetBinLowEdge(i)) > 1e-9)
                return false;
        }
        return true;
    };

    HistControl result(this->_varName);

    if (this->_templateHist && toAdd._templateHist)
    {
        if (!sameBinning(this->_templateHist, toAdd._templateHist))
        {
            throw std::runtime_error("Template histogram binning mismatch between the two HistControl objects.");
        }
        std::string templateName = this->_templateHist->GetName();
        result._templateHist = (TH1D *)this->_templateHist->Clone((templateName + "_sum").c_str());
        result._templateHist->SetDirectory(0);
    }
    else if (this->_templateHist)
    {
        std::string templateName = this->_templateHist->GetName();
        result._templateHist = (TH1D *)this->_templateHist->Clone((templateName + "_sum").c_str());
        result._templateHist->SetDirectory(0);
    }
    else if (toAdd._templateHist)
    {
        std::string templateName = toAdd._templateHist->GetName();
        result._templateHist = (TH1D *)toAdd._templateHist->Clone((templateName + "_sum").c_str());
        result._templateHist->SetDirectory(0);
    }

    for (const auto &[key, hist] : this->_histograms)
    {
        auto it = toAdd._histograms.find(key);
        if (it != toAdd._histograms.end())
        {
            if (!sameBinning(hist, it->second))
            {
                throw std::runtime_error("Binning mismatch for histogram key: " + key);
            }
            std::string sumName = hist->GetName();
            TH1D *sumHist = (TH1D *)hist->Clone((sumName + "_sum").c_str());
            sumHist->SetDirectory(0);
            sumHist->Add(it->second);
            result._histograms[key] = sumHist;
        }
        else
        {
            std::string clonedName = hist->GetName();
            TH1D *clonedHist = (TH1D *)hist->Clone((clonedName + "_cloned").c_str());
            clonedHist->SetDirectory(0);
            result._histograms[key] = clonedHist;
        }
    }

    for (const auto &[key, hist] : toAdd._histograms)
    {
        if (this->_histograms.find(key) == this->_histograms.end())
        {
            std::string clonedName = hist->GetName();
            TH1D *clonedHist = (TH1D *)hist->Clone((clonedName + "_cloned").c_str());
            clonedHist->SetDirectory(0);
            result._histograms[key] = clonedHist;
        }
    }

    return result;
}


void HistControl::absorbHistograms(const HistControl &toAdd)
{
    if (toAdd._histograms.empty())
        return;

    if (this->_histograms.empty())
    {
        *this = toAdd;
        return;
    }

    if (this->_varName != toAdd._varName)
    {
        throw std::runtime_error("Variable name mismatch between the two HistControl objects.");
    }

    auto sameBinning = [](const TH1D *lhs, const TH1D *rhs)
    {
        if (lhs == nullptr || rhs == nullptr)
            return lhs == rhs;
        if (lhs->GetNbinsX() != rhs->GetNbinsX())
            return false;
        for (int i = 1; i <= lhs->GetNbinsX() + 1; ++i)
        {
            if (std::abs(lhs->GetXaxis()->GetBinLowEdge(i) - rhs->GetXaxis()->GetBinLowEdge(i)) > 1e-9)
                return false;
        }
        return true;
    };

    if (this->_templateHist && toAdd._templateHist)
    {
        if (!sameBinning(this->_templateHist, toAdd._templateHist))
        {
            throw std::runtime_error("Template histogram binning mismatch between the two HistControl objects.");
        }
    }
    else if (!this->_templateHist && toAdd._templateHist)
    {
        std::string templateName = toAdd._templateHist->GetName();
        this->_templateHist = (TH1D *)toAdd._templateHist->Clone((templateName + "_cloned").c_str());
        this->_templateHist->SetDirectory(0);
    }

    for (const auto &[key, hist] : toAdd._histograms)
    {
        auto it = this->_histograms.find(key);
        if (it != this->_histograms.end())
        {
            if (!sameBinning(it->second, hist))
            {
                throw std::runtime_error("Binning mismatch for histogram key: " + key);
            }
            it->second->Add(hist);
        }
        else
        {
            std::string clonedName = hist->GetName();
            TH1D *clonedHist = (TH1D *)hist->Clone((clonedName + "_cloned").c_str());
            clonedHist->SetDirectory(0);
            this->_histograms[key] = clonedHist;
        }
    }
}

HistControl::HistControl(const HistControl &other)
{
    _varName = other._varName;
    if (other._templateHist)
    {
        std::string oriName = other._templateHist->GetName();
        _templateHist = (TH1D *)other._templateHist->Clone((oriName + "_copied").c_str());
        _templateHist->SetDirectory(0);
    }
    for (const auto &[key, hist] : other._histograms)
    {
        std::string oriName = hist->GetName();
        _histograms[key] = (TH1D *)hist->Clone((oriName + "_copied").c_str());
        _histograms[key]->SetDirectory(0);
    }
}

HistControl &HistControl::operator=(const HistControl &other)
{
    if (this != &other)
    {
        delete _templateHist;
        _templateHist = nullptr;
        for (auto &[key, hist] : _histograms)
        {
            delete hist;
        }
        _histograms.clear();

        _varName = other._varName;
        if (other._templateHist)
        {
            std::string oriName = other._templateHist->GetName();
            _templateHist = (TH1D *)other._templateHist->Clone((oriName + "_copied").c_str());
            _templateHist->SetDirectory(0);
        }
        for (const auto &[key, hist] : other._histograms)
        {
            std::string oriName = hist->GetName();
            _histograms[key] = (TH1D *)hist->Clone((oriName + "_copied").c_str());
            _histograms[key]->SetDirectory(0);
        }
    }
    return *this;
}

HistControl::HistControl(HistControl &&other) noexcept
{
    _varName = std::move(other._varName);
    _histograms = std::move(other._histograms);
    _templateHist = other._templateHist;
    other._templateHist = nullptr;
    other._histograms.clear();
}

HistControl &HistControl::operator=(HistControl &&other) noexcept
{
    if (this != &other)
    {
        delete _templateHist;
        for (auto &[key, hist] : _histograms)
        {
            delete hist;
        }
        _histograms.clear();

        _varName = std::move(other._varName);
        _histograms = std::move(other._histograms);
        _templateHist = other._templateHist;
        other._templateHist = nullptr;
        other._histograms.clear();
    }
    return *this;
}


std::map<std::string, TH1D*> HistControl::getHistInstance()
{
    std::cout << "You are directly accessing the histogram!" << std::endl;
    return this->_histograms;
}





