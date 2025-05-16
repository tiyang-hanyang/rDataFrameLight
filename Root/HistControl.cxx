#include "HistControl.h"
#include "Utility.h"

#include "TFile.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <iostream>

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
/// @return histogram extracted
TH1D *HistControl::createHistogram(ROOT::RDF::RNode &rnode, const std::string &datasetName, const HistBinning *binning, const std::string &varName, const std::string &weightName, const std::string &outDir)
{
    // make sure the variables are correct, one HistControl object should only deal with one variable
    if (this->_varName.empty())
    {
        if (varName.empty())
            rdfWS_utility::messageERROR("HistControl", "Variable name for histograms must be explicitly provided as no internal variable name is set.");
        this->_varName = varName;
    }
    else if (!varName.empty() && this->_varName != varName)
        rdfWS_utility::messageERROR("HistControl", "Variable name mismatch. Expected: " + this->_varName + ", but got: " + varName);

    // avoid double creating
    if (this->_histograms.find(datasetName) != this->_histograms.end())
        rdfWS_utility::messageERROR("HistControl", "Histogram for dataset " + datasetName + " already exists.");

    // Need to generate _templateHist, if there is no one. Using provided HistBinning or auto binning
    if (!this->_templateHist)
    {
        int tempBinning = 0;
        if (binning == nullptr)
        {
            rdfWS_utility::messageINFO("HistControl", "No template hist exists, no HistBinning provided, will get minimum and maximum from dataset");
            binning = new HistBinning;
            tempBinning = 1;
        }

        double min = binning->min;
        double max = binning->max;
        int bins = binning->nBins;
        // bins should only be -1 with default value of HistBinning
        if (bins == -1)
        {
            std::tie(min, max) = this->getMinMax(rnode, this->_varName, binning->stripeLow, binning->stripeHigh);
            bins = binning->defaultNBins;
        }
        if (tempBinning)
        {
            delete binning;
            binning = nullptr;
        }

        generateTemplateFromBinning(bins, min, max);
    }

    // template hist always with priority !
    int nBins = this->_templateHist->GetNbinsX();
    double min = this->_templateHist->GetXaxis()->GetXmin();
    double max = this->_templateHist->GetXaxis()->GetXmax();

    // extract histogram from rdataframe
    if (!rnode.HasColumn("one"))
    {
        rnode.Define("one", "1.0");
    }
    auto hist = rnode.Histo1D(ROOT::RDF::TH1DModel(datasetName.c_str(), datasetName.c_str(), nBins, min, max), this->_varName, weightName);
    std::string histName = datasetName + "_" + this->_varName + "_" + weightName;
    TH1D *histClone = (TH1D *)hist.GetPtr()->Clone(histName.c_str());
    histClone->SetDirectory(0);
    // store in the internal histogram map, save a copy to the disk, and return this hist if need to use
    this->_histograms.emplace(datasetName, histClone);
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
TH1D *HistControl::loadHistogram(const std::string &fileName, const std::string &histName, const std::string &histKey, const std::string &varName, const std::string &additionalName)
{
    // make sure the variables are correct, one HistControl object should only deal with one variable
    // plottng multiple histograms together would later need the other expantions
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

    TFile *f1 = TFile::Open(fileName.c_str(), "READ");
    if (!f1 || f1->IsZombie())
    {
        throw std::runtime_error("Failed to open hist TFile: " + fileName);
    }

    auto hist = (TH1D *)f1->Get(histName.c_str());
    if (!hist)
    {
        f1->Close();
        throw std::runtime_error("Histogram " + histName + " not found in file " + fileName);
    }
    hist->SetDirectory(0);

    // compare to template to make sure the same binning
    if (this->_templateHist == nullptr)
    {
        this->_templateHist = (TH1D *)hist->Clone(("template_" + this->_varName + additionalName).c_str());
        this->_templateHist->SetDirectory(0);
        this->_templateHist->Reset();
    }
    else
    {
        if (this->_templateHist->GetNbinsX() != hist->GetNbinsX() ||
            this->_templateHist->GetXaxis()->GetXmin() != hist->GetXaxis()->GetXmin() ||
            this->_templateHist->GetXaxis()->GetXmax() != hist->GetXaxis()->GetXmax())
        {
            throw std::runtime_error("Binning of histogram " + histName + " does not match the template histogram.");
        }
    }

    hist->SetDirectory(0);
    std::string histFullName = hist->GetName();
    histFullName += additionalName;
    auto histClone = (TH1D *)hist->Clone(histFullName.c_str());
    histClone->SetDirectory(0);
    f1->Close();
    delete hist;
    if (_histograms.find(histName) != _histograms.end())
    {
        throw std::runtime_error("Histogram with name " + histName + " already exists.");
    }
    this->_histograms.emplace(histKey, histClone);
    return histClone;
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
    auto saveFile = new TFile(outFilePath.c_str(), "RECREATE");
    hist->Write();
    saveFile->Close();
}

// take some modification to restrict extraction range
// so that the merged hist storage would not cause trouble
std::map<std::string, TH1D *> HistControl::getHists(std::vector<std::string> histKeys)
{
    if (this->_histograms.size() == 0)
    {
        throw std::runtime_error("No histograms prepared");
    }
    // no given keys means extract everything
    if (histKeys.size() == 0)
        return this->_histograms;
    // with given keys
    std::map<std::string, TH1D *> extractedHist;
    for (auto key : histKeys)
    {
        extractedHist.emplace(key, (TH1D *)this->_histograms[key]->Clone((key + "_copy").c_str()));
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

    // get a new object with the same varName but cropped histograms
    int newFirstBin = this->_templateHist->FindBin(newMin);
    int newLastBin = std::min(this->_templateHist->FindBin(newMax),  this->_templateHist->GetNbinsX());
    int nNewBins = newLastBin - newFirstBin + 1;

    double croppedMin = this->_templateHist->GetXaxis()->GetBinLowEdge(newFirstBin);
    double croppedMax = this->_templateHist->GetXaxis()->GetBinUpEdge(newLastBin);

    TH1D *newTemplate = new TH1D("croppedTemplate", "Cropped Template Histogram",
                                 nNewBins, croppedMin, croppedMax);
    newTemplate->SetDirectory(0);
    HistControl croppedControl(this->_varName, newTemplate);

    for (const auto &[key, hist] : this->_histograms)
    {
        std::string newName = std::string(hist->GetName()) + "_cropped_" + std::to_string(newMin) + "_" + std::to_string(newMax);

        TH1D *croppedHist = new TH1D(newName.c_str(), hist->GetTitle(),
                                     nNewBins, croppedMin, croppedMax);
        croppedHist->SetDirectory(0);

        for (int i = newFirstBin; i <= newLastBin; ++i)
        {
            double content = hist->GetBinContent(i);
            double error = hist->GetBinError(i);
            croppedHist->SetBinContent(i - newFirstBin + 1, content);
            croppedHist->SetBinError(i - newFirstBin + 1, error);
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
        if (doNormalize && numerInt>0) {
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
    // Check addbility
    // varName must match
    if (this->_varName != toAdd._varName)
    {
        throw std::runtime_error("Variable name mismatch between the two HistControl objects.");
    }

    // Create a new HistControl for the result
    HistControl result(this->_varName);

    // Check to make sure template binning match
    if (this->_templateHist && toAdd._templateHist)
    {
        if (this->_templateHist->GetNbinsX() != toAdd._templateHist->GetNbinsX() ||
            this->_templateHist->GetXaxis()->GetXmin() != toAdd._templateHist->GetXaxis()->GetXmin() ||
            this->_templateHist->GetXaxis()->GetXmax() != toAdd._templateHist->GetXaxis()->GetXmax())
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
        std::string templateName = this->_templateHist->GetName();
        result._templateHist = (TH1D *)toAdd._templateHist->Clone((templateName + "_sum").c_str());
        result._templateHist->SetDirectory(0);
    }

    // Add histograms with matching keys
    for (const auto &[key, hist] : this->_histograms)
    {
        auto it = toAdd._histograms.find(key);
        if (it != toAdd._histograms.end())
        {
            // Ensure binning matches
            if (hist->GetNbinsX() != it->second->GetNbinsX() ||
                hist->GetXaxis()->GetXmin() != it->second->GetXaxis()->GetXmin() ||
                hist->GetXaxis()->GetXmax() != it->second->GetXaxis()->GetXmax())
            {
                throw std::runtime_error("Binning mismatch for histogram key: " + key);
            }

            // Create a new histogram for the sum
            std::string sumName = hist->GetName();
            TH1D *sumHist = (TH1D *)hist->Clone((sumName + "_sum").c_str());
            sumHist->SetDirectory(0);
            sumHist->Add(it->second);
            result._histograms[key] = sumHist;
        }
        else
        {
            // If the key does not exist in the toAdd HistControl, clone the original histogram
            std::string clonedName = hist->GetName();
            TH1D *clonedHist = (TH1D *)hist->Clone((clonedName + "_cloned").c_str());
            clonedHist->SetDirectory(0);
            result._histograms[key] = clonedHist;
        }
    }

    // Add histograms that are only in the toAdd HistControl
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

// fix copy
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
        // Clean up existing resources
        delete _templateHist;
        for (auto &[key, hist] : _histograms)
        {
            delete hist;
        }
        _histograms.clear();

        // Copy resources
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
