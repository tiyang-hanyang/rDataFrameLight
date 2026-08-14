#include "PlotControl.h"
#include "Utility.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

#include "TLatex.h"
#include "TLine.h"
#include "TBox.h"
#include "TGaxis.h"
#include "TStyle.h"

namespace
{
bool hasSuffix(const std::string &value, const std::string &suffix)
{
    if (value.size() < suffix.size())
    {
        return false;
    }
    return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

double getTotalDataIntegral(const std::map<std::string, TH1D *> &hists, const std::map<std::string, int> &isData)
{
    double totalDataIntegral = 0.0;
    for (const auto &[name, hist] : hists)
    {
        if (!hist)
        {
            continue;
        }
        auto iter = isData.find(name);
        if (iter != isData.end() && iter->second)
        {
            totalDataIntegral += hist->Integral();
        }
    }
    return totalDataIntegral;
}

int countDataHistograms(const std::map<std::string, int> &isData)
{
    int count = 0;
    for (const auto &[name, flag] : isData)
    {
        if (flag)
        {
            count++;
        }
    }
    return count;
}

int getDatasetColorOrDefault(const std::map<std::string, int> &colorScheme, const std::string &dataset)
{
    auto iter = colorScheme.find(dataset);
    if (iter != colorScheme.end())
    {
        return iter->second;
    }
    rdfWS_utility::messageERROR("PlotControl", "No configured color for dataset " + dataset + ". Please define it in plot json colorMapping or json/general_config/dataset_color.json.");
}

int getNiceSignalScaleFactor(double referenceValue, double signalValue)
{
    if (referenceValue <= 0.0 || signalValue <= 0.0)
    {
        return 1;
    }

    const double ratio = referenceValue / signalValue;
    if (ratio < 1.0)
    {
        return 1;
    }

    const double exponent = std::floor(std::log10(ratio));
    const double base = std::pow(10.0, exponent);
    for (const double candidate : {5.0, 2.0, 1.0})
    {
        const double scaledCandidate = candidate * base;
        if (scaledCandidate <= ratio)
        {
            return static_cast<int>(scaledCandidate);
        }
    }

    return 1;
}

std::string formatSignalScaleFactor(int scaleFactor)
{
    if (scaleFactor <= 1)
    {
        return "";
    }

    int exponent = 0;
    int mantissa = scaleFactor;
    while (mantissa >= 10 && mantissa % 10 == 0)
    {
        mantissa /= 10;
        ++exponent;
    }

    if (mantissa == 1)
    {
        return " #times 10^{" + std::to_string(exponent) + "}";
    }
    if (mantissa == 2 || mantissa == 5)
    {
        return " #times " + std::to_string(mantissa) + "#times 10^{" + std::to_string(exponent) + "}";
    }
    return " #times " + std::to_string(scaleFactor);
}

double findSmallestPositiveBinContent(const std::map<std::string, TH1D *> &hists)
{
    double minPositive = std::numeric_limits<double>::infinity();
    for (const auto &[name, hist] : hists)
    {
        if (!hist)
        {
            continue;
        }
        for (int i = 1; i <= hist->GetNbinsX(); ++i)
        {
            const double content = hist->GetBinContent(i);
            if (content > 0.0 && content < minPositive)
            {
                minPositive = content;
            }
        }
    }

    if (!std::isfinite(minPositive))
    {
        return 0.0;
    }
    return minPositive;
}

double getNormalizedLogMinimum(double minPositive)
{
    if (minPositive <= 0.0)
    {
        return 1e-4;
    }

    const double exponent = std::floor(std::log10(minPositive)) - 1.0;
    const double axisMin = std::pow(10.0, exponent);
    if (!std::isfinite(axisMin) || axisMin <= 0.0)
    {
        return 1e-4;
    }
    return axisMin;
}

std::string formatEdgeValue(double value)
{
    std::ostringstream stream;
    stream.setf(std::ios::fmtflags(0), std::ios::floatfield);
    stream.precision(12);
    stream << value;
    return stream.str();
}

bool isBlindedValue(double xValue, const std::vector<std::pair<double, double>> &blindRanges)
{
    for (const auto &[rangeMin, rangeMax] : blindRanges)
    {
        if (xValue >= rangeMin && xValue <= rangeMax)
            return true;
    }
    return false;
}

double getOriginalBinCenterForBlindDecision(
    const TH1D *hist,
    int binIndex,
    bool displayEqualWidthBins,
    const std::vector<double> &displayBinEdges)
{
    if (hist == nullptr)
        return 0.0;

    if (displayEqualWidthBins && displayBinEdges.size() == static_cast<size_t>(hist->GetNbinsX()) + 1)
    {
        return 0.5 * (displayBinEdges[binIndex - 1] + displayBinEdges[binIndex]);
    }
    return hist->GetXaxis()->GetBinCenter(binIndex);
}

void maskBlindedBinsForDisplay(
    TH1D *hist,
    const std::vector<std::pair<double, double>> &blindRanges,
    bool displayEqualWidthBins = false,
    const std::vector<double> &displayBinEdges = {})
{
    if (hist == nullptr || blindRanges.empty())
        return;

    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {
        const double binCenter = getOriginalBinCenterForBlindDecision(
            hist,
            i,
            displayEqualWidthBins,
            displayBinEdges);
        if (!isBlindedValue(binCenter, blindRanges))
            continue;
        hist->SetBinContent(i, -1.0);
        hist->SetBinError(i, 0.0);
    }
}

double remapBlindCoordinateToEqualWidth(double xValue, const std::vector<double> &binEdges)
{
    if (binEdges.size() < 2)
        return xValue;

    const int nBins = static_cast<int>(binEdges.size()) - 1;
    if (xValue <= binEdges.front())
        return 0.0;
    if (xValue >= binEdges.back())
        return static_cast<double>(nBins);

    for (int i = 0; i < nBins; ++i)
    {
        const double lowEdge = binEdges[i];
        const double highEdge = binEdges[i + 1];
        if (xValue <= highEdge || i == nBins - 1)
        {
            const double width = highEdge - lowEdge;
            const double fraction = width > 0.0 ? (xValue - lowEdge) / width : 0.0;
            return static_cast<double>(i) + std::clamp(fraction, 0.0, 1.0);
        }
    }
    return static_cast<double>(nBins);
}

std::vector<std::pair<double, double>> buildDisplayBlindRanges(
    const std::vector<std::pair<double, double>> &blindRanges,
    bool displayEqualWidthBins,
    const std::vector<double> &binEdges,
    double axisMin,
    double axisMax)
{
    std::vector<std::pair<double, double>> displayRanges;
    for (const auto &[rangeMin, rangeMax] : blindRanges)
    {
        double xMin = rangeMin;
        double xMax = rangeMax;
        if (displayEqualWidthBins)
        {
            xMin = remapBlindCoordinateToEqualWidth(rangeMin, binEdges);
            xMax = remapBlindCoordinateToEqualWidth(rangeMax, binEdges);
        }
        if (xMax < axisMin || xMin > axisMax)
            continue;
        xMin = std::max(xMin, axisMin);
        xMax = std::min(xMax, axisMax);
        if (xMax <= xMin)
            continue;
        displayRanges.emplace_back(xMin, xMax);
    }
    return displayRanges;
}

std::vector<TBox *> drawBlindBoxes(
    const std::vector<std::pair<double, double>> &blindRanges,
    double yMin,
    double yMax,
    double alpha = 0.28,
    int fillColor = kGray + 1)
{
    std::vector<TBox *> boxes;
    for (const auto &[xMin, xMax] : blindRanges)
    {
        auto *box = new TBox(xMin, yMin, xMax, yMax);
        box->SetFillColorAlpha(fillColor, alpha);
        box->SetLineColor(0);
        box->SetFillStyle(1001);
        box->Draw("SAME");
        boxes.push_back(box);
    }
    return boxes;
}

TH1D *buildEqualWidthDisplayHistogram(const TH1D *sourceHist, const std::string &nameSuffix)
{
    if (sourceHist == nullptr)
    {
        return nullptr;
    }

    const int nBins = sourceHist->GetNbinsX();
    auto *displayHist = new TH1D(
        (std::string(sourceHist->GetName()) + nameSuffix).c_str(),
        sourceHist->GetTitle(),
        nBins,
        0.0,
        static_cast<double>(nBins));
    displayHist->SetDirectory(0);
    displayHist->Sumw2();

    for (int i = 0; i <= nBins + 1; ++i)
    {
        displayHist->SetBinContent(i, sourceHist->GetBinContent(i));
        displayHist->SetBinError(i, sourceHist->GetBinError(i));
    }
    displayHist->SetEntries(sourceHist->GetEntries());

    displayHist->SetMarkerColor(sourceHist->GetMarkerColor());
    displayHist->SetLineColor(sourceHist->GetLineColor());
    displayHist->SetFillColor(sourceHist->GetFillColor());
    displayHist->SetFillStyle(sourceHist->GetFillStyle());
    displayHist->SetLineStyle(sourceHist->GetLineStyle());
    displayHist->SetLineWidth(sourceHist->GetLineWidth());
    displayHist->SetMarkerStyle(sourceHist->GetMarkerStyle());
    displayHist->SetMarkerSize(sourceHist->GetMarkerSize());

    displayHist->GetXaxis()->SetTitle(sourceHist->GetXaxis()->GetTitle());
    displayHist->GetXaxis()->SetLabelFont(sourceHist->GetXaxis()->GetLabelFont());
    displayHist->GetXaxis()->SetLabelSize(sourceHist->GetXaxis()->GetLabelSize());
    displayHist->GetXaxis()->SetLabelOffset(sourceHist->GetXaxis()->GetLabelOffset());
    displayHist->GetXaxis()->SetTitleFont(sourceHist->GetXaxis()->GetTitleFont());
    displayHist->GetXaxis()->SetTitleSize(sourceHist->GetXaxis()->GetTitleSize());
    displayHist->GetXaxis()->SetTitleOffset(sourceHist->GetXaxis()->GetTitleOffset());

    displayHist->GetYaxis()->SetTitle(sourceHist->GetYaxis()->GetTitle());
    displayHist->GetYaxis()->SetLabelFont(sourceHist->GetYaxis()->GetLabelFont());
    displayHist->GetYaxis()->SetLabelSize(sourceHist->GetYaxis()->GetLabelSize());
    displayHist->GetYaxis()->SetLabelOffset(sourceHist->GetYaxis()->GetLabelOffset());
    displayHist->GetYaxis()->SetTitleFont(sourceHist->GetYaxis()->GetTitleFont());
    displayHist->GetYaxis()->SetTitleSize(sourceHist->GetYaxis()->GetTitleSize());
    displayHist->GetYaxis()->SetTitleOffset(sourceHist->GetYaxis()->GetTitleOffset());

    return displayHist;
}

std::map<std::string, TH1D *> buildEqualWidthDisplayHistMap(const std::map<std::string, TH1D *> &sourceHists)
{
    std::map<std::string, TH1D *> displayHists;
    for (const auto &[name, hist] : sourceHists)
    {
        displayHists.emplace(name, buildEqualWidthDisplayHistogram(hist, "_displayEqualWidth"));
    }
    return displayHists;
}

std::vector<std::pair<TH1D *, TH1D *>> buildEqualWidthDisplayPairs(
    const std::vector<std::pair<TH1D *, TH1D *>> &sourcePairs)
{
    std::vector<std::pair<TH1D *, TH1D *>> displayPairs;
    displayPairs.reserve(sourcePairs.size());
    for (size_t index = 0; index < sourcePairs.size(); ++index)
    {
        const auto &[upHist, downHist] = sourcePairs[index];
        TH1D *displayUp = buildEqualWidthDisplayHistogram(upHist, "_displayEqualWidth_up_" + std::to_string(index));
        TH1D *displayDown = buildEqualWidthDisplayHistogram(downHist, "_displayEqualWidth_down_" + std::to_string(index));
        displayPairs.emplace_back(displayUp, displayDown);
    }
    return displayPairs;
}

std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> buildEqualWidthDisplayPairMap(
    const std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> &sourcePairMap)
{
    std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> displayPairMap;
    for (const auto &[name, pairs] : sourcePairMap)
    {
        displayPairMap.emplace(name, buildEqualWidthDisplayPairs(pairs));
    }
    return displayPairMap;
}

void deleteHistogramMap(std::map<std::string, TH1D *> &hists)
{
    for (auto &[name, hist] : hists)
    {
        if (hist != nullptr)
        {
            hist->SetDirectory(0);
            delete hist;
            hist = nullptr;
        }
    }
    hists.clear();
}

void deleteHistogramPairs(std::vector<std::pair<TH1D *, TH1D *>> &pairs)
{
    for (auto &[upHist, downHist] : pairs)
    {
        if (upHist != nullptr)
        {
            upHist->SetDirectory(0);
            delete upHist;
            upHist = nullptr;
        }
        if (downHist != nullptr)
        {
            downHist->SetDirectory(0);
            delete downHist;
            downHist = nullptr;
        }
    }
    pairs.clear();
}

void deleteHistogramPairMap(std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> &pairMap)
{
    for (auto &[name, pairs] : pairMap)
    {
        deleteHistogramPairs(pairs);
    }
    pairMap.clear();
}

void configureEqualWidthAxisStyle(TH1D *hist, bool showTitle)
{
    if (hist == nullptr)
    {
        return;
    }

    hist->GetXaxis()->SetLabelSize(0.0);
    hist->GetXaxis()->SetNdivisions(-hist->GetNbinsX(), false);
    if (showTitle)
    {
        hist->GetXaxis()->SetTitleOffset(0.9);
    }
    else
    {
        hist->GetXaxis()->SetTitleSize(0.0);
        hist->GetXaxis()->SetTitleOffset(0.0);
    }
}

void drawEqualWidthAxisLabels(TVirtualPad *pad, const std::vector<double> &binEdges, double scale, bool isRatio = false)
{
    if (pad == nullptr || binEdges.size() < 2)
    {
        return;
    }

    pad->cd();
    const double leftMargin = pad->GetLeftMargin();
    const double rightMargin = pad->GetRightMargin();
    const double bottomMargin = pad->GetBottomMargin();
    const double usableWidth = 1.0 - leftMargin - rightMargin;
    const int nBins = static_cast<int>(binEdges.size()) - 1;

    const double textSize = (nBins >= 18 ? 0.024 : 0.03) / scale;
    const double labelY = isRatio
        ? std::max(0.06, bottomMargin - 0.85 * textSize)
        : std::max(0.01, bottomMargin - 0.35 * textSize);
    const double edgeInset = isRatio ? 0.015 : 0.008;

    for (int edgeIndex = 0; edgeIndex <= nBins; ++edgeIndex)
    {
        double xNdc = leftMargin + usableWidth * (static_cast<double>(edgeIndex) / nBins);
        if (edgeIndex == 0)
            xNdc += edgeInset;
        else if (edgeIndex == nBins)
            xNdc -= edgeInset;
        auto *label = new TLatex();
        label->SetNDC();
        label->SetTextFont(42);
        label->SetTextSize(textSize);
        label->SetTextAlign(21);
        label->DrawLatex(xNdc, labelY, formatEdgeValue(binEdges[edgeIndex]).c_str());
    }
}

TGraphAsymmErrors *buildHistogramUncertaintyBand(
    TH1D *nominalHist,
    const std::vector<std::pair<TH1D *, TH1D *>> &shapeUncerts,
    bool includeMCStat)
{
    if (nominalHist == nullptr)
    {
        return nullptr;
    }
    if (!includeMCStat && shapeUncerts.empty())
    {
        return nullptr;
    }

    TGraphAsymmErrors *systBand = new TGraphAsymmErrors(nominalHist);
    for (int i = 1; i <= nominalHist->GetNbinsX(); ++i)
    {
        const double binCenter = nominalHist->GetBinCenter(i);
        const double binWidth = nominalHist->GetBinWidth(i) / 2.0;
        const double nominalVal = nominalHist->GetBinContent(i);

        double errUp2 = 0.0;
        double errDown2 = 0.0;
        if (includeMCStat)
        {
            const double statErr = nominalHist->GetBinError(i);
            errUp2 += statErr * statErr;
            errDown2 += statErr * statErr;
        }

        auto addEnvelope = [&](double upVal, double downVal)
        {
            const double upDelta = upVal - nominalVal;
            const double downDelta = downVal - nominalVal;
            const double upErr = std::max({0.0, upDelta, downDelta});
            const double downErr = std::max({0.0, -upDelta, -downDelta});
            errUp2 += upErr * upErr;
            errDown2 += downErr * downErr;
        };

        for (const auto &[shapeUp, shapeDown] : shapeUncerts)
        {
            if (shapeUp == nullptr || shapeDown == nullptr)
            {
                continue;
            }
            addEnvelope(shapeUp->GetBinContent(i), shapeDown->GetBinContent(i));
        }

        systBand->SetPoint(i - 1, binCenter, nominalVal);
        systBand->SetPointError(i - 1, binWidth, binWidth, std::sqrt(errDown2), std::sqrt(errUp2));
    }

    systBand->SetFillColor(nominalHist->GetLineColor());
    systBand->SetFillStyle(3004);
    systBand->SetLineWidth(0);
    systBand->SetMarkerSize(0);
    return systBand;
}
}

PlotControl::PlotControl(const std::string &name) : _controllerName(name)
{
    this->_canvas = nullptr;
    this->_abovePad = nullptr;
    this->_belowPad = nullptr;
    this->_scale = {1.0, 1.0};
}

PlotControl::~PlotControl()
{
    // clear legends
    for (auto len : this->_legends)
    {
        if (len)
        {
            len->Clear();
            delete len;
        }
    }
    this->_legends.clear();

    // delete canvas would be enough, do not delete tpad manually, it crahsed
    delete this->_canvas;
}

// setup canvas and pad style, record the scale if do ratio plot
void PlotControl::setHanyangCanvas(double xSize, double ySize, int doLog, int doRatio)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetEndErrorSize(0);
    gStyle->SetErrorX(0);

    if (this->_canvas != nullptr)
    {
        rdfWS_utility::messageERROR("PlotControl", "Internal canvas already exist");
    }
    std::string canvasName = "c_" + this->_controllerName;
    this->_canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), xSize, ySize);

    // margin setup
    double tMargin = 0.08;
    double bMargin = 0.15;
    this->_topMargin = tMargin;
    this->_bottomMargin = bMargin;

    // ticks and grid (if log)
    this->_canvas->SetTickx(1);
    this->_canvas->SetTicky(1);
    if (doLog)
    {
        this->_canvas->SetLogy();
        // Not want grid even in log scale
        // this->_canvas->SetGrid();
    }

    // margin scheme
    this->_canvas->SetLeftMargin(0.15);
    this->_canvas->SetRightMargin(0.05);
    this->_canvas->SetTopMargin(tMargin);
    this->_canvas->SetBottomMargin(bMargin);

    // when do ratio
    std::vector<TBox *> ratioBlindBoxes;
    if (doRatio)
    {
        double aRatio = 0.75;
        double bRatio = 1.0 - aRatio;

        // record scale infor
        double splitingLine = bMargin + (1.0 - bMargin - tMargin) * bRatio;
        this->_scale = {1 - splitingLine, splitingLine};

        std::string apName = "pHist_" + this->_controllerName;
        this->_abovePad = new TPad(apName.c_str(), apName.c_str(), 0.0, splitingLine, 1.0, 1.0);
        this->_abovePad->SetTickx(1);
        this->_abovePad->SetTicky(1);
        this->_abovePad->SetLeftMargin(0.15);
        this->_abovePad->SetRightMargin(0.05);
        this->_abovePad->SetTopMargin(tMargin / (1.0 - splitingLine));
        this->_abovePad->SetBottomMargin(0.01 / (1.0 - splitingLine));
        if (doLog)
        {
            this->_abovePad->SetLogy();
            // this->_abovePad->SetGrid();
        }

        std::string bpName = "pRatio_" + this->_controllerName;
        this->_belowPad = new TPad(bpName.c_str(), bpName.c_str(), 0.0, 0.0, 1.0, splitingLine);
        this->_belowPad->SetTickx(1);
        this->_belowPad->SetTicky(1);
        this->_belowPad->SetLeftMargin(0.15);
        this->_belowPad->SetRightMargin(0.05);
        this->_belowPad->SetTopMargin(0.01 / splitingLine);
        this->_belowPad->SetBottomMargin(bMargin / splitingLine);
    }
}

// setup the style of histogram
void PlotControl::setHanyangHist(TH1D *hist, int color, int isData, const std::vector<std::string> &binLabels, double scale, std::string xTitle, std::string yTitle, bool isRatio)
{
    hist->SetTitle("");
    hist->SetStats(0);
    // color
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);

    // bin labels
    // auto nBins = hist->GetNbinsX();
    // if (binLabels.size() == nBins)
    // {
    //     for (int i = 0; i < nBins; ++i)
    //     {
    //         hist->GetXaxis()->SetBinLabel(i + 1, binLabels[i].c_str());
    //     }
    // }
    if (binLabels.size() > 0)
    {
        auto nBins = binLabels.size();
        for (int i = 0; i < nBins; ++i)
        {
            hist->GetXaxis()->SetBinLabel(i + 1, binLabels[i].c_str());
        }
    }


    // data has marker style while has MC line style
    if (isData)
    {
        // data draw option is ep
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.2);
        hist->SetLineStyle(0);
        hist->SetLineWidth(1);
    }
    else
    {
        // MC draw option is e hist
        hist->SetLineWidth(2);
        hist->SetLineStyle(1);
    }

    // axis style setting
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.03 / scale);
    hist->GetXaxis()->SetLabelOffset(0.007 / scale);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleSize(0.04 / scale);
    hist->GetXaxis()->SetTitleOffset(0.9);

    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelSize(0.03 / scale);
    hist->GetYaxis()->SetLabelOffset(0.007);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleSize(0.04 / scale);
    hist->GetYaxis()->SetTitleOffset(1.5 * scale);
    if (isRatio)
        hist->GetYaxis()->SetNdivisions(4, false);
}

TLegend *PlotControl::setHanyangLegend(int entries, double textSize, double scale,
                                        bool ratioLayout, bool ratioTopLayout)
{
    if (ratioLayout)
    {
        const int nColumns = std::max(1, static_cast<int>(std::ceil(entries / 6.0)));
        const int nRows = std::max(1, (entries + nColumns - 1) / nColumns);
        const double tunedTextSize = entries > 12 ? textSize * 0.78 : (entries > 6 ? textSize * 0.90 : textSize);
        const double x2 = 0.91;
        const double x1 = nColumns >= 3 ? 0.54 : (nColumns == 2 ? 0.64 : 0.72);
        const double y2 = 0.84;
        const double legendHeight = (0.018 + nRows * tunedTextSize * 0.92) / scale;
        const double y1 = std::max(0.52, y2 - legendHeight);

        auto len = new TLegend(x1, y1, x2, y2, NULL, "brNDC");
        len->SetNColumns(nColumns);
        len->SetBorderSize(0);
        len->SetTextFont(42);
        len->SetTextSize(tunedTextSize / scale);
        len->SetLineColor(0);
        len->SetLineStyle(1);
        len->SetLineWidth(1);
        len->SetFillColor(0);
        len->SetFillStyle(0);
        len->SetMargin(nColumns > 1 ? 0.16 : 0.22);
        len->SetColumnSeparation(nColumns >= 3 ? 0.035 : 0.02);
        return len;
    }

    if (ratioTopLayout)
    {
        // Keep the ratio upper-pad legend in its original coordinate system.
        // The non-ratio layout is intentionally tuned separately below.
        const int nColumns = std::max(1, static_cast<int>(std::ceil(entries / 6.0)));
        const int nRows = std::max(1, (entries + nColumns - 1) / nColumns);
        const double tunedTextSize = entries > 12 ? textSize * 0.78 : (entries > 6 ? textSize * 0.90 : textSize);
        const double x2 = 0.91;
        const double x1 = nColumns >= 3 ? 0.54 : (nColumns == 2 ? 0.64 : 0.72);
        const double y2 = 0.84;
        const double legendHeight = (0.018 + nRows * tunedTextSize * 0.92) / scale;
        const double y1 = std::max(0.52, y2 - legendHeight);

        auto len = new TLegend(x1, y1, x2, y2, NULL, "brNDC");
        len->SetNColumns(nColumns);
        len->SetBorderSize(0);
        len->SetTextFont(42);
        len->SetTextSize(tunedTextSize / scale);
        len->SetLineColor(0);
        len->SetLineStyle(1);
        len->SetLineWidth(1);
        len->SetFillColor(0);
        len->SetFillStyle(0);
        len->SetMargin(nColumns > 1 ? 0.16 : 0.22);
        len->SetColumnSeparation(nColumns >= 3 ? 0.035 : 0.02);
        return len;
    }

    // Prefer wider columns for long labels. Three narrow columns make labels
    // such as nonprompt DD and scaled signal names overlap.
    const int nColumns = entries > 16 ? 3 : (entries > 6 ? 2 : 1);
    const int nRows = std::max(1, (entries + nColumns - 1) / nColumns);
    const double tunedTextSize = entries > 12 ? textSize * 0.80 : (entries > 6 ? textSize * 0.90 : textSize);
    const double x2 = 0.96;
    const double x1 = nColumns >= 3 ? 0.42 : (nColumns == 2 ? 0.48 : 0.69);
    const double y2 = 0.86;
    const double legendHeight = (0.018 + nRows * tunedTextSize * 1.05) / scale;
    const double y1 = std::max(0.45, y2 - legendHeight);

    auto len = new TLegend(x1, y1, x2, y2, NULL, "brNDC");
    len->SetNColumns(nColumns);
    len->SetBorderSize(0);
    len->SetTextFont(42);
    len->SetTextSize(tunedTextSize / scale);
    len->SetLineColor(0);
    len->SetLineStyle(1);
    len->SetLineWidth(1);
    len->SetFillColor(0);
    len->SetFillStyle(0);
    len->SetMargin(nColumns > 1 ? 0.10 : 0.18);
    len->SetColumnSeparation(nColumns >= 3 ? 0.045 : 0.035);
    return len;
}

void PlotControl::setMax(std::map<std::string, TH1D *> hists, int doLog, int isRatio)
{
    double maxVal(0.0), minVal(1e9);
    if (isRatio)
    {
        for (const auto &[histName, hist] : hists)
        {
            if (!hist)
            {
                rdfWS_utility::messageERROR("PlotControl", "Histogram " + histName + " is null.");
            }

            if (hist->GetMaximum() > maxVal)
                maxVal = hist->GetMaximum();
            if (hist->GetMinimum() < minVal)
                minVal = hist->GetMinimum();
        }
        // maxVal *= 1.2;
        // minVal *= 0.8;

        // take the fixed default ratio range
        maxVal = 1.5;
        minVal = 0.5;
        
    }
    else
    {
        for (const auto &[histName, hist] : hists)
        {
            if (!hist)
            {
                rdfWS_utility::messageERROR("PlotControl", "Histogram " + histName + " is null.");
            }

            if (hist->GetMaximum() > maxVal)
                maxVal = hist->GetMaximum();
        }
        if (doLog)
        {
            maxVal *= 50000.0;
            minVal = 1.0;
        }
        else
        {
            maxVal *= 1.5;
            minVal = 0.0;
        }
    }

    for (const auto &[histName, hist] : hists)
    {
        hist->SetMaximum(maxVal);
        hist->SetMinimum(minVal);
    }
}

// Nov 24:
// should not set maximum here, as this is done before scaling 
// for ratio, can set maximum and minimum values here
std::map<std::string, TH1D *> PlotControl::setupHists(std::map<std::string, TH1D *> hists, PlotContext setup, const std::map<std::string, int> &colorScheme, const std::vector<std::string> &binLabels, int isRatio)
{
    // std::vector<int> colorScheme = {2, 3, 4, 5, 6, 8, 9, 11, 28, 42};
    int colorIndex(0);

    std::map<std::string, TH1D *> styledHists;
    const int nDataHists = countDataHistograms(setup.isData[isRatio]);
    // Clone and style each hist, not affect the original ones
    for (const auto &[histChannel, hist] : hists)
    {
        if (!hist)
            rdfWS_utility::messageERROR("PlotControl", "Histogram " + histChannel + " is null.");
        std::string newHistName = std::string("plot_") + hist->GetName();
        TH1D *clonedHist = dynamic_cast<TH1D *>(hist->Clone(newHistName.c_str()));
        clonedHist->SetDirectory(0);

        // Apply styles (e.g., isData, labels)
        // add check for if isData set
        int isData = setup.isData[isRatio].find(histChannel) != setup.isData[isRatio].end() && setup.isData[isRatio].at(histChannel);
        if (isData)
        {
            if (nDataHists > 1 && !isRatio && colorScheme.find(histChannel) != colorScheme.end())
                setHanyangHist(clonedHist, getDatasetColorOrDefault(colorScheme, histChannel), isData, binLabels, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio], isRatio);
            else
                setHanyangHist(clonedHist, 1, isData, binLabels, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio], isRatio);
        }
        else
        {
            rdfWS_utility::messageINFO("PlotControl", "Setup hist " + histChannel);
            if (!isRatio)
            {
                setHanyangHist(clonedHist, getDatasetColorOrDefault(colorScheme, histChannel), isData, binLabels, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio], false);
            }
            else
                setHanyangHist(clonedHist, 1, isData, binLabels, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio], true);
            // setHanyangHist(clonedHist, colorScheme[colorIndex], isData, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio]);
            // colorIndex++;
        }

        if (hasSuffix(histChannel, "__nominal"))
        {
            clonedHist->SetLineWidth(3);
            clonedHist->SetLineStyle(1);
        }
        else if (hasSuffix(histChannel, "__syst_a_up") || hasSuffix(histChannel, "__syst_b_down"))
        {
            clonedHist->SetLineWidth(2);
            clonedHist->SetLineStyle(1);
        }

        styledHists.emplace(histChannel, clonedHist);
    }

    // setup maximum value according to all hists
    if (isRatio)
        setMax(styledHists, setup.doLog, isRatio);

    return styledHists;
}

// function to containing all hist drawing for non-stack histograms
// the input hists should be already striped out stacked part
// same option is for if the first histogram need option same
void PlotControl::drawNonStackedHists(
    std::map<std::string, TH1D *> &hists,
    std::map<std::string, int> &isData,
    int same, int isRatio, float mcScaling)
{
    // reserve a place for the dataHist, so that I would always plot this last
    // std::vector<TH1D *> dataHists;
    std::vector<std::string> dataHists;
    for (auto &[histName, hist] : hists)
    {
        if (isData[histName])
            dataHists.push_back(histName);
        else
        {
            hist->Scale(mcScaling);
            if (isRatio)
            {
                //if (hist->GetMinimum() > 1.0) hist->SetMinimum(0.8);
                //if (hist->GetMaximum() < 1.0) hist->SetMaximum(1.2);

                // force setting a fixed min and max
                // hist->SetMinimum(0.8);
                // hist->SetMaximum(1.2);
                hist->SetMinimum(0.5);
                hist->SetMaximum(1.5);
            }
            else
            {
                hist->SetMaximum(hist->GetMaximum() * mcScaling);
            }
            std::string drawOption = "HIST E";
            // if (isRatio)
            // {
            //     drawOption += " E";
            // }
            if (same)
            {
                drawOption += " SAME";
            }

            hist->Draw(drawOption.c_str());
            std::cout << histName << " contribution: " << hist->Integral() << std::endl;
            same = true; // Ensure subsequent histograms use "SAME"
        }
    }
    for (auto histName : dataHists)
    {
        std::string drawOption = "E0P";
        if (same)
        {
            drawOption += " SAME";
        }

        hists[histName]->Draw(drawOption.c_str());
        std::cout << histName << " contribution: " << hists[histName]->Integral() << std::endl;
        same = true;
    }
}

// temp: not design every text drawing together yet
// Now, I would want to move the heading to the above
// \bf{CMS} \it{Preliminary} above the panel, to the left
// %1.0f fb^{-1} (13.6 TeV) above the panel, to the right
void PlotControl::drawHeader(const std::vector<std::string>& header,double scale, int doRatio)
{
    if (doRatio)
        this->_abovePad->cd();

    double CMSsize = 0.04;
    double defaultSize = 0.035;

    TLatex headerText;
    headerText.SetNDC();
    headerText.SetTextFont(42);

    // Draw header first
    if (doRatio) {
        headerText.SetTextSize(CMSsize / scale);
        headerText.DrawLatex(0.15, 0.9, "#bf{CMS}");
        headerText.SetTextSize(defaultSize / scale);
        std::string CMSheader = "#it{"+header[0]+"}";
        headerText.DrawLatex(0.25, 0.9, CMSheader.c_str());
        headerText.DrawLatex(0.65, 0.9, header[1].c_str());
    }
    else
    {
        headerText.SetTextSize(CMSsize / scale);
        headerText.DrawLatex(0.15, 0.95, "#bf{CMS}");
        headerText.SetTextSize(defaultSize / scale);
        std::string CMSheader = "#it{"+header[0]+"}";
        headerText.DrawLatex(0.25, 0.95, CMSheader.c_str());
        headerText.DrawLatex(0.68, 0.95, header[1].c_str());
    }
}

// latex
void PlotControl::drawTexts(const std::vector<std::string> &texts, double scale)
{
    double defaultSize = 0.028;

    TLatex padText;
    padText.SetNDC();
    padText.SetTextFont(42);
    padText.SetTextSize(defaultSize / scale);

    // Then draw the above text
    // double initHeight = 1 - (1 + defaultSize - 0.88) * scale;
    double initHeight = 1 - (1 - 0.87) / scale;
    for (const auto &t : texts)
    {
        padText.DrawLatex(0.2, initHeight, t.c_str());
        initHeight -= 0.047 / scale;
    }
}

// re-implement
void PlotControl::drawCanvasHists(
    std::map<std::string, TH1D *> &plotHists,
    double scale,
    std::map<std::string, int> &isData,
    int doLegend,
    int isRatio,
    std::vector<std::string> plotTexts,
    float mcScaling)
{
    // draw all histograms not in stack
    drawNonStackedHists(plotHists, isData, 0, isRatio, mcScaling);
    // non-stacked version of legend
    if (doLegend)
    {
        rdfWS_utility::messageINFO("PlotControl", "legend plot hist size: " + std::to_string(plotHists.size()));
        rdfWS_utility::messageINFO("PlotControl", "legend scale: " + std::to_string(scale));
        auto len = setHanyangLegend(plotHists.size(), 0.03, scale, true);
        const int nDataHists = countDataHistograms(isData);
        for (auto &[histName, hist] : plotHists)
        {
            if (isData.find(histName) != isData.end() && isData.at(histName))
            {
                if (nDataHists == 1)
                    len->AddEntry(hist, "Data");
                else
                    len->AddEntry(hist, histName.c_str());
            }
            else
            {
                len->AddEntry(hist, histName.c_str());
            }
        }

        len->Draw("same");
        this->_legends.push_back(len);
    }
    // Draw texts
    drawTexts(plotTexts, scale);
}

void PlotControl::saveCanvas(const std::string &fileName, int doLog)
{
    if (!this->_canvas)
    {
        rdfWS_utility::messageERROR("PlotControl", "Internal canvas not properly setup: " + this->_controllerName);
    }

    std::string fileDir = fileName.substr(0, fileName.rfind("/"));
    std::string outFileName = fileName.substr(fileName.rfind("/")+1);
    if (doLog) 
        outFileName = "log_"+outFileName;
    else
        outFileName = "linear_"+outFileName;
    outFileName = fileDir+"/"+outFileName;
    this->_canvas->SaveAs((outFileName + ".png").c_str());
    this->_canvas->SaveAs((outFileName + ".pdf").c_str());
    this->_canvas->SaveAs((outFileName + ".eps").c_str());
}

// prepare THStack
THStack *PlotControl::prepareStackHists(std::map<std::string, TH1D *> &hists, std::vector<std::string> &stackOrder, std::map<std::string, int> isData, int reOrder, int doNormalize, float mcScaling)
{
    THStack *stack = new THStack(("stack_" + _controllerName).c_str(), "Stacked Histograms");
    std::vector<int> colorScheme = {2, 3, 4, 5, 6, 8, 9, 12, 28, 42};
    int colorIndex = 0;

    std::vector<std::string> sortedOrder = {};
    std::map<std::string, float> histIntegrals = {};
    for (const auto &name : stackOrder)
    {
        if (hists.find(name) == hists.end())
            continue;
        hists[name]->SetMaximum(hists[name]->GetMaximum() * mcScaling);
        hists[name]->Scale(mcScaling);
        float histInt = hists[name]->Integral();
        if (histInt > 1e-6)
        {
            sortedOrder.push_back(name);
            histIntegrals.emplace(name, histInt);
        }
    }

    // auto ordering
    if (reOrder)
    {
        std::sort(sortedOrder.begin(), sortedOrder.end(), [&histIntegrals](const std::string &a, const std::string &b)
                  { return histIntegrals.at(a) < histIntegrals.at(b); });
    }
    stackOrder = sortedOrder;

    // normalization of stack hist
    if (doNormalize)
    {
        float total = std::accumulate(histIntegrals.begin(), histIntegrals.end(), 0.0f, [](float sum, const std::pair<const std::string, float> &p)
                                      { return sum + p.second; });
        for (const auto &name : sortedOrder)
        {
            auto hist = hists.at(name);
            hist->Scale(1.0 / total);
        }
    }

    std::vector<int> reverseColorScheme;
    auto nStack = stackOrder.size();
    // for (size_t i = 0; i < nStack; i++)
    // {
    //     reverseColorScheme.push_back(colorScheme[nStack - i - 1]);
    // }

    for (const auto &name : sortedOrder)
    {
        auto hist = hists.at(name);
        if (!isData.at(name))
        {
            // hist->SetFillColor(reverseColorScheme[colorIndex]);
            hist->SetFillColor(hist->GetLineColor());
            hist->SetLineColor(1);
            hist->SetLineWidth(0);
            hist->SetFillStyle(1001);
            stack->Add(hist);
            colorIndex++;
        }
    }

    std::cout << "inside the stack, debugging 5" << std::endl;

    return stack;
}

TGraphAsymmErrors *uncertaintyHist(
    THStack *stackHist,
    std::map<std::string, TH1D *> &hists,
    const std::vector<std::string> &stackOrder,
    std::map<std::string, double> stackUp,
    std::map<std::string, double> stackDown,
    const std::vector<std::pair<TH1D *, TH1D *>> &shapeUncerts)
{
    if (stackOrder.size() == 0 || !(stackHist))
    {
        return nullptr;
    }

    const bool hasGlobalSyst = stackUp.size() > 0 && stackDown.size() > 0;
    const bool hasShapeSyst = !shapeUncerts.empty();
    if (!hasGlobalSyst && !hasShapeSyst)
        return nullptr;

    TH1D *histUp = nullptr;
    TH1D *histDown = nullptr;
    if (hasGlobalSyst)
    {
        histUp = dynamic_cast<TH1D *>(hists[stackOrder[0]]->Clone("stack_Up"));
        histDown = dynamic_cast<TH1D *>(hists[stackOrder[0]]->Clone("stack_Down"));
        if (stackUp.find(stackOrder[0]) != stackUp.end())
            histUp->Scale(stackUp[stackOrder[0]]);
        if (stackDown.find(stackOrder[0]) != stackDown.end())
            histDown->Scale(stackDown[stackOrder[0]]);

        // These are global/channel-level normalisation factors. Shape/weight
        // variations are handled separately below and combined per bin.
        for (auto iter = stackOrder.begin() + 1; iter != stackOrder.end(); iter++)
        {
            double scaleUp = 1.0;
            double scaleDown = 1.0;
            if (stackUp.find(*iter) != stackUp.end())
                scaleUp = stackUp.at(*iter);
            if (stackDown.find(*iter) != stackDown.end())
                scaleDown = stackDown.at(*iter);
            histUp->Add(hists.at(*iter), scaleUp);
            histDown->Add(hists.at(*iter), scaleDown);
        }
    }

    TH1D *nominalHist = dynamic_cast<TH1D *>(stackHist->GetStack()->Last()->Clone("stackNom"));
    TGraphAsymmErrors *systBand = new TGraphAsymmErrors(nominalHist);
    for (int i = 1; i < nominalHist->GetNbinsX() + 1; i++)
    {
        auto binCenter = nominalHist->GetBinCenter(i);
        auto binWidth = nominalHist->GetBinWidth(i) / 2.0;
        auto nominalVal = nominalHist->GetBinContent(i);
        double errUp2 = 0.0;
        double errDown2 = 0.0;

        auto addEnvelope = [&](double upVal, double downVal)
        {
            const double upDelta = upVal - nominalVal;
            const double downDelta = downVal - nominalVal;
            const double upErr = std::max({0.0, upDelta, downDelta});
            const double downErr = std::max({0.0, -upDelta, -downDelta});
            errUp2 += upErr * upErr;
            errDown2 += downErr * downErr;
        };

        if (hasGlobalSyst)
            addEnvelope(histUp->GetBinContent(i), histDown->GetBinContent(i));
        for (const auto &[shapeUp, shapeDown] : shapeUncerts)
        {
            if (shapeUp == nullptr || shapeDown == nullptr)
                continue;
            addEnvelope(shapeUp->GetBinContent(i), shapeDown->GetBinContent(i));
        }

        const double errDown = std::sqrt(errDown2);
        const double errUp = std::sqrt(errUp2);
        systBand->SetPoint(i - 1, binCenter, nominalVal);
        systBand->SetPointError(i - 1, binWidth, binWidth, errDown, errUp);
    }

    systBand->SetFillColor(15);
    systBand->SetFillStyle(3004);
    systBand->SetLineWidth(0);
    systBand->SetMarkerSize(0);

    delete nominalHist;
    delete histUp;
    delete histDown;

    return systBand;
}

TGraphAsymmErrors *getRatioUncert(TGraphAsymmErrors *origUncert)
{
    if (origUncert == nullptr)
        return nullptr;

    int nPoints = origUncert->GetN();
    TGraphAsymmErrors *ratioUncert = new TGraphAsymmErrors(nPoints);
    for (int i = 0; i < nPoints; ++i)
    {
        double x, y;
        origUncert->GetPoint(i, x, y);
        double exL = origUncert->GetErrorXlow(i);
        double exH = origUncert->GetErrorXhigh(i);
        double eyL = origUncert->GetErrorYlow(i);
        double eyH = origUncert->GetErrorYhigh(i);

        if (y != 0)
        {
            ratioUncert->SetPoint(i, x, 1.0);
            ratioUncert->SetPointError(i, exL, exH, eyL / y, eyH / y);
        }
        else
        {
            ratioUncert->SetPoint(i, x, 1.0);
            ratioUncert->SetPointError(i, exL, exH, 0, 0);
        }
    }

    ratioUncert->SetFillColor(15);
    ratioUncert->SetFillStyle(3004);
    ratioUncert->SetLineWidth(0);
    ratioUncert->SetMarkerSize(0);

    return ratioUncert;
}

void PlotControl::drawStackHistWithRatio(
    const std::map<std::string, TH1D *> &hists,
    const std::vector<std::string> &stackOrder,
    const std::map<std::string, double> &stackUp,
    const std::map<std::string, double> &stackDown,
    int reOrder,
    const std::map<std::string, TH1D *> &ratioHists,
    PlotContext setup,
    float mcScaling,
    const std::map<std::string, int> &colorScheme,
    const std::map<std::string, std::string> &labels,
    const std::vector<std::string> &headerTexts,
    const std::vector<std::string> &aboveTexts,
    const std::vector<std::string> &belowTexts,
    const std::vector<std::pair<TH1D *, TH1D *>> &shapeUncerts,
    const std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> &nonStackShapeUncerts,
    bool includeNonStackMCStat,
    const std::vector<std::string> &binLabels)
{
    // checking input compatibility
    if (hists.empty())
    {
        rdfWS_utility::messageERROR("PlotControl", "No histograms provided for plotting.");
    }

    int doStack = 1;
    if (stackOrder.empty())
    {
        rdfWS_utility::messageINFO("PlotControl", "No stack order list provided.");
        doStack = 0;
    }

    int doRatio = 1;
    if (ratioHists.empty())
    {
        rdfWS_utility::messageINFO("PlotControl", "No ratio histogram provided, will only plot the above pad.");
        doRatio = 0;
    }

    for (const auto &name : stackOrder)
    {
        if (hists.find(name) == hists.end())
        {
            rdfWS_utility::messageERROR("PlotControl", "Histogram '" + name + "' specified in stackOrder does not exist in the provided histograms.");
        }
        std::cout << name << " contribution: " << hists.at(name)->Integral() << std::endl;
    }

    // need to setup canvas in advance for scales needed by hist
    setHanyangCanvas(setup.xSize, setup.ySize, setup.doLog, doRatio);

    // prepare new hists with proper style setup from cloned hists
    std::map<std::string, TH1D *> plotHists = setupHists(hists, setup, colorScheme, binLabels, 0);
    std::vector<std::pair<TH1D *, TH1D *>> displayShapeUncerts = shapeUncerts;
    std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> displayNonStackShapeUncerts = nonStackShapeUncerts;
    if (setup.displayEqualWidthBins)
    {
        auto remappedPlotHists = buildEqualWidthDisplayHistMap(plotHists);
        deleteHistogramMap(plotHists);
        plotHists = std::move(remappedPlotHists);
        displayShapeUncerts = buildEqualWidthDisplayPairs(shapeUncerts);
        displayNonStackShapeUncerts = buildEqualWidthDisplayPairMap(nonStackShapeUncerts);
    }
    std::map<std::string, int> isData = setup.isData[0];
    std::vector<std::string> isSignal = setup.isSignal;

    std::vector<std::string> reStackOrder = stackOrder;
    THStack *stack = prepareStackHists(plotHists, reStackOrder, isData, reOrder, setup.doNormalize, mcScaling);

    // Add stack histogram uncertainties
    auto systStack = uncertaintyHist(stack, plotHists, reStackOrder, stackUp, stackDown, displayShapeUncerts);
    const bool signalScaledForDisplay = !setup.doNormalize;
    const double totalDataIntegral = getTotalDataIntegral(plotHists, isData);
    const bool hasDataReference = totalDataIntegral > 0.0;

    // to make sure the stacked hist maximum also covers the non-stack hist range
    double nonStackMax = stack->GetMaximum();
    std::map<std::string, TH1D *> nonStackedHists;
    std::map<std::string, std::vector<std::pair<TH1D *, TH1D *>>> styledNonStackShapeUncerts;
    std::map<std::string, int> signalDisplayScaleFactors;
    for (auto &[name, hist] : plotHists)
    {
        if (std::find(stackOrder.begin(), stackOrder.end(), name) == stackOrder.end())
        {
            auto uncertIter = displayNonStackShapeUncerts.find(name);
            if (uncertIter != displayNonStackShapeUncerts.end())
            {
                styledNonStackShapeUncerts.emplace(name, uncertIter->second);
            }
            if (setup.doNormalize)
            {
                const double integral = hist->Integral();
                if (integral > 0.0)
                {
                    hist->Scale(1.0 / integral);
                }
                auto styledIter = styledNonStackShapeUncerts.find(name);
                if (styledIter != styledNonStackShapeUncerts.end())
                {
                    for (auto &[upHist, downHist] : styledIter->second)
                    {
                        if (upHist != nullptr && upHist->Integral() > 0.0)
                            upHist->Scale(1.0 / upHist->Integral());
                        if (downHist != nullptr && downHist->Integral() > 0.0)
                            downHist->Scale(1.0 / downHist->Integral());
                    }
                }
            }
            else if (std::find(isSignal.begin(), isSignal.end(), name) != isSignal.end())
            {
                const double referenceValue = hasDataReference ? totalDataIntegral : stack->GetMaximum();
                const double signalReferenceValue = hasDataReference ? hist->Integral() : hist->GetMaximum();
                const int signalScaleFactor = getNiceSignalScaleFactor(referenceValue, signalReferenceValue);
                hist->Scale(signalScaleFactor);
                hist->SetLineStyle(9);
                hist->SetLineWidth(2);
                signalDisplayScaleFactors.emplace(name, signalScaleFactor);
                auto styledIter = styledNonStackShapeUncerts.find(name);
                if (styledIter != styledNonStackShapeUncerts.end())
                {
                    for (auto &[upHist, downHist] : styledIter->second)
                    {
                        if (upHist != nullptr)
                            upHist->Scale(signalScaleFactor);
                        if (downHist != nullptr)
                            downHist->Scale(signalScaleFactor);
                    }
                }
            }
            nonStackedHists.emplace(name, hist);
            if (hist->GetMaximum() > nonStackMax) nonStackMax = hist->GetMaximum();
        }
    }

    if (!setup.blindRanges.empty())
    {
        for (const auto &[name, hist] : nonStackedHists)
        {
            auto isDataIter = isData.find(name);
            if (isDataIter != isData.end() && isDataIter->second)
                maskBlindedBinsForDisplay(
                    hist,
                    setup.blindRanges,
                    setup.displayEqualWidthBins,
                    setup.displayBinEdges);
        }
    }

    std::map<std::string, TGraphAsymmErrors *> nonStackSystBands;
    for (const auto &[name, variations] : styledNonStackShapeUncerts)
    {
        auto nominalIter = nonStackedHists.find(name);
        if (nominalIter == nonStackedHists.end())
        {
            continue;
        }
        auto *band = buildHistogramUncertaintyBand(nominalIter->second, variations, includeNonStackMCStat);
        if (band != nullptr)
        {
            band->SetFillStyle(3005);
            band->SetFillColor(nominalIter->second->GetLineColor());
            nonStackSystBands.emplace(name, band);
        }
    }

    // generate an empty histogram for the axis style
    TH1D *axisHist = dynamic_cast<TH1D *>(plotHists.begin()->second->Clone((this->_controllerName + "_axisHist").c_str()));
    axisHist->SetDirectory(0);
    axisHist->SetStats(0);
    if (!axisHist)
    {
        rdfWS_utility::messageERROR("PlotControl", "Failed to clone reference axis histogram.");
    }
    axisHist->Reset();
    if (setup.displayEqualWidthBins)
    {
        configureEqualWidthAxisStyle(axisHist, !doRatio);
    }
    if (setup.doLog)
    {
        if (setup.doNormalize)
        {
            const double minPositive = findSmallestPositiveBinContent(plotHists);
            const double currentMax = std::max(axisHist->GetMaximum(), nonStackMax);
            axisHist->SetMaximum(std::max(100.0, currentMax * 100.0));
            axisHist->SetMinimum(getNormalizedLogMinimum(minPositive));
        }
        else
        {
            axisHist->SetMaximum(std::max(axisHist->GetMaximum(), nonStackMax) * 5000.);
            axisHist->SetMinimum(1.);
        }
    }
    else
    {
        axisHist->SetMaximum(std::max(axisHist->GetMaximum(), nonStackMax) * 1.5);
        axisHist->SetMinimum(0.);
    }

    // plot
    this->_canvas->cd();
    if (doRatio)
    {
        this->_abovePad->Draw();
        this->_abovePad->cd();
        TGaxis::SetExponentOffset(-0.08, 0.02, "y");
    }
    else
    {
        TGaxis::SetExponentOffset(-0.08, 0.02, "y");
    }
    axisHist->Draw();
    if (doStack)
        stack->Draw("HIST SAME");
    if (systStack != nullptr)
    {
        rdfWS_utility::messageINFO("PlotControl", "Draw with systs in stack.");
        systStack->Draw("E2 SAME");
    }
    for (const auto &[name, band] : nonStackSystBands)
    {
        if (band != nullptr)
            band->Draw("E2 SAME");
    }
    std::vector<TBox *> topBlindBoxes;
    if (!setup.blindRanges.empty())
    {
        const auto displayBlindRanges = buildDisplayBlindRanges(
            setup.blindRanges,
            setup.displayEqualWidthBins,
            setup.displayBinEdges,
            axisHist->GetXaxis()->GetXmin(),
            axisHist->GetXaxis()->GetXmax());
        topBlindBoxes = drawBlindBoxes(displayBlindRanges, axisHist->GetMinimum(), axisHist->GetMaximum());
    }
    drawNonStackedHists(nonStackedHists, isData, 1, 0, 1.0);
    gPad->RedrawAxis();
    if (setup.displayEqualWidthBins && !doRatio)
    {
        drawEqualWidthAxisLabels(gPad, setup.displayBinEdges, this->_scale[0], false);
    }

    // legend and latex
    auto len = setHanyangLegend(hists.size(), 0.03, this->_scale[0], false, doRatio);
    const int nDataHists = countDataHistograms(isData);
    for (const auto &name : reStackOrder)
    {
        if ((labels.size() == 0) || (labels.find(name) == labels.end()))
            len->AddEntry(plotHists.at(name), name.c_str(), "f");
        else
            len->AddEntry(plotHists.at(name), labels.at(name).c_str(), "f");
    }
    for (const auto &[name, hist] : nonStackedHists)
    {
        if (isData[name])
        {
            if (nDataHists == 1)
                len->AddEntry(hist, "Data", "ep");
            else if ((labels.size() == 0) || (labels.find(name) == labels.end()))
                len->AddEntry(hist, name.c_str(), "ep");
            else
                len->AddEntry(hist, labels.at(name).c_str(), "ep");
        }
        else
        {
            std::string labelName = name;
            if ((labels.size() == 0) || (labels.find(name) == labels.end()))
            {
            }
            else
                labelName = labels.at(name);
            // if (colorScheme.find(name) != colorScheme.end())
            if (signalScaledForDisplay && std::find(isSignal.begin(), isSignal.end(), name) != isSignal.end())
            {
                const auto scaleIter = signalDisplayScaleFactors.find(name);
                if (scaleIter != signalDisplayScaleFactors.end() && scaleIter->second > 1)
                {
                    labelName += formatSignalScaleFactor(scaleIter->second);
                }
            }
            len->AddEntry(hist, labelName.c_str(), "l");
            // if ((labels.size() == 0) || (labels.find(name) == labels.end()))
            //     len->AddEntry(hist, name.c_str(), "l");
            // else
            //     len->AddEntry(hist, labels.at(name).c_str(), "l");
        }
    }

    len->Draw("SAME");
    this->_legends.push_back(len);
    drawHeader(headerTexts, this->_scale[0], doRatio);
    drawTexts(aboveTexts, this->_scale[0]);

    // draw ratio
    std::map<std::string, TH1D *> plotRatioHists;
    TLine *refLine;
    auto ratioSystStack = getRatioUncert(systStack);
    std::vector<TBox *> ratioBlindBoxes;
    if (doRatio)
    {
        this->_canvas->cd();
        this->_belowPad->Draw();
        this->_belowPad->cd();
        TGaxis::SetExponentOffset(-0.08, 0.02, "y");
        plotRatioHists = setupHists(ratioHists, setup, colorScheme, binLabels, 1);
        if (setup.displayEqualWidthBins)
        {
            auto remappedRatioHists = buildEqualWidthDisplayHistMap(plotRatioHists);
            deleteHistogramMap(plotRatioHists);
            plotRatioHists = std::move(remappedRatioHists);
            for (auto &[name, hist] : plotRatioHists)
            {
                configureEqualWidthAxisStyle(hist, true);
                hist->GetYaxis()->SetNdivisions(4, false);
            }
        }
        if (!setup.blindRanges.empty())
        {
            for (const auto &[name, hist] : plotRatioHists)
            {
                auto isDataIter = setup.isData[1].find(name);
                if (isDataIter != setup.isData[1].end() && isDataIter->second)
                    maskBlindedBinsForDisplay(
                        hist,
                        setup.blindRanges,
                        setup.displayEqualWidthBins,
                        setup.displayBinEdges);
            }
        }
        // in the ratio hist, mc in the denominator, need to scale by the inverse
        drawCanvasHists(plotRatioHists, this->_scale[1], setup.isData[1], 0, 1, belowTexts, 1.0/mcScaling);
        drawTexts(belowTexts, this->_scale[1]);
        refLine = new TLine(
            plotRatioHists.begin()->second->GetXaxis()->GetXmin(),
            1.0,
            plotRatioHists.begin()->second->GetXaxis()->GetXmax(),
            1.0);
        refLine->SetLineStyle(2);
        refLine->SetLineWidth(1);
        refLine->SetLineColor(1);
        if (ratioSystStack != nullptr)
        {
            ratioSystStack->Draw("E2 SAME");
        }
        if (!setup.blindRanges.empty())
        {
            const auto displayBlindRanges = buildDisplayBlindRanges(
                setup.blindRanges,
                setup.displayEqualWidthBins,
                setup.displayBinEdges,
                plotRatioHists.begin()->second->GetXaxis()->GetXmin(),
                plotRatioHists.begin()->second->GetXaxis()->GetXmax());
            // The ratio must be fully covered, including the reference line
            // and uncertainty band; keep the upper-pad blind overlay unchanged.
            ratioBlindBoxes = drawBlindBoxes(displayBlindRanges, 0.5, 1.5, 1.0, kGray);
        }
        // Keep the unit-ratio reference visible on top of the opaque blind box.
        refLine->Draw("same");
        gPad->RedrawAxis();
        if (setup.displayEqualWidthBins)
        {
            drawEqualWidthAxisLabels(gPad, setup.displayBinEdges, this->_scale[1], true);
        }
        this->_canvas->Update();
    }

    saveCanvas(this->_controllerName, setup.doLog);

    delete axisHist;
    delete stack;
    for (const auto &[name, hist] : plotHists)
    {
        if (hist)
        {
            hist->SetDirectory(0);
            delete hist;
        }
    }
    if (doRatio)
    {
        delete refLine;
        for (const auto &[name, hist] : plotRatioHists)
        {
            if (hist)
            {
                hist->SetDirectory(0);
                delete hist;
            }
        }
    }
    if (systStack != nullptr)
    {
        delete systStack;
    }
    if (ratioSystStack != nullptr)
    {
        delete ratioSystStack;
    }
    for (auto *box : topBlindBoxes)
        delete box;
    for (auto *box : ratioBlindBoxes)
        delete box;
    for (auto &[name, band] : nonStackSystBands)
    {
        delete band;
    }
    if (setup.displayEqualWidthBins)
    {
        deleteHistogramPairs(displayShapeUncerts);
        deleteHistogramPairMap(displayNonStackShapeUncerts);
    }
}
