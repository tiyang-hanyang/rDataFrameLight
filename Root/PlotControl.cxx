#include "PlotControl.h"
#include "Utility.h"

#include <exception>
#include <iostream>
#include <numeric>

#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"

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
void PlotControl::setHanyangHist(TH1D *hist, int color, int isData, const std::vector<std::string> &binLabels, double scale, std::string xTitle, std::string yTitle)
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
}

TLegend *PlotControl::setHanyangLegend(int entries, double textSize, double scale)
{
    // introduce the concept of effective entries to enable two columns
    int equEntries = entries;
    if (entries > 8)
    {
        equEntries = (equEntries + 1) / 2;
    }
    auto len = new TLegend(0.6, 1.0 - (1 - 0.9 + equEntries * textSize) / scale, 0.89, 1.0 - (1.0 - 0.89) / scale, NULL, "brNDC");
    if (entries > 8)
    {
        len->SetNColumns(2);
    }
    len->SetBorderSize(1);
    len->SetTextFont(42);
    len->SetTextSize(textSize / scale);
    len->SetLineColor(0);
    len->SetLineStyle(1);
    len->SetLineWidth(1);
    len->SetFillColor(0);
    len->SetFillStyle(1001);
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

        // take the fixed range option now
        maxVal = 1.2;
        minVal = 0.8;
        
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
            setHanyangHist(clonedHist, 1, isData, binLabels, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio]);
        else
        {
            rdfWS_utility::messageINFO("PlotControl", "Setup hist " + histChannel);
            if (!isRatio)
            {
                setHanyangHist(clonedHist, colorScheme.at(histChannel), isData, binLabels, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio]);
            }
            else
                setHanyangHist(clonedHist, 1, isData, binLabels, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio]);
            // setHanyangHist(clonedHist, colorScheme[colorIndex], isData, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio]);
            // colorIndex++;
        }
        if (isRatio)
            clonedHist->GetYaxis()->SetNdivisions(505);

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
                hist->SetMinimum(0.0);
                // hist->SetMinimum(0.5);
                hist->SetMaximum(2.0);
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
    // if (texts.empty())
    // {
    //     return;
    // }

    // the previous input style
    // "texts": [
    //     "CMS Preliminary",
    //     "#sqrt{s} = 13.6 TeV, L= %1.0f fb^{-1}",
    //     "S1 dimuon, tight ID jets, UParT"
    // ],

    double defaultSize = 0.03;

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
        initHeight -= 0.05 / scale;
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
        auto len = setHanyangLegend(plotHists.size(), 0.03, scale);
        for (auto &[histName, hist] : plotHists)
        {
            len->AddEntry(hist, histName.c_str());
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

TGraphAsymmErrors *uncertaintyHist(THStack *stackHist, std::map<std::string, TH1D *> &hists, const std::vector<std::string> &stackOrder, std::map<std::string, double> stackUp, std::map<std::string, double> stackDown)
{
    if (stackOrder.size() == 0 || !(stackHist))
    {
        return nullptr;
    }

    if (stackUp.size() == 0 || stackDown.size() ==0) return nullptr;

    TH1D *histUp = dynamic_cast<TH1D *>(hists[stackOrder[0]]->Clone("stack_Up"));
    TH1D *histDown = dynamic_cast<TH1D *>(hists[stackOrder[0]]->Clone("stack_Down"));
    if (stackUp.find(stackOrder[0]) != stackUp.end())
        histUp->Scale(stackUp[stackOrder[0]]);
    if (stackDown.find(stackOrder[0]) != stackDown.end())
        histDown->Scale(stackDown[stackOrder[0]]);

    // TODO
    // I found a problem here, the hists in the PlotControl are already the merged version instead of the single channels, thus the syst not applied here.
    // What I should do is to add structure inside the HistControl for holding all the syst uncertainties
    // Now just write simply the syst manually after merge
    for (auto iter = stackOrder.begin() + 1; iter != stackOrder.end(); iter++)
    {
        double scaleUp = 1;
        double scaleDown = 1.0;
        if (stackUp.find(*iter) != stackUp.end())
            scaleUp = stackUp.at(*iter);
        if (stackDown.find(*iter) != stackDown.end())
            scaleDown = stackDown.at(*iter);
        histUp->Add(hists.at(*iter), scaleUp);
        histDown->Add(hists.at(*iter), scaleDown);
    }

    TH1D *nominalHist = dynamic_cast<TH1D *>(stackHist->GetStack()->Last()->Clone("stackNom"));
    TGraphAsymmErrors *systBand = new TGraphAsymmErrors(nominalHist);
    for (int i = 1; i < nominalHist->GetNbinsX() + 1; i++)
    {
        auto binCenter = nominalHist->GetBinCenter(i);
        auto binWidth = nominalHist->GetBinWidth(i) / 2.0;
        auto nominalVal = nominalHist->GetBinContent(i);
        auto errDown = std::max(0.0, nominalVal - histDown->GetBinContent(i));
        auto errUp = std::max(0.0, histUp->GetBinContent(i) - nominalVal);
        systBand->SetPoint(i - 1, binCenter, nominalVal);
        systBand->SetPointError(i - 1, binWidth, binWidth, errDown, errUp);
    }

    systBand->SetFillColor(921);
    systBand->SetFillStyle(3004);
    systBand->SetLineWidth(0);
    systBand->SetMarkerSize(0);

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

    ratioUncert->SetFillColor(921);
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
    std::map<std::string, int> isData = setup.isData[0];
    std::vector<std::string> isSignal = setup.isSignal;

    std::vector<std::string> reStackOrder = stackOrder;
    THStack *stack = prepareStackHists(plotHists, reStackOrder, isData, reOrder, setup.doNormalize, mcScaling);

    // Add stack histogram uncertainties
    auto systStack = uncertaintyHist(stack, plotHists, reStackOrder, stackUp, stackDown);

    // to make sure the stacked hist maximum also covers the non-stack hist range
    double nonStackMax = stack->GetMaximum();
    std::map<std::string, TH1D *> nonStackedHists;
    for (auto &[name, hist] : plotHists)
    {
        if (std::find(stackOrder.begin(), stackOrder.end(), name) == stackOrder.end())
        {
            if (setup.doNormalize)
            {
                hist->Scale(1.0 / hist->Integral());
            }
            else if (std::find(isSignal.begin(), isSignal.end(), name) != isSignal.end())
            {
                hist->Scale(100000.);
                hist->SetLineStyle(9);
                hist->SetLineWidth(2);
            }
            nonStackedHists.emplace(name, hist);
            if (hist->GetMaximum() > nonStackMax) nonStackMax = hist->GetMaximum();
        }
    }

    // generate an empty histogram for the axis style
    TH1D *axisHist = dynamic_cast<TH1D *>(plotHists.begin()->second->Clone((this->_controllerName + "_axisHist").c_str()));
    axisHist->SetDirectory(0);
    if (!axisHist)
    {
        rdfWS_utility::messageERROR("PlotControl", "Failed to clone reference axis histogram.");
    }
    axisHist->Reset();
    if (setup.doLog)
    {
        axisHist->SetMaximum(std::max(axisHist->GetMaximum(), nonStackMax) * 5000.);
        axisHist->SetMinimum(1.);
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
    }
    axisHist->Draw();
    if (doStack)
        stack->Draw("HIST SAME");
    if (systStack != nullptr)
    {
        rdfWS_utility::messageINFO("PlotControl", "Draw with systs in stack.");
        systStack->Draw("E2 SAME");
    }
    drawNonStackedHists(nonStackedHists, isData, 1, doRatio, mcScaling);
    gPad->RedrawAxis();

    // legend and latex
    auto len = setHanyangLegend(hists.size(), 0.03, this->_scale[0]);
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
            // if ((labels.size() == 0) || (labels.find(name) == labels.end()))
            //     len->AddEntry(hist, name.c_str(), "ep");
            // else
            //     len->AddEntry(hist, labels.at(name).c_str(), "ep");
            len->AddEntry(hist, "Data", "ep");
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
            if (std::find(isSignal.begin(), isSignal.end(), name) != isSignal.end())
                labelName += " #times 1e5";
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
    if (doRatio)
    {
        this->_canvas->cd();
        this->_belowPad->Draw();
        this->_belowPad->cd();
        plotRatioHists = setupHists(ratioHists, setup, colorScheme, binLabels, 1);
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
        refLine->Draw("same");
        if (ratioSystStack != nullptr)
        {
            ratioSystStack->Draw("E2 SAME");
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
}
