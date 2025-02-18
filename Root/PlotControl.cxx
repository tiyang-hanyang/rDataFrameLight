#include "PlotControl.h"

#include <exception>
#include <iostream>

#include "TLatex.h"
#include "TLine.h"

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
        delete len;
    }
    this->_legends.clear();

    // clear canvas
    if (this->_canvas)
    {
        std::cout << "Removing all histograms from canvas before delete" << std::endl;
        this->_canvas->GetListOfPrimitives()->Clear();
        this->_canvas->Modified();
        this->_canvas->Update();
        delete this->_canvas;
        this->_canvas = nullptr;
    }

    // clear pads
    if (this->_abovePad)
    {
        // delete this->_abovePad;
        this->_abovePad = nullptr;
    }
    if (this->_belowPad)
    {
        // delete this->_belowPad;
        this->_belowPad = nullptr;
    }
}

// setup canvas and pad style, record the scale if do ratio plot
void PlotControl::setHanyangCanvas(double xSize, double ySize, int doLog, int doRatio)
{
    if (this->_canvas != nullptr)
    {
        throw std::runtime_error("internal canvas already exist");
    }
    std::string canvasName = "c_" + this->_controllerName;
    this->_canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), xSize, ySize);

    // margin setup
    double tMargin = 0.05;
    double bMargin = 0.15;
    this->_topMargin = tMargin;
    this->_bottomMargin = bMargin;

    // ticks and grid (if log)
    this->_canvas->SetTickx(1);
    this->_canvas->SetTicky(1);
    if (doLog)
    {
        this->_canvas->SetLogy();
        this->_canvas->SetGrid();
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
            this->_abovePad->SetGrid();
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
void PlotControl::setHanyangHist(TH1D *hist, int color, int isData, double scale, std::string xTitle, std::string yTitle)
{
    hist->SetTitle("");
    hist->SetStats(0);
    // color
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);

    // data has marker style while has MC line style
    if (isData)
    {
        // data draw option is ep
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(0.5);
        hist->SetLineStyle(0);
    }
    else
    {
        // MC draw option is e hist
        hist->SetLineWidth(2);
    }

    // axis style setting
    hist->GetXaxis()->SetTitle(xTitle.c_str());
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.04 / scale);
    hist->GetXaxis()->SetLabelOffset(0.007 / scale);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetTitleSize(0.055 / scale);
    hist->GetXaxis()->SetTitleOffset(0.9);

    hist->GetYaxis()->SetTitle(yTitle.c_str());
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelSize(0.04 / scale);
    hist->GetYaxis()->SetLabelOffset(0.007);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleSize(0.055 / scale);
    hist->GetYaxis()->SetTitleOffset(1.25 * scale);
}

TLegend *PlotControl::setHanyangLegend(int entries, double textSize, double scale)
{
    auto len = new TLegend(0.6, 1.0 - (1 - 0.92 + entries * textSize) / scale, 0.9, 1.0 - (1.0 - 0.92) / scale, NULL, "brNDC");
    len->SetBorderSize(1);
    len->SetTextFont(62);
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
    if (hists.empty())
    {
        throw std::runtime_error("No histograms provided for setup.");
    }

    double maxVal(0.0), minVal(1e9);
    if (isRatio)
    {
        for (const auto &[histName, hist] : hists)
        {
            if (!hist)
            {
                throw std::runtime_error("Histogram " + histName + " is null.");
            }

            if (hist->GetMaximum() > maxVal)
                maxVal = hist->GetMaximum();
            if (hist->GetMinimum() < minVal)
                minVal = hist->GetMinimum();
        }
        maxVal *= 1.2;
        minVal *= 0.8;
    }
    else
    {
        for (const auto &[histName, hist] : hists)
        {
            if (!hist)
            {
                throw std::runtime_error("Histogram " + histName + " is null.");
            }

            if (hist->GetMaximum() > maxVal)
                maxVal = hist->GetMaximum();
        }
        if (doLog)
        {
            maxVal *= 100.0;
            minVal = 1.0;
        }
        else
        {
            maxVal *= 1.4;
            minVal = 0.0;
        }
    }

    for (const auto &[histName, hist] : hists)
    {
        hist->SetMaximum(maxVal);
        hist->SetMinimum(minVal);
    }
}

std::map<std::string, TH1D *> PlotControl::setupHists(std::map<std::string, TH1D *> hists, PlotContext setup, int isRatio)
{
    if (hists.empty())
    {
        throw std::runtime_error("No histograms provided for setup.");
    }
    std::vector<int> colorScheme = {
        2, 3, 4, 5, 6, 7};
    int colorIndex(0);

    std::map<std::string, TH1D *> styledHists;
    // Clone and style each hist
    for (const auto &[histName, hist] : hists)
    {
        if (!hist)
        {
            throw std::runtime_error("Histogram " + histName + " is null.");
        }
        // Clone histogram
        std::string newHistName = std::string("plot_") + hist->GetName();
        TH1D *clonedHist = dynamic_cast<TH1D *>(hist->Clone(newHistName.c_str()));
        clonedHist->SetDirectory(0);
        // Apply styles (e.g., isData, labels)
        // add check for if isData set
        int isData = setup.isData[isRatio].find(histName) != setup.isData[isRatio].end() && setup.isData[isRatio].at(histName);
        if (isData)
            setHanyangHist(clonedHist, 1, isData, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio]);
        else
        {
            setHanyangHist(clonedHist, colorScheme[colorIndex], isData, this->_scale.at(isRatio), setup.xLabel, setup.yLabel[isRatio]);
            colorIndex++;
        }
        if (isRatio)
            clonedHist->GetYaxis()->SetNdivisions(505);

        styledHists.emplace(histName, clonedHist);
    }
    // setup maximum value according to all hists
    setMax(styledHists, setup.doLog, isRatio);

    return styledHists;
}

// function to containing all hist drawing for non-stack histograms
// the input hists should be already striped out stacked part
// same option is for if the first histogram need option same
void PlotControl::drawNonStackedHists(
    std::map<std::string, TH1D *> &hists,
    std::map<std::string, int> &isData,
    int same)
{
    // reserve a place for the dataHist, so that I would always plot this last
    std::vector<TH1D *> dataHists;
    for (auto &[histName, hist] : hists)
    {
        if (!hist)
        {
            throw std::runtime_error("Histogram " + histName + " is null.");
        }
        //     std::string drawOption = isData[histName] ? "ep" : "HIST";
        //     if (same)
        //     {
        //         drawOption += " SAME";
        //     }

        //     auto histCloned = (TH1D*)hist->Clone((std::string(hist->GetName())+"_copyDraw").c_str());
        //     histCloned->SetDirectory(0);
        //     histCloned->Draw(drawOption.c_str());
        //     same = true; // Ensure subsequent histograms use "SAME"
        // }
        if (isData[histName])
        {
            auto histCloned = (TH1D *)hist->Clone((std::string(hist->GetName()) + "_copyDraw").c_str());
            histCloned->SetDirectory(0);
            dataHists.push_back(histCloned);
        }
        else
        {
            std::string drawOption = "HIST";
            if (same)
            {
                drawOption += " SAME";
            }

            auto histCloned = (TH1D *)hist->Clone((std::string(hist->GetName()) + "_copyDraw").c_str());
            histCloned->SetDirectory(0);
            histCloned->Draw(drawOption.c_str());
            same = true; // Ensure subsequent histograms use "SAME"
        }
    }
    for (auto hist : dataHists)
    {
        std::string drawOption = "ep";
        if (same)
        {
            drawOption += " SAME";
        }

        auto histCloned = (TH1D *)hist->Clone((std::string(hist->GetName()) + "_copyDraw").c_str());
        histCloned->SetDirectory(0);
        histCloned->Draw(drawOption.c_str());
        same = true;
    }
}

// latex
void PlotControl::drawTexts(const std::vector<std::string> &texts, double scale)
{
    if (texts.empty())
    {
        return;
    }

    double defaultSize = 0.04;

    TLatex padText;
    padText.SetNDC();
    padText.SetTextSize(defaultSize / scale);

    // double initHeight = 1 - (1 + defaultSize - 0.88) * scale;
    double initHeight = 1 - (1 - 0.88) / scale;
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
    std::vector<std::string> plotTexts)
{
    // draw all histograms not in stack
    drawNonStackedHists(plotHists, isData, 0);
    // non-stacked version of legend
    if (doLegend)
    {
        std::cout << "legend plot hist size: " << plotHists.size() << std::endl;
        std::cout << "legend scale: " << scale << std::endl;
        auto len = setHanyangLegend(plotHists.size(), 0.04, scale);
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

void PlotControl::saveCanvas(const std::string &fileName)
{
    if (!this->_canvas)
    {
        throw std::runtime_error("Internl canvas not properly setup: " + this->_controllerName);
    }
    this->_canvas->SaveAs((fileName + ".png").c_str());
    this->_canvas->SaveAs((fileName + ".pdf").c_str());
    this->_canvas->SaveAs((fileName + ".eps").c_str());
}

void PlotControl::drawHist(std::map<std::string, TH1D *> hists, PlotContext setup, std::vector<std::string> text)
{
    if (hists.empty())
    {
        throw std::runtime_error("No histograms provided for plotting.");
    }

    // Check if isData[0] matches the number of histograms
    if (setup.isData[0].size() != hists.size())
    {
        throw std::runtime_error("Mismatch between the number of histograms and isData[0].");
    }

    // set style canvas
    setHanyangCanvas(setup.xSize, setup.ySize, setup.doLog, 0);
    if (!this->_canvas)
    {
        throw std::runtime_error("Internl canvas not properly setup: " + this->_controllerName);
    }

    // set style hist
    std::map<std::string, TH1D *> plotHists = setupHists(hists, setup, 0);

    // plot hists
    this->_canvas->cd();
    drawCanvasHists(plotHists, this->_scale[0], setup.isData[0], 1, text);
    // save
    saveCanvas(this->_controllerName);
}

void PlotControl::drawRatioHist(std::map<std::string, TH1D *> hists, std::map<std::string, TH1D *> ratioHists, PlotContext setup, std::vector<std::string> aboveText, std::vector<std::string> belowText)
{
    if (hists.empty())
    {
        throw std::runtime_error("No histograms provided for plotting.");
    }
    if (ratioHists.empty())
    {
        throw std::runtime_error("No ratio histograms provided for plotting.");
    }

    // Check if isData matches the number of histograms
    if (setup.isData.size() != 2)
    {
        throw std::runtime_error("No two bunches of isData infor provided for ratio plot.");
    }
    if (setup.isData[0].size() != hists.size())
    {
        throw std::runtime_error("Mismatch between the number of histograms and isData[0].");
    }
    if (setup.isData[1].size() != ratioHists.size())
    {
        throw std::runtime_error("Mismatch between the number of ratio histograms and isData[1].");
    }

    // set style canvas
    setHanyangCanvas(setup.xSize, setup.ySize, setup.doLog, 1);
    if (!this->_canvas)
    {
        throw std::runtime_error("Internl canvas not properly setup: " + this->_controllerName);
    }

    // set style hist
    std::map<std::string, TH1D *> plotHists = setupHists(hists, setup, 0);
    std::map<std::string, TH1D *> plotRatioHists = setupHists(ratioHists, setup, 1);

    // plot hists
    this->_canvas->cd();
    this->_abovePad->Draw();
    this->_abovePad->cd();
    drawCanvasHists(plotHists, this->_scale[0], setup.isData[0], 1, aboveText);
    // below pannel does not need legend
    this->_canvas->cd();
    this->_belowPad->Draw();
    this->_belowPad->cd();
    drawCanvasHists(plotRatioHists, this->_scale[1], setup.isData[1], 0, belowText);

    // save
    saveCanvas(this->_controllerName);
}

// prepare THStack
THStack *PlotControl::prepareStackHists(std::map<std::string, TH1D *> &hists, const std::vector<std::string> &stackOrder, std::map<std::string, int> isData)
{
    THStack *stack = new THStack(("stack_" + _controllerName).c_str(), "Stacked Histograms");
    std::vector<int> colorScheme = {2, 3, 4, 5, 6, 7, 9, 12};
    int colorIndex = 0;

    std::vector<int> reverseColorScheme;
    auto nStack = stackOrder.size();
    for (size_t i = 0; i < nStack; i++)
    {
        reverseColorScheme.push_back(colorScheme[nStack - i - 1]);
    }

    for (const auto &name : stackOrder)
    {
        if (hists.find(name) == hists.end())
        {
            continue;
        }
        auto hist = hists.at(name);
        if (!isData.at(name))
        {
            hist->SetFillColor(reverseColorScheme[colorIndex]);
            hist->SetLineColor(1);
            hist->SetLineWidth(1);
            hist->SetFillStyle(1001);
            stack->Add(hist);
            colorIndex++;
        }
    }
    return stack;
}

void PlotControl::drawStackHist(
    const std::map<std::string, TH1D *> &hists,
    const std::vector<std::string> &stackOrder,
    const PlotContext &setup,
    const std::vector<std::string> &plotTexts)
{
    // checking input compatibility
    if (hists.empty() || stackOrder.empty())
    {
        throw std::runtime_error("No histograms or stack order provided for plotting.");
    }
    for (const auto &name : stackOrder)
    {
        if (hists.find(name) == hists.end())
        {
            throw std::runtime_error("Histogram '" + name + "' specified in stackOrder does not exist in the provided histograms.");
        }
    }

    // need to setup canvas in advance for scales needed by hist
    setHanyangCanvas(setup.xSize, setup.ySize, setup.doLog, 1);
    if (!this->_canvas)
    {
        throw std::runtime_error("Internal canvas not properly setup: " + this->_controllerName);
    }

    // prepare hists
    std::map<std::string, TH1D *> plotHists = setupHists(hists, setup, 0);
    std::map<std::string, int> isData = setup.isData[0];

    THStack *stack = prepareStackHists(plotHists, stackOrder, isData);
    std::map<std::string, TH1D *> nonStackedHists;
    for (const auto &[name, hist] : plotHists)
    {
        if (std::find(stackOrder.begin(), stackOrder.end(), name) == stackOrder.end())
        {
            nonStackedHists.emplace(name, hist);
        }
    }

    // generate an empty histogram for the axis style
    TH1D *axisHist = dynamic_cast<TH1D *>(plotHists.begin()->second->Clone((this->_controllerName + "_axisHist").c_str()));
    axisHist->SetDirectory(0);
    if (!axisHist)
    {
        throw std::runtime_error("Failed to clone reference histogram.");
    }
    axisHist->Reset();
    if (setup.doLog)
    {
        axisHist->SetMaximum(std::max(axisHist->GetMaximum(), stack->GetMaximum()) * 100);
        axisHist->SetMinimum(1.);
    }
    else
    {
        axisHist->SetMaximum(std::max(axisHist->GetMaximum(), stack->GetMaximum()) * 1.4);
        axisHist->SetMinimum(0.);
    }

    // draw hists
    this->_canvas->cd();
    axisHist->Draw();
    stack->Draw("HIST SAME");
    drawNonStackedHists(nonStackedHists, isData, 1);
    this->_canvas->Update();

    // legend and text
    auto len = setHanyangLegend(hists.size(), 0.04, this->_scale[0]);
    // make sure stack is in order
    for (const auto &[name, hist] : nonStackedHists)
    {
        if (isData[name])
        {
            len->AddEntry(hist, name.c_str(), "ep");
        }
        else
        {
            len->AddEntry(hist, name.c_str(), "l");
        }
    }
    for (const auto &name : stackOrder)
    {
        len->AddEntry(plotHists.at(name), name.c_str(), "f");
    }
    len->Draw("SAME");
    this->_legends.push_back(len);
    drawTexts(plotTexts, 1);

    saveCanvas(this->_controllerName);
}

void PlotControl::drawStackHistWithRatio(
    const std::map<std::string, TH1D *> &hists,
    const std::vector<std::string> &stackOrder,
    const std::map<std::string, TH1D *> &ratioHists,
    PlotContext setup,
    const std::vector<std::string> &aboveTexts,
    const std::vector<std::string> &belowTexts)
{
    // checking input compatibility
    if (hists.empty() || stackOrder.empty() || ratioHists.empty())
    {
        throw std::runtime_error("No histograms, stack order, or ratio histograms provided.");
    }
    for (const auto &name : stackOrder)
    {
        if (hists.find(name) == hists.end())
        {
            throw std::runtime_error("Histogram '" + name + "' specified in stackOrder does not exist in the provided histograms.");
        }
    }

    // need to setup canvas in advance for scales needed by hist
    setHanyangCanvas(setup.xSize, setup.ySize, setup.doLog, 1);

    // prepare hists
    std::map<std::string, TH1D *> plotHists = setupHists(hists, setup, 0);
    std::map<std::string, TH1D *> plotRatioHists = setupHists(ratioHists, setup, 1);
    std::map<std::string, int> isData = setup.isData[0];

    THStack *stack = prepareStackHists(plotHists, stackOrder, isData);
    std::map<std::string, TH1D *> nonStackedHists;
    for (const auto &[name, hist] : plotHists)
    {
        if (std::find(stackOrder.begin(), stackOrder.end(), name) == stackOrder.end())
        {
            nonStackedHists.emplace(name, hist);
        }
    }

    // generate an empty histogram for the axis style
    TH1D *axisHist = dynamic_cast<TH1D *>(plotHists.begin()->second->Clone((this->_controllerName + "_axisHist").c_str()));
    axisHist->SetDirectory(0);
    if (!axisHist)
    {
        throw std::runtime_error("Failed to clone reference histogram.");
    }
    axisHist->Reset();
    if (setup.doLog)
    {
        axisHist->SetMaximum(std::max(axisHist->GetMaximum(), stack->GetMaximum()) * 100);
        axisHist->SetMinimum(1.);
    }
    else
    {
        axisHist->SetMaximum(std::max(axisHist->GetMaximum(), stack->GetMaximum()) * 1.4);
        axisHist->SetMinimum(0.);
    }

    // plot
    this->_canvas->cd();
    this->_abovePad->Draw();
    this->_abovePad->cd();
    axisHist->Draw();
    stack->Draw("HIST SAME");
    drawNonStackedHists(nonStackedHists, isData, 1);
    gPad->RedrawAxis();

    // legend and latex
    auto len = setHanyangLegend(hists.size(), 0.04, this->_scale[0]);
    for (const auto &name : stackOrder)
    {
        len->AddEntry(plotHists.at(name), name.c_str(), "f");
    }
    for (const auto &[name, hist] : nonStackedHists)
    {
        if (isData[name])
        {
            len->AddEntry(hist, name.c_str(), "ep");
        }
        else
        {
            len->AddEntry(hist, name.c_str(), "l");
        }
    }
    len->Draw("SAME");
    this->_legends.push_back(len);
    drawTexts(aboveTexts, this->_scale[0]);

    // draw ratio
    this->_canvas->cd();
    this->_belowPad->Draw();
    this->_belowPad->cd();
    drawCanvasHists(plotRatioHists, this->_scale[1], setup.isData[1], 0, belowTexts);
    drawTexts(belowTexts, this->_scale[1]);
    TLine *refLine = new TLine(
        plotRatioHists.begin()->second->GetXaxis()->GetXmin(),
        1.0,
        plotRatioHists.begin()->second->GetXaxis()->GetXmax(),
        1.0);
    refLine->SetLineStyle(2);
    refLine->SetLineWidth(1);
    refLine->SetLineColor(1);
    refLine->Draw("same");
    this->_canvas->Update();

    saveCanvas(this->_controllerName);

    delete axisHist;
    for (const auto &[name, hist] : plotHists)
    {
        if (hist)
        {
            hist->SetDirectory(0);
            delete hist;
        }
    }
    delete stack;
}
