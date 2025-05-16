#ifndef PLOT_CONTROLLER_H
#define PLOT_CONTROLLER_H

#include <string>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TPad.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"

// class to guiding processing he
struct PlotContext
{
    // Whether each histogram is data, data and MC has different styles
    // The first vector above is for yield; if ratio plot, the below is for ratio
    std::vector<std::map<std::string, int>> isData;
    std::vector<std::string> isSignal;

    // labels
    std::string xLabel = "m_{#mu#mu} / GeV";
    std::vector<std::string> yLabel = {"Yield", "Data/MC"};
    int doLog;
    int doNormalize=0;

    // canvas size
    double xSize;
    double ySize;
};

// class to keep an internal canvas
// each time input histograms to draw, do the setup and afterwards clear the canvas
class PlotControl
{
private:
    // ID name for safely manipulate canvas and pads
    std::string _controllerName;
    TCanvas *_canvas;
    TPad *_abovePad;
    TPad *_belowPad;
    // lengend & drawing
    std::vector<TLegend *> _legends;
    // scale for axis label tuning, by default {1.0, 1.0}
    std::vector<double> _scale;
    double _histXSize, _histYSize;

    // need to store the top & bottom margin
    double _topMargin;
    double _bottomMargin;

public:
    explicit PlotControl(const std::string &name);
    PlotControl() = delete;
    ~PlotControl();

    // implement proper function to clear the contents in the canvas?
    // void clearCanvas();

    // hanyang style canvas
    void setHanyangCanvas(double xSize, double ySize, int doLog, int doRatio);

    // hanyang style hists
    void setHanyangHist(TH1D *hist, int color, int isData, const std::vector<std::string>& binLabels={}, double scale = 1.0, std::string xTitle = "m_{#mu#mu}", std::string yTitle = "yield");
    void setMax(std::map<std::string, TH1D *> hists, int doLog, int isRatio);
    std::map<std::string, TH1D *> setupHists(std::map<std::string, TH1D *> hists, PlotContext setup, const std::map<std::string, int> &colorScheme, const std::vector<std::string>& binLabels = {}, int isRatio = 0);
    // get stack hist
    THStack *prepareStackHists(std::map<std::string, TH1D *> &hists, std::vector<std::string> &stackOrder, std::map<std::string, int> isData, int reOrder, int doNormalize=0, float mcScaling=1);

    // hanyang style legend
    TLegend *setHanyangLegend(int entries, double textSize = 0.04, double scale = 1.0);

    // function for drawing non-stacked histograms
    void drawNonStackedHists(std::map<std::string, TH1D *> &hists, std::map<std::string, int> &isData, int same = 0, int isRatio = 0, float mcScaling = 1);

    // latex
    void drawTexts(const std::vector<std::string> &texts, double scale);

    // function for plot hist
    void drawCanvasHists(std::map<std::string, TH1D *> &plotHists, double scale, std::map<std::string, int> &isData, int doLegend, int isRatio = 0, std::vector<std::string> plotTexts = {}, float mcScaling = 1);

    // save
    void saveCanvas(const std::string &fileName);

    // utilize function
    void drawStackHistWithRatio(const std::map<std::string, TH1D *> &hists, const std::vector<std::string> &stackOrder, const std::map<std::string, double> &stackUp, const std::map<std::string, double> &stackDown, int reOrder, const std::map<std::string, TH1D *> &ratioHists, PlotContext setup, float mcScaling, const std::map<std::string, int> &colorScheme, const std::map<std::string, std::string> &labels, const std::vector<std::string> &aboveTexts, const std::vector<std::string> &belowTexts = {}, const std::vector<std::string>& binLabels={});
};

#endif