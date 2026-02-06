#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath> 
#include <algorithm> 
#include <cstdint>

// look up table
struct Interval { 
    uint32_t ls_begin; 
    uint32_t ls_end; // should be inclusive bound
    float    w; 
};

struct RunBlock {
    uint32_t run;
    uint32_t offset;
    uint32_t n; 
};

struct LineCode {
    uint32_t run;
    uint32_t ls_start;
    uint32_t ls_end;
    float weight;
};

// allowing to re-register these
// should not be re-register, but increasing the registration to ensuring can do the snapshot
std::map<std::string, std::vector<RunBlock>> _prescaleRunInfo;
std::map<std::string, std::vector<Interval>> _prescaleLSWeight;
// std::vector<RunBlock> _prescaleRunInfo;
// std::vector<Interval> _prescaleLSWeight;

LineCode lineReader(const std::string& line)
{
    // run \t ls_init \t ls_end \t weight
    // using stringsteam to partition lines
    std::stringstream ss(line);
    std::string run_s, ls1_s, ls2_s, w_s;
    if (!std::getline(ss, run_s, '\t') ||
        !std::getline(ss, ls1_s, '\t') ||
        !std::getline(ss, ls2_s, '\t') ||
        !std::getline(ss, w_s,   '\t')) 
    {
        throw std::runtime_error(std::string("Bad TSV line (missing fields): ") + line);
    }

    LineCode out{
        static_cast<uint32_t>(std::stoul(run_s)),
        static_cast<uint32_t>(std::stoul(ls1_s)),
        static_cast<uint32_t>(std::stoul(ls2_s)),
        std::stof(w_s)
    };

    if (out.ls_end < out.ls_start) {
        throw std::runtime_error(std::string("Bad LS range (end < start) in line: ") + line);
    }
    if (!std::isfinite(out.weight)) {
        throw std::runtime_error(std::string("Bad weight (NaN/Inf) in line: ") + line);
    }
    
    return out;
}

void PreScale_init(const std::string& tpName , const std::string& prescaleTablePath)
{
    // // clear folder
    // _prescaleRunInfo.clear();
    // _prescaleLSWeight.clear();
    std::vector<RunBlock> new_prescaleRunInfo;
    std::vector<Interval> new_prescaleLSWeight;

    std::ifstream fin(prescaleTablePath);
    if (!fin.is_open())
    {
        throw std::runtime_error(std::string("The fin ") + prescaleTablePath + " cannot be opened");
    }
    std::string line;

    // skip the title line
    if (!std::getline(fin, line)) {
        throw std::runtime_error(std::string("TSV is empty: ") + prescaleTablePath);
    }

    // my run should by default be increasing according to the creation inside the python
    // skip the sort part
    /* sort func */

    // process first entry specially to get init values
    // find first data line (skip empty / comment)
    while (true) {
        if (!std::getline(fin, line)) 
        {
            throw std::runtime_error(std::string("No data lines in TSV: ") + prescaleTablePath);
        }
        if (!line.empty() && line[0] != '#') break;
    }
    auto firstLine = lineReader(line);
    uint32_t lastRun = firstLine.run;
    uint32_t lastOffset = 0; // offset of lastRun
    uint32_t lastN = 1; // number for entries in lastEntry
    Interval firstLS{firstLine.ls_start, firstLine.ls_end, firstLine.weight}; 
    new_prescaleLSWeight.push_back(firstLS);

    // processing other lines
    while (std::getline(fin, line))
    {
        // in case if any empty lines or command
        if (line.empty() || line[0] == '#') continue;
        auto newLine = lineReader(line);
        uint32_t newLineRun = newLine.run;
        // line changes
        if (newLineRun != lastRun)
        {
            if (newLineRun < lastRun) {
                throw std::runtime_error("Run is not strictly increasing; TSV not sorted by run");
            }

            // finish one run, now update to the global run table, and clear the previous run count
            RunBlock finishedRun{lastRun, lastOffset, lastN}; 
            new_prescaleRunInfo.push_back(finishedRun);
            lastRun = newLineRun;
            lastN=0;
            lastOffset = new_prescaleLSWeight.size();
        }
        // deal with the new LS
        Interval newLS{newLine.ls_start, newLine.ls_end, newLine.weight}; 
        new_prescaleLSWeight.push_back(newLS);
        lastN++;
    }

    // deal with the last run, no getline possible
    RunBlock finalRun{lastRun, lastOffset, lastN};
    new_prescaleRunInfo.push_back(finalRun);

    // push back into the original
    _prescaleRunInfo.emplace(tpName, new_prescaleRunInfo);
    _prescaleLSWeight.emplace(tpName, new_prescaleLSWeight);
}

float readTriggerPS(const std::string& tpName, unsigned int run, unsigned int luminosityBlock)
{
    // if (_prescaleRunInfo.empty() || _prescaleLSWeight.empty()) {
    //     return 0.0f; // not initialized or empty table
    // }

    if (_prescaleRunInfo.find(tpName) == _prescaleRunInfo.end() || _prescaleLSWeight.find(tpName) == _prescaleLSWeight.end()) {
        return 0.0f; // not initialized or empty table
    }

    auto prescaleRunInfo = _prescaleRunInfo[tpName];
    auto prescaleLSWeight = _prescaleLSWeight[tpName];

    // 1) find the run inside Runblock
    auto it = std::lower_bound(
        prescaleRunInfo.begin(), prescaleRunInfo.end(), static_cast<uint32_t>(run),
        [](const RunBlock& rb, unsigned int r){ return rb.run < r; }
    );
    // the run not mentioned by the prescale table are not triggered
    if (it == prescaleRunInfo.end() || it->run != run) {
        return 0.0;
    }

    // to find the corresponding luminosity blocks range
    const uint32_t offset = it->offset;
    const uint32_t n      = it->n;
    const Interval* base = prescaleLSWeight.data() + offset;
    const Interval* end  = base + n;

    // finding the corresponding luminosity block range by initial luminosity block
    const uint32_t lumi = static_cast<uint32_t>(luminosityBlock);
    auto it_ls = std::upper_bound(
        base, end, lumi,
        [](uint32_t ls, const Interval& lsi){ return ls < lsi.ls_begin; }
    );

    if (it_ls == base ) {
        return 0.0f; // first available LS range still have larger ls_begin than the requested events, not triggered.
    }

    // make sure really inside the range
    const Interval& ls_cand = *(it_ls - 1);
    if (lumi >= ls_cand.ls_begin && lumi <= ls_cand.ls_end) {
        return ls_cand.w;
    }
    return 0.0f;
}