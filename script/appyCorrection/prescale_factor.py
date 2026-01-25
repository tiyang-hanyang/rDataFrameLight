import json
import ROOT
import os
import correctionlib

def parse_prescale(jsonPath, outPath):
    with open(jsonPath) as jFile:
        scaleJson = json.load(jFile)
        # One thing: the HLT 20 uses HLT 17 name here so that goes wrong
        trigger_str = scaleJson['pathname']
        scale_data = scaleJson['prescales']

    # creating the TSV
    header = "run\tls_begin\tls_end\tprescale\n"
    tableEntries = []
    for (run, lsInfo) in scale_data.items():
        for lsEntry in lsInfo:
            line = run + "\t" + str(lsEntry['lumisec_start']) + "\t" + str(lsEntry['lumisec_end']) + "\t" + str(lsEntry['hltprescale']*lsEntry['l1prescale']) + "\n"
            tableEntries.append(line)

    # output
    #with open(trigger_str+".txt", "w") as oFile:
    with open(outPath, "w") as oFile:
        oFile.write(header)
        for line in tableEntries:
            oFile.write(line)


# PU-reweighting from the LUM POG
# In principal, I would need to let the HLT trigger path to be automatically specified but not direct input
# As the running would not need to do that in the same chain, so I now just have the hard-coding, later need to change this
def processing(rdf, recordedModules, branchArray, era, ds=""):
    jsonFolder = "/home/tiyang/public/rDataFrameLight_git/source/json/prescale/"
    jsonLists = [
        "HLT_Mu3_PFJet40_v",
        "HLT_Mu8_v",
        "HLT_Mu17_v",
        "HLT_Mu20_v",
        "HLT_Mu27_v",
    ]

    # to convert all the possible needed prescales
    for tp in jsonLists:
        jsonPath = jsonFolder + tp + "_prescales.json"  # e.g. HLT_Mu8_v_prescales.json
        prescaleTablePath = tp + ".txt"    # e.g. HLT_Mu8_v.txt
        if (not os.path.isfile(prescaleTablePath)):
            parse_prescale(jsonPath, prescaleTablePath)
    
    if "prescale_factor.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "prescale_factor.C"')
        recordedModules.append("prescale_factor.C")
        # ROOT.gInterpreter.ProcessLine('PreScale_init("'+prescaleTablePath+'")')

    if "HLT_Mu3_PFJet40_prescale" not in branchArray:
        toload = "HLT_Mu3_PFJet40_v.txt"   
        ROOT.gInterpreter.ProcessLine('PreScale_init("HLT_Mu3_PFJet40", "'+toload+'")')
        rdf = rdf.Define("HLT_Mu3_PFJet40_prescale", "readTriggerPS(\"HLT_Mu3_PFJet40\", run, luminosityBlock)")
        branchArray.append("HLT_Mu3_PFJet40_prescale")

    if "HLT_Mu8_prescale" not in branchArray:
        toload = "HLT_Mu8_v.txt"   
        ROOT.gInterpreter.ProcessLine('PreScale_init("HLT_Mu8", "'+toload+'")')
        rdf = rdf.Define("HLT_Mu8_prescale", "readTriggerPS(\"HLT_Mu8\", run, luminosityBlock)")
        branchArray.append("HLT_Mu8_prescale")

    if "HLT_Mu17_prescale" not in branchArray:
        toload = "HLT_Mu17_v.txt"   
        ROOT.gInterpreter.ProcessLine('PreScale_init("HLT_Mu17", "'+toload+'")')
        rdf = rdf.Define("HLT_Mu17_prescale", "readTriggerPS(\"HLT_Mu17\", run, luminosityBlock)")
        branchArray.append("HLT_Mu17_prescale")

    if "HLT_Mu20_prescale" not in branchArray:
        toload = "HLT_Mu20_v.txt"   
        ROOT.gInterpreter.ProcessLine('PreScale_init("HLT_Mu20", "'+toload+'")')
        rdf = rdf.Define("HLT_Mu20_prescale", "readTriggerPS(\"HLT_Mu20\", run, luminosityBlock)")
        branchArray.append("HLT_Mu20_prescale")

    if "HLT_Mu27_prescale" not in branchArray:
        toload = "HLT_Mu27_v.txt"   
        ROOT.gInterpreter.ProcessLine('PreScale_init("HLT_Mu27", "'+toload+'")')
        rdf = rdf.Define("HLT_Mu27_prescale", "readTriggerPS(\"HLT_Mu27\", run, luminosityBlock)")
        branchArray.append("HLT_Mu27_prescale")

    return rdf, recordedModules, branchArray