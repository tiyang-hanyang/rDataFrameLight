import ROOT
import os
import correctionlib

# PU-reweighting from the LUM POG
def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "PUWeight.C" not in recordedModules:
        correctionlib.register_pyroot_binding()
        this_dir = os.path.dirname(os.path.abspath(__file__))
        ROOT.gInterpreter.AddIncludePath(this_dir)
        ROOT.gInterpreter.ProcessLine('#include "PUWeight.C"')
        recordedModules.append("PUWeight.C")
    ROOT.gInterpreter.ProcessLine('PU_weight_init("'+era+'")')

    if "PUWeight" not in branchArray:
        rdf = rdf.Define("PUWeight", "PUReweightFunc(Pileup_nTrueInt, \"nominal\")")
        branchArray.append("PUWeight")
    if "PUWeight_up" not in branchArray:
        rdf = rdf.Define("PUWeight_up", "PUReweightFunc(Pileup_nTrueInt, \"up\")")
        branchArray.append("PUWeight_up")
    if "PUWeight_down" not in branchArray:
        rdf = rdf.Define("PUWeight_down", "PUReweightFunc(Pileup_nTrueInt, \"down\")")
        branchArray.append("PUWeight_down")

    return rdf, recordedModules, branchArray
