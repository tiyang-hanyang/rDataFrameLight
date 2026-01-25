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
        rdf = rdf.Define("PUWeight", "PUReweightFunc(Pileup_nTrueInt)")
        branchArray.append("PUWeight")

    return rdf, recordedModules, branchArray