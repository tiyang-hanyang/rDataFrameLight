import ROOT

# according to the tutorial https://gitlab.cern.ch/cms-analysis/jme/jerc-application-tutorial/-/blob/master/ApplyOnNanoAOD/
# 2024 data only has the L2 & L3 JES correction
# 2024 MC has L2 correction and JER correction
def processing(rdf, recordedModules, branchArray, era, ds=""):
    if "JetVetoMap.C" not in recordedModules:
        ROOT.gInterpreter.ProcessLine('#include "JetVetoMap.C"')
        recordedModules.append("JetVetoMap.C")
        ROOT.gInterpreter.ProcessLine('JVM_init("'+era+'")')

    # then this PassJetVeto branch is not really needed by this, but instead a sign already taken the JVM corrections
    if "JVMweight" not in branchArray:
        rdf = rdf.Define("PassJetVeto", "passJetVetoFunc(Jet_eta, Jet_phi)")
        # tight ID in certain eta range, common requirement in loose jet and selected jet
        if "2024" in era:
            rdf = rdf.Define("JetIdTight", "(abs(Jet_eta)<2.5) && (Jet_neHEF<0.99) && (Jet_neEmEF<0.90) && (Jet_nConstituents>1) && (Jet_muEF<0.80) && (Jet_chHEF>0.01) && (Jet_chMultiplicity>0) && (Jet_chEmEF<0.80)")
        else:
            rdf = rdf.Define("JetIdTight", "(Jet_jetId==6)")
        # loose jets is used for testing the JVM pass or fail
        rdf = rdf.Define("looseJetCond", "(Jet_pt_JEC > 15.0) && (Jet_rawFactor<0.9) && (JetIdTight) && (Jet_neEmEF + Jet_chEmEF < 0.9)")
        rdf = rdf.Define("isLooseJetFailJVM", "(looseJetCond && !(PassJetVeto))")
        rdf = rdf.Define("JVMweight", "float(Nonzero(isLooseJetFailJVM).size()==0)")
        # good jet selections
        rdf = rdf.Define("Jet_mediumPtTight", "(Jet_pt_JEC > 30.0) && (Jet_rawFactor<0.9) && (JetIdTight) && (Jet_neEmEF + Jet_chEmEF < 0.9)")
        # will return 0 if does not pass Jet_mediumPtTight, so Jet_drFromMuon>0.4 including the good selections.
        rdf = rdf.Define("Jet_drFromMuon", "minDistanceFromMuon(Jet_mediumPtTight, Jet_eta, Jet_phi, leadingMuonIdx, subleadingMuonIdx, Muon_eta, Muon_phi)")
        rdf = rdf.Define("GoodJetCond", "(Jet_drFromMuon>0.4)")
        # and also not overlapping with the muon
        branchArray.append("PassJetVeto")
        branchArray.append("JetIdTight")
        branchArray.append("looseJetCond")
        branchArray.append("isLooseJetFailJVM")
        branchArray.append("JVMweight")
        branchArray.append("Jet_mediumPtTight")
        branchArray.append("Jet_drFromMuon")
        branchArray.append("GoodJetCond")

    return rdf, recordedModules, branchArray