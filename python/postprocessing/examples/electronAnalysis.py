#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class ElectronAnalysis(Module):
    def __init__(self):
        self.writeHistFile = True

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.h_ele_pt = ROOT.TH1F("ele_pt",";p_{T} [GeV];entries",20,0,10)
        self.addObject(self.h_ele_pt)

    def analyze(self, event):
        all_leptons = Collection(event, "LepAll")
        electrons = filter( lambda lep: abs(lep.pdgId)==11, all_leptons)
        muons = filter( lambda lep: abs(lep.pdgId)==13, all_leptons)

        for ele in electrons:
            self.h_ele_pt.Fill( ele.pt )

        # return True to keep the event if we're saving another tree
        return True


preselection = "GenMET_pt>50"
files = ["file:test_files/Higgsino_N2N1_Chunk0.root"]
p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=[
                  ElectronAnalysis()], noOut=True, histFileName="histOut.root", histDirName="plots")
p.run()
