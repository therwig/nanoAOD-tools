#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

def equals(a,b,err=0.005): return abs(a-b)<err

class ElectronAnalysis(Module):
    def __init__(self):
        self.writeHistFile = True

        # helper function
        f = lambda p: (abs(p.pdgId) in [11,13] # lepton
                       and p.statusFlags & (1<<13) # last copy
                       and ((p.statusFlags & (1<<10)) # tau OR
                            or ((p.statusFlags & (1<<0)) and (p.statusFlags & (1<<8))))) # prompt+ hard proc
        self.genLeptonSelector = f

    # quick add a TH1
    def bookTH1(self, hist_name, hist_title, nbins, xlo, xhi):
        h = ROOT.TH1F(hist_name, hist_title, nbins, xlo, xhi)
        setattr(self, hist_name, h)
        self.addObject( getattr(self, hist_name) )
    def bookTH2(self, hist_name, hist_title, nxbins, xlo, xhi, nybins, ylo, yhi):
        h = ROOT.TH2F(hist_name, hist_title, nxbins, xlo, xhi, nybins, ylo, yhi)
        setattr(self, hist_name, h)
        self.addObject( getattr(self, hist_name) )

    # quick add a generic object
    def bookObject(self, obj_name, obj):
        setattr(self, obj_name, obj)
        self.addObject( getattr(self, obj_name) )

    def MakeRatio(self, ratio_name, num, den):
        r=ROOT.TGraphAsymmErrors()
        r.Divide(num, den, "b(1,1)mode")
        r.SetName(ratio_name)
        r.GetXaxis().SetTitle( num.GetXaxis().GetTitle() )
        r.GetYaxis().SetTitle( "Efficiency" )
        self.bookObject(ratio_name, r)
        
    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        #
        # reco electrons
        self.bookTH1("h_ele_pt",";Reco electron p_{T} [GeV];entries",20,0,10)
        self.bookObject("h_ele_pt_eta", ROOT.TH2F("h_ele_pt_eta",";Reco electron p_{T} [GeV];Reco electron #eta", 20,0,10, 10,-2.7,2.7))

        # reco electrons matching truth
        self.bookTH1("h_ele_pt_match",";Reco electron p_{T} [GeV];entries",20,0,10)

        # truth electrons
        self.bookTH1("h_gen_ele_pt",";Gen electron p_{T} [GeV];entries",20,0,10)

        # truth electrons matching reco
        self.bookTH1("h_gen_ele_pt_match",";Gen electron p_{T} [GeV];entries",20,0,10)

        #
        # reco muons
        self.bookTH1("h_mu_pt",";Reco muon p_{T} [GeV];entries",20,0,10)
        self.bookObject("h_mu_pt_eta", ROOT.TH2F("h_mu_pt_eta",";Reco muon p_{T} [GeV];Reco muon #eta", 20,0,10, 10,-2.7,2.7))

        # reco muons matching truth
        self.bookTH1("h_mu_pt_match",";Reco muon p_{T} [GeV];entries",20,0,10)

        # truth muons
        self.bookTH1("h_gen_mu_pt",";Gen muon p_{T} [GeV];entries",20,0,10)

        # truth muons matching reco
        self.bookTH1("h_gen_mu_pt_match",";Gen muon p_{T} [GeV];entries",20,0,10)

        #
        # event level info
        #
        self.bookTH1("h_gen_mu2_pt",";Sub-leading gen muon p_{T} [GeV];entries",20,0,10)
        self.bookTH1("h_gen_mu2_pt_match",";Sub-leading gen muon p_{T} [GeV];entries",20,0,10)
        self.bookTH1("h_gen_ele2_pt",";Sub-leading gen electron p_{T} [GeV];entries",20,0,10)
        self.bookTH1("h_gen_ele2_pt_match",";Sub-leading gen electron p_{T} [GeV];entries",20,0,10)

        self.bookTH1("h_dM","; dM;entries", 40,0,20)
        self.bookTH1("h_nGen_leptons","; Gen lepton multiplicity;entries", 6,-0.5,5.5)
        self.bookTH1("h_gen_pt1_over_dm","; Gen lepton p_{T,1}/dM; entries", 30,0,1.5)
        self.bookTH1("h_gen_pt2_over_dm","; Gen lepton p_{T,2}/dM; entries", 30,0,1.5)
        self.bookTH2("h_gen_pt1pt2_over_dm","; Gen lepton p_{T,1}/dM; Gen lepton p_{T,2}/dM;", 30,0,1.5, 30,0,1.5)

        self.bookTH1("h_mll","; mll;entries", 40,0,20)
        self.bookTH1("h_gen_pt1_over_mll","; Gen lepton p_{T,1}/mll; entries", 30,0,3)
        self.bookTH1("h_gen_pt2_over_mll","; Gen lepton p_{T,2}/mll; entries", 30,0,3)
        self.bookTH2("h_gen_pt1pt2_over_mll","; Gen lepton p_{T,1}/mll; Gen lepton p_{T,2}/mll;", 30,0,3, 30,0,3)

        
    def endJob(self):
        self.MakeRatio("gen_ele_efficiency", self.h_gen_ele_pt_match, self.h_gen_ele_pt )
        self.MakeRatio("gen_mu_efficiency", self.h_gen_mu_pt_match, self.h_gen_mu_pt )

        self.MakeRatio("two_reco_mu_efficiency", self.h_gen_mu2_pt_match, self.h_gen_mu2_pt )
        self.MakeRatio("two_reco_ele_efficiency", self.h_gen_ele2_pt_match, self.h_gen_ele2_pt )
        
        Module.endJob(self)
        
    def analyze(self, event):
        # all reco leptons
        all_leptons = Collection(event, "LepAll")
        electrons = filter( lambda lep: abs(lep.pdgId)==11, all_leptons)
        muons = filter( lambda lep: abs(lep.pdgId)==13, all_leptons)

        # truth / generator-level leptons
        gen_particles = Collection(event, "GenPart")
        for ip, p in enumerate(gen_particles): p.idx = ip
        gen_leptons = filter(self.genLeptonSelector, gen_particles)
        gen_electrons = filter( lambda lep: abs(lep.pdgId)==11, gen_leptons)
        gen_muons = filter( lambda lep: abs(lep.pdgId)==13, gen_leptons)
        gen_n2s = filter( lambda lep: abs(lep.pdgId)==1000023, gen_particles)
        gen_n1s = filter( lambda lep: abs(lep.pdgId)==1000022, gen_particles)
        dM=-1
        if len(gen_n2s) and len(gen_n1s): dM=gen_n2s[0].mass-gen_n1s[0].mass

        # for signal studies, record the two leptons
        # gen_leps = []
        reco_matches = []
        
        # record truth_to_reco mapping
        truth_to_reco = dict()

        #
        # object-level analysis
        #       
        
        # reco electron loop
        for ele_idx, ele in enumerate(electrons):
            self.h_ele_pt.Fill( ele.pt )
            self.h_ele_pt_eta.Fill( ele.pt, ele.eta )
            isTruthMatch = (ele.genPartIdx >= 0 and ele.genPartFlav==1)
            if isTruthMatch:
                truth_to_reco[ele.genPartIdx] = ele_idx
                self.h_ele_pt_match.Fill( ele.pt )
                #reco_leps.append(ele)
            
        # gen electron loop
        for gen_ele in gen_electrons:
            self.h_gen_ele_pt.Fill( gen_ele.pt )
            isRecoMatch = (gen_ele.idx in truth_to_reco)
            if isRecoMatch:
                self.h_gen_ele_pt_match.Fill( gen_ele.pt )
                reco_matches.append( electrons[truth_to_reco[gen_ele.idx]] )
                # gen_leps.append(gen_ele)

        # reco muon loop
        for mu_idx, mu in enumerate(muons):
            self.h_mu_pt.Fill( mu.pt )
            self.h_mu_pt_eta.Fill( mu.pt, mu.eta )
            isTruthMatch = (mu.genPartIdx >= 0 and mu.genPartFlav==1)
            if isTruthMatch:
                truth_to_reco[mu.genPartIdx] = mu_idx
                self.h_mu_pt_match.Fill( mu.pt )
                #reco_leps.append(mu)
            
        # gen muon loop
        for gen_mu in gen_muons:
            self.h_gen_mu_pt.Fill( gen_mu.pt )
            isRecoMatch = (gen_mu.idx in truth_to_reco)
            if isRecoMatch:
                self.h_gen_mu_pt_match.Fill( gen_mu.pt )
                reco_matches.append( muons[truth_to_reco[gen_mu.idx]] )
                # gen_leps.append(gen_mu)

        #
        # event-level analysis
        #
        # print "\nEvt", [p.pdgId for p in gen_leptons], [(p.pt,p.eta,p.phi) for p in gen_leptons]
        
        if len(gen_electrons)==2:
            pt2 = min(gen_electrons[0].pt, gen_electrons[1].pt)
            self.h_gen_ele2_pt.Fill( pt2 )
            if len(reco_matches)==2:
                self.h_gen_ele2_pt_match.Fill( pt2 )
        if len(gen_muons)==2:
            pt2 = min(gen_muons[0].pt, gen_muons[1].pt)
            self.h_gen_mu2_pt.Fill( pt2 )
            if len(reco_matches)==2:
                self.h_gen_mu2_pt_match.Fill( pt2 )

        self.h_dM.Fill(dM)
        self.h_nGen_leptons.Fill( len(gen_leptons) )
        if len(gen_leptons)==2:
            pt1 = max( gen_leptons[0].pt, gen_leptons[1].pt )
            pt2 = min( gen_leptons[0].pt, gen_leptons[1].pt )
            self.h_gen_pt1_over_dm.Fill(pt1/dM)
            self.h_gen_pt2_over_dm.Fill(pt2/dM)
            self.h_gen_pt1pt2_over_dm.Fill(pt1/dM,pt2/dM)

            mll = ( gen_leptons[0].p4() + gen_leptons[1].p4() ).M()
            self.h_mll.Fill(mll)
            self.h_gen_pt1_over_mll.Fill(pt1/mll)
            self.h_gen_pt2_over_mll.Fill(pt2/mll)
            self.h_gen_pt1pt2_over_mll.Fill(pt1/mll,pt2/mll)
            
        # self.bookTH1("h_nGen_leptons","; Gen lepton multiplicity;entries", 6,-0.5,5.5)
        # self.bookTH1("h_gen_pt1_over_dm","; Gen lepton p_{T,1}/dM; entries", 22,0,1.1)
        # self.bookTH1("h_gen_pt2_over_dm","; Gen lepton p_{T,2}/dM; entries", 22,0,1.1)
        # self.bookTH2("h_gen_pt1pt2_over_dm","; Gen lepton p_{T,1}/dM; Gen lepton p_{T,2}/dM;", 22,0,1.1, 22,0,1.1)
        
        # return True to keep the event if we're saving another tree
        return True


preselection = "GenMET_pt>100"
files = ["file:test_files/Higgsino_N2N1.root"]
p = PostProcessor(".", files, cut=preselection, branchsel=None, modules=[
                  ElectronAnalysis()], noOut=True, histFileName="histOut.root", histDirName="plots")
#p.maxEntries = 100
p.run()
