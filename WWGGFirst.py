from bamboo.analysismodules import AnalysisModule, HistogramsModule
import logging

logger = logging.getLogger(__name__)

class CMSPhase2SimRTBModule(AnalysisModule):
    """ Base module for processing Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimRTBModule, self).__init__(args)
        self._h_genwcount = {}
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import decorateCMSPhase2SimTree
        from bamboo.dataframebackend import DataframeBackend
        t = decorateCMSPhase2SimTree(tree, isMC=True)
        be, noSel = DataframeBackend.create(t)
        from bamboo.root import gbl
        self._h_genwcount[sample] = be.rootDF.Histo1D(
                gbl.ROOT.RDF.TH1DModel("h_count_genweight", "genweight sum", 1, 0., 1.),
                "_zero_for_stats",
                "genweight"
                )
        return t, noSel, be, tuple()
    def mergeCounters(self, outF, infileNames, sample=None):
        outF.cd()
        self._h_genwcount[sample].Write("h_count_genweight")
    def readCounters(self, resultsFile):
        return {"sumgenweight": resultsFile.Get("h_count_genweight").GetBinContent(1)}

class CMSPhase2SimRTBHistoModule(CMSPhase2SimRTBModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimRTBHistoModule, self).__init__(args)

################################
## An analysis module example ##
################################

class SnowmassExample(CMSPhase2SimRTBHistoModule):
    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("--mvaSkim", action="store_true", help="Produce MVA training skims")
        parser.add_argument("--datacards", action="store_true", help="Produce histograms for datacards")
        parser.add_argument("--mvaEval", action="store_true", help="Import MVA model and evaluate it on the dataframe")

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, SummedPlot
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op
        
        #count no of events here 
        noSel = noSel.refine("withgenweight", weight=t.genweight)
        plots = []
        #yields
        yields = CutFlowReport("yields", recursive=True, printInLog=True)
        plots.append(yields)
        yields.add(noSel, title= 'noSel')

        #H->gg 

        #selection of photons with eta in the detector acceptance
        photons = op.select(t.gamma, lambda ph : op.AND(op.abs(ph.eta)<2.5, ph.pt >25.)) 
        #sort photons by pT 
        sort_ph = op.sort(photons, lambda ph : -ph.pt)

        #selection of photons with loose ID        
        idPhotons = op.select(sort_ph, lambda ph : ph.idpass & (1<<0))
        #H->WW->2q1l1nu
       
        electrons = op.select(t.elec, lambda el : op.AND(
        el.pt > 10., op.abs(el.eta) < 2.5
        ))
        
        clElectrons = op.select(electrons, lambda el : op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(el.p4, ph.p4) < 0.4))) #apply a deltaR 
        sort_el = op.sort(clElectrons, lambda el : -el.pt)        
        idElectrons = op.select(sort_el, lambda el : el.idpass & (1<<0))  #apply loose ID   
           

        muons = op.select(t.muon, lambda mu : op.AND(
        mu.pt > 10., op.abs(mu.eta) < 2.5
        ))

        clMuons = op.select(muons, lambda mu : op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(mu.p4, ph.p4) < 0.4 )))
        sort_mu = op.sort(clMuons, lambda mu : -mu.pt)
        idMuons = op.select(sort_mu, lambda mu : mu.idpass & (1<<2)) #apply loose ID  
        #isoMuons = op.select(idMuons, lambda mu : mu.isopass & (1<<2)) #apply tight isolation 
     
        #combine leptons
        #lepton = op.combine((idElectrons,idMuons))

        #select jets with pt>25 GeV end eta in the detector acceptance
        jets = op.select(t.jetpuppi, lambda jet : op.AND(jet.pt > 30., op.abs(jet.eta) < 2.5))
         
        clJets = op.select(jets, lambda j : op.AND(
            op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(ph.p4, j.p4) < 0.4) ),
            op.NOT(op.rng_any(idElectrons, lambda el : op.deltaR(el.p4, j.p4) < 0.4) ),  
            op.NOT(op.rng_any(idMuons, lambda mu : op.deltaR(mu.p4, j.p4) < 0.4) )
        ))
        sort_jets = op.sort(clJets, lambda jet : -jet.pt)  
        idJets = op.select(sort_jets, lambda j : j.idpass & (1<<2))
        
        mGG = op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4)
        hGG = op.sum(idPhotons[0].p4, idPhotons[1].p4)
        mJets= op.invariant_mass(idJets[0].p4, idJets[1].p4)
        mJets_SL= op.invariant_mass(idJets[1].p4, idJets[2].p4)
        hJets = op.sum(idJets[0].p4, idJets[1].p4)
       
        #missing transverse energy
        met = op.select(t.metpuppi)      

        #define more variables for ease of use
        nElec = op.rng_len(idElectrons)
        nMuon = op.rng_len(idMuons)
        nJet = op.rng_len(idJets)
        nPhoton = op.rng_len(idPhotons)

        #defining more DNN variables
        pT_mGGL = op.product(idPhotons[0].pt, op.pow(mGG, -1)) 
        pT_mGGSL = op.product(idPhotons[1].pt, op.pow(mGG, -1)) 
        E_mGGL = op.product(idPhotons[0].p4.energy(), op.pow(mGG, -1))
        E_mGGSL = op.product(idPhotons[1].p4.energy(), op.pow(mGG, -1))


        #selections for efficiency check

        sel1_p = noSel.refine("2Photon", cut = op.AND((op.rng_len(sort_ph) >= 2), (sort_ph[0].pt > 35.)))

        sel2_p = sel1_p.refine("idPhoton", cut = op.AND((op.rng_len(idPhotons) >= 2), (idPhotons[0].pt > 35.)))

        sel1_e = noSel.refine("OneE", cut = op.rng_len(sort_el) >= 1)
        
        sel2_e = sel1_e.refine("idElectron", cut = op.rng_len(idElectrons) >= 1)

        sel1_m = noSel.refine("OneM", cut = op.rng_len(sort_mu) >= 1)
        
        sel2_m = sel1_m.refine("idMuon", cut = op.rng_len(idMuons) >= 1)
         
        #sel4 = sel3.refine("TwoPhLNuTwoJ", cut = op.AND((op.rng_len(cleanedJets) >= 2),(met[0].pt > 30)))

        #selection: 2 photons (at least) in an event 
        hasTwoPh = sel2_p.refine("hasTwoPh", cut= op.rng_len(idPhotons) >= 2)

        yields.add(hasTwoPh, title='hasTwoPh')

        #selections for the event inv mass of photons within the 100-180 window
        hasInvM = hasTwoPh.refine("hasInvM", cut= op.AND(
            (op.in_range(100, op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4), 180)) 
        ))
        yields.add(hasInvM, title='hasInvM')

        #selections for semileptonic final state
        hasOneL = hasInvM.refine("hasOneL", cut = op.OR(nElec == 1, nMuon == 1))
        yields.add(hasOneL, title='hasOneL')

        hasOneEl = hasInvM.refine("hasOneEl", cut = op.AND(nElec == 1, nMuon == 0))
        yields.add(hasOneEl, title='hasOneEl')

        hasOneMu = hasInvM.refine("hasOneMu", cut = op.AND(nElec == 0, nMuon == 1))
        yields.add(hasOneMu, title='hasOneMu')

        #adding jets on the semileptonic final state
        hasOneJ = hasOneL.refine("hasOneJ", cut = nJet >= 1)
        yields.add(hasOneJ, title='hasOneJ')

        hasTwoJ = hasOneJ.refine("hasTwoJ", cut = nJet >= 2)
        yields.add(hasTwoJ, title='hasTwoJ')
      
        hasThreeJ = hasTwoJ.refine("hasThreeJ", cut = nJet >= 3)
        yields.add(hasThreeJ, title='hasThreeJ')

        #plots
        
        #sel1_p
        plots.append(Plot.make1D("LeadingPhotonPTNoID", sort_ph[0].pt, sel1_p, EqB(30, 0., 300.), title="Leading Photon pT"))        
        plots.append(Plot.make1D("SubLeadingPhotonPTNoID", sort_ph[1].pt, sel1_p, EqB(30, 0., 300.), title="SubLeading Photon pT"))
       
        #sel2_p  
        plots.append(Plot.make1D("LeadingPhotonPTID", idPhotons[0].pt, sel2_p, EqB(30, 0., 300.), title="Leading Photon pT"))    
        plots.append(Plot.make1D("SubLeadingPhotonPTID", idPhotons[0].pt, sel2_p, EqB(30, 0., 300.), title="SubLeading Photon pT")) 
       
        #sel1_e
        plots.append(Plot.make1D("LeadingElectronNoID", sort_el[0].pt, sel1_e, EqB(30, 0., 300.), title="Leading Electron pT"))
        #sel2_e
        plots.append(Plot.make1D("LeadingElectronID", idElectrons[0].pt, sel2_e, EqB(30, 0., 300.), title="Leading Electron pT"))
        #sel3_e
        #plots.append(Plot.make1D("LeadingElectronNoZee", slElectrons[0].pt, sel3_e, EqB(30, 0., 300.), title="Leading Electron pT"))

        #sel1_m
        plots.append(Plot.make1D("LeadingMuonNoID", sort_mu[0].pt, sel1_m, EqB(30, 0., 300.), title="Leading Muon pT"))       
        #sel2_m
        plots.append(Plot.make1D("LeadingMuonID", idMuons[0].pt, sel2_m, EqB(30, 0., 300.), title="Leading Muon pT"))
        #sel3_m
        #plots.append(Plot.make1D("LeadingMuonIso", isoMuons[0].pt, sel3_m, EqB(30, 0., 300.), title="Leading Muon pT"))

        #hasTwoPh
        plots.append(Plot.make1D("LeadingPhotonPtTwoPh", idPhotons[0].pt, hasTwoPh, EqB(30, 0., 300.), title="Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPtTwoPh", idPhotons[1].pt, hasTwoPh, EqB(30, 0., 300.), title="SubLeading Photon pT"))
        plots.append(Plot.make1D("nElectronsTwoPh", nElec, hasTwoPh, EqB(10, 0., 10.), title="Number of electrons"))
        plots.append(Plot.make1D("nMuonsTwoPh", nMuon, hasTwoPh, EqB(10, 0., 10.), title="Number of Muons"))
        plots.append(Plot.make1D("nJetsTwoPh", nJet, hasTwoPh, EqB(10, 0., 10.), title="Number of Jets"))
        plots.append(Plot.make1D("nPhotonsTwoPh", nPhoton, hasTwoPh, EqB(10, 0., 10.), title="Number of Photons"))
        plots.append(Plot.make1D("Inv_mass_gghasTwoPh",mGG,hasTwoPh,EqB(50, 100.,150.), title = "m_{\gamma\gamma}"))
        #plots.append(Plot.make1D("LeadingJetPtTwoPh", idJets[0].pt, hasTwoPh, EqB(10, 0., 10.), title = 'Leading Jet pT'))

        #hasInvM
        plots.append(Plot.make1D("LeadingPhotonPtInvM", idPhotons[0].pt, hasInvM, EqB(30, 0., 300.), title="Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPtInvM", idPhotons[1].pt, hasInvM, EqB(30, 0., 300.), title="SubLeading Photon pT"))
        plots.append(Plot.make1D("nElectronsInvM", nElec, hasInvM, EqB(10, 0., 10.), title="Number of electrons"))
        plots.append(Plot.make1D("nMuonsInvM", nMuon, hasInvM, EqB(10, 0., 10.), title="Number of Muons"))
        plots.append(Plot.make1D("nJetsInvM", nJet, hasInvM, EqB(10, 0., 10.), title="Number of Jets"))
        plots.append(Plot.make1D("Inv_mass_gghasInvM",mGG,hasInvM,EqB(50, 100.,150.), title = "m_{\gamma\gamma}"))
        #plots.append(Plot.make1D("LeadingJetPtInvM", idJets[0].pt, hasInvM, EqB(10, 0., 10.), title = 'Leading Jet pT'))

        #hasOneL
        plots.append(Plot.make1D("LeadingPhotonPtOneL", idPhotons[0].pt, hasOneL, EqB(30, 0., 300.), title="Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPtOneL", idPhotons[1].pt, hasOneL, EqB(30, 0., 300.), title="SubLeading Photon pT"))
        plots.append(Plot.make1D("LeadingPhotonEtaOneL", idPhotons[0].eta, hasOneL, EqB(80, -4., 4.), title="Leading Photon eta"))
        plots.append(Plot.make1D("SubLeadingPhotonEtaOneL", idPhotons[1].eta, hasOneL, EqB(80, -4., 4.), title="SubLeading Photon eta"))
        plots.append(Plot.make1D("LeadingPhotonPhiOneL", idPhotons[0].phi, hasOneL, EqB(100, -3.5, 3.5), title="Leading Photon phi"))
        plots.append(Plot.make1D("SubLeadingPhotonPhiOneL", idPhotons[1].phi, hasOneL, EqB(100, -3.5, 3.5), title="SubLeading Photon phi"))
        plots.append(Plot.make1D("nElectronsOneL", nElec, hasOneL, EqB(10, 0., 10.), title="Number of electrons"))
        plots.append(Plot.make1D("nMuonsOneL", nMuon, hasOneL, EqB(10, 0., 10.), title="Number of Muons"))
        plots.append(Plot.make1D("nJetsOneL", nJet, hasOneL, EqB(10, 0., 10.), title="Number of Jets"))
        plots.append(Plot.make1D("Inv_mass_gghasOneL",mGG , hasOneL, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
        plots.append(Plot.make1D("LeadingPhotonpT_mGGLhasOneL", pT_mGGL, hasOneL,EqB(100, 0., 5.) ,title = "Leading Photon p_{T}/m_{\gamma\gamma}"))  
        plots.append(Plot.make1D("SubLeadingPhotonpT_mGGLhasOneL", pT_mGGSL, hasOneL,EqB(100, 0., 5.) ,title = "SubLeading Photon p_{T}/m_{\gamma\gamma}"))
        plots.append(Plot.make1D("LeadingPhotonE_mGGLhasOneL", E_mGGL, hasOneL,EqB(100, 0., 5.) ,title = "Leading Photon E/m_{\gamma\gamma}"))
        plots.append(Plot.make1D("SubLeadingPhotonE_mGGLhasOneL", E_mGGSL, hasOneL,EqB(100, 0., 5.) ,title = "SubLeading Photon E/m_{\gamma\gamma}")) 
        plots.append(Plot.make1D("MET", met, hasOneL,EqB(100, 0., 50.) ,title="MET"))
        
        #Lepton Plots
        ElectronpT = Plot.make1D("ElectronpT", idElectrons[0].pt, hasOneEl, EqB(30, 0., 300.), title = 'Leading Electron pT')
        MuonpT = Plot.make1D("MuonpT", idMuons[0].pt, hasOneMu, EqB(30, 0., 300.), title = 'Leading Muon pT')
        LeptonpT = SummedPlot('LeptonpT', 
                                [ElectronpT,MuonpT],
                                xTitle = 'Leading Lepton pT')
        plots.append(ElectronpT)
        plots.append(MuonpT)
        plots.append(LeptonpT)

        ElectronE = Plot.make1D("ElectronE", idElectrons[0].p4.E(), hasOneEl, EqB(50, 0., 500.), title = 'Leading Electron E')
        MuonE = Plot.make1D("MuonE", idMuons[0].p4.E(), hasOneMu, EqB(50, 0., 500.), title = 'Leading Muon E')
        LeptonE = SummedPlot('LeptonE', 
                                [ElectronE,MuonE],
                                xTitle = 'Leading Lepton E')
        plots.append(ElectronE)
        plots.append(MuonE)
        plots.append(LeptonE)
    
        ElectronEta = Plot.make1D("ElectronEta", idElectrons[0].eta, hasOneEl, EqB(80, -4., 4.), title = 'Leading Electron eta')
        MuonEta = Plot.make1D("MuonEta", idMuons[0].eta, hasOneMu, EqB(80, -4., 4.), title = 'Leading Muon eta')
        LeptonEta = SummedPlot('LeptonEta', 
                                [ElectronEta,MuonEta],
                                xTitle = 'Leading Lepton Eta')
        plots.append(ElectronEta)
        plots.append(MuonEta)
        plots.append(LeptonEta)

        ElectronPhi = Plot.make1D("ElectronPhi", idElectrons[0].phi, hasOneEl, EqB(100, -3.5, 3.5), title = 'Leading Electron phi')
        MuonPhi = Plot.make1D("MuonPhi", idMuons[0].phi, hasOneMu, EqB(100, -3.5, 3.5), title = 'Leading Muon phi')
        LeptonPhi = SummedPlot('LeptonPhi', 
                                [ElectronPhi,MuonPhi],
                                xTitle = 'Leading Lepton Phi')
        plots.append(ElectronPhi)
        plots.append(MuonPhi)
        plots.append(LeptonPhi)

        #hasOneJ
        plots.append(Plot.make1D("LeadingPhotonPtOneJ", idPhotons[0].pt, hasOneJ, EqB(30, 0., 300.), title="Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPtOneJ", idPhotons[1].pt, hasOneJ, EqB(30, 0., 300.), title="SubLeading Photon pT"))
        plots.append(Plot.make1D("nElectronsOneJ", nElec, hasOneJ, EqB(10, 0., 10.), title="Number of electrons"))
        plots.append(Plot.make1D("nMuonsOneJ", nMuon, hasOneJ, EqB(10, 0., 10.), title="Number of Muons"))
        plots.append(Plot.make1D("nJetsOneJ", nJet, hasOneJ, EqB(10, 0., 10.), title="Number of Jets"))
        plots.append(Plot.make1D("Inv_mass_ggOneJ",mGG , hasOneJ, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
        plots.append(Plot.make1D("LeadingJetPtOneJ", idJets[0].pt, hasOneJ, EqB(30, 0., 300.), title = 'Leading Jet pT'))
        plots.append(Plot.make1D("LeadingJetEtaOneJ", idJets[0].eta, hasOneJ, EqB(80, -4., 4.), title="Leading Jet eta"))
        plots.append(Plot.make1D("LeadingJetPhiOneJ", idJets[0].phi, hasOneJ, EqB(100, -3.5, 3.5), title="Leading Jet phi"))
        plots.append(Plot.make1D("LeadingJetEOnej", idJets[0].p4.energy(), hasOneJ, EqB(50, 0.,500.), title = 'Leading Jet E'))
        
        #hasTwoJ
        plots.append(Plot.make1D("LeadingPhotonPtTwoJ", idPhotons[0].pt, hasTwoJ, EqB(30, 0., 300.), title="Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPtTwoJ", idPhotons[1].pt, hasTwoJ, EqB(30, 0., 300.), title="SubLeading Photon pT"))
        plots.append(Plot.make1D("nElectronsTwoJ", nElec, hasTwoJ, EqB(10, 0., 10.), title="Number of electrons"))
        plots.append(Plot.make1D("nMuonsOneTwoJ", nMuon, hasTwoJ, EqB(10, 0., 10.), title="Number of Muons"))
        plots.append(Plot.make1D("nJetsOneTwoJ", nJet, hasTwoJ, EqB(10, 0., 10.), title="Number of Jets"))
        plots.append(Plot.make1D("LeadingJetPthasTwoJ", idJets[0].pt, hasTwoJ, EqB(30, 0., 300.), title = 'Leading Jet pT'))
        plots.append(Plot.make1D("SubLeadingJetPtTwoJ", idJets[1].pt, hasTwoJ, EqB(30, 0., 300.), title = 'SubLeading Jet pT'))
        plots.append(Plot.make1D("Inv_mass_jjTwoJ",mJets,hasTwoJ,EqB(80, 100.,180.), title = "m_{jets}"))
        plots.append(Plot.make1D("LeadingJetEtaTwoJ", idJets[0].eta, hasTwoJ, EqB(80, -4., 4.), title="Leading Jet eta"))
        plots.append(Plot.make1D("SubLeadingJetEtaTwoJ", idJets[1].eta, hasTwoJ, EqB(80, -4., 4.), title="SubLeading Jet eta"))
        plots.append(Plot.make1D("LeadingJetPhiTwoJ", idJets[0].phi, hasTwoJ, EqB(100, -3.5, 3.5), title="Leading Jet phi"))
        plots.append(Plot.make1D("SubLeadingJetPhiTwoJ", idJets[1].phi, hasTwoJ, EqB(100, -3.5, 3.5), title="SubLeading Jet phi"))
        plots.append(Plot.make1D("LeadingJetETwoJ", idJets[0].p4.energy(), hasTwoJ, EqB(50, 0.,500.), title = 'Leading Jet E'))
        plots.append(Plot.make1D("SubLeadingJetETwoJ", idJets[1].p4.energy(), hasTwoJ, EqB(50, 0.,500.), title = 'SubLeading Jet E'))

        #hasThreeJ
        plots.append(Plot.make1D("Inv_mass_jjThreeJ",mJets_SL,hasThreeJ,EqB(80, 100.,180.), title = "m_{jets}"))

        mvaVariables = {
                "weight": noSel.weight,
                "njets": nJet,
                "Eta_ph1": idPhotons[0].eta,
                "Phi_ph1": idPhotons[0].phi,
                "E_mGG_ph1": E_mGGL,
                "pT_mGG_ph1": pT_mGGL,
                "Eta_ph2": idPhotons[1].eta,
                "Phi_ph2": idPhotons[1].phi,
                "E_mGG_ph2": E_mGGSL,
                "pT_mGG_ph2": pT_mGGSL,
                "E_jet1": idJets[0].p4.E(),
                "pt_jet1": idJets[0].pt,
                "Eta_jet1": idJets[0].eta,
                "Phi_jet1": idJets[0].phi,
                "E_jet2": idJets[1].p4.E(),
                "pt_jet2": idJets[1].pt,
                "Eta_jet2": idJets[1].eta,
                "Phi_jet2": idJets[1].phi,
                #"Lepton_E": LeptonE,
                #"Lepton_Phi": LeptonPhi,
                #"Lepton_Eta": LeptonEta,
                #"Lepton_pT": LeptonpT,
                "InvM_jets1": mJets,
                "InvM_jets2": mJets_SL
                }

        #save mvaVariables to be retrieved later in the postprocessor and saved in a parquet file
        if self.args.mvaSkim or self.args.mvaEval:
            from bamboo.plots import Skim
            plots.append(Skim("allevts", mvaVariables, hasThreeJ))

        #evaluate dnn model on data
        #if self.args.mvaEval:
        #   mvaVariables.pop("weight", None)
        #   dnn = op.mvaEvaluator("./model.onnx", mvaType = "ONNXRuntime", otherArgs = "predictions")
        #   plots.append(Plot.make1D("dnn_score", dnn(*mvaVariables.values()),hasTwoPhTwoB,EqB(20, 1, 1.)))
        
        return plots

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(SnowmassExample, self).postProcess(taskList, config=config, workdir=workdir, resultsdir=resultsdir)
        from bamboo.plots import Plot, DerivedPlot
        plotList = [ ap for ap in self.plotList if ( isinstance(ap, Plot) or isinstance(ap, DerivedPlot) ) ]
        from bamboo.analysisutils import loadPlotIt
        p_config, samples, plots, systematics, legend = loadPlotIt(config, plotList, eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
        
        #mvaSkim 
        import os.path 
        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]
        if self.args.mvaSkim and skims:
            from bamboo.analysisutils import loadPlotIt
            p_config, samples, _, systematics, legend = loadPlotIt(config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
            try:
                from bamboo.root import gbl
                import pandas as pd
                for skim in skims:
                    frames = []
                    for smp in samples:
                        for cb in (smp.files if hasattr(smp, "files") else [smp]):  # could be a helper in plotit
                            # Take specific columns
                            tree = cb.tFile.Get(skim.treeName)
                            if not tree:
                                print( f"KEY TTree {skim.treeName} does not exist, we are gonna skip this {smp}\n")
                            else:
                                N = tree.GetEntries()
                            cols = gbl.ROOT.RDataFrame(cb.tFile.Get(skim.treeName)).AsNumpy()
                            cols["weight"] *= cb.scale
                            cols["process"] = [smp.name]*len(cols["weight"])
                            frames.append(pd.DataFrame(cols))
                    df = pd.concat(frames)
                    df["process"] = pd.Categorical(df["process"], categories=pd.unique(df["process"]), ordered=False)
                    pqoutname = os.path.join(resultsdir, f"{skim.name}.parquet")
                    df.to_parquet(pqoutname)
                    logger.info(f"Dataframe for skim {skim.name} saved to {pqoutname}")
            except ImportError as ex:
                logger.error("Could not import pandas, no dataframes will be saved")
        
        #produce histograms "with datacard conventions"
        if self.args.datacards:
            datacardPlots = [ap for ap in self.plotList if ap.name == "Empty_histo" or ap.name =="Inv_mass_gg" or ap.name =="Inv_mass_bb" or ap.name =="Inv_mass_HH" or (self.args.mvaEval and ap.name =="dnn_score")]
            p_config, samples, plots_dc, systematics, legend = loadPlotIt(
                config, datacardPlots, eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
            dcdir = os.path.join(workdir, "datacard_histograms")
            import os
            import numpy as np
            os.makedirs(dcdir, exist_ok=True)
            def _saveHist(obj, name, tdir=None):
                if tdir:
                    tdir.cd()
                obj.Write(name)
            from functools import partial
            import plotit.systematics
            from bamboo.root import gbl
            
            for era in (self.args.eras[1] or config["eras"].keys()):
                f_dch = gbl.TFile.Open(os.path.join(dcdir, f"histo_for_combine_{era}.root"), "RECREATE")
                saveHist = partial(_saveHist, tdir=f_dch)
                smp = next(smp for smp in samples if smp.cfg.type == "SIGNAL")
                plot =  next(plot for plot in plots_dc if plot.name == "Empty_histo")
                h = smp.getHist(plot, eras=era)
                saveHist(h.obj, f"data_obs")
        
                for plot in plots_dc:   
                    if plot.name != "Empty_histo":
                       for smp in samples:
                           smpName = smp.name
                           if smpName.endswith(".root"):
                               smpName = smpName[:-5]
                           h = smp.getHist(plot, eras=era)
                           saveHist(h.obj, f"h_{plot.name}_{smpName}")
            
            f_dch.Close()    