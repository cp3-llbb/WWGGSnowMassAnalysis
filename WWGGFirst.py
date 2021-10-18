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
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op
        
        #count no of events here 

        noSel = noSel.refine("withgenweight", weight=t.genweight)

        plots = []

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
        idMuons = op.select(sort_mu, lambda mu : mu.idpass & (1<<2)) #apply tight ID  
        isoMuons = op.select(idMuons, lambda mu : mu.isopass & (1<<2)) #apply tight isolation 
     

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
        hJets = op.sum(idJets[0].p4, idJets[1].p4)
       
        #missing transverse energy
        met = op.select(t.metpuppi)      

        #define more variables for ease of use
        nElec = op.rng_len(idElectrons)
        nMuon = op.rng_len(idMuons)
        nJet = op.rng_len(idJets)
        nPhoton = op.rng_len(idPhotons)

        #selections for efficiency check

        sel1_p = noSel.refine("2Photon", cut = op.AND((op.rng_len(sort_ph) >= 2), (sort_ph[0].pt > 35.)))

        sel2_p = sel1_p.refine("idPhoton", cut = op.AND((op.rng_len(idPhotons) >= 2), (idPhotons[0].pt > 35.)))

        sel1_e = noSel.refine("OneE", cut = op.AND(op.rng_len(sort_el) >= 1))
        
        sel2_e = sel1_e.refine("idElectron", cut = op.AND(op.rng_len(idElectrons) >= 1))

        sel1_m = noSel.refine("OneM", cut = op.AND(op.rng_len(sort_mu) >= 1))
        
        sel2_m = sel1_e.refine("idMuon", cut = op.AND(op.rng_len(idMuons) >= 1))
         
        #sel4 = sel3.refine("TwoPhLNuTwoJ", cut = op.AND((op.rng_len(cleanedJets) >= 2),(met[0].pt > 30)))


        #selections for final desired final state
      
        #sel3 = sel2_p.refine("TwoPhOneL", cut = op.OR(op.rng_len(clElectrons) == 1, op.rng_len(clMuons) == 1))

        #sel4 = sel3.refine("TwoPhOneLTwoJ", cut = op.rng_len(clJets) >= 2)

        #selection: 2 photons (at least) in an event 
        hasTwoPh = sel2_p.refine("hasTwoPh", cut= op.AND(
            (op.rng_len(idPhotons) >= 2)
            #(op.in_range(100, op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4), 180)) 
        ))

        #selections for the event inv mass of photons within the 100-180 window
        hasInvM = hasTwoPh.refine("hasInvM", cut= op.AND(
            (op.in_range(100, op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4), 180)) 
        ))

        #selections for semileptonic final state
        hasOneL = hasInvM.refine("hasOneL", cut = op.AND(op.OR(nElec == 1, nMuon == 1)))

        #adding two jets on the semileptonic final state
        #hasTwoJ = hasOneL.refine("hasTwoJ", cut = op.rng_len(idJets) >= 2)

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

       #sel1_m
        plots.append(Plot.make1D("LeadingMuonNoID", sort_mu[0].pt, sel1_m, EqB(30, 0., 300.), title="Leading Muon pT"))
       
       #sel2_m
        plots.append(Plot.make1D("LeadingMuonID", idMuons[0].pt, sel2_m, EqB(30, 0., 300.), title="Leading Muon pT"))

       #sel3
        #plots.append(Plot.make1D("LeadingPhotonPtSel3", idPhotons[0].pt, sel3, EqB(30, 0., 250.), title="Leading Photon pT"))
        #plots.append(Plot.make1D("SubLeadingPhotonPtSel3", idPhotons[1].pt, sel3, EqB(30, 0., 250.), title="SubLeading Photon pT"))
    

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
        plots.append(Plot.make1D("nElectronsOneL", nElec, hasOneL, EqB(10, 0., 10.), title="Number of electrons"))
        plots.append(Plot.make1D("nMuonsOneL", nMuon, hasOneL, EqB(10, 0., 10.), title="Number of Muons"))
        plots.append(Plot.make1D("nJetsOneL", nJet, hasOneL, EqB(10, 0., 10.), title="Number of Jets"))
        plots.append(Plot.make1D("Inv_mass_gghasOneL",mGG , hasOneL, EqB(50, 100.,180.), title = "m_{\gamma\gamma}"))
        #plots.append(Plot.make1D("LeadingJetPtOneL", idJets[0].pt, hasOneL, EqB(10, 0., 10.), title = 'Leading Jet pT'))

        #hasTwoJ
        #plots.append(Plot.make1D("LeadingPhotonPtTwoJ", idPhotons[0].pt, hasTwoJ, EqB(30, 0., 300.), title="Leading Photon pT"))
        #plots.append(Plot.make1D("SubLeadingPhotonPtTwoJ", idPhotons[1].pt, hasTwoJ, EqB(30, 0., 300.), title="SubLeading Photon pT"))
        #plots.append(Plot.make1D("nElectronsTwoJ", nElec, hasTwoJ, EqB(10, 0., 10.), title="Number of electrons"))
        #plots.append(Plot.make1D("nMuonsOneTwoJ", nMuon, hasTwoJ, EqB(10, 0., 10.), title="Number of Muons"))
        #plots.append(Plot.make1D("nJetsOneTwoJ", nJet, hasTwoJ, EqB(10, 0., 10.), title="Number of Jets"))
        #plots.append(Plot.make1D("Inv_mass_gghasTwoJ",mGG , hasTwoJ, EqB(50, 100.,180.), title = "m_{\gamma\gamma}"))
        #plots.append(Plot.make1D("LeadingJetPtTwoJ", idJets[0].pt, hasTwoJ, EqB(10, 0., 10.), title = 'Leading Jet pT'))

       #HH invariant mass  
       # plots.append(Plot.make1D("Inv_mass_HH",mHH,sel4,EqB(50, 200.,1000.), title = "m_{HH}"))


        #yields
        yields = CutFlowReport("yields", recursive=True, printInLog=True)
        plots.append(yields)
        yields.add(noSel, title= 'noSel')
        #yields.add(sel1_p, title='sel1_p')
        #yields.add(sel2_p, title='sel2_p')
        #yields.add(sel1_e, title='sel1_e')
        #yields.add(sel2_e, title='sel2_e')
        #yields.add(sel1_m, title='sel1_m')
        #yields.add(sel2_m, title='sel2_m')
        yields.add(hasTwoPh, title='hasTwoPh')
        yields.add(hasInvM, title='hasInvM')
        yields.add(hasOneL, title='hasOneL')
        #yields.add(hasTwoJ, title='hasTwoJ')

        return plots

