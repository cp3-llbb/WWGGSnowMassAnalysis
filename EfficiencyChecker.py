from bamboo.analysismodules import AnalysisModule, HistogramsModule

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
 
        #selection: 2 photons (at least) in an event with invariant mass within [100,150]
        hasTwoPh = noSel.refine("hasMassPhPh", cut= op.AND(
            (op.rng_len(idPhotons) >= 2), 
            (op.in_range(100, op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4), 180)) 
        ))
        
        #mGG = op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4)
        #hGG = op.sum(sorted_ph[0].p4, sorted_ph[1].p4)
        
        #H->WW->2q1l1nu
        
        electrons = op.select(t.elec, lambda el : op.AND(
            el.pt > 10., op.abs(el.eta) < 2.5
        ))
        sort_el = op.sort(electrons, lambda el : -el.pt)        
        idElectrons = op.select(sort_el, lambda el : el.idpass & (1<<0) )  #apply loose ID  
        
        clElectrons = op.select(idElectrons, lambda el : op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(el.p4, ph.p4) < 0.4))) #apply a deltaR 
        
        muons = op.select(t.muon, lambda mu : op.AND(
            mu.pt > 10., op.abs(mu.eta) < 2.5
        ))
        sort_mu = op.sort(muons, lambda mu : -mu.pt)
        idMuons = op.select(sort_mu, lambda mu : mu.idpass & (1<<0) ) #apply loose ID  
        clMuons = op.select(idMuons, lambda mu : op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(mu.p4, ph.p4) < 0.4 )))
              
        ##select jets with pt>25 GeV end eta in the detector acceptance
        jets = op.select(t.jetpuppi, lambda jet : op.AND(jet.pt > 30., op.abs(jet.eta) < 2.5))
        sort_jets = op.sort(jets, lambda jet : -jet.pt)        
        
        idJets = op.select(sort_jets, lambda j : j.idpass & (1<<2))
        clJets = op.select(idJets, lambda j : op.AND(
            op.NOT(op.rng_any(clElectrons, lambda el : op.deltaR(el.p4, j.p4) < 0.4) ),  
            op.NOT(op.rng_any(clMuons, lambda mu : op.deltaR(mu.p4, j.p4) < 0.4) )
        ))
         
        #mJets= op.invariant_mass(cleanedJets[0].p4, cleanedJets[1].p4)
        #hJets = op.sum(cleanedJets[0].p4, cleanedJets[1].p4)
       
        #missing transverse energy
        met = op.select(t.metpuppi)      

        #selections
         
        sel1_p = noSel.refine("2Photon", cut = op.AND((op.rng_len(sort_ph) >= 2), (sort_ph[0].pt > 35.)))

        sel2_p = sel1_p.refine("idPhoton", cut = op.AND((op.rng_len(idPhotons) >= 2), (idPhotons[0].pt > 35.)))

        sel1_e = noSel.refine("OneE", cut = op.AND(op.rng_len(sort_el) >= 1))
        sel2_e = sel1_e.refine("idElectron", cut = op.AND(op.rng_len(idElectrons) >= 1))

        sel1_m = noSel.refine("OneM", cut = op.AND(op.rng_len(sort_mu) >= 1))
        sel2_m = sel1_e.refine("idMuon", cut = op.AND(op.rng_len(idMuons) >= 1))

        
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

    
       ##sel4
       # plots.append(Plot.make1D("LeadingPhotonPtSel4", sorted_ph[0].pt, sel4, EqB(30, 0., 250.), title="Leading Photon pT"))
       #
       # plots.append(Plot.make1D("SubLeadingPhotonPtSel4", sorted_ph[1].pt, sel4, EqB(30, 0., 250.), title="SubLeading Photon pT"))    
       #
       ##sel5
       # plots.append(Plot.make1D("LeadingElectronPtNOID", sort_el[0].pt, noSel, EqB(30, 0., 250.), title="Leading Electron pT")) 
       # plots.append(Plot.make1D("LeadingElectronPtID", sorted_el[0].pt, noSel, EqB(30, 0., 250.), title="Leading Electron pT"))

   

       #diphoton invariant mass plot
        #plots.append(Plot.make1D("Inv_mass_ggSel1",mGG,sel1,EqB(50, 100.,150.), title = "m_{\gamma\gamma}")) #segmentation error? how?
        #plots.append(Plot.make1D("Inv_mass_ggSel2",mGG,sel2,EqB(50, 100.,150.), title = "m_{\gamma\gamma}"))
        #plots.append(Plot.make1D("Inv_mass_ggSel3",mGG,sel3,EqB(50, 100.,150.), title = "m_{\gamma\gamma}"))
        #plots.append(Plot.make1D("Inv_mass_ggSel4",mGG,sel4,EqB(50, 100.,150.), title = "m_{\gamma\gamma}"))
         
       #HH invariant mass  
       # plots.append(Plot.make1D("Inv_mass_HH",mHH,sel4,EqB(50, 200.,1000.), title = "m_{HH}"))


       #yields
        yields = CutFlowReport("yields", recursive=True, printInLog=True)
        plots.append(yields)
        yields.add(noSel, title= 'noSel')
        yields.add(sel1_p, title='sel1')
        yields.add(sel2_p, title='sel2')

        #yields.add(sel3, title='sel3')
        #yields.add(sel4, title='sel4')
        #yields.add(sel5, title='sel5')

        return plots
