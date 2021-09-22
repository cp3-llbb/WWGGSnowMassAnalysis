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
     
        cleanedPhotons = op.select(photons, lambda ph : ph.idpass & (1<<0))
 
        #sort photons by pT 
        sort_ph = op.sort(photons, lambda ph : -ph.pt)
        #apply pt cut on the leading photon
        #sorted_ph = op.select(sort_ph, lambda ph : op.AND(ph[0].pt > 35.))  

        #sortcleanphotons
        sorted_ph = op.sort(cleanedPhotons, lambda ph : -ph.pt)

        #selection: 2 photons (at least) in an event with invariant mass within [100,150]
        hasTwoPh = noSel.refine("hasMassPhPh", cut= op.AND(
           (op.rng_len(sort_ph) >= 2), 
           (op.in_range(100, op.invariant_mass(sort_ph[0].p4, sort_ph[1].p4), 180)) 
           ))

        mGG = op.invariant_mass(sort_ph[0].p4, sort_ph[1].p4)
        hGG = op.sum(sort_ph[0].p4, sort_ph[1].p4)


       #selections

        sel1 = noSel.refine("LPhoton", cut = op.AND((op.rng_len(sort_ph) >= 2), (sort_ph[0].pt > 35.)))

        sel2 = noSel.refine("Lphoton", cut = op.AND((op.rng_len(sorted_ph) >= 2), (sorted_ph[0].pt > 35.)))

    #plots
       
       #sel1
        plots.append(Plot.make1D("LeadingPhotonPTSel1", sort_ph[0].pt, sel1, EqB(30, 0., 250.), title="Leading Photon pT"))
            
        plots.append(Plot.make1D("SubLeadingPhotonPTSel1", sort_ph[1].pt, sel1, EqB(30, 0., 250.), title="SubLeading Photon pT"))
       
       #sel2 
        
        plots.append(Plot.make1D("LeadingPhotonPTSel2", sorted_ph[0].pt, sel2, EqB(30, 0., 250.), title="Leading Photon pT"))
            
        plots.append(Plot.make1D("SubLeadingPhotonPTSel2", sorted_ph[1].pt, sel2, EqB(30, 0., 250.), title="SubLeading Photon pT")) 
       
    #yields
       # yields = CutFlowReport("yields", recursive=True, printInLog=True)
       # plots.append(yields)
       # yields.add(noSel, title= 'noSel')
       # yields.add(sel1, title='sel1')
       # yields.add(sel2, title='sel2')
       #yields.add(hasTwoPhTwoB, title='2_photons_2_b')


        return plots
