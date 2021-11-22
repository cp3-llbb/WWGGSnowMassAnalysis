import logging
from bamboo.analysisutils import loadPlotIt
import os.path
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
            gbl.ROOT.RDF.TH1DModel("h_count_genweight",
                                   "genweight sum", 1, 0., 1.),
            "_zero_for_stats",
            "genweight"
        )
        return t, noSel, be, tuple()

    def mergeCounters(self, outF, infileNames, sample=None):
        outF.cd()
        self._h_genwcount[sample].Write("h_count_genweight")

    def readCounters(self, resultsFile):
        return {"sumgenweight": resultsFile.Get("h_count_genweight").GetBinContent(1)}

# BEGIN cutflow reports, adapted from bamboo.analysisutils


logger = logging.getLogger(__name__)

_yieldsTexPreface = "\n".join(f"{ln}" for ln in
                              r"""\documentclass[12pt, landscape]{article}
\usepackage[margin=0.2in, a3paper]{geometry}
\begin{document}
""".split("\n"))


def _texProcName(procName):
    if ">" in procName:
        procName = procName.replace(">", "$ > $")
    if "=" in procName:
        procName = procName.replace("=", "$ = $")
    if "_" in procName:
        procName = procName.replace("_", "\_")
    return procName


def _makeYieldsTexTable(MCevents, report, samples, entryPlots, stretch=1.5, orientation="v", align="c", yieldPrecision=1, ratioPrecision=2):
    if orientation not in ("v", "h"):
        raise RuntimeError(
            f"Unsupported table orientation: {orientation} (valid: 'h' and 'v')")
    import plotit.plotit
    from plotit.plotit import Stack
    import numpy as np
    from itertools import repeat, count

    def getHist(smp, plot):
        try:
            h = smp.getHist(plot)
            h.contents  # check
            return h
        except KeyError:
            return None

    def colEntriesFromCFREntryHists(report, entryHists, precision=1, showUncert=True):
        stacks_t = []
        colEntries = []
        for entries in report.titles.values():
            s_entries = []
            for eName in entries:
                eh = entryHists[eName]
                if eh is not None:
                    if (not isinstance(eh, Stack)) or eh.entries:
                        s_entries.append(eh)
            st_t = Stack(entries=s_entries)
            if s_entries:
                uncert = " \pm {{:.{}f}}".format(precision).format(
                    np.sqrt(st_t.sumw2+st_t.syst2)[1]) if showUncert else ""
                colEntries.append("${{0:.2e}}$".format(
                    precision).format(st_t.contents[1]))
                stacks_t.append(st_t)
            else:
                colEntries.append("---")
                stacks_t.append(None)
        return stacks_t, colEntries

    def colEntriesFromCFREntryHists_forEff(report, entryHists, precision=1, showUncert=True):
        stacks_t = []
        colEntries = []
        for entries in report.titles.values():  # selection names
            s_entries = []
            for eName in entries:
                eh = entryHists[eName]
                if eh is not None:
                    if (not isinstance(eh, Stack)) or eh.entries:
                        s_entries.append(eh)
            st_t = Stack(entries=s_entries)
            if s_entries:
                uncert = " \pm {{:.{}f}}".format(precision).format(
                    np.sqrt(st_t.sumw2+st_t.syst2)[1]) if showUncert else ""
                colEntries.append("{{0}}".format(
                    precision).format(st_t.contents[1]))
                stacks_t.append(st_t)
            else:
                colEntries.append("---")
                stacks_t.append(None)
        return stacks_t, colEntries

    smp_signal = [smp for smp in samples if smp.cfg.type == "SIGNAL"]
    smp_mc = [smp for smp in samples if smp.cfg.type == "MC"]
    smp_data = [smp for smp in samples if smp.cfg.type == "DATA"]
    sepStr = "|l|"
    smpHdrs = []
    titles = list(report.titles.keys())  # titles are selections
    entries_smp = []
    stTotSig, stTotMC, stTotData = None, None, None
    if smp_signal:
        sepStr += "|"
        sel_list = []
        for sigSmp in smp_signal:
            _, colEntries = colEntriesFromCFREntryHists(report,
                                                        {eName: getHist(sigSmp, p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            sepStr += f"{align}|"
            smpHdrs.append(
                f"${_texProcName(sigSmp.cfg.yields_group)}$")  # sigSmp.cfg.yields_group is the name in the legend
            _, colEntries_forEff = colEntriesFromCFREntryHists_forEff(report, {eName: sigSmp.getHist(
                p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            colEntries_matrix = np.array(colEntries_forEff)
            sel_eff = np.array([100])
            for i in range(1, len(report.titles)):
                sel_eff = np.append(sel_eff, [float(
                    colEntries_matrix[i]) / float(colEntries_matrix[0]) * 100]).tolist()
            for i in range(len(report.titles)):
                sel_eff[i] = str(f"({sel_eff[i]:.2f}\%)")
            colEntries_withEff = []
            for i, entry in enumerate(colEntries):
                colEntries_withEff.append("{0} {1} {2}".format(
                    entry, sel_eff[i], MCevents[sigSmp.cfg.pretty_name.rstrip(".root")][0][i]))
            entries_smp.append(colEntries_withEff)
        if len(smp_signal) > 1:
            sepStr += f"|{align}|"
            smpHdrs.append("Signal")
            stTotSig, colEntries = colEntriesFromCFREntryHists(report, {eName: Stack(entries=[h for h in (getHist(
                smp, p) for smp in smp_signal) if h]) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            entries_smp.append(colEntries)
    if smp_mc:
        sepStr += "|"
        for mcSmp in smp_mc:
            stTotMC, colEntries = colEntriesFromCFREntryHists(report,
                                                              {eName: getHist(mcSmp, p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            sepStr += f"{align}|"
            if isinstance(mcSmp, plotit.plotit.Group):
                smpHdrs.append(_texProcName(mcSmp.name))
            else:
                smpHdrs.append(f"${_texProcName(mcSmp.cfg.yields_group)}$")
            _, colEntries_forEff = colEntriesFromCFREntryHists_forEff(report, {eName: mcSmp.getHist(
                p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            colEntries_matrix = np.array(colEntries_forEff)
            sel_eff = np.array([100])
            for i in range(1, len(report.titles)):
                sel_eff = np.append(sel_eff, [float(
                    colEntries_matrix[i]) / float(colEntries_matrix[0]) * 100]).tolist()
            for i in range(len(report.titles)):
                sel_eff[i] = str(f"({sel_eff[i]:.2f}\%)")
            colEntries_withEff = []
            for i, entry in enumerate(colEntries):
                colEntries_withEff.append("{0} {1}".format(entry, sel_eff[i]))
            entries_smp.append(colEntries_withEff)
        if len(smp_mc) > 1:
            sepStr += f"|{align}|"
            smpHdrs.append("Background")
            stTotMC, colEntries = colEntriesFromCFREntryHists(report, {eName: Stack(entries=[h for h in (getHist(
                smp, p) for smp in smp_mc) if h]) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            entries_smp.append(colEntries)
    if smp_data:
        sepStr += f"|{align}|"
        smpHdrs.append("Data")
        stTotData, colEntries = colEntriesFromCFREntryHists(report, {eName: Stack(entries=[h for h in (getHist(
            smp, p) for smp in smp_data) if h]) for eName, p in entryPlots.items()}, precision=0, showUncert=False)
        entries_smp.append(colEntries)
    if smp_data and smp_mc:
        sepStr += f"|{align}|"
        smpHdrs.append("Data/MC")
        colEntries = []
        import numpy.ma as ma
        for stData, stMC in zip(stTotData, stTotMC):
            if stData is not None and stMC is not None:
                dtCont = stData.contents
                mcCont = ma.array(stMC.contents)
                ratio = dtCont/mcCont
                ratioErr = np.sqrt(mcCont**2*stData.sumw2 +
                                   dtCont**2*(stMC.sumw2+stMC.syst2))/mcCont**2
                if mcCont[1] != 0.:
                    colEntries.append("${{0:.{0}f}}$".format(
                        ratioPrecision).format(ratio[1]))
                else:
                    colEntries.append("---")
            else:
                colEntries.append("---")
        entries_smp.append(colEntries)
    c_bySmp = entries_smp
    c_byHdr = [[smpEntries[i] for smpEntries in entries_smp]
               for i in range(len(titles))]
    if orientation == "v":
        rowHdrs = titles  # selections
        colHdrs = ["Selections"]+smpHdrs  # samples
        c_byRow = c_byHdr
        c_byCol = c_bySmp
    else:  # horizontal
        sepStr = "|l|{0}|".format("|".join(repeat(align, len(titles))))
        rowHdrs = smpHdrs  # samples
        colHdrs = ["Samples"]+titles  # selections
        c_byRow = c_bySmp
        c_byCol = c_byHdr
    if entries_smp:
        colWidths = [max(len(rh) for rh in rowHdrs)+1]+[max(len(hdr), max(len(c)
                                                                          for c in col))+1 for hdr, col in zip(colHdrs[1:], c_byCol)]
        return "\n".join([
            f"\\renewcommand{{\\arraystretch}}{{{stretch}}}",
            f"\\begin{{tabular}}{{ {sepStr} }}",
            "    \\hline",
            "    {0} \\\\".format(" & ".join(h.ljust(cw)
                                  for cw, h in zip(colWidths, colHdrs))),
            "    \\hline"]+[
                "    {0} \\\\".format(" & ".join(en.rjust(cw)
                                      for cw, en in zip(colWidths, [rh]+rowEntries)))
                for rh, rowEntries in zip(rowHdrs, c_byRow)
        ]+[
            "    \\hline",
            "\\end{tabular}"
            "\\end{document}"
        ])


def printCutFlowReports(config, reportList, workdir=".", resultsdir=".", suffix=None, readCounters=lambda f: -1., eras=("all", None), verbose=False):
    """
    Print yields to the log file, and write a LaTeX yields table for each

    Samples can be grouped (only for the LaTeX table) by specifying the
    ``yields-group`` key (overriding the regular ``groups`` used for plots).
    The sample (or group) name to use in this table should be specified
    through the ``yields-title`` sample key.

    In addition, the following options in the ``plotIt`` section of
    the YAML configuration file influence the layout of the LaTeX yields table:

    - ``yields-table-stretch``: ``\\arraystretch`` value, 1.15 by default
    - ``yields-table-align``: orientation, ``h`` (default), samples in rows, or ``v``, samples in columns
    - ``yields-table-text-align``: alignment of text in table cells (default: ``c``)
    - ``yields-table-numerical-precision-yields``: number of digits after the decimal point for yields (default: 1)
    - ``yields-table-numerical-precision-ratio``: number of digits after the decimal point for ratios (default: 2)
    """
    eraMode, eras = eras
    if not eras:  # from config if not specified
        eras = list(config["eras"].keys())
    ## helper: print one bamboo.plots.CutFlowReport.Entry

    def printEntry(entry, printFun=logger.info, recursive=True, genEvents=None):
        if entry.nominal is not None:
            effMsg = ""
            if entry.parent:
                sumPass = entry.nominal.GetBinContent(1)
                sumTotal = (entry.parent.nominal.GetBinContent(
                    1) if entry.parent.nominal is not None else 0.)
                if sumTotal != 0.:
                    effMsg = f", Eff={sumPass/sumTotal:.2%}"
                    if genEvents:
                        effMsg += f", TotalEff={sumPass/genEvents:.2%}"
            printFun(
                f"Selection {entry.name}: N={entry.nominal.GetEntries()}), SumW={entry.nominal.GetBinContent(1)}{effMsg}")
            printFun(f"Selection {entry.name}: N={entry.nominal.GetEntries()}")
        if recursive:
            for c in entry.children:
                printEntry(c, printFun=printFun,
                           recursive=recursive, genEvents=genEvents)

    def unwMCevents(entry, smp, mcevents, genEvents=None):
        mcevents.append(entry.nominal.GetEntries())
        for c in entry.children:
            unwMCevents(c, smp, mcevents, genEvents=genEvents)
        return mcevents

    ## retrieve results files, get generated events for each sample
    from bamboo.root import gbl
    resultsFiles = dict()
    generated_events = dict()
    for smp, smpCfg in config["samples"].items():
        if "era" not in smpCfg or smpCfg["era"] in eras:
            resF = gbl.TFile.Open(os.path.join(resultsdir, f"{smp}.root"))
            resultsFiles[smp] = resF
            genEvts = None
            if "generated-events" in smpCfg:
                if isinstance(smpCfg["generated-events"], str):
                    genEvts = readCounters(resF)[smpCfg["generated-events"]]
                else:
                    genEvts = smpCfg["generated-events"]
            generated_events[smp] = genEvts
    has_plotit = None
    try:
        import plotit.plotit
        has_plotit = True
    except ImportError:
        has_plotit = False
    from bamboo.plots import EquidistantBinning as EqB

    class YieldPlot:
        def __init__(self, name):
            self.name = name
            self.plotopts = dict()
            self.axisTitles = ("Yield",)
            self.binnings = [EqB(1, 0., 1.)]
    for report in reportList:
        smpReports = {smp: report.readFromResults(
            resF) for smp, resF in resultsFiles.items()}
        ## debug print
        MCevents = {}
        for smp, smpRep in smpReports.items():
            #if smpRep.printInLog:
            logger.info(f"Cutflow report {report.name} for sample {smp}")
            MCevents[smp] = []
            for root in smpRep.rootEntries():
                printEntry(root, genEvents=generated_events[smp])
                mcevents = []
                MCevents[smp].append(unwMCevents(
                    root, smp, mcevents, genEvents=generated_events[smp]))
        ## save yields.tex (if needed)
        if any(len(cb) > 1 or tt != cb[0] for tt, cb in report.titles.items()):
            if not has_plotit:
                logger.error(
                    f"Could not load plotit python library, no TeX yields tables for {report.name}")
            else:
                yield_plots = [YieldPlot(f"{report.name}_{eName}")
                               for tEntries in report.titles.values() for eName in tEntries]
                out_eras = []
                if len(eras) > 1 and eraMode in ("all", "combined"):
                    nParts = [report.name]
                    if suffix:
                        nParts.append(suffix)
                    out_eras.append(("{0}.tex".format("_".join(nParts)), eras))
                if len(eras) == 1 or eraMode in ("split", "all"):
                    for era in eras:
                        nParts = [report.name]
                        if suffix:
                            nParts.append(suffix)
                        nParts.append(era)
                        out_eras.append(
                            ("{0}.tex".format("_".join(nParts)), [era]))
                for outName, iEras in out_eras:
                    pConfig, samples, plots, _, _ = loadPlotIt(
                        config, yield_plots, eras=iEras, workdir=workdir, resultsdir=resultsdir, readCounters=readCounters)
                    tabBlock = _makeYieldsTexTable(MCevents, report, samples,
                                                   {p.name[len(
                                                       report.name)+1:]: p for p in plots},
                                                   stretch=pConfig.yields_table_stretch,
                                                   orientation=pConfig.yields_table_align,
                                                   align=pConfig.yields_table_text_align,
                                                   yieldPrecision=pConfig.yields_table_numerical_precision_yields,
                                                   ratioPrecision=pConfig.yields_table_numerical_precision_ratio)
                    if tabBlock:
                        with open(os.path.join(workdir, outName), "w") as ytf:
                            ytf.write("\n".join((_yieldsTexPreface, tabBlock)))
                        logger.info("Yields table for era(s) {0} was written to {1}".format(
                            ",".join(iEras), os.path.join(workdir, outName)))
                    else:
                        logger.warning(
                            f"No samples for era(s) {','.join(iEras)}, so no yields.tex")

# END cutflow reports, adapted from bamboo.analysisutils


class CMSPhase2SimHistoModule(CMSPhase2SimRTBModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """

    def __init__(self, args):
        super(CMSPhase2SimHistoModule, self).__init__(args)

    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        """ Customised cutflow reports and plots """
        if not self.plotList:
            self.plotList = self.getPlotList(resultsdir=resultsdir)
        from bamboo.plots import Plot, DerivedPlot, CutFlowReport
        plotList_cutflowreport = [
            ap for ap in self.plotList if isinstance(ap, CutFlowReport)]
        plotList_plotIt = [ap for ap in self.plotList if (isinstance(
            ap, Plot) or isinstance(ap, DerivedPlot)) and len(ap.binnings) == 1]
        eraMode, eras = self.args.eras
        if eras is None:
            eras = list(config["eras"].keys())
        if plotList_cutflowreport:
            printCutFlowReports(config, plotList_cutflowreport, workdir=workdir, resultsdir=resultsdir,
                                readCounters=self.readCounters, eras=(eraMode, eras), verbose=self.args.verbose)
        if plotList_plotIt:
            from bamboo.analysisutils import writePlotIt, runPlotIt
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                        readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
            runPlotIt(cfgName, workdir=workdir, plotIt=self.args.plotIt,
                      eras=(eraMode, eras), verbose=self.args.verbose)

################################
  ## Actual analysis module ##
################################


class CMSPhase2Sim(CMSPhase2SimHistoModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        # count no of events here

        genweightsel = noSel.refine("withgenweight", weight=t.genweight)

        plots = []

        # select photons
        photons = op.select(t.gamma, lambda ph: op.AND(op.abs(ph.eta) < 3, op.NOT(
            op.in_range(1.442, op.abs(ph.eta), 1.566)), ph.pt > 25))

        # select loose ID photon
        looseIDPhotons = op.select(
            photons, lambda ph: ph.idpass & (1 << 0))  # looseID

        # sortIDphotons
        sortedIDphotons = op.sort(looseIDPhotons, lambda ph: -ph.pt)

        mgg = op.invariant_mass(sortedIDphotons[0].p4, sortedIDphotons[1].p4)

        # selection: at least 2 photons
        twoPhotonsSel = noSel.refine(
            "hasInvMassPhPh", cut=op.AND(op.rng_len(sortedIDphotons) >= 2, sortedIDphotons[0].pt > 35, sortedIDphotons[1].pt > 25))  # sel1

        # pT/InvM(gg) > 0.33 selection for leading photon
        pTmggRatio_sel = twoPhotonsSel.refine(
            "ptMggLeading", cut=(op.AND(op.product(sortedIDphotons[0].pt, op.pow(mgg, -1)) > 0.33), op.product(sortedIDphotons[1].pt, op.pow(mgg, -1)) > 0.25))
        # pTmggRatio_sel = pTmggRatioLeading_sel.refine(
        #     "ptMggLead_Subleading", cut=op.product(sortedIDphotons[1].pt, op.pow(mgg, -1)) > 0.25)  # sel2

        mgg_sel = pTmggRatio_sel.refine("mgg_sel", cut=[mgg > 100])  # sel3

        # electrons

        electrons = op.select(t.elec, lambda el: op.AND(op.abs(el.eta) < 3, op.NOT(
            op.in_range(1.442, op.abs(el.eta), 1.566)), el.pt > 10.))

        isolatedElectrons = op.select(
            electrons, lambda el: el.isopass & (1 << 2))

        IDelectrons = op.select(
            isolatedElectrons, lambda el: el.idpass & (1 << 0))  # loose ID

        cleanedElectrons = op.select(IDelectrons, lambda el: op.NOT(
            op.rng_any(sortedIDphotons, lambda ph: op.deltaR(el.p4, ph.p4) < 0.2)))

        # muons

        muons = op.select(t.muon, lambda mu: op.AND(
            mu.pt > 10., op.abs(mu.eta) < 3))

        isolatedMuons = op.select(muons, lambda mu: mu.isopass & (1 << 2))

        IDmuons = op.select(
            isolatedMuons, lambda mu: mu.idpass & (1 << 0))  # loose ID

        cleanedMuons = op.select(IDmuons, lambda mu: op.NOT(
            op.rng_any(sortedIDphotons, lambda ph: op.deltaR(mu.p4, ph.p4) < 0.2)))

        # taus

        taus = op.select(t.tau, lambda tau: op.AND(
            tau.pt > 20., op.abs(tau.eta) < 3))

        isolatedTaus = op.select(taus, lambda tau: tau.isopass & (1 << 2))

        cleanedTaus = op.select(isolatedTaus, lambda tau: op.AND(
            op.NOT(op.rng_any(sortedIDphotons,
                   lambda ph: op.deltaR(tau.p4, ph.p4) < 0.2)),
            op.NOT(op.rng_any(cleanedElectrons,
                   lambda el: op.deltaR(tau.p4, el.p4) < 0.2)),
            op.NOT(op.rng_any(cleanedMuons,
                   lambda mu: op.deltaR(tau.p4, mu.p4) < 0.2))
        ))

        twoTausSel = mgg_sel.refine(
            "twotausel", cut=[op.rng_len(cleanedTaus) >= 2])  # sel4

        # def nDaughters(gen):
        #     """Return the number of daughters of a given object. """
        #     return gen.d2() - gen.d1()

        # genTaus = op.select(t.genpart, lambda g: op.abs(g.pid) == 15)

        # oneGenTauSel = mgg_sel.refine("onegentau", cut = [op.rng_len(genTaus) >= 1])

        # jets

        jets = op.select(t.jetpuppi, lambda jet: op.AND(
            jet.pt > 30., op.abs(jet.eta) < 3))

        IDJets = op.select(jets, lambda j: j.idpass & (1 << 2))  # tight ID

        cleanedJets = op.select(IDJets, lambda j: op.AND(
            op.NOT(op.rng_any(cleanedElectrons,
                   lambda el: op.deltaR(j.p4, el.p4) < 0.4)),
            op.NOT(op.rng_any(cleanedMuons, lambda mu: op.deltaR(j.p4, mu.p4) < 0.4)),
            op.NOT(op.rng_any(cleanedTaus, lambda tau: op.deltaR(j.p4, tau.p4) < 0.4)),
            op.NOT(op.rng_any(sortedIDphotons,
                   lambda ph: op.deltaR(j.p4, ph.p4) < 0.4))
        ))

        btaggedJets = op.select(
            cleanedJets, lambda j: j.btag & (1 << 1))  # medium  WP

        # mJets = op.invariant_mass(cleanedJets[0].p4, cleanedJets[1].p4)
        # hJets = op.sum(cleanedJets[0].p4, cleanedJets[1].p4)

        # met = op.select(t.metpuppi)

      # selections

        # sel1 = noSel.refine("DiPhoton", cut=op.AND(
        # (op.rng_len(looseIDPhotons) >= 2), (looseIDPhotons[0].pt > 35.)))

        OneJetSel = twoTausSel.refine(
            "twojetsel", cut=op.rng_len(cleanedJets) >= 1)

        btaggedJetSel = OneJetSel.refine(
            "btaggedjet", cut=op.rng_len(btaggedJets) >= 1)

       # plots

       # sel1: twoPhotonsSel
       # sel2: pTmggRatio_sel
       # sel3: mgg_sel
       # sel4: twoTausSel

        plots.append(Plot.make1D("LeadingPhotonPTSel1", sortedIDphotons[0].pt, twoPhotonsSel, EqB(
            30, 0., 250.), title="Leading Photon pT"))
        plots.append(Plot.make1D("LeadingPhotonPTSel2", sortedIDphotons[0].pt, pTmggRatio_sel, EqB(
            30, 0., 250.), title="Leading Photon pT"))
        plots.append(Plot.make1D("LeadingPhotonPTSel3", sortedIDphotons[0].pt, mgg_sel, EqB(
            30, 0., 250.), title="Leading Photon pT"))
        plots.append(Plot.make1D("LeadingPhotonPTSel4", sortedIDphotons[0].pt, twoTausSel, EqB(
            30, 0., 250.), title="Leading Photon pT"))

        plots.append(Plot.make1D("SubLeadingPhotonPTSel1", sortedIDphotons[1].pt, twoPhotonsSel, EqB(
            30, 0., 250.), title="Sub-Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPTSel2", sortedIDphotons[1].pt, pTmggRatio_sel, EqB(
            30, 0., 250.), title="Sub-Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPTSel3", sortedIDphotons[1].pt, mgg_sel, EqB(
            30, 0., 250.), title="Sub-Leading Photon pT"))
        plots.append(Plot.make1D("SubLeadingPhotonPTSel4", sortedIDphotons[1].pt, twoTausSel, EqB(
            30, 0., 250.), title="Sub-Leading Photon pT"))

        plots.append(Plot.make1D("leadingTau_ptSel4", cleanedTaus[0].pt, twoTausSel, EqB(
            30, 0., 250.), title="Leading Tau p_{T}"))

        plots.append(Plot.make1D("MggSel1", mgg, twoPhotonsSel, EqB(
            30, 100., 180.), title="M_{\gamma\gamma}"))
        plots.append(Plot.make1D("MggSel2", mgg, pTmggRatio_sel, EqB(
            30, 100., 180.), title="M_{\gamma\gamma}"))
        plots.append(Plot.make1D("MggSel3", mgg, mgg_sel, EqB(
            30, 100., 180.), title="M_{\gamma\gamma}"))
        plots.append(Plot.make1D("MggSel4", mgg, twoTausSel, EqB(
            30, 100., 180.), title="M_{\gamma\gamma}"))

        cfr = CutFlowReport("yields", recursive=True, printInLog=False)
        plots.append(cfr)
        cfr.add(noSel, "No selection")
        cfr.add(twoPhotonsSel, "Two Photons")
        cfr.add(pTmggRatio_sel, "pT/mgg Ratio")
        cfr.add(mgg_sel, "Inv(M) Sel")
        cfr.add(twoTausSel, "Two Taus")

        return plots
