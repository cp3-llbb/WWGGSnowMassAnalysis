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
                f"${sigSmp.cfg.yields_group}$")  # sigSmp.cfg.yields_group is the name in the legend
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
            smpHdrs.append("\\textbf{Signal}")
            stTotSig, colEntries = colEntriesFromCFREntryHists(report, {eName: Stack(entries=[h for h in (getHist(
                smp, p) for smp in smp_signal) if h]) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            stTotSig, colEntries_forEff = colEntriesFromCFREntryHists_forEff(report, {eName: Stack(entries=[h for h in (getHist(
                smp, p) for smp in smp_signal) if h]) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            colEntries_matrix = np.array(colEntries_forEff)
            sel_eff = np.array([100])
            for i in range(1, len(report.titles)):
                sel_eff = np.append(sel_eff, [float(
                    colEntries_matrix[i]) / float(colEntries_matrix[0]) * 100]).tolist()
            for i in range(len(report.titles)):
                sel_eff[i] = str(f"({sel_eff[i]:.2f}\%)")
            colEntries_withEff = []
            for i, entry in enumerate(colEntries):
                colEntries_withEff.append("{0} {1}".format(
                    entry, sel_eff[i]))
            entries_smp.append(colEntries_withEff)
    if smp_mc:
        sepStr += "|"
        for mcSmp in smp_mc:
            stTotMC, colEntries = colEntriesFromCFREntryHists(report,
                                                              {eName: getHist(mcSmp, p) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            sepStr += f"{align}|"
            if isinstance(mcSmp, plotit.plotit.Group):
                smpHdrs.append(f"${mcSmp.name}$")
            else:
                smpHdrs.append(f"${mcSmp.cfg.yields_group}$")
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
                colEntries_withEff.append("{0} {1} {2}".format(
                    entry, sel_eff[i], MCevents[mcSmp.cfg.pretty_name.rstrip(".root")][0][i]))
            entries_smp.append(colEntries_withEff)
        if len(smp_mc) > 1:
            sepStr += f"|{align}|"
            smpHdrs.append("\\textbf{Background}")
            stTotMC, colEntries = colEntriesFromCFREntryHists(report, {eName: Stack(entries=[h for h in (getHist(
                smp, p) for smp in smp_mc) if h]) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            stTotMC, colEntries_forEff = colEntriesFromCFREntryHists_forEff(report, {eName: Stack(entries=[h for h in (getHist(
                smp, p) for smp in smp_mc) if h]) for eName, p in entryPlots.items()}, precision=yieldPrecision)
            colEntries_matrix = np.array(colEntries_forEff)
            sel_eff = np.array([100])
            for i in range(1, len(report.titles)):
                sel_eff = np.append(sel_eff, [float(
                    colEntries_matrix[i]) / float(colEntries_matrix[0]) * 100]).tolist()
            for i in range(len(report.titles)):
                sel_eff[i] = str(f"({sel_eff[i]:.2f}\%)")
            colEntries_withEff = []
            for i, entry in enumerate(colEntries):
                colEntries_withEff.append("{0} {1}".format(
                    entry, sel_eff[i]))
            entries_smp.append(colEntries_withEff)
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
    # helper: print one bamboo.plots.CutFlowReport.Entry

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
                f"Selection {entry.name}: N={entry.nominal.GetEntries()}, SumW={entry.nominal.GetBinContent(1)}{effMsg}")
            printFun(f"Selection {entry.name}: N={entry.nominal.GetEntries()}")
        if recursive:
            for c in entry.children:
                printEntry(c, printFun=printFun,
                           recursive=recursive, genEvents=genEvents)

    def unwMCevents(entry, smp, mcevents, genEvents=None):
        if entry.nominal is not None:
            mcevents.append(entry.nominal.GetEntries())
        for c in entry.children:
            unwMCevents(c, smp, mcevents, genEvents=genEvents)
        return mcevents

    # retrieve results files, get generated events for each sample
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
        # debug print
        MCevents = {}
        for smp, smpRep in smpReports.items():
            # if smpRep.printInLog:
            logger.info(f"Cutflow report {report.name} for sample {smp}")
            MCevents[smp] = []
            for root in smpRep.rootEntries():
                printEntry(root, genEvents=generated_events[smp])
                mcevents = []
                MCevents[smp].append(unwMCevents(
                    root, smp, mcevents, genEvents=generated_events[smp]))
        # save yields.tex (if needed)
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


        plots = []

        # Photons
        photons = op.select(t.gamma, lambda ph: op.AND(op.abs(ph.eta) < 3, ph.pt > 30))

        ISOphotons = op.select(photons, lambda ph: ph.isopass & (1 << 0)) # loose working point

        IDphotons = op.select(
            ISOphotons, lambda ph: ph.idpass & (1 << 0))

        # di-Photon mass
        mgg = op.invariant_mass(IDphotons[0].p4, IDphotons[1].p4)

        # Photon selection 1: at least 2 photons with leading photon p_T > 35 and sub-leading photon p_T > 25
        twoPhotonsSel = noSel.refine(
            "twoPhotons", cut = op.AND(op.rng_len(IDphotons) >= 2, IDphotons[0].pt > 35, IDphotons[1].pt > 25))

        # Photon selection 2: pT/mgg > 0.33 for leading photon and 0.25 for sub-leading photon
        pTmggRatio_sel = twoPhotonsSel.refine(
            "ptMggRatio", cut = op.AND(IDphotons[0].pt / mgg > 0.33, IDphotons[1].pt / mgg > 0.25))

        # Electrons
        electrons = op.select(t.elec, lambda el: op.AND(op.abs(el.eta) < 3, el.pt > 30.))

        ISOelectrons = op.select(electrons, lambda el: el.isopass & (1 << 0))

        IDelectrons = op.select(
            ISOelectrons, lambda el: el.idpass & (1 << 0))

        cleanedElectrons = op.select(IDelectrons, lambda el: op.NOT(
            op.rng_any(IDphotons, lambda ph: op.deltaR(el.p4, ph.p4) < 0.2)))

        # Muons
        muons = op.select(t.muon, lambda mu: op.AND(
            mu.pt > 30., op.abs(mu.eta) < 2.8))

        ISOmuons = op.select(muons, lambda mu: mu.isopass & (1 << 0))

        IDmuons = op.select(
            ISOmuons, lambda mu: mu.idpass & (1 << 0))

        cleanedMuons = op.select(IDmuons, lambda mu: op.NOT(
            op.rng_any(IDphotons, lambda ph: op.deltaR(mu.p4, ph.p4) < 0.2)))

        # taus

        taus = op.select(t.tau, lambda tau: op.AND(
            tau.pt > 20., op.abs(tau.eta) < 3))

        isolatedTaus = op.select(taus, lambda tau: tau.isopass & (1 << 2)) # tight working point

        cleanedTaus = op.select(isolatedTaus, lambda tau: op.AND(
            op.NOT(op.rng_any(IDphotons,
                   lambda ph: op.deltaR(tau.p4, ph.p4) < 0.2)),
            op.NOT(op.rng_any(cleanedElectrons,
                   lambda el: op.deltaR(tau.p4, el.p4) < 0.2)),
            op.NOT(op.rng_any(cleanedMuons,
                   lambda mu: op.deltaR(tau.p4, mu.p4) < 0.2))
        ))
        
        # diTaus = op.combine(cleanedTaus, N=2, pred=lambda tau1, tau2: tau1.charge != tau2.charge)
        
        # diTauMass = op.invariant_mass(diTaus[0][0].p4, diTaus[0][1].p4)

        # jets

        jets = op.select(t.jetpuppi, lambda jet: op.AND(
            jet.pt > 25., op.abs(jet.eta) < 3))

        cleanedJets = op.select(jets, lambda j: op.AND(
            op.NOT(op.rng_any(cleanedElectrons,
                   lambda el: op.deltaR(j.p4, el.p4) < 0.4)),
            op.NOT(op.rng_any(cleanedMuons, lambda mu: op.deltaR(j.p4, mu.p4) < 0.4)),
            op.NOT(op.rng_any(cleanedTaus, lambda tau: op.deltaR(j.p4, tau.p4) < 0.4)),
            op.NOT(op.rng_any(IDphotons,
                   lambda ph: op.deltaR(j.p4, ph.p4) < 0.4))
        ))

        IDJets = op.select(cleanedJets, lambda j: j.idpass & (1 << 2))

        btaggedJets = op.select(
            IDJets, lambda j: j.btag & (1 << 1))  # medium working point

        hJets = op.sum(cleanedJets[0].p4, cleanedJets[1].p4)



        hGG = op.sum(IDphotons[0].p4, IDphotons[1].p4)
        mJets= op.invariant_mass(IDJets[0].p4, IDJets[1].p4)
        mJets_SL= op.invariant_mass(IDJets[1].p4, IDJets[2].p4)
        hJets = op.sum(IDJets[0].p4, IDJets[1].p4)
       
        #missing transverse energy
        met = op.select(t.metpuppi)      

        #define more variables for ease of use
        nElec = op.rng_len(cleanedElectrons)
        nMuon = op.rng_len(cleanedMuons)
        nJet = op.rng_len(IDJets)
        nPhoton = op.rng_len(IDphotons)
        
        # varibles for DNN
        pT_mGGL = IDphotons[0].pt/ mgg
        pT_mGGSL = IDphotons[1].pt/ mgg
        E_mGGL = IDphotons[0].p4.energy() / mgg
        E_mGGSL = IDphotons[1].p4.energy() / mgg

        # Categories
        
        hasOneTauOneElectron = pTmggRatio_sel.refine(
            "hasOneTauOneElectron", cut = op.AND(op.rng_len(cleanedTaus) == 1, op.rng_len(cleanedElectrons) == 1, cleanedTaus[0].charge != cleanedElectrons[0].charge))
        
        hasOneTauOneMuon = pTmggRatio_sel.refine(
            "hasOneTauOneMuon", cut = op.AND(op.rng_len(cleanedTaus) == 1, op.rng_len(cleanedMuons) == 1, cleanedTaus[0].charge != cleanedMuons[0].charge))
        
        hasTwoTaus = pTmggRatio_sel.refine(
            "hasTwoTaus", cut = op.AND(op.rng_len(cleanedTaus) == 2, cleanedTaus[0].charge != cleanedTaus[1].charge)) # chose the pair having InvM close to 125 GeV
        
        hasOneTauNoLept = pTmggRatio_sel.refine(
            "hasOneTauNoLept", cut = op.AND(op.rng_len(cleanedTaus) == 1, op.rng_len(cleanedElectrons) == 0, op.rng_len(cleanedMuons) == 0))
        
        mTauElec = op.invariant_mass(cleanedTaus[0].p4, cleanedElectrons[0].p4)
        
        mTauMuon = op.invariant_mass(cleanedTaus[0].p4, cleanedMuons[0].p4)
        
        mTauTau = op.invariant_mass(cleanedTaus[0].p4, cleanedTaus[1].p4)
        
        hasOneTauOneElectron_Zveto = hasOneTauOneElectron.refine("hasOneTauOneElectron_Zveto", cut = op.NOT(op.in_range(80, mTauElec, 100)))

        hasOneTauOneMuon_Zveto = hasOneTauOneMuon.refine("hasOneTauOneMuon_Zveto", cut = op.NOT(op.in_range(80, mTauMuon, 100)))
        
        hasTwoTaus_Zveto = hasTwoTaus.refine("hasTwoTaus_Zveto", cut = op.NOT(op.in_range(80, mTauTau, 100)))
        

       # plots

        plots.append(Plot.make1D("LeadingPhotonPTSel1", IDphotons[0].pt, twoPhotonsSel, EqB(
            30, 30., 1000.), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("LeadingPhotonPTSel2", IDphotons[0].pt, pTmggRatio_sel, EqB(
            30, 30., 1000.), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))

        plots.append(Plot.make1D("SubLeadingPhotonPTSel1", IDphotons[1].pt, twoPhotonsSel, EqB(
            20, 20., 700.), title="Sub-Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubLeadingPhotonPTSel2", IDphotons[1].pt, pTmggRatio_sel, EqB(
            20, 20., 700.), title="Sub-Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))

        plots.append(Plot.make1D("leadingTau_ptSel3", cleanedTaus[0].pt, hasTwoTaus, EqB(
            20, 0., 500.), title="Leading Tau p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("leadingTau_ptSel4", cleanedTaus[0].pt, hasTwoTaus_Zveto, EqB(
            20, 0., 500.), title="Leading Tau p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubleadingTau_ptSel3", cleanedTaus[1].pt, hasTwoTaus, EqB(
            20, 0., 500.), title="Subleading Tau p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubleadingTau_ptSel4", cleanedTaus[1].pt, hasTwoTaus_Zveto, EqB(
            20, 0., 500.), title="Subleading Tau p_{T} [GeV]", plotopts={"log-y": True}))

        plots.append(Plot.make1D("MggSel1", mgg, twoPhotonsSel, EqB(
            30, 0., 1000.), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("MggSel2", mgg, pTmggRatio_sel, EqB(
            30, 0., 1000.), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        
        plots.append(Plot.make1D("MttSel3", mTauTau, hasTwoTaus, EqB(
            30, 0, 200.), title="M_{\tau\tau}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("MttSel4", mTauTau, hasTwoTaus_Zveto, EqB(
            30, 0, 200.), title="M_{\tau\tau}", plotopts={"log-y": True}))

        cfr = CutFlowReport("yields", recursive=True, printInLog=False)
        plots.append(cfr)
        
        cfr.add(noSel, "No selection")
        cfr.add(hasOneTauNoLept, "One Tau No Lept")
        cfr.add(hasOneTauOneMuon_Zveto, "One Tau One Muon")
        cfr.add(hasTwoTaus_Zveto, "Two Taus")
        cfr.add(hasOneTauOneElectron_Zveto, "One Tau One Electron")

        return plots
