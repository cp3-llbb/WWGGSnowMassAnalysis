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
                sel_eff[i] = str(f"({sel_eff[i]:.3f}\%)")
            colEntries_withEff = []
            for i, entry in enumerate(colEntries):
                colEntries_withEff.append("{0} {1}".format(
                    entry, sel_eff[i]))
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
                sel_eff[i] = str(f"({sel_eff[i]:.3f}\%)")
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
                sel_eff[i] = str(f"({sel_eff[i]:.3f}\%)")
            colEntries_withEff = []
            for i, entry in enumerate(colEntries):
                colEntries_withEff.append("{0} {1}".format(
                    entry, sel_eff[i]))
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
                sel_eff[i] = str(f"({sel_eff[i]:.3f}\%)")
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
            import os.path
            from bamboo.analysisutils import writePlotIt, runPlotIt
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                        readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
            runPlotIt(cfgName, workdir=workdir, plotIt=self.args.plotIt,
                      eras=(eraMode, eras), verbose=self.args.verbose)

        # mvaSkim
        #import os.path
        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]
        if self.args.mvaSkim and skims:
            from bamboo.analysisutils import loadPlotIt
            p_config, samples, _, systematics, legend = loadPlotIt(
                config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
            # try:
            from bamboo.root import gbl
            import pandas as pd
            import os.path
            # except ImportError as ex:
            #logger.error("Could not import pandas, no dataframes will be saved")
            for skim in skims:
                frames = []
                for smp in samples:
                    # could be a helper in plotit
                    for cb in (smp.files if hasattr(smp, "files") else [smp]):
                        # Take specific columns
                        tree = cb.tFile.Get(skim.treeName)
                        if not tree:
                            print(
                                f"KEY TTree {skim.treeName} does not exist, we are gonna skip this {smp}\n")
                        else:
                            N = tree.GetEntries()
                            cols = gbl.ROOT.RDataFrame(tree).AsNumpy()
                            cols["weight"] *= cb.scale
                            cols["process"] = [smp.name]*len(cols["weight"])
                            frames.append(pd.DataFrame(cols))
                df = pd.concat(frames)
                df["process"] = pd.Categorical(
                    df["process"], categories=pd.unique(df["process"]), ordered=False)
                pqoutname = os.path.join(resultsdir, f"{skim.name}.parquet")
                df.to_parquet(pqoutname)
                logger.info(
                    f"Dataframe for skim {skim.name} saved to {pqoutname}")

        # produce histograms "with datacard conventions"
        if self.args.datacards:
            from bamboo.analysisutils import loadPlotIt
            datacardPlots = [ap for ap in self.plotList if ap.name == "Mgg_c3" or ap.name == "Mgg_c1_Zveto" or ap.name ==
                             "Mgg_c2_Zveto" or ap.name == "Mgg_c4_Zveto" or (self.args.mvaEval and ap.name == "dnn_score")]
            p_config, samples, plots_dc, systematics, legend = loadPlotIt(
                config, datacardPlots, eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir,
                readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
            dcdir = os.path.join(workdir, "datacard_histograms")
            import os
            os.makedirs(dcdir, exist_ok=True)

            def _saveHist(obj, name, tdir=None):
                if tdir:
                    tdir.cd()
                obj.Write(name)
            from functools import partial
            from bamboo.root import gbl

            for era in (self.args.eras[1] or config["eras"].keys()):
                f_dch = gbl.TFile.Open(os.path.join(
                    dcdir, f"histo_for_combine_{era}.root"), "RECREATE")
                saveHist = partial(_saveHist, tdir=f_dch)
                smp = next(smp for smp in samples if smp.cfg.type == "SIGNAL")
                plot = next(
                    plot for plot in plots_dc if plot.name == "Empty_histo")
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

################################
  ## Actual analysis module ##
################################


class CMSPhase2Sim(CMSPhase2SimHistoModule):
    def addArgs(self, parser):
        super().addArgs(parser)
        parser.add_argument("--mvaSkim", action="store_true",
                            help="Produce MVA training skims")
        parser.add_argument("--datacards", action="store_true",
                            help="Produce histograms for datacards")
        parser.add_argument("--mvaEval", action="store_true",
                            help="Import MVA model and evaluate it on the dataframe")

    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport
        from bamboo.plots import EquidistantBinning as EqB
        from bamboo import treefunctions as op

        plots = []

        noSel = noSel.refine("withgenweight", weight=t.genweight)

        # Photons
        photons = op.select(t.gamma, lambda ph: op.AND(
            op.abs(ph.eta) < 3.0, ph.pt > 25))

        # sort photons by pT
        sort_ph = op.sort(photons, lambda ph: -ph.pt)

        # select photons with loose ISO & ID
        isoPhotons = op.select(sort_ph, lambda ph: ph.isopass & (1 << 0))

        IDphotons = op.select(isoPhotons, lambda ph: ph.idpass & (1 << 0))

        # di-Photon mass
        mgg = op.invariant_mass(IDphotons[0].p4, IDphotons[1].p4)

        # di-Photon preselection 1: at least 2 photons with leading photon p_T > 35 and sub-leading photon p_T > 25
        twoPhotonsSel = noSel.refine(
            "twoPhotons", cut=op.AND(op.rng_len(IDphotons) >= 2, IDphotons[0].pt > 35, IDphotons[1].pt > 25))

        # di-Photon preselection 2: pT/mgg > 0.33 for leading photon and 0.25 for sub-leading photon
        pTmggRatio_sel = twoPhotonsSel.refine(
            "ptMggRatio", cut=op.AND(IDphotons[0].pt / mgg > 0.33, IDphotons[1].pt / mgg > 0.25))

        # di-Photon preselection 3: Invarient mass cut
        mgg_sel = pTmggRatio_sel.refine("mgg", cut=op.in_range(100, mgg, 180))

        # Electrons

        electrons = op.sort(op.select(t.elec, lambda el: op.AND(
            op.abs(el.eta) < 3, el.pt > 30)), lambda el: -el.pt)

        cleanedElectrons = op.select(electrons, lambda el: op.NOT(
            op.rng_any(IDphotons, lambda ph: op.deltaR(el.p4, ph.p4) < 0.4)))

        ISOelectrons = op.select(
            cleanedElectrons, lambda el: el.isopass & (1 << 0))

        IDelectrons = op.select(
            ISOelectrons, lambda el: el.idpass & (1 << 0))

        # Jets
        jets = op.sort(op.select(t.jetpuppi, lambda j: op.AND(
            j.pt > 30, op.abs(j.eta) < 5)), lambda j: -j.pt)

        # Muons

        muons = op.sort(op.select(t.muon, lambda mu: op.AND(
            op.abs(mu.eta) < 2.8, mu.pt > 30)), lambda mu: -mu.pt)

        cleanedMuons = op.select(muons, lambda mu: op.AND(op.NOT(
            op.rng_any(IDphotons, lambda ph: op.deltaR(mu.p4, ph.p4) < 0.4)),
            op.NOT(op.rng_any(jets, lambda j: op.deltaR(mu.p4, j.p4) < 0.4))))

        ISOmuons = op.select(cleanedMuons, lambda mu: mu.isopass & (1 << 0))

        IDmuons = op.select(
            ISOmuons, lambda mu: mu.idpass & (1 << 0))

        # taus

        taus = op.sort(op.select(t.tau, lambda tau: op.AND(
            tau.pt > 20, op.abs(tau.eta) < 3)), lambda tau: -tau.pt)

        isolatedTaus = op.select(taus, lambda tau: tau.isopass & (1 << 0))

        cleanedTaus = op.select(isolatedTaus, lambda tau: op.AND(
            op.NOT(op.rng_any(IDphotons,
                   lambda ph: op.deltaR(tau.p4, ph.p4) < 0.2)),
            op.NOT(op.rng_any(cleanedElectrons,
                   lambda el: op.deltaR(tau.p4, el.p4) < 0.2)),
            op.NOT(op.rng_any(cleanedMuons,
                   lambda mu: op.deltaR(tau.p4, mu.p4) < 0.2))
        ))

        # Higgs mass
        mH = 125

        # All tau pairs
        allTauPairs = op.combine(
            cleanedTaus, N=2, pred=lambda t1, t2: t1.charge != t2.charge)

        # Best tau pair with invariant mass closest to Higgs mass
        bestTauPair = op.rng_min_element_by(
            allTauPairs, lambda tt: op.abs(op.invariant_mass(tt[0].p4, tt[1].p4)-mH))

        # Cleaned Jets and B-Jets

        cleanedJets = op.select(jets, lambda j: op.AND(
            op.NOT(op.rng_any(cleanedElectrons,
                   lambda el: op.deltaR(j.p4, el.p4) < 0.4)),
            op.NOT(op.rng_any(cleanedMuons, lambda mu: op.deltaR(j.p4, mu.p4) < 0.4)),
            op.NOT(op.rng_any(cleanedTaus, lambda tau: op.deltaR(j.p4, tau.p4) < 0.4)),
            op.NOT(op.rng_any(IDphotons, lambda ph: op.deltaR(j.p4, ph.p4) < 0.4))
        ))

        IDJets = op.select(cleanedJets, lambda j: j.idpass &
                           (1 << 2))  # tight working point

        bJets = op.select(
            IDJets, lambda j: j.btag & (1 << 1))  # medium working point

        # missing transverse energy
        met = op.select(t.metpuppi)

        ########## Variables for ease of use ##########

        nTaus = op.rng_len(cleanedTaus)
        nElec = op.rng_len(IDelectrons)
        nMuon = op.rng_len(IDmuons)
        nJets = op.rng_len(cleanedJets)

        ########## End of Variables ##########

        sel1_p = noSel.refine("OnePhoton", cut=op.rng_len(isoPhotons) >= 1)

        sel2_p = sel1_p.refine("IDPhoton", cut=op.rng_len(IDphotons) >= 1)

        sel1_e = noSel.refine("OneElec", cut=op.rng_len(ISOelectrons) >= 1)

        sel2_e = sel1_e.refine("IDElec", cut=nElec >= 1)

        sel1_m = noSel.refine("OneMuon", cut=op.rng_len(ISOmuons) >= 1)

        sel2_m = sel1_m.refine("IDMuon", cut=nMuon >= 1)

        plots.append(Plot.make1D("LeadingPhotonISO", op.map(
            isoPhotons, lambda p: p.pt), sel1_p, EqB(30, 0, 300), title="Leading Photon pT"))

        plots.append(Plot.make1D("LeadingPhotonIDISO", op.map(
            IDphotons, lambda p: p.pt), sel2_p, EqB(30, 0, 300), title="Leading Photon pT"))

        plots.append(Plot.make1D("LeadingElectronISO", ISOelectrons[0].pt, sel1_e, EqB(
            30, 0, 300), title="Leading Electron pT"))

        plots.append(Plot.make1D("LeadingElectronIDISO", IDelectrons[0].pt, sel2_e, EqB(
            30, 0, 300), title="Leading Electron pT"))

        plots.append(Plot.make1D("LeadingMuonISO", ISOmuons[0].pt, sel1_m, EqB(
            30, 0, 300), title="Leading Muon pT"))

        plots.append(Plot.make1D("LeadingMuonIDISO", IDmuons[0].pt, sel2_m, EqB(
            30, 0, 300), title="Leading Muon pT"))

        ## Categories ##

        c1 = mgg_sel.refine("hasOneTauOneElec", cut=op.AND(
            nTaus == 1,
            op.rng_len(cleanedElectrons) == 1,
            op.rng_len(cleanedMuons) == 0,
            cleanedTaus[0].charge != cleanedElectrons[0].charge
        ))

        c2 = mgg_sel.refine("hasOneTauOneMuon", cut=op.AND(
            nTaus == 1,
            op.rng_len(cleanedMuons) == 1,
            op.rng_len(cleanedElectrons) == 0,
            cleanedTaus[0].charge != cleanedMuons[0].charge
        ))

        c3 = mgg_sel.refine("hasOneTauNoLept", cut=op.AND(
            nTaus == 1,
            op.rng_len(cleanedElectrons) == 0,
            op.rng_len(cleanedMuons) == 0
        ))

        c4 = mgg_sel.refine("hasTwoTaus", cut=op.AND(
            nTaus >= 2,
            op.rng_len(cleanedElectrons) == 0,
            op.rng_len(cleanedMuons) == 0,
            #op.deltaR(bestTauPair[0].p4, bestTauPair[1].p4) > 0.2
        ))

        ## End of Categories ##

        ########## Z veto ##########

        mTauElec = op.invariant_mass(cleanedTaus[0].p4, cleanedElectrons[0].p4)

        mTauMuon = op.invariant_mass(cleanedTaus[0].p4, cleanedMuons[0].p4)

        mTauTau = op.invariant_mass(cleanedTaus[0].p4, cleanedTaus[1].p4)

        c1_Zveto = c1.refine(
            "hasOneTauOneElec_Zveto", cut=op.NOT(op.in_range(80, mTauElec, 100)))

        c2_Zveto = c2.refine(
            "hasOneTauOneMuon_Zveto", cut=op.NOT(op.in_range(80, mTauMuon, 100)))

        c4_Zveto = c4.refine(
            "hasTwoTaus_Zveto", cut=op.NOT(op.in_range(80, mTauTau, 100)))

        ########## End of Z veto ############

        # plots

        # Leading Photon p_T plots
        # plots.append(Plot.make1D("LeadingPhotonPTtwoPhotonsSel", IDphotons[0].pt, twoPhotonsSel, EqB(
        #     100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("LeadingPhotonPTpTmggRatio_sel", IDphotons[0].pt, pTmggRatio_sel, EqB(
        #     100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("LeadingPhotonPTmgg_sel", IDphotons[0].pt, mgg_sel, EqB(
        #     100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("LeadingPhotonPTc1", IDphotons[0].pt, c1, EqB(
        #     100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("LeadingPhotonPTc2", IDphotons[0].pt, c2, EqB(
        #     100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("LeadingPhotonPTc3", IDphotons[0].pt, c3, EqB(
            100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("LeadingPhotonPTc4", IDphotons[0].pt, c4, EqB(
        #     100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("LeadingPhotonPTc1_Zveto", IDphotons[0].pt, c1_Zveto, EqB(
            100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("LeadingPhotonPTc2_Zveto", IDphotons[0].pt, c2_Zveto, EqB(
            100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("LeadingPhotonPTc4_Zveto", IDphotons[0].pt, c4_Zveto, EqB(
            100, 30, 1000), title="Leading Photon p_{T} [GeV]", plotopts={"log-y": True}))

        # Sub-leading Photon p_T plots
        # plots.append(Plot.make1D("SubleadingPhotonPTtwoPhotonsSel", IDphotons[1].pt, twoPhotonsSel, EqB(
        #     100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("SubleadingPhotonPTpTmggRatio_sel", IDphotons[1].pt, pTmggRatio_sel, EqB(
        #     100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("SubleadingPhotonPTmgg_sel", IDphotons[1].pt, mgg_sel, EqB(
        #     100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("SubleadingPhotonPTc1", IDphotons[1].pt, c1, EqB(
        #     100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("SubleadingPhotonPTc2", IDphotons[1].pt, c2, EqB(
        #     100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubleadingPhotonPTc3", IDphotons[1].pt, c3, EqB(
            100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("SubleadingPhotonPTc4", IDphotons[1].pt, c4, EqB(
        #     100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubleadingPhotonPTc1_Zveto", IDphotons[1].pt, c1_Zveto, EqB(
            100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubleadingPhotonPTc2_Zveto", IDphotons[1].pt, c2_Zveto, EqB(
            100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubleadingPhotonPTc4_Zveto", IDphotons[1].pt, c4_Zveto, EqB(
            100, 25., 1000), title="Subleading Photon p_{T} [GeV]", plotopts={"log-y": True}))

        # Leading Tau p_T plots
        plots.append(Plot.make1D("leadingTauPT_c4", cleanedTaus[0].pt, c4, EqB(
            100, 25, 500), title="Leading Tau p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("leadingTauPT_c4Zveto", cleanedTaus[0].pt, c4_Zveto, EqB(
            100, 25, 500), title="Leading Tau p_{T} [GeV]", plotopts={"log-y": True}))

        # Sub-leading Tau p_T plots
        plots.append(Plot.make1D("SubleadingTauPTc4", cleanedTaus[1].pt, c4, EqB(
            100, 25, 500), title="Sub-leading Tau p_{T} [GeV]", plotopts={"log-y": True}))
        plots.append(Plot.make1D("SubleadingTauPTc4_Zveto", cleanedTaus[1].pt, c4_Zveto, EqB(
            100, 25, 500), title="Subleading Tau p_{T} [GeV]", plotopts={"log-y": True}))

        # di-Photon mass plots
        # plots.append(Plot.make1D("Mgg_twoPhotonsSel", mgg, twoPhotonsSel, EqB(
        #     100, 0, 1000), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_pTmggRatio_sel", mgg, pTmggRatio_sel, EqB(
        #     100, 0, 1000), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_mggsel", mgg, mgg_sel, EqB(
        #     80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c1", mgg, c1, EqB(
        #     80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c2", mgg, c2, EqB(
        #     80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c4", mgg, c4, EqB(
        #     80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c1_Zveto_180", mgg, c1_Zveto, EqB(
        #     80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c2_Zveto_180", mgg, c2_Zveto, EqB(
        #     80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c3_180", mgg, c3, EqB(
            80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c4_Zveto_180", mgg, c4_Zveto, EqB(
            80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))

        # plots.append(Plot.make1D("Mgg_c1_Zveto_140", mgg, c1_Zveto, EqB(
        #     40, 100, 140), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c2_Zveto_140", mgg, c2_Zveto, EqB(
        #     40, 100, 140), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c3_140", mgg, c3, EqB(
            40, 100, 140), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c4_Zveto_140", mgg, c4_Zveto, EqB(
            40, 100, 140), title="M_{\gamma\gamma}", plotopts={"log-y": True}))

        # plots.append(Plot.make1D("Mgg_c1_Zveto_105_145", mgg, c1_Zveto, EqB(
        #     40, 105, 145), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c2_Zveto_105_145", mgg, c2_Zveto, EqB(
        #     40, 105, 145), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c3_105_145", mgg, c3, EqB(
            40, 105, 145), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c4_Zveto_105_145", mgg, c4_Zveto, EqB(
            40, 105, 145), title="M_{\gamma\gamma}", plotopts={"log-y": True}))

        # plots.append(Plot.make1D("Mgg_c1_Zveto_150", mgg, c1_Zveto, EqB(
        #     50, 100, 150), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        # plots.append(Plot.make1D("Mgg_c2_Zveto_150", mgg, c2_Zveto, EqB(
        #     50, 100, 150), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c3_150", mgg, c3, EqB(
            50, 100, 150), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c4_Zveto_150", mgg, c4_Zveto, EqB(
            50, 100, 150), title="M_{\gamma\gamma}", plotopts={"log-y": True}))

        plots.append(Plot.make1D("Mgg_c3_135", mgg, c3, EqB(
            20, 115, 135), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("Mgg_c4_Zveto_135", mgg, c4_Zveto, EqB(
            20, 115, 135), title="M_{\gamma\gamma}", plotopts={"log-y": True}))

        # Cutflow report
        cfr = CutFlowReport("yields", recursive=True, printInLog=False)
        plots.append(cfr)
        cfr.add(noSel, "No selection")
        # cfr.add(c1_Zveto, "One Tau One Electron")
        # cfr.add(c2_Zveto, "One Tau One Muon")
        cfr.add(c4_Zveto, "Two Taus")
        cfr.add(c3, "One Tau No Lept")

        mvaVariables_c4_zveto = {
            "weight": noSel.weight,

            # Event level variables
            "nTaus": op.rng_len(cleanedTaus),
            "nJets": op.rng_len(cleanedJets),
            "nBJets": op.rng_len(bJets),
            "metPt": met[0].pt,

            # Photon and di-Photon variables
            "L_pt_mgg": IDphotons[0].pt / mgg,
            "L_photon_eta": IDphotons[0].eta,
            "L_photon_phi": IDphotons[0].phi,
            "L_photon_ID": IDphotons[0].idpass,
            "SL_pt_mgg": IDphotons[1].pt / mgg,
            "SL_photon_eta": IDphotons[1].eta,
            "SL_photon_phi": IDphotons[1].phi,
            "SL_photon_ID": IDphotons[1].idpass,
            "diP_pt_mgg": (IDphotons[0].pt + IDphotons[1].pt) / mgg,
            "diP_eta": op.sum(IDphotons[0].eta, IDphotons[1].eta) / 2,
            "diP_DR": op.deltaR(IDphotons[0].p4, IDphotons[1].p4),
            "diP_Phi": op.deltaPhi(IDphotons[0].p4, IDphotons[1].p4),

            "LtauPt":  cleanedTaus[0].pt,
            "LtauEta": cleanedTaus[0].eta,
            "SLtauPt": cleanedTaus[1].pt,
            "SLtauEta": cleanedTaus[1].eta,

            "DRtautau": op.deltaR(
                cleanedTaus[0].p4, cleanedTaus[1].p4),
            "DPhitautau": op.deltaPhi(
                cleanedTaus[0].p4, cleanedTaus[1].p4),
            "Mtautau": op.invariant_mass(
                cleanedTaus[0].p4, cleanedTaus[1].p4),
            "pTtautau": op.sum(cleanedTaus[0].pt +
                               cleanedTaus[1].pt),
            "eta_tautau": op.sum(cleanedTaus[0].eta + cleanedTaus[1].eta) / 2,

            "Ljet_Pt": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].pt),
            "Ljet_Eta": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].eta),
            "Ljet_ID": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].idpass),
            "Ljet_btag": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].btag),
            "SLjet_Pt": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].pt),
            "SLjet_Eta": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].eta),
            "SLjet_ID": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].idpass),
            "SLjet_btag": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].btag)
        }

        mvaVariables_c3 = {
            "weight": noSel.weight,

            # Event level variables
            "nJets": op.rng_len(cleanedJets),
            "nBJets": op.rng_len(bJets),
            "metPt": met[0].pt,

            # Photon and di-Photon variables
            "L_pt_mgg": IDphotons[0].pt / mgg,
            "L_photon_eta": IDphotons[0].eta,
            "L_photon_phi": IDphotons[0].phi,
            "L_photon_ID": IDphotons[0].idpass,
            "SL_pt_mgg": IDphotons[1].pt / mgg,
            "SL_photon_eta": IDphotons[1].eta,
            "SL_photon_phi": IDphotons[1].phi,
            "SL_photon_ID": IDphotons[1].idpass,
            "diP_pt_mgg": (IDphotons[0].pt + IDphotons[1].pt) / mgg,
            "diP_eta": op.sum(IDphotons[0].eta, IDphotons[1].eta) / 2,
            "diP_DR": op.deltaR(IDphotons[0].p4, IDphotons[1].p4),
            "diP_Phi": op.deltaPhi(IDphotons[0].p4, IDphotons[1].p4),

            "LtauPt": cleanedTaus[0].pt,
            "LtauEta": cleanedTaus[0].eta,

            "Ljet_Pt": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].pt),
            "Ljet_Eta": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].eta),
            "Ljet_ID": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].idpass),
            "Ljet_btag": op.switch(nJets == 0, op.c_float(0.), cleanedJets[0].btag),
            "SLjet_Pt": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].pt),
            "SLjet_Eta": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].eta),
            "SLjet_ID": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].idpass),
            "SLjet_btag": op.switch(nJets < 2, op.c_float(0.), cleanedJets[1].btag)
        }

        # save mvaVariables to be retrieved later in the postprocessor and save in a parquet file
        if self.args.mvaSkim or self.args.mvaEval:
            from bamboo.plots import Skim
            plots.append(Skim("c4_Zveto", mvaVariables_c4_zveto, c4_Zveto))
            plots.append(Skim("c3", mvaVariables_c3, c3))

        # evaluate dnn model on data
        if self.args.mvaEval:
            #from IPython import embed
            DNNmodel_path = "DNN_HHWWGG/model.onnx"
            mvaVariables_c4_zveto.pop("weight", None)
            mvaVariables_c3.pop("weight", None)
            dnn = op.mvaEvaluator(
                DNNmodel_path, mvaType="ONNXRuntime", otherArgs="predictions")
            inputs = op.array(
                'float', *[op.static_cast('float', val) for val in mvaVariables_c4_zveto.values()])
            output = dnn(inputs)

            plots.append(Plot.make1D(
                "dnn_score", output, c4_Zveto, EqB(50, 0, 1.)))
            hasDNNscore = c4_Zveto.refine("hasDNNscore", cut=output[0] > 0.58)
            cfr.add(hasDNNscore, title='hasDNNscore')
            plots.append(Plot.make1D("Mgg_c4Zveto_hasDNNscore", mgg, hasDNNscore, EqB(
                80, 100., 180.), title="m_{\gamma\gamma}"))
            plots.append(Plot.make1D("DNN_output", op.rng_len(
                output), hasDNNscore, EqB(20, 0, 10), title="dnn_output"))
            # embed()

        return plots
