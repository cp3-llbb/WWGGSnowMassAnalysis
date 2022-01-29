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
                              r"""\documentclass{report}
\usepackage{graphicx}     
\usepackage[a4paper,bindingoffset=0.5cm,left=0cm,right=1cm,top=2cm,bottom=2cm,footskip=0.25cm]{geometry}                         
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
    if smp_signal and smp_mc:
        sepStr += f"|{align}|"
        smpHdrs.append("$S/\sqrt{B}$")
        colEntries = []
        import numpy.ma as ma
        for stSig, stMC in zip(stTotSig, stTotMC):
            if stSig is not None and stMC is not None:
                dtCont = stSig.contents
                mcCont = ma.array(stMC.contents)
                ratio = dtCont/np.sqrt(mcCont)
                ratioErr = np.sqrt(mcCont**2*stSig.sumw2 +
                                   dtCont**2*(stMC.sumw2+stMC.syst2))/mcCont**2
                if mcCont[1] != 0.:
                    colEntries.append("${0:.5f}$".format(ratio[1]))
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
                f"\\resizebox{{\\textwidth}}{{!}}{{",
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
                "}"
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


class CMSPhase2SimRTBHistoModule(CMSPhase2SimRTBModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimRTBHistoModule, self).__init__(args)
    
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(CMSPhase2SimRTBHistoModule, self).postProcess(taskList, config=config, workdir=workdir, resultsdir=resultsdir)
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
            import os.path
            cfgName = os.path.join(workdir, "plots.yml")
            writePlotIt(config, plotList_plotIt, cfgName, eras=eras, workdir=workdir, resultsdir=resultsdir,
                        readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes, plotDefaults=self.plotDefaults)
            runPlotIt(cfgName, workdir=workdir, plotIt=self.args.plotIt,
                      eras=(eraMode, eras), verbose=self.args.verbose)
   
        #mvaSkim 
        #import os.path 
        from bamboo.plots import Skim
        skims = [ap for ap in self.plotList if isinstance(ap, Skim)]
        if self.args.mvaSkim and skims:
            from bamboo.analysisutils import loadPlotIt
            p_config, samples, _, systematics, legend = loadPlotIt(config, [], eras=self.args.eras[1], workdir=workdir, resultsdir=resultsdir, readCounters=self.readCounters, vetoFileAttributes=self.__class__.CustomSampleAttributes)
            #try:
            from bamboo.root import gbl
            import pandas as pd
            import os.path
            #except ImportError as ex:
                #logger.error("Could not import pandas, no dataframes will be saved")
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
                            cols = gbl.ROOT.RDataFrame(tree).AsNumpy()
                            cols["weight"] *= cb.scale
                            cols["process"] = [smp.name]*len(cols["weight"])
                            frames.append(pd.DataFrame(cols))
                df = pd.concat(frames)
                df["process"] = pd.Categorical(df["process"], categories=pd.unique(df["process"]), ordered=False)
                pqoutname = os.path.join(resultsdir, f"{skim.name}.parquet")
                df.to_parquet(pqoutname)
                logger.info(f"Dataframe for skim {skim.name} saved to {pqoutname}")    
        
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

        #selection of photons with eta in the detector acceptance
        photons = op.select(t.gamma, lambda ph : op.AND(op.abs(ph.eta)<2.5, ph.pt >25.)) 
        #sort photons by pT 
        sort_ph = op.sort(photons, lambda ph : -ph.pt)
        #selection of photons with loose ID        
        isoPhotons = op.select(sort_ph, lambda ph : ph.isopass & (1<<0)) #switched to tight ID on 26/11
        idPhotons = op.select(isoPhotons, lambda ph : ph.idpass & (1<<0))

        electrons = op.select(t.elec, lambda el : op.AND(el.pt > 10., op.abs(el.eta) < 2.5))
        
        jets = op.select(t.jetpuppi, lambda jet : op.AND(jet.pt > 30., op.abs(jet.eta) < 5))
        clElectrons = op.select(electrons, lambda el : op.AND(
            op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(el.p4, ph.p4) < 0.4 )),
            #op.NOT(op.rng_any(jets, lambda j : op.deltaR(el.p4, j.p4) < 0.4 ))
            ))
        sort_el = op.sort(clElectrons, lambda el : -el.pt)        
        isoElectrons = op.select(sort_el, lambda el : el.isopass & (1<<0))
        idElectrons = op.select(isoElectrons, lambda el : el.idpass & (1<<0))     
        #slElectrons = op.select(idElectrons, lambda el : op.NOT(op.in_range(86.187, op.rng_any(idPhotons,lambda ph:op.invariant_mass(el.p4, ph.p4)), 90.187000))) #apply the removal of rmZee peak   
        
        muons = op.select(t.muon, lambda mu : op.AND(mu.pt > 10., op.abs(mu.eta) < 2.5))
        clMuons = op.select(muons, lambda mu : op.AND(op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(mu.p4, ph.p4) < 0.4 )),
                                                      op.NOT(op.rng_any(jets, lambda j : op.deltaR(mu.p4, j.p4) < 0.4 ))))
        sort_mu = op.sort(clMuons, lambda mu : -mu.pt)
        idMuons = op.select(sort_mu, lambda mu : mu.idpass & (1<<2)) #apply tight ID  
        isoMuons = op.select(idMuons, lambda mu : mu.isopass & (1<<2)) #apply tight isolation 
        
        taus = op.sort(op.select(t.tau, lambda tau: op.AND(tau.pt > 20., op.abs(tau.eta) < 3)), lambda tau: -tau.pt)
        cleanedTaus = op.select(taus, lambda tau: op.AND(op.NOT(op.rng_any(idPhotons, lambda ph: op.deltaR(tau.p4, ph.p4) < 0.2)),
                                                         op.NOT(op.rng_any(idElectrons, lambda el: op.deltaR(tau.p4, el.p4) < 0.2)),
                                                         op.NOT(op.rng_any(isoMuons,lambda mu: op.deltaR(tau.p4, mu.p4) < 0.2))
                                                     ))
        isolatedTaus = op.select(cleanedTaus, lambda tau: tau.isopass & (1 << 2)) # tight working point Oguz is using loose ISO

        # All tau pairs
        allTauPairs = op.combine(isolatedTaus, N=2, pred=lambda t1, t2: t1.charge != t2.charge)

        # Best tau pair with invariant mass closest to Higgs mass
        bestTauPair = op.rng_min_element_by(allTauPairs, lambda tt: op.abs(op.invariant_mass(tt[0].p4, tt[1].p4)-125))

        clJets = op.select(jets, lambda j : op.AND(
            op.NOT(op.rng_any(idPhotons, lambda ph : op.deltaR(ph.p4, j.p4) < 0.4) ),
            op.NOT(op.rng_any(idElectrons, lambda el : op.deltaR(el.p4, j.p4) < 0.4) ),  
            op.NOT(op.rng_any(isoMuons, lambda mu : op.deltaR(mu.p4, j.p4) < 0.4) ),
            op.NOT(op.rng_any(isolatedTaus, lambda tau: op.deltaR(j.p4, tau.p4) < 0.4))
        ))
        sort_jets = op.sort(clJets, lambda jet : -jet.pt)  
        idJets = op.select(sort_jets, lambda j : j.idpass & (1<<2))
        bJets = op.select(idJets, lambda j: j.btag & (1 << 1))  

        mGG     = op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4)
        mTauTau = op.invariant_mass(isolatedTaus[0].p4, isolatedTaus[1].p4)
        pTGG    = op.sum(idPhotons[0].pt, idPhotons[1].pt)
        mJets   = op.invariant_mass(idJets[0].p4, idJets[1].p4)
        mJets_SL= op.invariant_mass(idJets[1].p4, idJets[2].p4)
       
        #Fully leptonic FL invmasses
        mE   = op.invariant_mass(idElectrons[0].p4, idElectrons[1].p4)
        mMu  = op.invariant_mass(idMuons[0].p4, idMuons[1].p4)
        mEMu = op.invariant_mass(idElectrons[0].p4, idMuons[0].p4)
        #missing transverse energy
        met  = op.select(t.metpuppi)  
        metPt= met[0].pt

        #define more variables for ease of use
        nElec = op.rng_len(idElectrons)
        nMuon = op.rng_len(isoMuons)
        nJet = op.rng_len(idJets)
        nPhoton = op.rng_len(idPhotons)
        nTau = op.rng_len(isolatedTaus) 
        
        #defining more DNN variables
        pT_mGGL = op.product(idPhotons[0].pt, op.pow(mGG, -1)) 
        pT_mGGSL = op.product(idPhotons[1].pt, op.pow(mGG, -1)) 
        E_mGGL = op.product(idPhotons[0].p4.energy(), op.pow(mGG, -1))
        E_mGGSL = op.product(idPhotons[1].p4.energy(), op.pow(mGG, -1))

        #selections for efficiency check

        sel1_p = noSel.refine("2Photon", cut = op.AND((op.rng_len(sort_ph) >= 2), (sort_ph[0].pt > 35.)))
        sel2_p = sel1_p.refine("idPhoton", cut = op.AND((op.rng_len(idPhotons) >= 2), (idPhotons[0].pt > 35.)))
        #selections for the event inv mass of photons within the 100-180 window
        hasInvM = sel2_p.refine("hasInvM", cut= op.AND(
            (op.in_range(100, op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4), 180)) 
        ))

        sel2_p_80 = sel1_p.refine("LP80", cut = op.AND((op.rng_len(idPhotons) >= 2), (idPhotons[0].pt > 80.)))
        hasInvM80 = sel2_p_80.refine("hasInvM80", cut= op.AND(
            (op.in_range(115, op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4), 135)) 
        ))

        sel2_p_100 = sel1_p.refine("LP100", cut = op.AND((op.rng_len(idPhotons) >= 2), (idPhotons[0].pt > 100.)))
        hasInvM100 = sel2_p_100.refine("hasInvM100", cut= op.AND(
            (op.in_range(115, op.invariant_mass(idPhotons[0].p4, idPhotons[1].p4), 135)) 
        ))


        sel1_e = noSel.refine("OneE", cut = op.rng_len(sort_el) >= 1)
        sel2_e = sel1_e.refine("idElectron", cut = op.rng_len(idElectrons) >= 1)
        sel3_e = sel2_e.refine("slElectron", cut = op.AND(op.rng_len(idElectrons) >= 1))

        sel1_m = noSel.refine("OneM", cut = op.rng_len(sort_mu) >= 1)
        sel2_m = sel1_m.refine("idMuon", cut = op.rng_len(idMuons) >= 1)
        sel3_m = sel2_m.refine("isoMuon", cut = op.AND(op.rng_len(isoMuons) >= 1))

       
        ## Categories ##
        #selections for semileptonic final state
        hasOneL = hasInvM.refine("hasOneL", cut = op.OR(op.AND(nElec == 1, nMuon == 0), op.AND(nElec == 0, nMuon == 1)))
        yields.add(hasOneL, title='hasOneL')
        hasOneL80 = hasInvM80.refine("hasOneL80", cut = op.AND(op.OR(op.AND(nElec == 1, nMuon == 0), op.AND(nElec == 0, nMuon == 1)), met[0].pt > 80))
        yields.add(hasOneL80, title='hasOneL80')
        hasOneL100 = hasInvM100.refine("hasOneL100", cut = op.AND(op.OR(op.AND(nElec == 1, nMuon == 0), op.AND(nElec == 0, nMuon == 1)), met[0].pt > 100))
        yields.add(hasOneL100, title='hasOneL100')
        hasOneEl = hasInvM.refine("hasOneEl", cut = op.AND(nElec == 1, nMuon == 0))
        hasOneMu = hasInvM.refine("hasOneMu", cut = op.AND(nElec == 0, nMuon == 1))
        #adding jets on the semileptonic final state
        hasOneJ = hasOneL.refine("hasOneJ", cut = nJet >= 1)
        hasTwoJ = hasOneJ.refine("hasTwoJ", cut = nJet >= 2)
        hasThreeJ = hasTwoJ.refine("hasThreeJ", cut = nJet >= 3)
        
        hasTwoL = hasInvM.refine('hasTwoL', cut = op.AND(
            op.OR(
            op.AND(op.AND(nElec >= 2, nMuon == 0), idElectrons[0].charge != idElectrons[1].charge, op.NOT(op.deltaR(idElectrons[0].p4, idElectrons[1].p4) < 0.4), op.OR(mE < 80, mE >100)),
            op.AND(op.AND(nElec >= 1, nMuon == 1), idElectrons[0].charge != idMuons[0].charge, op.NOT(op.deltaR(idElectrons[0].p4, idMuons[0].p4) < 0.4), op.OR(mEMu < 80, mEMu >100)),
            op.AND(op.AND(nElec == 1, nMuon >= 1), idElectrons[0].charge != idMuons[0].charge, op.NOT(op.deltaR(idElectrons[0].p4, idMuons[0].p4) < 0.4), op.OR(mEMu < 80, mEMu >100)),
            op.AND(op.AND(nMuon >= 2, nElec == 0), idMuons[0].charge != idMuons[1].charge, op.NOT(op.deltaR(idMuons[0].p4, idMuons[1].p4) < 0.4), op.OR(mMu < 80, mMu >100))),
            pTGG > 91,
            #op.AND(idElectrons[2].pt > 10, idMuons[2].pt > 10),
            #bJets.pt < 20,
            met[0].pt > 20   
            ))

        yields.add(hasTwoL, title='hasTwoL')
        
        #hasZeroL = hasInvM.refine('hasZeroL', cut = op.AND(nJet >= 4, nElec == 0, nMuon == 0, nTau == 0))
        #yields.add(hasZeroL, title='hasZeroL')

        c3 = hasInvM.refine("hasOneTauNoLept", cut=op.AND( nTau == 1, op.rng_len(idElectrons) == 0, op.rng_len(isoMuons) == 0 ))
        yields.add(c3, "One Tau No Lept")

        c4 = hasInvM.refine("hasTwoTaus", cut=op.AND(nTau >= 2, op.rng_len(idElectrons) == 0, op.rng_len(isoMuons) == 0 ))
        ########## Z veto ##########
        c4_Zveto = c4.refine( "hasTwoTaus_Zveto", cut=op.NOT(op.in_range(80, mTauTau, 100)))
        yields.add(c4_Zveto, "Two Taus")

        ## End of Categories ##
        #plots       
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
        plots.append(Plot.make1D("LeadingPhotonpT_mGGLhasOneL", pT_mGGL, hasOneL,EqB(100, 0., 5.) ,title = "Leading Photon p_{T}/m_{\gamma\gamma}"))  
        plots.append(Plot.make1D("SubLeadingPhotonpT_mGGLhasOneL", pT_mGGSL, hasOneL,EqB(100, 0., 5.) ,title = "SubLeading Photon p_{T}/m_{\gamma\gamma}"))
        plots.append(Plot.make1D("LeadingPhotonE_mGGLhasOneL", E_mGGL, hasOneL,EqB(100, 0., 5.) ,title = "Leading Photon E/m_{\gamma\gamma}"))
        plots.append(Plot.make1D("SubLeadingPhotonE_mGGLhasOneL", E_mGGSL, hasOneL,EqB(100, 0., 5.) ,title = "SubLeading Photon E/m_{\gamma\gamma}")) 
        plots.append(Plot.make1D("MET", metPt, hasOneL,EqB(80, 0., 800.) ,title="MET"))
        plots.append(Plot.make1D("Inv_mass_gghasOneL",mGG , hasOneL, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
        plots.append(Plot.make1D("Inv_mass_gghasOneL_135",mGG , hasOneL, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))
        plots.append(Plot.make1D("Inv_mass_gghasOneL80",mGG , hasOneL80, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))
        plots.append(Plot.make1D("Inv_mass_gghasOneL100",mGG , hasOneL100, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))
 
        #Leading electron Plots
        ElectronpT = Plot.make1D("ElectronpT", idElectrons[0].pt, hasOneEl, EqB(30, 0., 300.), title = 'Leading Electron pT')
        ElectronE = Plot.make1D("ElectronE", idElectrons[0].p4.E(), hasOneEl, EqB(50, 0., 500.), title = 'Leading Electron E')
        ElectronEta = Plot.make1D("ElectronEta", idElectrons[0].eta, hasOneEl, EqB(80, -4., 4.), title = 'Leading Electron eta')
        ElectronPhi = Plot.make1D("ElectronPhi", idElectrons[0].phi, hasOneEl, EqB(100, -3.5, 3.5), title = 'Leading Electron phi')

        #Leading muon Plots
        MuonpT = Plot.make1D("MuonpT", idMuons[0].pt, hasOneMu, EqB(30, 0., 100.), title = 'Leading Muon pT')
        MuonE = Plot.make1D("MuonE", idMuons[0].p4.E(), hasOneMu, EqB(50, 0., 500.), title = 'Leading Muon E')
        MuonEta = Plot.make1D("MuonEta", idMuons[0].eta, hasOneMu, EqB(80, -4., 4.), title = 'Leading Muon eta')
        MuonPhi = Plot.make1D("MuonPhi", idMuons[0].phi, hasOneMu, EqB(100, -3.5, 3.5), title = 'Leading Muon phi')
        
        #Lepton Plots
        LeptonpT = SummedPlot('LeptonpT', [ElectronpT,MuonpT], xTitle = 'Leading Lepton pT')
        plots.append(ElectronpT)
        plots.append(MuonpT)
        plots.append(LeptonpT)
        LeptonE = SummedPlot('LeptonE', [ElectronE,MuonE], xTitle = 'Leading Lepton E')
        plots.append(ElectronE)
        plots.append(MuonE)
        plots.append(LeptonE)
        LeptonEta = SummedPlot('LeptonEta', [ElectronEta,MuonEta], xTitle = 'Leading Lepton Eta')
        plots.append(ElectronEta)
        plots.append(MuonEta)
        plots.append(LeptonEta)
        LeptonPhi = SummedPlot('LeptonPhi', [ElectronPhi,MuonPhi], xTitle = 'Leading Lepton Phi')
        plots.append(ElectronPhi)
        plots.append(MuonPhi)
        plots.append(LeptonPhi)

        #hasTwoL
        plots.append(Plot.make1D("Inv_mass_gghasTwoL",mGG , hasTwoL, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
        plots.append(Plot.make1D("Inv_mass_gghasTwoL_135",mGG , hasTwoL, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))

        #hasZeroL
        #plots.append(Plot.make1D("Inv_mass_gghasZeroL",mGG , hasZeroL, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
        
        # tau category plots
        plots.append(Plot.make1D("mGG_c3", mGG, c3, EqB(80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("mGG_c3_135", mGG, c3, EqB(20, 115, 135), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("mGG_c4_Zveto", mGG, c4_Zveto, EqB(80, 100, 180), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        plots.append(Plot.make1D("mGG_c4_Zveto_135", mGG, c4_Zveto, EqB(20, 115, 135), title="M_{\gamma\gamma}", plotopts={"log-y": True}))
        
        #hasOneJ
        plots.append(Plot.make1D("LeadingJetPtOneJ", idJets[0].pt, hasOneJ, EqB(30, 0., 300.), title = 'Leading Jet pT'))
        plots.append(Plot.make1D("LeadingJetEtaOneJ", idJets[0].eta, hasOneJ, EqB(80, -4., 4.), title="Leading Jet eta"))
        plots.append(Plot.make1D("LeadingJetPhiOneJ", idJets[0].phi, hasOneJ, EqB(100, -3.5, 3.5), title="Leading Jet phi"))
        plots.append(Plot.make1D("LeadingJetEOnej", idJets[0].p4.energy(), hasOneJ, EqB(50, 0.,500.), title = 'Leading Jet E'))
        
        #hasTwoJ
        plots.append(Plot.make1D("SubLeadingJetPtTwoJ", idJets[1].pt, hasTwoJ, EqB(30, 0., 300.), title = 'SubLeading Jet pT'))
        plots.append(Plot.make1D("SubLeadingJetEtaTwoJ", idJets[1].eta, hasTwoJ, EqB(80, -4., 4.), title="SubLeading Jet eta"))
        plots.append(Plot.make1D("SubLeadingJetPhiTwoJ", idJets[1].phi, hasTwoJ, EqB(100, -3.5, 3.5), title="SubLeading Jet phi"))
        plots.append(Plot.make1D("SubLeadingJetETwoJ", idJets[1].p4.energy(), hasTwoJ, EqB(50, 0.,500.), title = 'SubLeading Jet E'))
        plots.append(Plot.make1D("Inv_mass_jjTwoJ", mJets, hasTwoJ, EqB(80, 20.,220.), title = "m_{jets}"))

        #hasThreeJ
        plots.append(Plot.make1D("Inv_mass_jjThreeJ", mJets_SL, hasThreeJ, EqB(80, 100.,180.), title = "m_{jets}"))

        mvaVariables = {
            "weight": noSel.weight,
            "Eta_ph1": idPhotons[0].eta,
            "Phi_ph1": idPhotons[0].phi,
            "E_mGG_ph1": E_mGGL,
            "pT_mGG_ph1": pT_mGGL,
            "Eta_ph2": idPhotons[1].eta,
            "Phi_ph2": idPhotons[1].phi,
            "E_mGG_ph2": E_mGGSL,
            "pT_mGG_ph2": pT_mGGSL,
            "Electron_E": op.switch(op.rng_len(idElectrons)==0,op.c_float(0.),idElectrons[0].p4.E()), 
            "Electron_pT": op.switch(op.rng_len(idElectrons)==0,op.c_float(0.),idElectrons[0].pt),
            "Electron_Eta": op.switch(op.rng_len(idElectrons)==0,op.c_float(0.),idElectrons[0].eta),
            "Electron_Phi": op.switch(op.rng_len(idElectrons)==0,op.c_float(0.),idElectrons[0].phi),
            "Muon_E": op.switch(op.rng_len(idMuons)==0,op.c_float(0.),idMuons[0].p4.E()), 
            "Muon_pT": op.switch(op.rng_len(idMuons)==0,op.c_float(0.),idMuons[0].pt),
            "Muon_Eta": op.switch(op.rng_len(idMuons)==0,op.c_float(0.),idMuons[0].eta),
            "Muon_Phi": op.switch(op.rng_len(idMuons)==0,op.c_float(0.),idMuons[0].phi),
            "nJets": nJet,
            "E_jet1": op.switch(op.rng_len(idJets)==0,op.c_float(0.),idJets[0].p4.E()),   
            "pT_jet1": op.switch(op.rng_len(idJets)==0,op.c_float(0.),idJets[0].pt),
            "Eta_jet1": op.switch(op.rng_len(idJets)==0,op.c_float(0.),idJets[0].eta),
            "Phi_jet1": op.switch(op.rng_len(idJets)==0,op.c_float(0.),idJets[0].phi), 
            "E_jet2": op.switch(op.rng_len(idJets)<2,op.c_float(0.),idJets[1].p4.E()),   
            "pT_jet2": op.switch(op.rng_len(idJets)<2,op.c_float(0.),idJets[1].pt),
            "Eta_jet2": op.switch(op.rng_len(idJets)<2,op.c_float(0.),idJets[1].eta),
            "Phi_jet2": op.switch(op.rng_len(idJets)<2,op.c_float(0.),idJets[1].phi),  
            "InvM_jet": op.switch(op.rng_len(idJets)<2,op.c_float(0.),mJets),
            "InvM_jet2": op.switch(op.rng_len(idJets)<3,op.c_float(0.),mJets_SL),
            "met":metPt
        } 
              
        mvaVariables_c3 = {
            "weight": noSel.weight,
            # Event level variables
            "nJets": nJet,
            "nBJets": op.rng_len(bJets),
            "metPt": metPt,
            # Photon and di-Photon variables
            "L_pt_mGG": pT_mGGL,
            "L_photon_eta": idPhotons[0].eta,
            "L_photon_phi": idPhotons[0].phi,
            "E_mGG_ph1": E_mGGL,
            "E_mGG_ph2": E_mGGSL,
            "SL_pt_mGG": pT_mGGSL,
            "SL_photon_eta": idPhotons[1].eta,
            "SL_photon_phi": idPhotons[1].phi,
            "LTauE": isolatedTaus[0].p4.E(),
            "LtauPt": isolatedTaus[0].pt,
            "LtauEta": isolatedTaus[0].eta,
            "LtauPhi": isolatedTaus[0].phi,
            "Ljet_Pt": op.switch(nJet == 0, op.c_float(0.), idJets[0].pt),
            "Ljet_Eta": op.switch(nJet == 0, op.c_float(0.), idJets[0].eta),
            "SLjet_Pt": op.switch(nJet < 2, op.c_float(0.), idJets[1].pt),
            "SLjet_Eta": op.switch(nJet < 2, op.c_float(0.), idJets[1].eta),
        }
        mvaVariables_c4_zveto = {
            "weight": noSel.weight,
            # Event level variables
            "nTaus": nTau,
            "nJets": nJet,
            "nBJets": op.rng_len(bJets),
            "metPt": metPt,
            # Photon and di-Photon variables
            "L_pt_mGG": pT_mGGL,
            "L_photon_eta": idPhotons[0].eta,
            "L_photon_phi": idPhotons[0].phi,
            "E_mGG_ph1": E_mGGL,
            "E_mGG_ph2": E_mGGSL,
            "SL_pt_mGG": pT_mGGSL,
            "SL_photon_eta": idPhotons[1].eta,
            "SL_photon_phi": idPhotons[1].phi,
            "LtauPt":  isolatedTaus[0].pt,
            "LtauEta": isolatedTaus[0].eta,
            "SLtauPt": isolatedTaus[1].pt,
            "SLtauEta": isolatedTaus[1].eta,
            "LtauPhi": isolatedTaus[0].phi,
            "SLtauPhi": isolatedTaus[1].phi,
            "DRtautau": op.deltaR(isolatedTaus[0].p4, isolatedTaus[1].p4),
            "DPhitautau": op.deltaPhi(isolatedTaus[0].p4, isolatedTaus[1].p4),
            "Mtautau": op.invariant_mass(isolatedTaus[0].p4, isolatedTaus[1].p4),
            "LTauE": isolatedTaus[0].p4.E(),
            "SLTauE": isolatedTaus[1].p4.E(),
            "Ljet_Pt": op.switch(nJet == 0, op.c_float(0.), idJets[0].pt),
            "Ljet_Eta": op.switch(nJet == 0, op.c_float(0.), idJets[0].eta),
            "SLjet_Pt": op.switch(nJet < 2, op.c_float(0.), idJets[1].pt),
            "SLjet_Eta": op.switch(nJet < 2, op.c_float(0.), idJets[1].eta),
        }


        # save mvaVariables to be retrieved later in the postprocessor and save in a parquet file
        if self.args.mvaSkim or self.args.mvaEval:
            from bamboo.plots import Skim
            plots.append(Skim("Skim", mvaVariables,hasOneL))
            plots.append(Skim("c3", mvaVariables_c3, c3))
            plots.append(Skim("c4_Zveto", mvaVariables_c4_zveto, c4_Zveto))
            #plots.append(Skim("Skim_FH", mvaVariables_FH, hasZeroL))


        # evaluate dnn model on data
        if self.args.mvaEval:
            #from IPython import embed
            tt_DNNmodel_path_even  = "/home/ucl/cp3/sdonerta/DNN_Tautau/even_model.onnx" 
            tt_DNNmodel_path_odd  = "/home/ucl/cp3/sdonerta/DNN_Tautau/odd_model.onnx"     
            WW_DNNmodel_path_even  = "/home/ucl/cp3/sdonerta/DNN_new/even_model_tth.onnx" 
            WW_DNNmodel_path_odd  = "/home/ucl/cp3/sdonerta/DNN_new/odd_model_tth.onnx" 
            mvaVariables.pop("weight", None)
            mvaVariables_c3.pop("weight", None)
            from bamboo.root import loadHeader
            loadHeader("/home/ucl/cp3/sdonerta/bamboodev/WWGG/header_split.h") 

            split_evaluator = op.extMethod('split::Ph1_phi')
            split = split_evaluator(idPhotons[0].phi)

            if split == 0:
                tt_model = tt_DNNmodel_path_even      
                WW_model = WW_DNNmodel_path_even      
            else:
                tt_model = tt_DNNmodel_path_odd
                WW_model = WW_DNNmodel_path_odd

            dnn_ww = op.mvaEvaluator(WW_model, mvaType = "ONNXRuntime", otherArgs = "predictions")
            inputs_ww = op.array('float',*[op.static_cast('float',val) for val in mvaVariables.values()])
            output_ww = dnn_ww(inputs_ww)

            dnn_tt = op.mvaEvaluator(tt_model, mvaType = "ONNXRuntime", otherArgs = "predictions")
            inputs_tt = op.array('float',*[op.static_cast('float',val) for val in mvaVariables_c3.values()])
            output_tt = dnn_tt(inputs_tt)

            #hasDNNscore = hasOneL.refine("hasDNNscore", cut = output[0] < 0.6)
            hasDNNscore = hasOneL.refine("hasDNNscore", cut = op.in_range(0.1, output_ww[0], 0.6))
            yields.add(hasDNNscore, title='hasDNNscore')
            hasDNNscore2 = hasOneL.refine("hasDNNscore2", cut = op.in_range(0.6 ,output_ww[0], 0.8))
            yields.add(hasDNNscore2, title='hasDNNscore2')
            hasDNNscore3 = hasOneL.refine("hasDNNscore3", cut = op.in_range(0.8 ,output_ww[0], 0.92))
            yields.add(hasDNNscore3, title='hasDNNscore3')
            hasDNNscore4 = hasOneL.refine("hasDNNscore4", cut = output_ww[0] > 0.92)
            yields.add(hasDNNscore4, title='hasDNNscore4')
            

            hasDNNscore_tt = c3.refine("hasDNNscore_tt", cut=op.in_range(0.1, output_tt[0], 0.75))
            yields.add(hasDNNscore_tt, title='hasDNNscore_tt')
            hasDNNscore2_tt = c3.refine("hasDNNscore2_tt", cut=output_tt[0] > 0.75)
            yields.add(hasDNNscore2_tt, title='hasDNNscore2_tt')

            plots.append(Plot.make1D("dnn_score_ww", output_ww[0],hasOneL, EqB(50, 0, 1.)))
            plots.append(Plot.make1D("dnn_score_tt", output_tt[0],c3, EqB(50, 0, 1.)))

            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN"  ,mGG, hasDNNscore, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN_2",mGG, hasDNNscore2, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN_3",mGG, hasDNNscore3, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN_4",mGG, hasDNNscore4, EqB(80, 100.,180.), title = "m_{\gamma\gamma}"))
            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN_135"  ,mGG, hasDNNscore, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))
            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN_2_135",mGG, hasDNNscore2, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))
            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN_3_135",mGG, hasDNNscore3, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))
            plots.append(Plot.make1D("Inv_mass_gghasOneL_DNN_4_135",mGG, hasDNNscore4, EqB(20, 115.,135.), title = "m_{\gamma\gamma}"))

            plots.append(Plot.make1D("mGG_c3_hasDNNscore", mGG, hasDNNscore_tt, EqB(80, 100., 180.), title="m_{\gamma\gamma}"))
            plots.append(Plot.make1D("mGG_c3_hasDNNscore2", mGG, hasDNNscore2_tt, EqB(80, 100., 180.), title="m_{\gamma\gamma}"))
            plots.append(Plot.make1D("mGG_c3_hasDNNscore_135", mGG, hasDNNscore_tt, EqB(20, 115., 135.), title="m_{\gamma\gamma}"))
            plots.append(Plot.make1D("mGG_c3_hasDNNscore2_135", mGG, hasDNNscore2_tt, EqB(20, 115., 135.), title="m_{\gamma\gamma}"))

        return plots



    
