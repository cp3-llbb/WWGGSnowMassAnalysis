import ROOT
def doFitAndCreateWorkspace(fin,channel_name,channel_name2, nameout):

    fout = ROOT.TFile(nameout,"recreate")
    # output as a RooWorkspace
    rws = ROOT.RooWorkspace("HHWWgg","HHWWgg")

    # observable variable
    mass = ROOT.RooRealVar("mgg_%s"%channel_name2,"invariant mass",125,100,180)
    getattr(rws,"import")(mass)

    # background model
    alpha = ROOT.RooRealVar("alpha_%s"%channel_name2,"#alpha",-0.05,-0.2,0.01);
    #alpha.setConstant(1)
    expo = ROOT.RooExponential("exp_%s"%channel_name,"exponential function",mass,alpha);

    mch = fin.Get("Continuum_Bkg")
    dhbkg = ROOT.RooDataHist("mc_%s"%(channel_name),"Continuum_Bkg",ROOT.RooArgList(mass),mch)
    getattr(rws,"import")(dhbkg)

    expo.fitTo(dhbkg)
    # make a nice plot just to check 
    c = ROOT.TCanvas()
    pl = mass.frame()
    dhbkg.plotOn(pl)
    expo.plotOn(pl)
    pl.Draw()
    c.SaveAs("background_fit_%s.pdf"%channel_name2)
    alpha.setConstant(1)
    getattr(rws,"import")(expo)


    # fit gaussians in a certain range
    mass.setRange("r1",123,127)

    # single-Higgs backgrounds & signal we'll assume each is a Gaussian
    #for s in ["GGHH","GGH","VBFH","VH","ttH","tHq"]:
    #    mh    = ROOT.RooRealVar("mh_%s_%s"%(s,channel_name),"mean in GeV",125,123,127)
    #    sigma1 = ROOT.RooRealVar("sigma1_%s_%s"%(s,channel_name),"sigma1 in GeV",2.0,0.5,6)
    #    sigma2 = ROOT.RooRealVar("sigma2_%s_%s"%(s,channel_name),"sigma2 in GeV",4.0,0.5,9)
    #    gauss1 = ROOT.RooGaussian("gauss1_%s_%s"%(s,channel_name),"f1(m|M_{H},#sigma)",mass,mh,sigma1);
    #    gauss2 = ROOT.RooGaussian("gauss2_%s_%s"%(s,channel_name),"f2(m|M_{H},#sigma)",mass,mh,sigma2);
    #    frac   = ROOT.RooRealVar("frac_%s_%s"%(s,channel_name),"sigma1 in GeV",0.6,0.,1)
    #    gauss  = ROOT.RooAddPdf("gauss_%s_%s"%(s,channel_name),"",ROOT.RooArgList(gauss1,gauss2),ROOT.RooArgList(frac))
    #    th1   = fin.Get(s)
    #    dh1   = ROOT.RooDataHist("rdh_%s"%(s),"rdh_%s"%(s),ROOT.RooArgList(mass),th1)
    #    gauss.fitTo(dh1,ROOT.RooFit.Range("fitrange"))
    #
    #    print("Results peak -> %s, channel -> %s"%(s,channel_name))
    #    mh.Print()
    #    sigma1.Print()
    #    sigma2.Print()
    #    print("---------------------------------")
    #    
    #    getattr(rws,"import")(gauss)
    #
    #    # plot to check
    #    c = ROOT.TCanvas()
    #    pla = mass.frame()
    #    dh1.plotOn(pla)
    #    gauss.plotOn(pla)
    #    pla.Draw()
    #    c.SaveAs("signal_mc_fit_%s_%s.pdf"%(s,channel_name))
    #
    #for s in ["GGHH","GGH","VBFH","VH","ttH","tHq"]:
    #   # must set these parameters to constant for the limit setting
    #   rws.var("mh_%s_%s"%(s,channel_name)).setConstant(1)
    #   rws.var("sigma1_%s_%s"%(s,channel_name)).setConstant(1)
    #   rws.var("sigma2_%s_%s"%(s,channel_name)).setConstant(1)
    #   rws.var("frac_%s_%s"%(s,channel_name)).setConstant(1)
        


    #  lets import the signal and single H
    for s in ["GGHH","GGH","VBFH","VH","ttH","tHq"]:
         mh    = fin.Get(s)
         sig   = ROOT.RooDataHist("mh_%s_%s"%(s,channel_name),"mh_%s_%s"%(s,channel_name),ROOT.RooArgList(mass),mh)
         getattr(rws,"import")(sig)

    
    # and lets import the data obs too!
    datah = fin.Get("data_obs")
    dh1   = ROOT.RooDataHist("data_obs_%s"%(channel_name),"obs data",ROOT.RooArgList(mass),datah)
    getattr(rws,"import")(dh1)


    rws.Print()
    fout.WriteTObject(rws)
    fout.Close()
    



file_dict = {
    "Inv_mass_gghasOneL_DNN_1_HL.root": ["", "Cat1", 1000],
    "Inv_mass_gghasOneL_DNN_2_HL.root": ["", "Cat2", 100],
    "Inv_mass_gghasOneL_DNN_3_HL.root": ["", "Cat3", 50],
    "Inv_mass_gghasOneL_DNN_4_HL.root": ["", "Cat4", 10],
    "Inv_mass_gghasTwoL_HL.root": ["", "2L", 100],
    "Mgg_c3_DNN_1_HL.root": ["", "Cat1Tau", 5000],
    "Mgg_c3_DNN_2_HL.root": ["", "Cat2Tau", 100],
    "Mgg_c4_Zveto_HL.root": ["", "2T", 1000],
}

for file in file_dict:
    fIn = ROOT.TFile.Open(file)
    doFitAndCreateWorkspace(fIn,file_dict[file][0],file_dict[file][1], "roofit_"+file)
