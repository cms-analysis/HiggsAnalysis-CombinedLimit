import ROOT
import optparse

def makePostFitPlot(preFitHisto, preFitErrHisto, postFitNP):
    h = preFitHisto.Clone("PostFit")
    for i in range(1, h.GetNbinsX()+1):
        val = preFitHisto.GetBinContent(i)*(preFitErrHisto.GetBinContent(i)**postFitNP[i-1])
        h.SetBinContent(i,val)
    return h

def getMap(inFile, l):
    m = {}
    for bin in l: 
        h = inFile.Get(bin+"_qcdCR")
        hErr = inFile.Get(bin+"_qcdCRErr")
        for i in range(1, h.GetNbinsX()+1):
            r = h.GetBinContent(i)
            rPrime = hErr.GetBinContent(i)
            h.SetBinError(i, r*(rPrime-1))
        m[bin] = (h, hErr)
    return m

def getNPList(entry, bin, year):
    l = []
    if(year=="2016"):
        if  (bin=="D1"): l = [entry.np_tt_qcdCRErrD1Bin1_2016, entry.np_tt_qcdCRErrD1Bin2_2016, entry.np_tt_qcdCRErrD1Bin3_2016, entry.np_tt_qcdCRErrD1Bin4_2016, entry.np_tt_qcdCRErrD1Bin5_2016, entry.np_tt_qcdCRErrD1Bin6_2016]
        elif(bin=="D2"): l = [entry.np_tt_qcdCRErrD2Bin1_2016, entry.np_tt_qcdCRErrD2Bin2_2016, entry.np_tt_qcdCRErrD2Bin3_2016, entry.np_tt_qcdCRErrD2Bin4_2016, entry.np_tt_qcdCRErrD2Bin5_2016, entry.np_tt_qcdCRErrD2Bin6_2016]
        elif(bin=="D3"): l = [entry.np_tt_qcdCRErrD3Bin1_2016, entry.np_tt_qcdCRErrD3Bin2_2016, entry.np_tt_qcdCRErrD3Bin3_2016, entry.np_tt_qcdCRErrD3Bin4_2016, entry.np_tt_qcdCRErrD3Bin5_2016, entry.np_tt_qcdCRErrD3Bin6_2016]
        elif(bin=="D4"): l = [entry.np_tt_qcdCRErrD4Bin1_2016, entry.np_tt_qcdCRErrD4Bin2_2016, entry.np_tt_qcdCRErrD4Bin3_2016, entry.np_tt_qcdCRErrD4Bin4_2016, entry.np_tt_qcdCRErrD4Bin5_2016, entry.np_tt_qcdCRErrD4Bin6_2016]
    elif(year=="2017"):
        if  (bin=="D1"): l = [entry.np_tt_qcdCRErrD1Bin1_2017, entry.np_tt_qcdCRErrD1Bin2_2017, entry.np_tt_qcdCRErrD1Bin3_2017, entry.np_tt_qcdCRErrD1Bin4_2017, entry.np_tt_qcdCRErrD1Bin5_2017, entry.np_tt_qcdCRErrD1Bin6_2017]
        elif(bin=="D2"): l = [entry.np_tt_qcdCRErrD2Bin1_2017, entry.np_tt_qcdCRErrD2Bin2_2017, entry.np_tt_qcdCRErrD2Bin3_2017, entry.np_tt_qcdCRErrD2Bin4_2017, entry.np_tt_qcdCRErrD2Bin5_2017, entry.np_tt_qcdCRErrD2Bin6_2017]
        elif(bin=="D3"): l = [entry.np_tt_qcdCRErrD3Bin1_2017, entry.np_tt_qcdCRErrD3Bin2_2017, entry.np_tt_qcdCRErrD3Bin3_2017, entry.np_tt_qcdCRErrD3Bin4_2017, entry.np_tt_qcdCRErrD3Bin5_2017, entry.np_tt_qcdCRErrD3Bin6_2017]
        elif(bin=="D4"): l = [entry.np_tt_qcdCRErrD4Bin1_2017, entry.np_tt_qcdCRErrD4Bin2_2017, entry.np_tt_qcdCRErrD4Bin3_2017, entry.np_tt_qcdCRErrD4Bin4_2017, entry.np_tt_qcdCRErrD4Bin5_2017, entry.np_tt_qcdCRErrD4Bin6_2017]
    return l

if __name__ == '__main__':
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option ('-t', dest='dataType', type='string', default = 'data', help="Specify if running over data or pseudo data")
    parser.add_option ('-s', dest='sysFile',  type='string', default = 'ttbar_systematics.root', help="Specify if running over data or pseudo data")
    parser.add_option ('-f', dest='fitFile',  type='string', default = 'fitDiagnosticsComboRPV350.root', help="Specify if running over data or pseudo data")
    parser.add_option ('-y', dest='year',     type='string', default = '2016', help="year")
    parser.add_option ('-F', dest='fitType',  type='string', default = 'b',    help="year")
    options, args = parser.parse_args()

    inFile = ROOT.TFile.Open(options.sysFile)
    qcdSysMap = getMap(inFile, ["D1","D2","D3","D4"])
    fitFileCombo = ROOT.TFile.Open(options.fitFile)
    treeIter = "default"
    if options.fitType == "b": treeIter = iter(fitFileCombo.tree_fit_b)
    elif options.fitType == "sb": treeIter = iter(fitFileCombo.tree_fit_sb)    
    entry = treeIter.next()

    c = ROOT.TCanvas( "c", "c", 0, 0, 800, 800)
    c.Divide(2,2)
    i=1
    objects = []
    for bin in ["D1","D2","D3","D4"]:
        c.cd(i)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.15)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetTicks(1,1)
        ROOT.gStyle.SetOptStat(0); 

        leg = ROOT.TLegend(0.15, 0.65, 0.45, 0.88)
        objects.append(leg)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetLineWidth(1)
        leg.SetNColumns(1)
        leg.SetTextFont(42)

        qcdSysHist = qcdSysMap[bin][0]
        objects.append(qcdSysHist)
        qcdSysHist.SetMaximum(1.3)
        qcdSysHist.SetMinimum(0.7)
        qcdSysHist.SetLineColor(ROOT.kBlack)
        qcdSysHist.SetTitle("Pre and Post Fit qcdCR NP: "+bin+" "+options.dataType+" "+options.year+" "+options.fitType)
        qcdSysHist.Draw("E")
        leg.AddEntry(qcdSysHist, "Pre-Fit", "l")

        qcdSysHistPostFit = makePostFitPlot(qcdSysMap[bin][0], qcdSysMap[bin][1], getNPList(entry,bin,options.year))
        objects.append(qcdSysHistPostFit)
        qcdSysHistPostFit.SetLineColor(ROOT.kRed)
        leg.AddEntry(qcdSysHistPostFit, "Post-Fit", "l")
        qcdSysHistPostFit.Draw("hist same")

        line = ROOT.TF1("1" ,"1" ,-2000,20000)
        objects.append(line)
        line.SetLineColor(ROOT.kBlack)
        line.Draw("same")

        leg.Draw()
        i+=1
    c.SaveAs("qcdCR_rPrimePlotCompare_"+options.dataType+"_"+options.year+"_"+options.fitType+".pdf")
