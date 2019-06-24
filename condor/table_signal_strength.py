# make a table of the signal strength and associated signifances. 
# Do this for every signal point

# First do it in Data, for 2016, 2017, and Combo
import optparse
import ROOT
from array import array

def makePValuePlot(dataSet):
    # CL observed pvalues
    pvalue_2016 = dataSet["data"]["2016"]["pList"]
    pvalue_2017 = dataSet["data"]["2017"]["pList"]
    pvalue_Combo = dataSet["data"]["Combo"]["pList"]
    xpoints = dataSet["data"]["2016"]["mList"]
    npoints = len(pvalue_2016)
    Xmin = 300
    Xmax = 900
    Ymin = 0.0000005
    Ymax = 1

    c1 = ROOT.TCanvas("c1","PValues",1000,1000)
    c1.Divide(1, 2)    
    c1.SetFillColor(0)
    c1.cd(1)
    ROOT.gPad.SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, ROOT.kWhite, 0, 0)
    ROOT.gPad.SetBottomMargin(0.01)
    ROOT.gPad.SetLeftMargin(0.11)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.06 * (8.0 / 6.5))
    ROOT.gPad.SetLogy()
    ROOT.gPad.SetTicks(1,1)

    h = ROOT.TH1F("dummy","dummy",1, 300, 900)
    h.SetMaximum(Ymax)
    h.SetMinimum(Ymin)
    h.SetTitle("")
    h.SetStats(0)
    h.GetXaxis().SetLimits(Xmin,Xmax)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetTitleOffset(1.0)
    #h.GetYaxis().SetRangeUser(0, 2.2)
    #h.GetYaxis().SetNdivisions(4, 2, 0)
    h.GetXaxis().SetTitle("m_{#tilde{t}} [GeV]")
    h.GetYaxis().SetTitle("Local p-value")
    h.Draw()

    gr_Combo = ROOT.TGraph(npoints, array('d', xpoints), array('d', pvalue_Combo))
    gr_Combo.SetLineColor(ROOT.kBlack)
    gr_Combo.SetLineWidth(3)
    gr_Combo.Draw("same")

    gr_2016 = ROOT.TGraph(npoints, array('d', xpoints), array('d', pvalue_2016))
    gr_2016.SetLineColor(ROOT.kRed+1)
    gr_2016.SetLineStyle(7)
    gr_2016.SetLineWidth(3)

    gr_2017 = ROOT.TGraph(npoints, array('d', xpoints), array('d', pvalue_2017))
    gr_2017.SetLineColor(ROOT.kBlue+1)
    gr_2017.SetLineStyle(3)
    gr_2017.SetLineWidth(3)

    # Draw the 1sigma, 2sigma, and 3sigma lines
    # For 1 sigma: s = 0.68
    #   1 - (0.5 + s/2) = 0.5 - s/2
    entries = []
    numSigma = 4
    for s in range(1, numSigma+1):
        sigma = 0.5 - ROOT.TMath.Erf(s/ROOT.TMath.Sqrt(2))/2
        L = ROOT.TLine(Xmin, sigma, Xmax, sigma)
        L.SetLineColor(2)
        L.SetLineWidth(2)
        L.Draw("same")

        S = ROOT.TPaveText(900,sigma-0.25*sigma,920,sigma+0.5*sigma,"")
        S.SetBorderSize(0)
        S.SetFillStyle(0)
        S.SetTextColor(2)
        S.AddText( str(s)+"#sigma" )
        S.Draw("same")
        entries.append((L,S))

    gr_2016.Draw("L,same")
    gr_2017.Draw("L,same")
    gr_Combo.Draw("L,same")

    legend1 = ROOT.TLegend(0.44397,0.1,0.897487,0.25,"")
    if dataSet["model"]=="RPV":
        legend1.SetHeader("pp #rightarrow #tilde{t}#bar{#tilde{t}}, #tilde{t} #rightarrow t #tilde{#chi}^{0}_{1},  #tilde{#chi}^{0}_{1} #rightarrow jjj");
    elif dataSet["model"]=="SYY":
        legend1.SetHeader("pp #rightarrow #tilde{t}#bar{#tilde{t}}, SYY coupling");
    elif dataSet["model"]=="SHH":
        legend1.SetHeader("pp #rightarrow #tilde{t}#bar{#tilde{t}}, SHH coupling");
    legend1.AddEntry(gr_Combo, "Combined Observed, L_{Int}=77.4 fb^{-1}", "l")
    legend1.AddEntry(gr_2016, "2016 Observed, L_{Int}=35.9 fb^{-1}", "l")
    legend1.AddEntry(gr_2017, "2017 Observed, L_{Int}=41.5 fb^{-1}", "l")
    legend1.SetBorderSize(0)
    legend1.SetFillStyle(0)
    legend1.Draw("same")

    cmstext = ROOT.TPaveText(0.325, 0.92, 0.905, 0.98, "ndc")
    cmstext.SetBorderSize(0)
    cmstext.SetFillStyle(0)
    cmstext.SetTextAlign(12)
    cmstext.AddText("CMS Preliminary, #sqrt{s}=13 TeV")
    cmstext.Draw("same")

    c1.Update()
    c1.cd(2)
    ROOT.gPad.SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, ROOT.kWhite, 0, 0)
    ROOT.gPad.SetLeftMargin(0.11)
    ROOT.gPad.SetRightMargin(0.04)
    ROOT.gPad.SetTopMargin(0.01)
    ROOT.gPad.SetBottomMargin(0.37)
    ROOT.gPad.SetTicks(1,1)

    # Make ratio of Data over total stack
    hr = ROOT.TH1F("dummyr","dummyr",1, 300, 900)
    hr.SetStats(0)
    hr.SetTitle("")
    hr.GetXaxis().SetTitle("m_{#tilde{t}} [GeV]")
    hr.GetYaxis().SetTitle("#sigma_{meas.}/#sigma_{pred.}")
    hr.GetXaxis().SetLimits(Xmin,Xmax)
    hr.GetXaxis().SetLabelSize(0.14)
    hr.GetXaxis().SetTitleSize(0.15)
    hr.GetYaxis().SetLabelSize(0.13)
    hr.GetYaxis().SetTitleSize(0.15)
    hr.GetYaxis().SetTitleOffset(0.3)
    hr.GetYaxis().SetRangeUser(-0.2, 1.2)
    hr.GetYaxis().SetNdivisions(4, 2, 0)
    hr.Draw()

    rvalue_Combo = array('d', dataSet["data"]["Combo"]["rList"])
    rpvalue_Combo = array('d',dataSet["data"]["Combo"]["rpList"])
    rmvalue_Combo = array('d', dataSet["data"]["Combo"]["rmList"])
    zero = array('d', dataSet["data"]["2016"]["zero"])
    rband = ROOT.TGraphAsymmErrors(npoints, array('d', xpoints), rvalue_Combo, zero, zero, rmvalue_Combo, rpvalue_Combo)
    rband.SetFillColor(ROOT.kGreen+1)
    rband.Draw("3 same")
    r = ROOT.TGraph(npoints, array('d', xpoints), rvalue_Combo)
    r.SetLineColor(ROOT.kBlack)
    r.SetLineStyle(ROOT.kDashed)
    r.SetLineWidth(3)
    r.Draw("PL same")
    c1.Update()

    line = ROOT.TF1("line", "1", 300, 900)
    line.SetLineColor(ROOT.kRed)
    line.Draw("same")

    line2 = ROOT.TF1("line", "1", 300, 900)
    line2.SetLineColor(ROOT.kBlack)
    line2.Draw("same")

    c1.Print(dataSet["runtype"]+"_"+dataSet["model"]+".pdf")
    del c1

def makeSigTex(name, l):    
    f = open(name, 'w')

    f.write( "\\documentclass[12pt]{article}\n" )
    f.write( "\n" ) 
    f.write( "\\begin{document}\n" )
    f.write( "\n" )

    for dic in l:
        caption = "Best fit signal strengths for %s model in data" % dic["model"]
        if dic["runtype"] == "pseudoDataS": 
            caption = "Best fit signal strengths for %s model in MC with signal injection" % dic["model"] 
        f.write( "\\begin{table}[p]\n" )
        f.write( "\\centering\n" )
        f.write( "\\caption{%s}\n" % caption )
        f.write( "\\input{%s}\n" % dic["outFileName"] )
        f.write( "\\end{table}\n" )
        f.write( "\n" )

    f.write( "\\end{document}\n" )
    f.close()

def main():
    parser = optparse.OptionParser("usage: %prog [options]\n")
    parser.add_option ('--basedir', dest='basedir',  type='string', default = '.', help="Path to output files")
    options, args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)

    pre_tabular = """\\begin{tabular}{l l l l}
    Mass & Best fit signal strength & Observed Significance & p-value\\\\ \hline
    """    
    path = options.basedir
    runtypes = ["Data", "pseudoDataS"]
    models = ["RPV","SYY"]
    years = ["2016","2017","Combo"]
    masses = ["300","350","400","450","500","550","600","650","700","750","800","850","900"]

    # Loop over all jobs in get the info needed
    l=[]
    d=[]
    for runtype in runtypes:        
        for model in models:
            outFileName = "table_signal_strength_%s_%s_%s" % (model, runtype, path)
            file_table = open(outFileName,'w')
            file_table.write(pre_tabular)
            dataSet={"runtype":runtype,"model":model,"data":{}}
            for year in years:
                data={"mList":[],"rList":[],"rmList":[],"rpList":[],"sList":[],"pList":[],"zero":[]}
                file_table.write("\\multicolumn{4}{c}{%s} \\\\ \\hline \n"%year)
                for mass in masses:
                    print "Year %s, Model %s, Mass %s"%(year, model, mass)
                    filename_r = "%s/Fit_%s_%s/output-files/%s_%s_%s/log_%s%s%s_FitDiag.txt" % (path,runtype, year, model, mass, year, year, model, mass)
                    info_r = ["0","0","0"]
    
                    # Get r from fit jobs
                    if not ((model=="RPV" and year=="Combo" and mass=="0") or (model=="SYY" and year=="Combo" and mass=="0") ):
                        file_r=-1
                        try:
                            file_r = open(filename_r)
                        except:
                            print "File not found:",filename_r 
                            continue
                        line_r = ""
                        for line in file_r:
                            if "Best fit r" in line:        
                                if "Fit failed" in line:
                                    info_r = ["Fit failed"]
                                else:
                                    line_r = line.replace("Best fit r: ","").replace("(68% CL)","").strip().replace("/", " ")
                                    info_r = line_r.split() # best fit r, -error, +error
        
                    # Get sigma and p-value from fit jobs
                    line_sig = ""
                    line_pvalue = ""
                    if not ((model=="RPV" and year=="Combo" and mass=="0") or (model=="SYY" and year=="2017" and mass in ["0"])):
                        filename_sig = "%s/Fit_%s_%s/output-files/%s_%s_%s/log_%s%s%s_Sign_noSig.txt" % (path,runtype, year, model, mass, year, year, model, mass)
                        file_sig=-1
                        try:
                            file_sig = open(filename_sig)
                        except:
                            print "File not found:",filename_sig
                            continue
                        for line in file_sig:
                            if "Significance:" in line:
                                line_sig = line.replace("Significance:", "").strip()
                            elif "p-value" in line:
                                line_pvalue = line.replace("       (p-value =", "").strip()
                                line_pvalue = line_pvalue.replace(")","").strip()
    
                    # Fill lists of data
                    data["mList"].append(abs(float(mass)))
                    data["rList"].append(abs(float(info_r[0])))
                    data["rmList"].append(abs(float(info_r[1])))
                    data["rpList"].append(abs(float(info_r[2])))
                    data["sList"].append(abs(float(line_sig)))
                    data["pList"].append(abs(float(line_pvalue)))
                    data["zero"].append(0.0)
    
                    # Write out r, sigma, and p-value to file
                    if len(info_r) < 3:
                        file_table.write("%s & %s & %s & %s\\\\ \n" % (mass, info_r[0], line_sig, line_pvalue))
                    elif "#" in info_r[0]:
                        file_table.write("%s & Fit failed &  \\\\ \n" % (mass))
                    elif line_sig == "":
                        file_table.write("%s & $%.2f_{%.2f}^{%.2f}$ & %s & %s\\\\ \n" % (mass, float(info_r[0]), float(info_r[1]), float(info_r[2].replace("+-","+")), line_sig, line_pvalue))
                    else:
                        file_table.write("%s & $%.2f_{%.2f}^{%.2f}$ & %.2f & %s\\\\ \n" % (mass, float(info_r[0]), float(info_r[1]), float(info_r[2].replace("+-","+")), float(line_sig), line_pvalue))
                file_table.write("\\hline \n")
                dataSet["data"][year]=data
            file_table.write("\\end{tabular}\n")
            file_table.close()
            l.append({"model":model, "runtype":runtype, "outFileName":outFileName})
            d.append(dataSet)
    
    # Make tex file with all signal strengths
    makeSigTex("table_signal_strength.tex", l)
    
    for dataSet in d:        
        if dataSet["runtype"] == "Data":
            print "------------------------------------------"
            print dataSet["runtype"], dataSet["model"]        
            makePValuePlot(dataSet)

if __name__ == '__main__':
    main()
