# make a table of the signal strength and associated signifances. 
# Do this for every signal point

# First do it in Data, for 2016, 2017, and Combo
import optparse
import ROOT
from array import array

def makePValuePlot(dataSet):
    sohaStyle= ROOT.TStyle("SOHA","soha plot style")
    
    # use plain black on white colors
    sohaStyle.SetFrameBorderMode(0)
    sohaStyle.SetCanvasBorderMode(0)
    sohaStyle.SetPadBorderMode(0)
    sohaStyle.SetPadColor(0)
    sohaStyle.SetCanvasColor(0)
    sohaStyle.SetTitleColor(1)
    sohaStyle.SetStatColor(0)

    # set the paper & margin sizes
    sohaStyle.SetPaperSize(20,26)
    sohaStyle.SetPadTopMargin(0.05)
    sohaStyle.SetPadRightMargin(0.05)
    sohaStyle.SetPadBottomMargin(0.16)
    sohaStyle.SetPadLeftMargin(0.12)

    # use large Times-Roman fonts
    sohaStyle.SetTextFont(132)
    sohaStyle.SetTextSize(0.08)
    sohaStyle.SetLabelFont(132,"x")
    sohaStyle.SetLabelFont(132,"y")
    sohaStyle.SetLabelFont(132,"z")
    sohaStyle.SetLabelSize(0.05,"x")
    sohaStyle.SetTitleSize(0.06,"x")
    sohaStyle.SetLabelSize(0.05,"y")
    sohaStyle.SetTitleSize(0.06,"y")
    sohaStyle.SetLabelSize(0.05,"z")
    sohaStyle.SetTitleSize(0.06,"z")

    # use bold lines and markers
    sohaStyle.SetMarkerStyle(8)
    #sohaStyle.SetHistLineWidth(1.85)
    sohaStyle.SetLineStyleString(2,"[12 12]") # postscript dashes

    # do not display any of the standard histogram decorations
    sohaStyle.SetOptTitle(0)
    sohaStyle.SetOptStat(0)
    sohaStyle.SetOptFit(0)

    # put tick marks on top and RHS of plots
    sohaStyle.SetPadTickX(1)
    sohaStyle.SetPadTickY(1)
    sohaStyle.SetPadTopMargin(0.10)
    sohaStyle.SetTitleOffset(-0.1,"x")
    ROOT.gROOT.SetStyle("SOHA")

    # CL observed pvalues
    pvalue_2016 = dataSet["data"]["2016"]["pList"]
    pvalue_2017 = dataSet["data"]["2017"]["pList"]
    pvalue_Combo = dataSet["data"]["Combo"]["pList"]
    xpoints = dataSet["data"]["2016"]["mList"]
    npoints = len(pvalue_2016)
    Xmin = 300
    Xmax = 900
    Ymin = 0.000000001
    Ymax = 1

    c1 = ROOT.TCanvas("c1","Color Octet PValues 8 TeV",800,600)
    c1.SetFillColor(0)
    c1.cd()
    ROOT.gPad.SetLogy()

    gr_Combo = ROOT.TGraph(npoints, array('d', xpoints), array('d', pvalue_Combo))
    gr_Combo.SetLineColor(1)
    gr_Combo.SetLineWidth(3)
    gr_Combo.SetMaximum(Ymax)
    gr_Combo.SetMinimum(Ymin)
    gr_Combo.GetXaxis().SetLimits(Xmin,Xmax)
    gr_Combo.GetXaxis().SetTitle("m_{#Theta^{0}} [GeV]")
    #gr_Combo.GetXaxis().SetTitleSize(0.06)
    gr_Combo.GetXaxis().SetTitleOffset(1.0)
    gr_Combo.GetYaxis().SetTitle("local p-value")
    gr_Combo.Draw("AL")

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
    numSigma = 5
    for s in range(1, numSigma+1):
        sigma = 0.5 - ROOT.TMath.Erf(s/ROOT.TMath.Sqrt(2))/2
        L = ROOT.TLine(Xmin, sigma, Xmax, sigma)
        L.SetLineColor(2)
        L.SetLineWidth(2)
        L.Draw("same")

        S = ROOT.TPaveText(920,sigma-0.25*sigma,955,sigma+0.5*sigma,"")
        S.SetBorderSize(0)
        S.SetFillStyle(0)
        S.SetTextColor(2)
        S.AddText( str(s)+"#sigma" )
        S.Draw("same")
        entries.append((L,S))

    gr_2016.Draw("L,same")
    gr_2017.Draw("L,same")
    gr_Combo.Draw("L,same")

    legend1 = ROOT.TLegend(0.54397,0.183566,0.997487,0.363636,"")
    legend1.SetHeader("pp #rightarrow #Theta^{0}#Theta^{0} #rightarrow Zgbb (m_{G'} = 2.3 m_{#Theta^{0}})")
    legend1.AddEntry(gr_Combo, "Combined Observed", "l")
    legend1.AddEntry(gr_2016, "2016 Observed", "l")
    legend1.AddEntry(gr_2017, "2017 Observed", "l")
    legend1.SetBorderSize(0)
    legend1.SetFillStyle(0)
    legend1.Draw("same")

    cmstext = ROOT.TPaveText(0.233688, 0.910839, 0.895729, 0.963287, "ndc")
    cmstext.SetBorderSize(0)
    cmstext.SetFillStyle(0)
    cmstext.SetTextAlign(12)
    cmstext.AddText("CMS Preliminary, #sqrt{s}=8 TeV, L_{Int}=19.7 fb^{-1}")
    cmstext.Draw("same")

    c1.Update()
    c1.Print("test.pdf")

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
                data={"mList":[],"rList":[],"rmList":[],"rpList":[],"sList":[],"pList":[]}
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
                    data["mList"].append(float(mass))
                    data["rList"].append(float(info_r[0]))
                    data["rmList"].append(float(info_r[1]))
                    data["rpList"].append(float(info_r[2]))
                    data["sList"].append(float(line_sig))
                    data["pList"].append(float(line_pvalue))
    
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
        print "------------------------------------------"
        print dataSet["runtype"], dataSet["model"]        
        for year in dataSet["data"]:
            print year, dataSet["data"][year]["pList"]

    makePValuePlot(d[0])

if __name__ == '__main__':
    main()
