# make a table of the signal strength and associated signifances. 
# Do this for every signal point

# First do it in Data, for 2016, 2017, and Combo
import optparse

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
    
    path = options.basedir
    
    pre_tabular = """\\begin{tabular}{l l l l}
    Mass & Best fit signal strength & Observed Significance & p-value\\\\ \hline
    """    
    runtypes = ["Data", "pseudoDataS"]
    models = ["RPV","SYY"]
    years = ["2016","2017","Combo"]
    masses = ["300","350","400","450","500","550","600","650","700","750","800","850","900"]

    l = []
    for runtype in runtypes:
        for model in models:
            outFileName = "table_signal_strength_%s_%s_%s" % (model, runtype, path)
            file_table = open(outFileName,'w')
            file_table.write(pre_tabular)
            for year in years:
                file_table.write("\\multicolumn{4}{c}{%s} \\\\ \\hline \n"%year)
                for mass in masses:
                    print "Year %s, Model %s, Mass %s"%(year, model, mass)
                    filename_r = "%s/Fit_%s_%s/output-files/%s_%s_%s/log_%s%s%s_FitDiag.txt" % (path,runtype, year, model, mass, year, year, model, mass)
                    info_r = ["0","0","0"]
                    #if year in years:
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
                            #line_sig = file_sig.readlines()[-2]
                            if "Significance:" in line:
                                line_sig = line.replace("Significance:", "").strip()
                            elif "p-value" in line:
                                line_pvalue = line.replace("       (p-value =", "").strip()
                                line_pvalue = line_pvalue.replace(")","").strip()
                    if len(info_r) < 3:
                        file_table.write("%s & %s & %s & %s\\\\ \n" % (mass, info_r[0], line_sig, line_pvalue))
                    elif "#" in info_r[0]:
                        file_table.write("%s & Fit failed &  \\\\ \n" % (mass))
                    elif line_sig == "":
                        file_table.write("%s & $%.2f_{%.2f}^{%.2f}$ & %s & %s\\\\ \n" % (mass, float(info_r[0]), float(info_r[1]), float(info_r[2].replace("+-","+")), line_sig, line_pvalue))
                    else:
                        file_table.write("%s & $%.2f_{%.2f}^{%.2f}$ & %.2f & %s\\\\ \n" % (mass, float(info_r[0]), float(info_r[1]), float(info_r[2].replace("+-","+")), float(line_sig), line_pvalue))            
                file_table.write("\\hline \n")
            file_table.write("\\end{tabular}\n")
            file_table.close()
            l.append({"model":model, "runtype":runtype, "outFileName":outFileName})

    makeSigTex("table_signal_strength.txt", l)

if __name__ == '__main__':
    main()
