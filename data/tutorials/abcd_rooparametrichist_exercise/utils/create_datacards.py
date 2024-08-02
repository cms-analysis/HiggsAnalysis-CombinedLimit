
## Create workspace storing information with RooParametricHist

import ROOT
import os
import optparse

#add arguments
parser = optparse.OptionParser(description="Option parser")
parser.add_option('-m','--mass', dest='mass',help='Input mass point', default=1500,type=int)
parser.add_option('--deplete_crs_from_signal', dest='deplete_crs_from_signal', help='Deplete the control regions from signal', action='store_true', default=False)
(opt, args) = parser.parse_args()


#main code starting here
def main():

    print ("Creating datacards for the analysis with mass point: ", opt.mass)

    #depletion of control regions from signal
    deplete_str = ""
    if opt.deplete_crs_from_signal:
        print ("Depleting control regions from signal")
        deplete_str = "_depletedCRs"
    else:
        print ("No depletion of control regions from signal")
        deplete_str = ""

    #get current directory
    current_directory = os.getcwd()
    card_output_directory = current_directory + "/example_analysis%s/" % deplete_str + "datacards/" + "mPhi%s" % int(opt.mass) + "/"
    #check if output directory exists, if not exit with error
    if not os.path.exists(card_output_directory):
       print ("Error: Directory does not exist. Please run create_workspace.py first.")
       exit()



    signal = "mPhi%s" % int(opt.mass)

    #OOpen workspace
    workfile_name = "param_ws.root"
    wspace_name = "wspace"

    #systematics applied
    lumi_sys = 1.016
    non_closure_sys = 1.05


    #Adding header datacards
    cardSR_A  = "imax * number of bins \n"
    cardSR_A  += "jmax * number of processes minus 1 \n"
    cardSR_A  += "kmax * number of nuisance parameters\n"
    cardSR_A += "-----------------------------------------------------------------------------------\n"
    
    cardCR_B  = "imax * number of bins \n"
    cardCR_B  += "jmax * number of processes minus 1 \n"
    cardCR_B  += "kmax * number of nuisance parameters\n"
    cardCR_B += "-----------------------------------------------------------------------------------\n"

    cardCR_C  = "imax * number of bins \n"
    cardCR_C  += "jmax * number of processes minus 1 \n"
    cardCR_C  += "kmax * number of nuisance parameters\n"
    cardCR_C += "-----------------------------------------------------------------------------------\n"

    cardCR_D  = "imax * number of bins \n"
    cardCR_D  += "jmax * number of processes minus 1 \n"
    cardCR_D  += "kmax * number of nuisance parameters\n"
    cardCR_D += "-----------------------------------------------------------------------------------\n"


    #Adding shapes section
    cardSR_A  += "shapes   %s  %s    %s    %s\n" % ("data_obs", "A", workfile_name, wspace_name +":data_obs_A")
    cardSR_A  += "shapes   %s  %s    %s    %s\n" % ("Bkg","A", workfile_name, wspace_name +":bkg_A")
    cardSR_A  += "shapes   %s  %s    %s    %s\n" % (signal,"A", workfile_name, wspace_name +":" + signal+"_A")
    cardSR_A  += "-----------------------------------------------------------------------------------\n"
    
    cardCR_B  += "shapes   %s  %s    %s    %s\n" % ("data_obs", "B", workfile_name, wspace_name +":data_obs_B")
    cardCR_B  += "shapes   %s  %s    %s    %s\n" % ("Bkg","B", workfile_name, wspace_name +":bkg_B")
    cardCR_B  += "shapes   %s  %s    %s    %s\n" % (signal,"B", workfile_name, wspace_name +":" + signal+"_B")
    cardCR_B  += "-----------------------------------------------------------------------------------\n"

    cardCR_C  += "shapes   %s  %s    %s    %s\n" % ("data_obs", "C", workfile_name, wspace_name +":data_obs_C")
    cardCR_C  += "shapes   %s  %s    %s    %s\n" % ("Bkg","C", workfile_name, wspace_name +":bkg_C")
    cardCR_C  += "shapes   %s  %s    %s    %s\n" % (signal,"C", workfile_name, wspace_name +":" + signal+"_C")
    cardCR_C  += "-----------------------------------------------------------------------------------\n"

    cardCR_D  += "shapes   %s  %s    %s    %s\n" % ("data_obs", "D", workfile_name, wspace_name +":data_obs_D")
    cardCR_D  += "shapes   %s  %s    %s    %s\n" % ("Bkg","D", workfile_name, wspace_name +":bkg_D")
    cardCR_D  += "shapes   %s  %s    %s    %s\n" % (signal,"D", workfile_name, wspace_name +":" + signal+"_D")
    cardCR_D  += "-----------------------------------------------------------------------------------\n"

    #Adding bin section
    cardSR_A += "bin               %s\n" % "A"
    cardSR_A += "observation       %s\n" % ("-1")
    cardSR_A += "-----------------------------------------------------------------------------------\n"
    cardSR_A += "bin                                     %-43s  %-43s\n" % ("A","A")
    cardSR_A += "process                                 %-43s  %-43s\n" % ("Bkg", signal)
    cardSR_A += "process                                 %-43s %-43s\n" % ("1", "0")
    cardSR_A += "rate                                    %-43s %-43s\n" % ("1", "-1")
    cardSR_A += "-----------------------------------------------------------------------------------\n"

    cardCR_B += "bin               %s\n" % "B"
    cardCR_B += "observation       %s\n" % ("-1")
    cardCR_B += "-----------------------------------------------------------------------------------\n"
    cardCR_B += "bin                                     %-43s %-43s\n" % ("B", "B")
    cardCR_B += "process                                 %-43s %-43s\n" % ("Bkg", signal)
    cardCR_B += "process                                 %-43s %-43s\n" % ("1", "0")
    cardCR_B += "rate                                    %-43s %-43s\n" % ("1", "-1")
    cardCR_B += "-----------------------------------------------------------------------------------\n"

    cardCR_C += "bin               %s\n" % "C"
    cardCR_C += "observation       %s\n" % ("-1")
    cardCR_C += "-----------------------------------------------------------------------------------\n"
    cardCR_C += "bin                                     %-43s %-43s\n" % ("C", "C")
    cardCR_C += "process                                 %-43s %-43s\n" % ("Bkg", signal)
    cardCR_C += "process                                 %-43s %-43s\n" % ("1", "0")
    cardCR_C += "rate                                    %-43s %-43s\n" % ("1", "-1")
    cardCR_C += "-----------------------------------------------------------------------------------\n"

    cardCR_D += "bin               %s\n" % "D"
    cardCR_D += "observation       %s\n" % ("-1")
    cardCR_D += "-----------------------------------------------------------------------------------\n"
    cardCR_D += "bin                                     %-43s %-43s\n" % ("D", "D")
    cardCR_D += "process                                 %-43s %-43s\n" % ("Bkg", signal)
    cardCR_D += "process                                 %-43s %-43s\n" % ("1", "0")
    cardCR_D += "rate                                    %-43s %-43s\n" % ("1", "-1")
    cardCR_D += "-----------------------------------------------------------------------------------\n"

    #add systematics section
    cardSR_A += "lumi      lnN                          %-43s %-43s\n" % ("--", lumi_sys)
    cardSR_A += "BkgRate      lnN                       %-43s %-43s\n" % ("--", non_closure_sys)

    cardCR_B += "lumi      lnN                          %-43s %-43s\n" % ("--", lumi_sys)
    
    cardCR_C += "lumi       lnN                          %-43s %-43s\n" % ("--", lumi_sys)

    cardCR_D += "lumi      lnN                          %-43s %-43s\n" % ("--", lumi_sys)

    #write datacards in output directory
    cardSR_A_file = open(card_output_directory + signal + "_Catany_2018_SR.txt", "w")
    cardSR_A_file.write(cardSR_A)
    cardSR_A_file.close()

    cardCR_B_file = open(card_output_directory + signal + "_Catany_2018_CR_B.txt", "w")
    cardCR_B_file.write(cardCR_B)
    cardCR_B_file.close()

    cardCR_C_file = open(card_output_directory +  signal + "_Catany_2018_CR_C.txt", "w")
    cardCR_C_file.write(cardCR_C)
    cardCR_C_file.close()

    cardCR_D_file = open(card_output_directory + signal + "_Catany_2018_CR_D.txt", "w")
    cardCR_D_file.write(cardCR_D)
    cardCR_D_file.close()

    print ("Datacards created in directory: " + card_output_directory)


if __name__ == "__main__":

    main()






