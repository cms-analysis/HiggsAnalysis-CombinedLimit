#!/usr/bin/env python
import ROOT as r 
r.gROOT.SetBatch(True)
import os, sys
import optparse
import fnmatch
import pickle
import aggregateCFG
from array import array
from collections import OrderedDict as odict


def parse_args():

    usage = ('usage: %prog getCorrelationMatrix.py [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)

    parser.add_option('-i', '--inFile',      help='Input mlfit.root file')
    parser.add_option('-o', '--outFileName',default = "simplified_input.root",      help='Output file')
    parser.add_option('-w', '--whichFits',default = "fit_b,prefit,fit_s",  help='Which fit(s) to use')
    parser.add_option('--filterStrings',default = "*",help='Take only bins with this in name (can give comma separated list)')
    parser.add_option('--config',action="store_true",help='choose bins to combine for aggregate regions (NB this overwrites the filter)')
    parser.add_option('--threshold',default = 1.E-5, type=float,help='Only take bins with yield higher than threshold')

    options,args = parser.parse_args()
    return options 
def makeAggregate(aggregateDict,covarianceInput,totalBackground,totalSignal,total,totalData):
    binLabels = [covarianceInput.GetXaxis().GetBinLabel(iBin) for iBin in range(1,covarianceInput.GetNbinsX()+1)]
    # binLabelsAggregate = {}
    binLabelsAggregate = odict()
    for aggregateBinLabel,aggregateList in aggregateDict.iteritems():
        binLabelsAggregate[aggregateBinLabel] = []
        for binIndexMinusOne,binLabel in enumerate(binLabels):
            binIndex = binIndexMinusOne + 1
            if any([fnmatch.fnmatch(binLabel,aggregateFilter) for aggregateFilter in aggregateList]):
                binLabelsAggregate[aggregateBinLabel].append(binIndex)
    consistencyChecks(binLabelsAggregate)
    
    aggregateCovariance = r.TH2D("total_covar","covariance",len(aggregateDict),0,len(aggregateDict),len(aggregateDict),0,len(aggregateDict))
    aggregateBackground = r.TH1D("total_background","background",len(aggregateDict),0,len(aggregateDict))
    aggregateTotal = r.TH1D("total","total",len(aggregateDict),0,len(aggregateDict))
    aggregateData = r.TH1D("total_data","data",len(aggregateDict),0,len(aggregateDict))
    aggregateSignal = r.TH1D("total_signal","signal",len(aggregateDict),0,len(aggregateDict))

    binNumAgg1 = 1
    for label,binNums in binLabelsAggregate.iteritems():
        #set bin labels
        aggregateTotal.GetXaxis().SetBinLabel(binNumAgg1,label)
        aggregateBackground.GetXaxis().SetBinLabel(binNumAgg1,label)
        aggregateSignal.GetXaxis().SetBinLabel(binNumAgg1,label)
        aggregateData.GetXaxis().SetBinLabel(binNumAgg1,label)
        aggregateCovariance.GetXaxis().SetBinLabel(binNumAgg1,label)
        aggregateCovariance.GetYaxis().SetBinLabel(binNumAgg1,label)

        binNumAgg2 = 0
        for label2,binNums2 in binLabelsAggregate.iteritems():
            #covar of X,Y where X=sum(x_j), Y= sum(y_k) is sum_j(sum_k(x_j*y_k))
            totalCov = 0.
            for binNum in binNums:
                for binNum2 in binNums2:
                    totalCov += covarianceInput.GetBinContent(binNum,binNum2)
            aggregateCovariance.Fill(label,label2,totalCov)
            if label != label2:
                aggregateCovariance.Fill(label2,label,totalCov)

            if label == label2:
                break
            binNumAgg2 += 1

        dataX = 0.  
        dataY = 0.
        for binNum in binNums:
            #totalData.GetPoint(binNum-1,dataX[0],dataY[0])
	    dataX = totalData.GetX()[binNum-1]
	    dataY = totalData.GetY()[binNum-1]
            aggregateBackground.AddBinContent(binNumAgg1,totalBackground.GetBinContent(binNum))
            aggregateTotal.AddBinContent(binNumAgg1,total.GetBinContent(binNum))
            aggregateData.AddBinContent(binNumAgg1,dataY)
            aggregateSignal.AddBinContent(binNumAgg1,totalSignal.GetBinContent(binNum))

        aggregateTotal.SetBinError(binNumAgg1,(aggregateCovariance.GetBinContent(binNumAgg1,binNumAgg1))**0.5)
        binNumAgg1 += 1
    return aggregateCovariance,aggregateBackground,aggregateSignal,aggregateTotal,aggregateData


def consistencyChecks(binLabelsAggregate):
    emptyBins = []
    overallList = []
    for binLabelAggregate,binIndices in binLabelsAggregate.iteritems():
        overallList.extend(binIndices)
        if binIndices == []:
            emptyBins.append(binLabelAggregate)
    
    duplicates = set([x for x in overallList if overallList.count(x) > 1])

    if len(emptyBins) != 0:
        raise ValueError, "No bins found for aggregate bins "+str(emptyBins)

    if len(duplicates) != 0:
        raise ValueError, "Duplicates found for bins "+str(duplicates)
    return


def main(filterStrings,inFile,outFileName,whichFits,threshold,config):
    filters = ["*"+i.strip()+"*" for i in filterStrings.split(",")] if filterStrings else ["*"]
    whichFitList = [i.strip() for i in whichFits.split(",")] 
    outFile = r.TFile(outFileName,"RECREATE")
    outFile.cd()
    for whichFit in whichFitList:
        iFile = r.TFile(inFile)
        outDir = outFile.mkdir("shapes_{0}".format(whichFit))
        outDir.cd()
        # get inputs
        covarianceInput = iFile.Get("shapes_{0}/overall_total_covar".format(whichFit))
        totalBackground = iFile.Get("shapes_{0}/total_background".format(whichFit))
        totalSignal = iFile.Get("shapes_{0}/total_signal".format(whichFit))
        total = iFile.Get("shapes_{0}/total_overall".format(whichFit))
        totalData = iFile.Get("shapes_{0}/total_data".format(whichFit))
        totalWidth = iFile.Get("shapes_{0}/total_bin_width".format(whichFit))
        if totalSignal:
            totalSignal.SetName("signalTemp")

        #make restricted set of bins based on filters + thresh
        if covarianceInput:
            binLabels = [covarianceInput.GetXaxis().GetBinLabel(iBin) for iBin in range(1,covarianceInput.GetNbinsX()+1)]
        else:
            binLabels =[]
        binLabelsFiltered = []
        binDict = {}
        if config:
            outCovar,outBackground,outSignal,outTotal,outData = makeAggregate(aggregateCFG.aggregateDict,covarianceInput,totalBackground,totalSignal,total,totalData)
        else:
            for iBinMinusOne,binLabel in enumerate(binLabels):
                #check filters
                if any([fnmatch.fnmatch(binLabel,filterString) for filterString in filters]):
                    print("Found matching bin: " + binLabel)
                    #check threshold
                    if totalBackground.GetBinContent(iBinMinusOne+1) > threshold:
                        binLabelsFiltered.append(binLabel)
                binDict[binLabel] = iBinMinusOne+1
        
            #define outputs
            outBackground = r.TH1D("total_background","total_background",len(binLabelsFiltered),0,len(binLabelsFiltered))
            outSignal = r.TH1D("total_signal","total_signal",len(binLabelsFiltered),0,len(binLabelsFiltered))
            outTotal = r.TH1D("total","total",len(binLabelsFiltered),0,len(binLabelsFiltered))
            outData = r.TH1D("total_data","total_data",len(binLabelsFiltered),0,len(binLabelsFiltered))
            outCovar = r.TH2D("total_covar","total_covar",len(binLabelsFiltered),0,len(binLabelsFiltered),\
                    len(binLabelsFiltered),0,len(binLabelsFiltered))

            #set output contents/errors + labels
            for iBinMinusOne,binLabel in enumerate(binLabelsFiltered):
                #totalData.GetPoint(binDict[binLabel]-1,dataX[0],dataY[0])
	        dataX = totalData.GetX()[binDict[binLabel]-1]
	        dataY = totalData.GetY()[binDict[binLabel]-1]

                # All bin contents are normalized to bin width
                # -> We invert this here so that bin contents are absolute yields
                if totalWidth:
                    bin_width = totalWidth.GetBinContent(binDict[binLabel])
                else:
                    bin_width = 1.0  # No histogram total_bin_width in fitDiagnosticsTest.root
                outData.SetBinContent(iBinMinusOne+1, dataY*bin_width)
                outTotal.SetBinError(iBinMinusOne+1,total.GetBinError(binDict[binLabel])*bin_width)
                outTotal.SetBinContent(iBinMinusOne+1,total.GetBinContent(binDict[binLabel])*bin_width)
                outBackground.SetBinContent(iBinMinusOne+1,totalBackground.GetBinContent(binDict[binLabel])*bin_width)
                outBackground.SetBinError(iBinMinusOne+1,totalBackground.GetBinError(binDict[binLabel])*bin_width)
                outSignal.SetBinError(iBinMinusOne+1,totalSignal.GetBinError(binDict[binLabel])*bin_width)
                outSignal.SetBinContent(iBinMinusOne+1,totalSignal.GetBinContent(binDict[binLabel])*bin_width)

                outBackground.GetXaxis().SetBinLabel(iBinMinusOne+1,binLabel)
                outSignal.GetXaxis().SetBinLabel(iBinMinusOne+1,binLabel)
                outTotal.GetXaxis().SetBinLabel(iBinMinusOne+1,binLabel)
                outData.GetXaxis().SetBinLabel(iBinMinusOne+1,binLabel)
                outCovar.GetXaxis().SetBinLabel(iBinMinusOne+1,binLabel)
                outCovar.GetYaxis().SetBinLabel(iBinMinusOne+1,binLabel)
                for jBinMinusOne,binLabel2 in enumerate(binLabelsFiltered):
                    # The covariance is also bin width normalized
                    if totalWidth:
                        bin_width_2 =  totalWidth.GetBinContent(binDict[binLabel2])
                    else:
                        bin_width_2 = 1.0  # No histogram total_bin_width in fitDiagnosticsTest.root
                    if covarianceInput:
                        cov = covarianceInput.GetBinContent(binDict[binLabel],binDict[binLabel2]) * bin_width * bin_width_2
                    outCovar.SetBinContent(iBinMinusOne+1,jBinMinusOne+1,cov)
        iFile.Close()
        #write it!
        outTotal.Write()
        outBackground.Write()
        outSignal.Write()
        outCovar.Write()
        outData.Write()

    outFile.Close()



if __name__ == "__main__":
    main(**vars(parse_args()))

