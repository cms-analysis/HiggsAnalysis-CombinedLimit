import os

masses = ["350","450","550","650","750","850"]
signalTypes = ["RPV","SYY","SHH"]
years = ["2016","2017"]
dataTypes = ["Data","pseudoData"]
path = "DoubleTTStuff"

for mass in masses:
    for st in signalTypes:
        for year in years:
            for dt in dataTypes:
                command = "root -l -q -b 'fit_report_ESM.C(\"condor/{4}/Fit_{3}_{0}/output-files/{2}_{1}_{0}/fitDiagnostics{0}{2}{1}.root\",\"condor/{4}/Fit_{3}_{0}/output-files/{2}_{1}_{0}/\")'".format(year,mass,st,dt,path)
                print command
                os.system(command)
