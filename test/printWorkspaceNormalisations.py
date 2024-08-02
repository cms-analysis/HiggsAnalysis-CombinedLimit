from __future__ import absolute_import, print_function

import datetime
from collections import OrderedDict as od
from optparse import OptionParser
from sys import argv, exit, stderr, stdout
import io
import sys

from six.moves import range

import ROOT

hasHelp = False
for X in ("-h", "-?", "--help"):
    if X in argv:
        hasHelp = True
        argv.remove(X)
argv.append("-b-")

ROOT.gROOT.SetBatch(True)
# ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
argv.remove("-b-")
if hasHelp:
    argv.append("-h")

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option(
    "",
    "--printValueOnly",
    dest="printValueOnly",
    default=False,
    action="store_true",
    help="Just print the default value of the normalisation.",
)
parser.add_option(
    "",
    "--min_threshold",
    dest="min_threshold",
    default=-1.0,
    type="float",
    help="Only print values if yield is greater than this threshold.",
)
parser.add_option(
    "",
    "--max_threshold",
    dest="max_threshold",
    default=-1.0,
    type="float",
    help="Only print values if yield is less than this threshold.",
)
# parser.add_option(
#    "",
#    "--use-cms-histsum",
#    dest="use_cms_histsum",
#    default=False,
#    action="store_true",
#    help="Set to true if workspace built with --use-cms-histsum option.",
# )
parser.add_option(
    "",
    "--procFilter",
    dest="process_filter_list",
    default="",
    type="string",
    help="Filter (keep) only processes containing these names. Enter option as comma separated string",
)
parser.add_option(
    "-m",
    "--mass",
    dest="massVal",
    default=125,
    type=float,
    help="Set mass value in workspace (default=125).",
)
parser.add_option(
    "-f",
    "--format",
    dest="format",
    default="text",
    help="Choose output format [text,html]",
)
parser.add_option(
    "",
    "--output-json",
    dest="output_json",
    default="",
    help="Also output norms in json file with provided name. Note, filters/thresholds are ignored for this output",
)


(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

if options.format not in ["text", "html"]:
    exit("Error, --format must be either 'text' or 'html'")
if options.output_json != "" and options.output_html != "":
    exit("Error, only provide one output option (JSON or HTML)")

if options.max_threshold < options.min_threshold:
    exit("Error - require that --max_threshold is larger than --min_threshold!")

file_in = ROOT.TFile(args[0])
ws = file_in.Get("w")

if ws.var("MH"):
    ws.var("MH").setVal(options.massVal)


def find_chan_proc(name):
    chan = norm_name[norm_name.find("_bin") + len("_bin") : norm_name.find("_proc")]
    if "proc" not in name:
        return chan, ""
    proc = norm_name[norm_name.find("_proc_") + len("_proc_") :]
    return chan, proc


process_filter_list = options.process_filter_list.split(",")

chan_procs = {}

all_norms = ws.allFunctions().selectByName("n_exp_final*")
norm_it = all_norms.createIterator()
for i in range(all_norms.getSize()):
    norm = norm_it.Next()
    norm_name = norm.GetName()
    chan, proc = find_chan_proc(norm_name)
    if len(proc.strip()) == 0:
        continue  # ignore summations

    if chan in chan_procs:
        chan_procs[chan].append([proc, norm, 1])
    else:
        chan_procs[chan] = [[proc, norm, 1]]

# now look for cases where there wasn't a "final"
all_norms = ws.allFunctions().selectByName("n_exp*")
norm_it = all_norms.createIterator()
for i in range(all_norms.getSize()):
    norm = norm_it.Next()
    norm_name = norm.GetName()
    chan, proc = find_chan_proc(norm_name)
    if len(proc.strip()) == 0:
        continue  # ignore summations

    if chan in chan_procs:
        if proc in [chan_procs[chan][i][0] for i in range(len(chan_procs[chan]))]:
            continue
        chan_procs[chan].append([proc, norm, 0, 1])
    else:
        chan_procs[chan] = [[proc, norm, 0, 1]]

# Finally look for the simplest case where the normalisation is a constant fixed number
all_norms = ws.allVars().selectByName("n_exp*")
norm_it = all_norms.createIterator()
for i in range(all_norms.getSize()):
    norm = norm_it.Next()
    norm_name = norm.GetName()
    chan, proc = find_chan_proc(norm_name)
    if len(proc.strip()) == 0:
        continue  # ignore summations

    if chan in chan_procs:
        if proc in [chan_procs[chan][i][0] for i in range(len(chan_procs[chan]))]:
            continue
        chan_procs[chan].append([proc, norm, 0, 0])
    else:
        chan_procs[chan] = [[proc, norm, 0, 0]]

use_cms_histsum = False

all_props = ws.allFunctions().selectByName("prop_bin*")
if all_props.getSize() > 0:
    use_cms_histsum = True

# Look for cases where chan stored as CMSHistSum, set flag
if use_cms_histsum:
    chan_CMSHistSum_norms = {}
    # all_props = ws.allFunctions().selectByName("prop_bin*")
    for chan in chan_procs.keys():
        prop_it = all_props.createIterator()
        for i in range(all_props.getSize()):
            prop = prop_it.Next()
            prop_name = prop.GetName()
            if chan == prop_name.split("_bin")[-1]:
                types = [ROOT.CMSHistSum, ROOT.CMSHistErrorPropagator]
                if type(prop) in types:
                    chan_CMSHistSum_norms[chan] = dict(prop.getProcessNorms())


def checkFilter(proc):
    if not options.process_filter_list:
        return True
    for pc in process_filter_list:
        if pc in proc:
            return True
    return False


def printExpand(proc):
    # print("""<button type="button" class="collapsible">+</button>""")
    # print("""<div class="content">""")
    print('<details><summary style="font-size: 14px; color: blue;">expand %s</summary><p>' % proc)
    print("<pre><code>")


def printEndExpand():
    print("</code></pre>")
    print("</p></details>")


# Now print out information
default_norms = od()

# Save to html formal
if options.format == "html":
    print(
        """
    <html>
    <body>
    <style type="text/css">
    body { font-family: 'Consolas', 'Courier New', courier, monospace; font-size: normal; }
    td, th { border-bottom: 1px solid black; padding: 1px 1em; vertical-align: top; }
    td.channDetails { font-size: x-small; }
    </style>
    <title>Process Normalizations</title>
    </head><body>
    <h1>Process Normalizations</h1>
    """
    )

    print("Normalisation Values Evaluated at MH =", options.massVal)
    for chan in chan_procs.keys():
        default_norms[chan] = od()
        print("<hr>")
        print("<hr>")
        print("<h2>Channel - {chan}</h2>".format(chan=chan))
        chanInfo = chan_procs[chan]
        for proc in chanInfo:
            skipProc = False
            if not checkFilter(proc):
                skipProc = True
            if options.min_threshold > 0:
                skipProc = proc[1].getVal() < options.min_threshold
            if options.max_threshold > 0:
                skipProc = proc[1].getVal() > options.max_threshold
            if skipProc:
                continue
            print("<hr>")
            print("Top-level normalization for process {proc0} -> {proc1_name}<br>".format(proc0=proc[0], proc1_name=proc[1].GetName()))
            if use_cms_histsum:
                if chan in chan_CMSHistSum_norms:
                    default_val = chan_CMSHistSum_norms[chan][proc[1].GetName()]
                else:
                    default_val = proc[1].getVal()
            else:
                default_val = proc[1].getVal()
            default_norms[chan][proc[1].GetName()] = default_val

            if options.printValueOnly:
                print("default value = {default_val:.5f}<br>".format(default_val=default_val))
            else:
                if proc[2]:
                    proc_norm_var = ws.function("n_exp_bin%s_proc_%s" % (chan, proc[0]))
                    printExpand(proc[1].GetName())
                    proc[1].Print()
                    printEndExpand()
                    if proc_norm_var.Class().GetName() == ROOT.ProcessNormalization().Class().GetName():
                        print("... is a product, which contains {proc_norm_var_name}<br>".format(proc_norm_var_name=proc_norm_var.GetName()))
                        printExpand(proc_norm_var.GetName())
                        proc_norm_var.dump()
                        printEndExpand()
                    else:
                        print("... is a product, which contains {proc_norm_var_name}<br>".format(proc_norm_var_name=proc_norm_var.GetName()))
                        printExpand(proc_norm_var.GetName())
                        proc_norm_var.Print()
                        printEndExpand()
                else:
                    if proc[3]:
                        if proc[1].Class().GetName() == ROOT.ProcessNormalization().Class().GetName():
                            printExpand(proc[1].GetName())
                            proc[1].dump()
                            printEndExpand()

                        else:
                            printExpand(proc[1].GetName())
                            proc[1].Print()
                            print(" ... is a constant (formula)")
                            printEndExpand()

                    else:
                        printExpand(proc[1].GetName())
                        proc[1].Print()
                        print(" ... is a constant (RooRealVar)")
                        printEndExpand()

            print("  default value = %.5f " % default_val, "<br>")
            print("</tr>")

    print(
        """
    </body>
    </html>
    """
    )

if options.format == "text":
    print("Normalisation Values Evaluated at MH =", options.massVal)
    for chan in chan_procs.keys():
        default_norms[chan] = od()
        print("---------------------------------------------------------------------------")
        print("---------------------------------------------------------------------------")
        print("Channel - %s " % chan)
        chanInfo = chan_procs[chan]
        for proc in chanInfo:
            skipProc = False
            if not checkFilter(proc):
                skipProc = True
            if options.min_threshold > 0:
                skipProc = proc[1].getVal() < options.min_threshold
            if options.max_threshold > 0:
                skipProc = proc[1].getVal() > options.max_threshold
            if skipProc:
                continue
            print("---------------------------------------------------------------------------")
            print("  Top-level normalization for process %s -> %s" % (proc[0], proc[1].GetName()))
            print("---------------------------------------------------------------------------")
            if use_cms_histsum:
                if chan in chan_CMSHistSum_norms:
                    default_val = chan_CMSHistSum_norms[chan][proc[1].GetName()]
                else:
                    default_val = proc[1].getVal()
            else:
                default_val = proc[1].getVal()
            default_norms[chan][proc[1].GetName()] = default_val

            if options.printValueOnly:
                print("  default value = %.5f " % default_val)
            # if options.printValueOnly: print " --xcp %s:%s "%(chan,proc[0]),
            else:
                if proc[2]:
                    proc_norm_var = ws.function("n_exp_bin%s_proc_%s" % (chan, proc[0]))
                    proc[1].Print()
                    if proc_norm_var.Class().GetName() == ROOT.ProcessNormalization().Class().GetName():
                        print(" ... is a product, which contains ", proc_norm_var.GetName())
                        proc_norm_var.dump()
                    else:
                        print(" ... is a product, which contains ", proc_norm_var.GetName())
                        # proc_norm_var = ws.var("n_exp_bin%s_proc_%s"%(chan,proc[0]))
                        proc_norm_var.Print()
                else:
                    if proc[3]:
                        if proc[1].Class().GetName() == ROOT.ProcessNormalization().Class().GetName():
                            proc[1].dump()
                        else:
                            proc[1].Print()
                            print(" ... is a constant (formula)")
                    else:
                        proc[1].Print()
                        print(" ... is a constant (RooRealVar)")
                print("  -------------------------------------------------------------------------")
                print("  default value = %.5f " % default_val)

# Save norms to json file
if options.output_json != "":
    if ".json" not in options.output_json:
        options.output_json += ".json"
    with open(options.output_json, "w") as jf:
        jf.write("{\n")
        chans = list(default_norms.keys())
        for chan in chans:
            jf.write('    "%s":{\n' % chan)
            procs = list(default_norms[chan].keys())
            for proc in default_norms[chan]:
                if proc == procs[-1]:
                    jf.write('        "%s":%.8f\n' % (proc, default_norms[chan][proc]))
                else:
                    jf.write('        "%s":%.8f,\n' % (proc, default_norms[chan][proc]))
            if chan == chans[-1]:
                jf.write("    }\n")
            else:
                jf.write("    },\n")
        jf.write("}")
