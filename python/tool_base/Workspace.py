#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import ROOT

import HiggsAnalysis.CombinedLimit.tool_base.utils as utils
from HiggsAnalysis.CombinedLimit.tool_base.opts import OPTS

from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase


class PrintWorkspace(CombineToolBase):
    description = "Load a Workspace and call Print()"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("input", help="The input specified as FILE:WORKSPACE")

    def run_method(self):
        ws_in = self.args.input.split(":")
        f = ROOT.TFile(ws_in[0])
        ws = f.Get(ws_in[1])
        ws.Print()


class ModifyDataSet(CombineToolBase):
    description = "Change the name of a dataset in an existing workspace"
    requires_root = True

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument("input", help="The input specified as FILE:WORKSPACE:DATASET or FILE:WORKSPACE")
        group.add_argument("output", help="The output specified as FILE:WORKSPACE:DATASET or FILE:WORKSPACE")
        group.add_argument("-d", "--data", help="Source data from other file, either FILE:WORKSPACE:DATA or FILE:DATA")

    def run_method(self):
        ws_in = self.args.input.split(":")
        print(">> Input:  " + str(ws_in))
        ws_out = self.args.output.split(":")
        print(">> Output: " + str(ws_out))
        f = ROOT.TFile(ws_in[0])
        ws = f.Get(ws_in[1])
        if len(ws_in) == 3:
            data = ws.data(ws_in[2])
            if len(ws_out) == 3:
                data.SetName(ws_out[2])
        else:
            ws_d = self.args.data.split(":")
            print(">> Data: " + str(ws_d))
            f_d = ROOT.TFile(ws_d[0])
            if len(ws_d) == 2:
                data = f_d.Get(ws_d[1])
            else:
                data = f_d.Get(ws_d[1]).data(ws_d[2])
            if len(ws_out) == 3:
                data.SetName(ws_out[2])
            getattr(ws, "import")(data, ROOT.RooCmdArg())
        ws.SetName(ws_out[1])
        ws.writeToFile(ws_out[0])
