import json
from math import *


def isort(map):
    keys = list(map.keys())
    keys.sort()
    return [(k, map[k]) for k in keys]


def textReport(report):
    tlength = max([10] + [len(t) for t in report.keys()])
    slength = max([max([0] + [len(r) for r in v["results"]]) for v in report.values() if not v and "results" in v])
    tfmt = "%-" + str(tlength) + "s"
    sfmt = "%-" + str(slength) + "s"
    hasref = "has_ref" in report
    if hasref:
        print((tfmt + "  " + sfmt + "  %-40s  %-40s") % ("test", "subtest", "  result", "  reference"))
        print((tfmt + "  " + sfmt + "  %-40s  %-40s") % ("-" * tlength, "-" * slength, "-" * 40, "-" * 40))
    else:
        print((tfmt + "  " + sfmt + "  results") % ("test", "subtest"))
        print((tfmt + "  " + sfmt + "  %s") % ("-" * tlength, "-" * slength, "-" * 40))
    for tn, tv in isort(report):
        if tn == "has_ref":
            continue
        if "results" not in tv or len(tv["results"]) > 1:
            print((tfmt + "  %s: %s") % (tn, tv["status"].upper(), tv["comment"]))
        if "results" not in tv:
            continue
        one = len(tv["results"]) == 1
        for rn, rv in isort(tv["results"]):
            msg = "%-8s " % ((tv["status"] if one else rv["status"]) + ":")
            if rv["status"] != "aborted":
                if rv["limitErr"] != 0:
                    msg += "%-40s" % ("%6.3f +/- %5.3f  (%5.1f min)" % (rv["limit"], rv["limitErr"], rv["t_real"]))
                else:
                    msg += "%-40s" % ("%6.3f            (%5.1f min)" % (rv["limit"], rv["t_real"]))
                if "ref" in rv:
                    if "limit" in rv["ref"]:
                        if rv["ref"]["comment"] != "":
                            msg += "%-40s" % (
                                "%6.3f +/- %5.3f  (%s)"
                                % (
                                    rv["ref"]["limit"],
                                    rv["ref"]["limitErr"],
                                    rv["ref"]["comment"],
                                )
                            )
                        elif rv["ref"]["limitErr"] != 0:
                            msg += "%-40s" % (
                                "%6.3f +/- %5.3f  (%5.1f min)"
                                % (
                                    rv["ref"]["limit"],
                                    rv["ref"]["limitErr"],
                                    rv["ref"]["t_real"],
                                )
                            )
                        else:
                            msg += "%-40s" % ("%6.3f            (%5.1f min)" % (rv["ref"]["limit"], rv["ref"]["t_real"]))
                    else:
                        msg += "%-40s" % rv["ref"]["comment"]
            else:
                comm = tv["comment"] if one else rv["comment"]
                if comm != "":
                    msg += " (" + comm + ")"
            print((tfmt + "  " + sfmt + "  %s") % (tn if one else " ^", rn, msg))


def twikiReport(report):
    statusToIconMap = {
        "aborted": "stop",
        "partial": "processing",
        "mixed": "question",
        "done": "flag",
        "ok": "choice-yes",
        "failed": "choice-no",
        "warning": "warning",
    }
    hasref = "has_ref" in report
    if hasref:
        print("| *Test* | *Subtest* | *%ICON{info}%* | *Result* || *Comments* | *Reference* || *Ref. comments* |")
    else:
        print("| *Test* | *Subtest* | *%ICON{info}%* | *Result* || *Comments* |")
    for tn, tv in isort(report):
        if tn == "has_ref":
            continue
        header = "results" not in tv or (len(tv["results"]) != 1)
        if header:
            if hasref:
                cref = tv["ref_comment"] if "ref_comment" in tv else ""
                print("| <b>%s</b>  |  | %%ICON{%s}%% | <i>%s</i> ||| <i>%s</i>  |||" % (tn, statusToIconMap[tv["status"]], tv["comment"], cref))
            else:
                print("| <b>%s</b>  |  | %%ICON{%s}%% | <i>%s</i> ||  |" % (tn, statusToIconMap[tv["status"]], tv["comment"]))
        if "results" not in tv:
            continue
        for rn, rv in isort(tv["results"]):
            msg = "| %s | %s   | %%ICON{%s}%% | " % (
                "^" if header else "<b>" + tn,
                rn,
                statusToIconMap[rv["status"]],
            )
            if rv["status"] != "aborted":
                if rv["limitErr"] != 0:
                    msg += "%.3f &plusmn; %5.3f   | %.1f min   | %s  |" % (
                        rv["limit"],
                        rv["limitErr"],
                        rv["t_real"],
                        rv["comment"],
                    )
                else:
                    msg += "%.3f             | %.1f min   | %s  | " % (
                        rv["limit"],
                        rv["t_real"],
                        rv["comment"],
                    )
            else:
                msg += "<i>%s</i>  |||" % rv["comment"]
            if "ref" in rv:
                if "limit" in rv["ref"]:
                    if rv["limitErr"] != 0:
                        msg += "%.3f &plusmn; %5.3f   | %.1f min   | %s  |" % (
                            rv["ref"]["limit"],
                            rv["ref"]["limitErr"],
                            rv["ref"]["t_real"],
                            rv["ref"]["comment"],
                        )
                    else:
                        msg += "%.3f                  | %.1f min   | %s  |" % (
                            rv["ref"]["limit"],
                            rv["ref"]["t_real"],
                            rv["ref"]["comment"],
                        )
                else:
                    msg += "%s   |||" % rv["ref"]["comment"]
            elif hasref:
                msg += "  |||"
            print(msg)
