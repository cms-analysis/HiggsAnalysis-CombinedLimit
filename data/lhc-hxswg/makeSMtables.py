import logging
import sys

import xlrd

stdHeading = ("mH_GeV", "XS_pb", "Sca_Hi", "Sca_Lo", "Pdf_alpha_s", "Total_pos", "Total_neg")
xsecGroups = {
    "ggH": {"col": "A", "heading": ("mH_GeV", "XS_pb", "Sca_Hi", "Sca_Lo", "Pdf_minus_TH","Gauss","Pdf_alpha_s", "Total_pos", "Total_neg")},
    "vbfH": {"col": "L", "heading": stdHeading},
    "WH": {"col": "T", "heading": stdHeading + ("XS_WplusH_pb", "XS_WminusH_pb")},
    "ZH": {"col": "AD", "heading": stdHeading + ("XS_ggZH_pb",)},
    "ttH": {"col": "AM", "heading": stdHeading},
    "bbH": {
        "col": "AU",
        "heading": ("mH_GeV", "XS_pb", "Sca_Hi", "Sca_Lo", "Pdf_as_mb", "Total_pos", "Total_neg"),
    },
    "tH_tchan": {"col": "BC", "heading": stdHeading + ("XS_tH_pb", "XS_tbarH_pb")},
    "tH_schan": {"col": "BM", "heading": stdHeading + ("XS_tH_pb", "XS_tbarH_pb")},
    "ggZH": {"col": "AD", "heading": ("mH_GeV", "SKIP", "SKIP", "SKIP", "SKIP", "SKIP", "SKIP", "XS_pb")},
    # 'total':  {'col':'BT', 'heading':('XS_pb',)},
    #"WminusH_lv": {"col": "BX", "heading": stdHeading + ("XS_gamma_pb",)},
    #"WplusH_lv": {"col": "CG", "heading": stdHeading + ("XS_gamma_pb",)},
    #"ZH_ll": {"col": "CP", "heading": stdHeading + ("XS_ggZH_pb", "XS_gamma_pb")},
    #"ZH_vv": {"col": "CZ", "heading": stdHeading + ("XS_ggZH_pb", "XS_gamma_pb")},
    #"VBF_qqH_schan": {"col": "DJ", "heading": ("mH_GeV", "XS_pb")},
}

reducedHeading = ("mH_GeV", "XS_pb", "Sca_Hi", "Sca_Lo", "Pdf_alpha_s")
xsecGroupsBSM = {
    "ggH": {"col": "A", "heading": stdHeading + ("1_plus_dEW",)},
    "vbfH": {"col": "J", "heading": stdHeading},
    "WH": {"col": "S", "heading": reducedHeading},
    "ZH": {"col": "AD", "heading": reducedHeading},
    "bbH": {
        "col": "AW",
        "heading": ("mH_GeV", "XS_pb", "Sca_Pdf_mb_mub_Hi", "Sca_Pdf_mb_mub_Lo"),
    },
    "WminusH": {"col": "CF", "heading": reducedHeading},
    "WplusH": {"col": "CO", "heading": reducedHeading},
}

specs = {
    "7 TeV": {
        "rows": (6, 43),
        "groups": xsecGroups,
    },
    "8 TeV": {
        "rows": (6, 43),
        "groups": xsecGroups,
    },
    "13 TeV": {
        "rows": (6, 43),
        "groups": xsecGroups,
    },
    "13.6 TeV": {
        "rows": (6, 43),
        "groups": xsecGroups,
    },
    "14 TeV": {
        "rows": (6, 43),
        "groups": xsecGroups,
    },
    "BSM 7TeV": {
        "rows": (6, 119),
        "groups": xsecGroupsBSM,
    },
    "BSM 8TeV": {
        "rows": (6, 119),
        "groups": xsecGroupsBSM,
    },
    "BSM 13TeV": {
        "rows": (6, 119),
        "groups": xsecGroupsBSM,
    },
    "BSM 14TeV": {
        "rows": (6, 119),
        "groups": xsecGroupsBSM,
    },
    "SM BR": {
        "rows": (7, 44),
        "groups": {
            "BR1": {
                "col": "A",
                "heading": (
                    "mH_GeV",
                    "H_bb",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "H_tautau",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "H_mumu",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "H_ccbar",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                ),
            },
            "BR": {
                "col": "AS",
                "heading": (
                    "mH_GeV",
                    "H_gg",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "H_gamgam",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "H_Zgam",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "H_WW",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "H_ZZ",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                    "Total_Width_GeV",
                    "THU_Hi",
                    "THU_Lo",
                    "PU_mq_Hi",
                    "PU_mq_Lo",
                    "PU_as_Hi",
                    "PU_as_Lo",
                ),
            },
            "BR2": {
                "col": "CK",
                "heading": (
                    "mH_GeV",
                    "H_llll_emt",
                    "H_llll_em",
                    "H_eeee",
                    "H_eemm",
                    "H_llvv_emt",
                    "H_evev",
                    "H_llqq_emt",
                    "H_llqq_em",
                    "H_lvqq_em",
                    "H_vvqq",
                    "H_qqqq",
                    "H_ffff",
                    "DBR",
                ),
            },
        },
    },
}

morespecs = {
    #'sm/xs/7TeV/7TeV-ggH.txt'
}


# import prettytable
def print_table(table, save_file=None):
    col_width = [max(len(x) for x in col) for col in zip(*table)]
    for line in table:
        print("  ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)))

    # if given option, save_file, then also save to the file given too 
    if save_file is not None:
        with open(save_file, "w") as f:
            for line in table:
                f.write("  ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)) + "\n")
        print(f"Saved table to {save_file}")


# Based on http://stackoverflow.com/a/12640614/665025
def col2num(col_str):
    """Convert base26 column string to number."""
    expn = 0
    col_num = 0
    for char in reversed(col_str):
        col_num += (ord(char) - ord("A") + 1) * (26**expn)
        expn += 1
    return col_num


def find_filename(sheet_name, group_name):
    # given the sheet name and group name, find which of the dictionary keys match 
    # and return the filename in the format sm/xs/{sheet_name}/{group_name}.txt
    for key in specs.keys():
        if key in sheet_name:
            for group_key in specs[key]["groups"].keys():
                if group_key in group_name:
                    if "BR" in group_key:
                        filename = f"sm/br/{group_name}.txt"
                    else:
                        xsfoldername = sheet_name.replace(" ", "")
                        filename = f"sm/xs/{xsfoldername}/{xsfoldername}-{group_name}.txt"
    return filename

def formatval(v):
    try:
        # This is to remove pesky unicode symbols like \pm
        v = float(v.encode("ascii", "ignore"))
    except AttributeError:
        pass
    # want to have at most 4 decimal digits.
    v = round(v, 4)
    return "{:+}".format(v)


def main(o):
    f = xlrd.open_workbook(o.input)
    for s in f.sheets():
        try:
            spec = specs[s.name]
        except KeyError:
            logging.info("Skipping sheet [%s]: I do not have parsing rules for it.", s.name)
            continue
        logging.info("Processing sheet [" + s.name + "]")
        for group, props in spec["groups"].items():
            table = []
            logging.info("Processing [" + group + "] in [" + s.name + "]")
            # open output
            # dump heading
            heading = props["heading"]
            real_heading = [h for h in heading if h != "SKIP"]
            table.append(real_heading)
            startRow, endRow = spec["rows"]
            for r in range(startRow - 1, endRow):
                offset = col2num(props["col"]) - 1
                vals = s.row_values(r)[offset : offset + len(heading)]
                # remove headers that say "SKIP"
                vals_keep = []
                for i, v in enumerate(vals):
                    if heading[i] == "SKIP":
                        continue
                    vals_keep.append(v)
                if set(vals[1:]) == set(("",)):
                    continue
                try:
                    table.append(list(map(formatval, vals_keep)))
                except ValueError:
                    print("Could not parse the followig tuple: ")
                    print(vals_keep)
                    raise
            file_name=find_filename(s.name,group)
            print_table(table, save_file=file_name)


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser(usage="%prog -i FILE.xls[x]", version="%prog 3.141")

    parser.add_option("-i", "--input", type="string", help="HXSWG XSBR file", metavar="FILE")
    parser.add_option(
        "-l",
        "--log",
        default="INFO",
        metavar="LEVEL",
        help="Set the minimum logging level.",
    )
    parser.add_option(
        "-s",
        "--savefiles",
        action="store_true",
        help="Whether to save the tables to files in the sm directory.",
    )

    o, args = parser.parse_args()

    if not o.input:
        parser.error("Please specify an input Excel file from the LHC HXSWG.")

    logging.basicConfig(level=getattr(logging, o.log.upper()))

    logging.debug("%s" % str(o))

    sys.exit(main(o))
