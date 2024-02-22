from __future__ import absolute_import
from __future__ import print_function
import itertools
import HiggsAnalysis.CombinedLimit.tool_base.utils as utils
import json
import os
from HiggsAnalysis.CombinedLimit.tool_base.opts import OPTS

from HiggsAnalysis.CombinedLimit.tool_base.CombineToolBase import CombineToolBase


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


class T2W(CombineToolBase):
    """
    combineTool.py -M T2W [-m MASS] [--cc [card.txt]] [card1.txt some/dir/125/card2.txt some/dir some/dir/125 ...]

    Algorithm:
        1) Each argument is either a single datacard file or a directory
        2) If a datacard, and if the --cc option is not specified, go to the
           directory where the card is located, do text2workspace.py on that
           card. If -m is not set explicity and the enclosing directory name
           is convertible to float set the -m argument to this, otherwise the
           -m argument will not be used.
        3) If a directory, cd to it and combine all .txt files, then do
           text2workspace. If the --cc option is specified, use the given name
           for the combined card that is created, otherwise a default name will
           be used. Then do text2workspace.py on the combined card, following the
           same rule with the -m option as 2)
        4) If individual datacards are included in the list and the --cc
           option is used then combine all of these cards first. As these
           cards could be in different directories the combined card and
           workspace will be created in the current directory.
    """

    description = "Run text2workspace.py on multiple cards or directories"
    requires_root = False
    default_card = "combined.txt"

    def __init__(self):
        CombineToolBase.__init__(self)

    def attach_intercept_args(self, group):
        CombineToolBase.attach_intercept_args(self, group)
        group.add_argument(
            "-m",
            "--mass",
            help="""
            The mass value to set in the text2workspace.py call""",
        )

    def attach_args(self, group):
        CombineToolBase.attach_args(self, group)
        group.add_argument(
            "-i",
            "--input",
            nargs="+",
            help=""" A list of
            input datacards and directories. For the latter, all .txt files
            within the directory will be combined. If the -m option has not
            been specified and the enclosing directory is a number, this will
            be taken as the mass value to set. """,
        )
        group.add_argument(
            "--cc",
            nargs="?",
            const=self.default_card,
            default=None,
            help=""" Create a combined datacard
            with a specified name from the individual cards given by the -i
            option. Note that if this option is used without an argument a
            default card name will be used. For directory arguments in -i, the
            cards will be combined regardless of whether --cc is specified,
            but can still be used to set the name of the combined card that is
            created. """,
        )

    def set_args(self, known, unknown):
        CombineToolBase.set_args(self, known, unknown)

    def run_method(self):
        # The basic structure of each command - we'll fill in the blanks later
        proto = "pushd %(DIR)s; %(CC)stext2workspace.py %(PASSTHRU)s %(CARD)s; popd"
        proto_cc = "combineCards.py %(CARDS)s &> %(COMBINED)s"
        cc_cards_post = []
        for arg in self.args.input:
            passthru = [x for x in self.passthru]
            # Deal with the directory case first (3)
            if os.path.isdir(arg):
                print(">> Directory %s, looking for datacards" % arg)
                files = sorted([file for file in os.listdir(arg) if file.endswith(".txt")])
                if len(files) == 0:
                    print(">> No .txt files found, skipping this directory")
                    continue
                # else:
                # print '>> Will combine %i cards: %s' % (len(files), ' '.join(files))
                cc_cards = [os.path.splitext(file)[0] + "=" + file for file in files]
                cardname = self.args.cc if self.args.cc is not None else self.default_card
                # put an extra extension to avoid accidentally reusing this combined card
                # in a subsequent combination
                cardname += ".cmb"
                cc_cmd = proto_cc % ({"CARDS": " ".join(cc_cards), "COMBINED": cardname})
                base = os.path.basename(arg)
                if self.args.mass is None and isfloat(base):
                    print(">> Enclosing directory will be treated as mass value %s" % base)
                    passthru.extend(["-m", base])
                elif self.args.mass is not None:
                    passthru.extend(["-m", self.args.mass])
                cmd = proto % ({"DIR": arg, "PASSTHRU": " ".join(passthru), "CARD": cardname, "CC": cc_cmd + "; "})
                self.job_queue.append(cmd)
            # Now do case (2) of a single datacard and --cc isn't specified
            elif self.args.cc is None:
                dirname = os.path.dirname(arg)
                if dirname == "":
                    dirname = "."
                base = os.path.split(dirname)[-1]
                if self.args.mass is None and isfloat(base):
                    print(">> Enclosing directory will be treated as mass value %s" % base)
                    passthru.extend(["-m", base])
                elif self.args.mass is not None:
                    passthru.extend(["-m", self.args.mass])
                cmd = proto % ({"DIR": dirname, "PASSTHRU": " ".join(passthru), "CARD": os.path.basename(arg), "CC": ""})
                self.job_queue.append(cmd)
            # Case (2) where --cc is specified
            else:
                cc_cards_post.append(os.path.splitext(os.path.basename(arg))[0] + "=" + arg)

        # Check if we need to combine some individual cards
        if len(cc_cards_post) > 0:
            passthru = [x for x in self.passthru]
            if self.args.mass is not None:
                passthru.extend(["-m", self.args.mass])
            cc_cmd = proto_cc % ({"CARDS": " ".join(cc_cards_post), "COMBINED": self.args.cc})
            cmd = proto % ({"DIR": ".", "PASSTHRU": " ".join(passthru), "CARD": self.args.cc, "CC": cc_cmd + "; "})
            self.job_queue.append(cmd)
        self.flush_queue()

        # self.put_back_arg('name', '-n')
        # proto = 'text2workspace.py ' + (' '.join(self.passthru))
        # for it in itertools.product(*subbed_vars.values()):
        #     keys = subbed_vars.keys()
        #     dict = {}
        #     for i, k in enumerate(keys):
        #         for tuple_i, tuple_ele in enumerate(k):
        #             dict[tuple_ele] = it[i][tuple_i]
        #     self.job_queue.append(proto % dict)
        # self.flush_queue()
