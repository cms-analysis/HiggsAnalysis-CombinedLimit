"""
Performs rounding of values with uncertainties and produces output that can be used in ROOT or LaTeX

Written by andre.david@cern.ch
"""


from __future__ import absolute_import
from __future__ import print_function
from math import *
from decimal import *
from six.moves import range

###
def roundUnc(unc, method="Publication"):
    """By default, rounds uncertainty 'unc' according to the PDG rules plus one significant digit ("Publication").

    Optionally it rounds according with 'method':
        - "PDG" applies the PDG algorithm
        - "Publication" is like "PDG" with an extra significant digit (for results that need to be combined later)
        - "OneDigit" forces one single significant digit (useful when there are multiple uncertainties that vary by more than a factor 10 among themselves)

    Returns a tuple with (uncString, uncMagnitude), where magnitude is the power of 10 that applies to the string to recover the uncertainty.

    """

    # PDG rules (from the Introduction, Section 5.3)
    #
    # Uncertainty leading digits in range:
    #  100 to 354 -> keep 2 digits
    #  355 to 949 -> keep 1 digit
    #  950 to 999 -> keep 2 digits, rounding up to 1000 (e.g. 0.099 -> 0.10, not 0.1)

    uncDigs, uncMagnitude = getDigsMag(unc)

    prec = 1
    unc3Digs = int(round(100*uncDigs))

    if method=="SingleDigit":
        pass
    elif method=="PDG" or method=="Publication":
        if method=="Publication":
            prec += 1
        if 100 <= unc3Digs <= 354:
            prec += 1
    else:
        raise TypeError('Unknown precision method ("%s")'% method)

    uncStr = matchPrec(uncDigs, str(10 ** int(1-prec)))

    # put String in integer form
    uncString = str((Decimal(uncStr)*(10 ** int(prec-1))).quantize(Decimal("1")))
    uncMagnitude -= prec-1

    return (uncString, uncMagnitude)

###
def PDGRoundUnc(unc):
    """Rounds uncertainty unc according to the PDG rules."""

    return roundUnc(unc, "PDG")

###
def matchPrec(val, refStr):
    """Returns a string version of val matching refStr in terms of significant digits."""

    valDec = Decimal(str(val))
    refDec = Decimal(refStr)
    return str(valDec.quantize(refDec))

###
def getDigsMag(val):
    """Extracts the mantissa and exponent of val.

    Returns (valDigs, valMag)

    """
    try:
        valMag = int(floor(log10(val)))
        valDigs = val/pow(10,valMag)
    except:
        print(val)
        valDigs=1
        valMag=1

    return (valDigs, valMag)

###
def PDGRoundSym(val, unc):
    """Rounds a value with a single symmetric uncertainty according to the PDG rules and calculates the order of magnitude of both.    

    Returns (valStr, [uncStr], uncMag)

    """

    assert unc > 0
    uncStr, uncMag = PDGRoundUnc(unc)
    valStr = matchPrec(val/pow(10,uncMag), uncStr)
    return (valStr, [uncStr], uncMag)

###
def PDGRoundAsym(val, uncP, uncM):
    """Rounds a value with a single asymmetric uncertainty according to the PDG rules and calculates the order of magnitude of both.    

    Returns (valStr, [[uncPStr, uncMStr]], uncMag)

    """

    assert uncP > 0
    assert uncM > 0

    uncRef = min(uncP, uncM)
    uncRefStr, uncRefMag = PDGRoundUnc(uncRef)

    uncPStr = matchPrec(uncP/pow(10,uncRefMag), uncRefStr)
    uncMStr = matchPrec(uncM/pow(10,uncRefMag), uncRefStr)
    valStr = matchPrec(val/pow(10,uncRefMag), uncRefStr)

    return (valStr, [[uncPStr, uncMStr]], uncRefMag)

###
def roundMultiple(vals, uncs, method="PDG"):
    """Rounds value with multiple symmetric or asymmetric uncertainties, ignoring the PDG rule when the uncertainty values are too disparate.

    Uncertainties should be a tuple or list of the form
        uncs = (symunc1,(asymP2,asymM2),sym3,etc)

    Returns (valStr, [symunc1,[asymP2,asymM2],sym3,etc], order of magnitude)

    """

    uncList = list()

    if not isinstance(uncs, (list, tuple)):
        uncs = [uncs]

    for unc in uncs:
        try:
            uncList.append(unc[0])
            uncList.append(unc[1])
        except:
            uncList.append(unc)

    uncMin = min(uncList)
    uncMax = max(uncList)

    # If the discrepancy in the uncertainties is too big, downgrade the number of precision digits.
    if uncMax > 10*uncMin:
        if method=="Publication": method = "PDG"
        elif method=="PDG": method = "SingleDigit"

    uncRefStr, uncRefMag = roundUnc(uncMin, method)

    try:
        valsStr = [ matchPrec(val/pow(10,uncRefMag), uncRefStr) for val in vals]
        print("foo")
    except:
        valsStr = matchPrec(vals/pow(10,uncRefMag), uncRefStr)

    uncsStr = list()

    for unc in uncs:
        if isinstance(unc, (list, tuple)):
            elt = [matchPrec(x/pow(10,uncRefMag), uncRefStr) for x in unc]
        else:
            elt = matchPrec(unc/pow(10,uncRefMag), uncRefStr)
        uncsStr.append(elt)

#    print valsStr, uncsStr, uncRefMag

    return (valsStr, uncsStr, uncRefMag)

###
def downgradePrec(valStr, valMag):
    """Returns a string with valStr multiplied by the exponent valMag."""

    # assert valMag<=0
    mag = 10 ** int(valMag)
    return matchPrec(float(valStr) * mag, str(mag))

###
def toROOTRounded(vals, uncs, uncLbls=None, units=None):

    valStr, uncsStr, mag = roundMultiple(vals, uncs)
    return toROOTorLatex(valStr, uncsStr, mag, uncLbls, units, mode="ROOT")

###
def toLatexRounded(vals, uncs, uncLbls=None, units=None):

    valStr, uncsStr, mag = roundMultiple(vals, uncs)
    return toROOTorLatex(valStr, uncsStr, mag, uncLbls, units, mode="Latex")

commonSIPrefixes={
    12 : "T",
     9 : "G",
     6 : "M",
     3 : "k",
    -3 : "m",
    -6 : "\mu ",
    -9 : "n",
   -12 : "p",
   -15 : "f"
   }

###
def toROOT(valStr, uncsStr, mag, uncLbls=None, units=None):
    return toROOTorLatex(valStr, uncsStr, mag, uncLbls, units, mode="ROOT")

###
def toLatex(valStr, uncsStr, mag, uncLbls=None, units=None):
    return toROOTorLatex(valStr, uncsStr, mag, uncLbls, units, mode="Latex")

###
def toROOTorLatex(valStr, uncsStr, mag, uncLbls=None, units=None, mode=None):

    # Used http://www.codecogs.com/latex/eqneditor.php to check results

    if uncLbls:
        assert len(uncsStr) == len(uncLbls)

    salt = -1 if mag>=0 else 0
    magTen = 3 * int( (mag+salt)/3 + 1 )
    #if magTen==-3: magTen=0
    magTgt = mag - magTen

    if mode=="Latex":
        t={
           "sep": "\\",
           "space": "\;",
           "times": "\\times",
           "left": "\\left",
           "right": "\\right",
           }
    elif mode=="ROOT":
        t={
           "sep": "#",
           "space": "",
           "times": "#times",
           "left": "#left",
           "right": "#right",
           }
    else:
        raise TypeError('Unknown mode ("%s")'% mode)


    symUncStr = t["sep"]+"pm%s "
    asymUncStr = "^{+%s}_{-%s} "
    if units and magTen in list(commonSIPrefixes.keys()):
        pwrStr = t["space"]+t["sep"]+"mathrm{"+commonSIPrefixes[magTen]+units+"} "
    else:
        pwrStr = t["times"]+"10^{%d} " % magTen
    lblStr = t["sep"]+"mathrm{(%s)} "

    transform = lambda x: downgradePrec(x, magTgt)

    # Build the string
    outStr = ""

    if mode=="Latex":
        outStr += "$ "

    if magTen and not units: outStr += "[ " #"t["left"]+"( "

    outStr += transform(valStr) + " "

    for i, unc in enumerate(uncsStr):
        if isinstance(unc, (list, tuple)):
            outStr += asymUncStr % (transform(unc[0]), transform(unc[1]))
        else:
            outStr += symUncStr % transform(unc)
        if uncLbls:
            outStr += lblStr % uncLbls[i]


    if magTen:
        if not units: outStr += "] " #"t["right"]+") "
        outStr += pwrStr

    if mode=="Latex":
        outStr += "$ "

    return outStr




import unittest

class PDGRoundingTests(unittest.TestCase):

    knownValues = (
        (0.119,("12",-2)),
        (0.367,("4",-1)),
        (9.99,("10",0)),
        (35002,("35",3)),
        (10.54,("11",0)),
        (0.099,("10",-2)),
    )

    def testPDGRoundUnc(self):
        """Uncertainty roundings according to the PDG rules"""
        for toround, rounded in self.knownValues:
            result = PDGRoundUnc(toround)
            self.assertEquals( result, rounded)

class RoundSymUncTests(unittest.TestCase):

    knownValues = (
        ((0.827,0.119),("83",["12"],-2)),
        ((0.827,0.367),("8",["4"],-1)),
        ((0.827,0.99),("8",["10"],-1)),
        ((100.32,0.843),("1003",["8"],-1)),
        ((10032.34332,8.6234),("10032",["9"],0)),
        ((10583,984),("106",["10"],2)),
        ((10.543e5,73.42e4),("11",["7"],5)),
        ((1.030,0.032),("1030",["32"],-3)),
    )

    def testSymmErrors(self):
        """PDG rules: symmetric errors"""
        for toround, rounded in self.knownValues:
            result = PDGRoundSym(toround[0], toround[1])
            self.assertEquals( result, rounded)

class RoundAsymUncTests(unittest.TestCase):

    knownValues = (
        ((0.827,0.119,0.020),("827",[["119","20"]],-3)),
        ((0.827,0.260,0.025),("827",[["260","25"]],-3)),
    )

    def testAsymmErrors(self):
        """PDG rules: asymmetric errors"""
        for toround, rounded in self.knownValues:
            result = PDGRoundAsym(toround[0], toround[1], toround[2])
            self.assertEquals( result, rounded)

class RoundMultipleTests(unittest.TestCase):

    knownValues = (
        ((0.827,(0.119,(0.020,0.04))),("827",["119",["20","40"]],-3)),
        ((5.234,(0.035,0.361)),("523",["4","36"],-2)),
        ((0.827,[[0.260,0.025]]),("83",[["26","2"]],-2)),
        ((0.827,0.119),("83",["12"],-2)),
        ((1.030,0.032),("1030",["32"],-3)),
    )

    def testAsymmErrors(self):
        """Human rules: multiple symmetric and/or asymmetric errors"""
        for toround, rounded in self.knownValues:
            result = roundMultiple(toround[0], toround[1])
            self.assertEquals( result, rounded)




if __name__ == "__main__":

    # run latex trials here
    #
    print()
    for i in range(-6,6 +1):
        print(toLatexRounded(5.234*pow(10,i),(0.045*pow(10,i),0.361*pow(10,i)),None,"W"))
    print()
    for i in range(-6,6 +1):
        print(toLatexRounded(5.746*pow(10,i),(0.023*pow(10,i),0.954*pow(10,i))))
    print()
    print(toLatexRounded(0.8274e-18,(0.1191e-18,(0.0202e-18,0.0432e-18),0.0582e-18),("stat.","syst.","theo.")))
    print(toLatexRounded(0.8274e-4,(0.1191e-4,(0.0202e-6,0.0432e-4),0.0582e-4),("stat.","syst.","theo."),"b"))
    print()
    for i in range(-6,6 +1):
        print(toLatexRounded(1.030*pow(10,i),(0.032*pow(10,i))))
    print()
    for i in range(-6,6 +1):
        print(toLatexRounded(0.549*pow(10,i),(0.019*pow(10,i),0.063*pow(10,i),0.060*pow(10,i))))
    print()

    print(toROOTRounded(2850e9,(2850e9*0.11)))

    # unit tests exit after running
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner = runner)
