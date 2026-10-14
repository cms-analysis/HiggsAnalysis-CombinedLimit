#!/bin/env python

# This class implements the pdg rounding rules indicated in
# section 5.3 of doi:10.1088/0954-3899/33/1/001
#
# Note: because it uses round (and in general floats), it is affected
# by the limitations described in
# http://docs.python.org/library/functions.html#round
#  and
# http://dx.doi.org/10.1145/103162.103163
#
# davide.gerbaudo@gmail.com
# September 2012

from __future__ import print_function
def pdgRound(value, error) :
    "Given a value and an error, round and format them according to the PDG rules for significant digits"
    def threeDigits(value) :
        "extract the three most significant digits and return them as an int"
        return int(("%.2e"%float(error)).split('e')[0].replace('.','').replace('+','').replace('-',''))
    def nSignificantDigits(threeDigits) :
        assert threeDigits<1000,"three digits (%d) cannot be larger than 10^3"%threeDigits
        if threeDigits<101 : return 2 # not sure
        elif threeDigits<356 : return 2
        elif threeDigits<950 : return 1
        else : return 2
    def frexp10(value) :
        "convert to mantissa+exp representation (same as frex, but in base 10)"
        valueStr = ("%e"%float(value)).split('e')
        return float(valueStr[0]), int(valueStr[1])
    def nDigitsValue(expVal, expErr, nDigitsErr) :
        "compute the number of digits we want for the value, assuming we keep nDigitsErr for the error"
        return expVal-expErr+nDigitsErr
    def formatValue(value, exponent, nDigits, extraRound=0) :
        "Format the value; extraRound is meant for the special case of threeDigits>950"
        roundAt = nDigits-1-exponent - extraRound
        nDec = roundAt if exponent<nDigits else 0
        nDec = max([nDec, 0])
        return ('%.'+str(nDec)+'f')%round(value,roundAt)
    tD = threeDigits(error)
    nD = nSignificantDigits(tD)
    expVal, expErr = frexp10(value)[1], frexp10(error)[1]
    extraRound = 1 if tD>=950 else 0
    return (formatValue(value, expVal, nDigitsValue(expVal, expErr, nD), extraRound),
            formatValue(error,expErr, nD, extraRound))

def test(valueError=(0., 0.)) :
    val, err = valueError
    print(val,' +/- ',err,' --> ', end=' ')
    val, err = pdgRound(val, err)
    print(' ',val,' +/- ',err)

if __name__=='__main__' :
    for x in [(26710, 177)
                ] : test(x)
