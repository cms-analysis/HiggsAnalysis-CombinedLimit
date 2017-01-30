import xlrd
import sys
import logging

stdHeading = ('mH_GeV','XS_pb','Sca_Hi','Sca_Lo','Pdf_alpha_s','Pdf','alpha_s')
xsecGroups = {
    'ggH': {'col':'A',  'heading':stdHeading},
    'VBF': {'col':'I',  'heading':stdHeading},
    'WH':  {'col':'Q',  'heading':stdHeading+('XS_WplusH_pb','XS_WminusH_pb')},
    'ZH':  {'col':'AA', 'heading':stdHeading+('XS_ggZH_pb',)},
    'ttH': {'col':'AJ', 'heading':stdHeading},
    'bbH': {'col':'AR', 'heading':('mH_GeV','XS_pb','Sca_Pdf_mb_mub_Hi','Sca_Pdf_mb_mub_Lo')},
    'tH_tchan':  {'col':'AZ', 'heading':stdHeading+('XS_tH_pb','XS_tbarH_pb')},
    'tH_schan':  {'col':'BJ', 'heading':stdHeading+('XS_tH_pb','XS_tbarH_pb')},
    # 'total':  {'col':'BT', 'heading':('XS_pb',)},
    'WminusH_lv':  {'col':'BX',  'heading':stdHeading+('XS_gamma_pb',)},
    'WplusH_lv':  {'col':'CG',  'heading':stdHeading+('XS_gamma_pb',)},
    'ZH_ll':  {'col':'CP',  'heading':stdHeading+('XS_ggZH_pb','XS_gamma_pb')},
    'ZH_vv':  {'col':'CZ',  'heading':stdHeading+('XS_ggZH_pb','XS_gamma_pb')},
    'VBF_qqH_schan':  {'col':'DJ',  'heading':('mH_GeV','XS_pb')},
    }

reducedHeading = ('mH_GeV','XS_pb','Sca_Hi','Sca_Lo','Pdf_alpha_s')
xsecGroupsBSM = {
    'ggH': {'col':'A',  'heading':stdHeading+('1_plus_dEW',)},
    'VBF': {'col':'J',  'heading':stdHeading},
    'WH':  {'col':'S',  'heading':reducedHeading},
    'ZH':  {'col':'AD', 'heading':reducedHeading},
    'bbH': {'col':'AW', 'heading':('mH_GeV','XS_pb','Sca_Pdf_mb_mub_Hi','Sca_Pdf_mb_mub_Lo')},
    'WminusH':  {'col':'CF',  'heading':reducedHeading},
    'WplusH':  {'col':'CO',  'heading':reducedHeading},
    }

specs = {
'YR4 SM 7TeV' : {
    'rows' : (6,43),
    'groups' : xsecGroups,
    },

'YR4 SM 8TeV' : {
    'rows' : (6,43),
    'groups' : xsecGroups,
    },

'YR4 SM 13TeV' : {
    'rows' : (6,43),
    'groups' : xsecGroups,
    },

'YR4 SM 14TeV' : {
    'rows' : (6,43),
    'groups' : xsecGroups,
    },

'YR4 BSM 7TeV' : {
    'rows' : (6,119),
    'groups' : xsecGroupsBSM,
    },

'YR4 BSM 8TeV' : {
    'rows' : (6,119),
    'groups' : xsecGroupsBSM,
    },

'YR4 BSM 13TeV' : {
    'rows' : (6,119),
    'groups' : xsecGroupsBSM,
    },

'YR4 BSM 14TeV' : {
    'rows' : (6,119),
    'groups' : xsecGroupsBSM,
    },

'YR4 SM BR' : {
    'rows' : (7,44),
    'groups' : {
        'BR1': {'col':'A', 'heading':('mH_GeV',
            'H_bb', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'H_tautau', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'H_mumu', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'H_ccbar', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            )},
        'BR': {'col':'AS','heading':('mH_GeV',
            'H_gg', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'H_gamgam', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'H_Zgam', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'H_WW', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'H_ZZ', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            'Total_Width_GeV', 'THU_Hi','THU_Lo','PU_mq_Hi','PU_mq_Lo','PU_as_Hi','PU_as_Lo',
            )},
        'BR2': {'col':'CK','heading':('mH_GeV',
            'H_llll_emt', 'H_llll_em', 'H_eeee', 'H_eemm', 'H_llvv_emt',
            'H_evev', 'H_llqq_emt', 'H_llqq_em', 'H_lvqq_em', 'H_vvqq', 'H_qqqq', 'H_ffff', 'DBR',
            )},
        },
    },
}

morespecs = {
#'sm/xs/7TeV/7TeV-ggH.txt'

}

#import prettytable
def print_table(table):
    col_width = [max(len(x) for x in col) for col in zip(*table)]
    for line in table:
        print '  '.join('{:{}}'.format(x, col_width[i]) for i, x in enumerate(line))

# Based on http://stackoverflow.com/a/12640614/665025
def col2num(col_str):
    """ Convert base26 column string to number. """
    expn = 0
    col_num = 0
    for char in reversed(col_str):
        col_num += (ord(char) - ord('A') + 1) * (26 ** expn)
        expn += 1
    return col_num

#import urllib2
#response = urllib2.urlopen('https://twiki.cern.ch/twiki/pub/LHCPhysics/CERNYellowReportPageAt8TeV/Higgs_XSBR_YR3.xlsx')
#f = xlrd.open_workbook(file_contents=response.read())

def formatval(v):
    try:
        #This is to remove pesky unicode symbols like \pm
        v = float(v.encode('ascii', 'ignore'))
    except AttributeError:
        pass
    return '{:+}'.format(v)

def main(o):
    f = xlrd.open_workbook(o.input)
    for s in f.sheets():
        try:
            spec = specs[s.name]
        except KeyError:
            logging.info('Skipping sheet [%s]: I do not have parsing rules for it.', s.name)
            continue
        logging.info('Processing sheet ['+ s.name +']')
        for group,props in spec['groups'].iteritems():
            table = []
            logging.info('Processing ['+group+'] in ['+ s.name +']')
            # open output
            #dump heading
            heading = props['heading']
            table.append(heading)
            (startRow, endRow) = spec['rows']
            for r in xrange(startRow - 1, endRow):
                offset = col2num(props['col']) - 1
                vals = s.row_values(r)[offset:offset+len(heading)]
                if set(vals[1:]) == set(('',)):
                    continue
                try:
                    table.append(map(formatval,vals))
                except ValueError:
                    print 'Could not parse the followig tuple: '
                    print vals
                    raise


            print_table(table)

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog -i FILE.xls[x]',
                          version='%prog 3.141')

    parser.add_option('-i', '--input', type='string',
                      help='HXSWG XSBR file', metavar='FILE')
    parser.add_option('-l', '--log', default='INFO', metavar='LEVEL',
                      help='Set the minimum logging level.')

    (o, args) = parser.parse_args()

    if not o.input:
        parser.error("Please specify an input Excel file from the LHC HXSWG.")

    logging.basicConfig(level=getattr(logging, o.log.upper()))

    logging.debug('%s'% str(o))

    sys.exit(main(o))
