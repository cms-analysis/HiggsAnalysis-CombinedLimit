# Validating datacards

This section covers the main features of the datacard validation tool which helps you spot potential problems with your datacards at an early stage. The tool is implemented
in the [`CombineHarvester/CombineTools`](https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools) subpackage. See the [`combineTool`](/part1/gettingstarted#combine-tool) 
section of the documentation for checkout instructions.

The datacard validation tool contains a number of checks. It is possible to call sub-sets of these checks when creating datacards within CombineHarvester. However, for now we will only
describe the usage of the validation tool on already existing datacards. If you create your datacards with CombineHarvester and would like to include the checks at the datacard creation
stage, please contact us via [https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination.html](https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination.html).


## How to use the tool

The basic syntax is:

```sh
ValidateDatacards.py datacard.txt
```

This will write the results of the checks to a json file (default: `validation.json`), and will print a summary to the screen, for example:

```nohighlight
================================
=======Validation results=======
================================
>>>There were  7800 warnings of type  'up/down templates vary the yield in the same direction'
>>>There were  5323 warnings of type  'up/down templates are identical'
>>>There were no warnings of type  'At least one of the up/down systematic uncertainty templates is empty'
>>>There were  4406 warnings of type  'Uncertainty has normalisation effect of more than 10.0%'
>>>There were  8371 warnings of type  'Uncertainty probably has no genuine shape effect'
>>>There were no warnings of type 'Empty process'
>>>There were no warnings of type 'Bins of the template empty in background'
>>>INFO: there were  169  alerts of type  'Small signal process'
```

The meaning of each of these warnings/alerts is discussed [below](/part3/validation#details-on-checks).

The following arguments are possible:
```nohighlight
usage: ValidateDatacards.py [-h] [--printLevel PRINTLEVEL] [--readOnly]
                            [--checkUncertOver CHECKUNCERTOVER]
                            [--reportSigUnder REPORTSIGUNDER]
                            [--jsonFile JSONFILE] [--mass MASS]
                            cards

positional arguments:
  cards                 Specifies the full path to the datacards to check

optional arguments:
  -h, --help            show this help message and exit
  --printLevel PRINTLEVEL, -p PRINTLEVEL
                        Specify the level of info printing (0-3, default:1)
  --readOnly            If this is enabled, skip validation and only read the
                        output json
  --checkUncertOver CHECKUNCERTOVER, -c CHECKUNCERTOVER
                        Report uncertainties which have a normalisation effect
                        larger than this fraction (default:0.1)
  --reportSigUnder REPORTSIGUNDER, -s REPORTSIGUNDER
                        Report signals contributing less than this fraction of
                        the total in a channel (default:0.001)
  --jsonFile JSONFILE   Path to the json file to read/write results from
                        (default:validation.json)
  --mass MASS           Signal mass to use (default:*)
```
`printLevel` adjusts how much information is printed to the screen. When set to 0, the results are only written to the json file, but not to the screen. When set to 1 (default), the number of warnings/alerts
of a given type is printed to the screen. Setting this option to 2 prints the same information as level 1, and additionally which uncertainties are affected (if the check is related to uncertainties) or which processes are affected (if the check is related only to processes). When `printLevel` is set to 3, the information from level 2 is printed, and additionaly for checks related to uncertainties prints which processes are affected.

What this is doing is parsing the json file which contains the results of the validation checks, so if you have already run the validation tool and produced this json file, you can simply change the `printLevel` by re-running the tool with `printLevel` set to a different value, and enabling the `--readOnly` option.

The options `--checkUncertOver` and `--reportSigUnder` will be described in more detail in the section that discusses the checks for which they are relevant.

Note: the `--mass` argument should only be set if you normally use it when running Combine, otherwise you can leave it at the default.


## Details on checks 

### Uncertainties with large normalisation effect

This check highlights nuisance parameters which have a normalisation effect larger than the fraction set by the setting `--checkUncertOver`. The default value is 0.1, meaning that any uncertainties with a normalisation
effect larger than 10% are flagged up.

### At least one of the Up/Down systematic templates is empty

For shape uncertainties, this check reports all cases where the up and/or down template(s) are empty, when the nominal template is not.

### Identical Up/Down templates

This check applies to shape uncertainties only, and will highlight cases where the shape uncertainties have identical Up and Down templates (identical in shape and in normalisation).

### Up and Down templates vary the yield in the same direction

Again this check only applies to shape uncertainties - it highlights cases where the 'Up' template and the 'Down' template both have the effect of increasing or decreasing the normalisation of a process.

### Uncertainty probably has no genuine shape effect

Experimental!

### Empty process

### Bins which have signal but no background 

### Small signal process

## What to do in case of a warning
