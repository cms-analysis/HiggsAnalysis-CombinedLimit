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

```sh
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
