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

To print information to screen, the script parses the json file which contains the results of the validation checks, so if you have already run the validation tool and produced this json file, you can simply change the `printLevel` by re-running the tool with `printLevel` set to a different value, and enabling the `--readOnly` option.

The options `--checkUncertOver` and `--reportSigUnder` will be described in more detail in the section that discusses the checks for which they are relevant.

Note: the `--mass` argument should only be set if you normally use it when running Combine, otherwise you can leave it at the default.

The datacard validation tool is primarily intended for shape (histogram)-based analyses. However, when running on a parametric model or counting experiment the checks for small signal processes, empty processes and uncertainties with large normalisation effects will still be performed. 

## Details on checks 

### Uncertainties with large normalisation effect

This check highlights nuisance parameters which have a normalisation effect larger than the fraction set by the setting `--checkUncertOver`. The default value is 0.1, meaning that any uncertainties with a normalisation
effect larger than 10% are flagged up.

The output file contains the following information for this check:

```nohighlight
largeNormEff: {
  <Uncertainty name>: {
    <analysis category>: {
      <process>: {
        "value_d":<value>
        "value_u":<value>
      } 
    }
  }
}
```

Where `value_u` and `value_d` are the values of the 'up' and 'down' normalisation effects.


### At least one of the Up/Down systematic templates is empty

For shape uncertainties, this check reports all cases where the up and/or down template(s) are empty, when the nominal template is not.

The output file contains the following information for this check:

```nohighlight
emptySystematicShape: {
  <Uncertainty name>: {
    <analysis category>: {
      <process>: {
        "value_d":<value>
        "value_u":<value>
      } 
    }
  }
}
```
Where `value_u` and `value_d` are the values of the 'up' and 'down' normalisation effects.

### Identical Up/Down templates

This check applies to shape uncertainties only, and will highlight cases where the shape uncertainties have identical Up and Down templates (identical in shape and in normalisation).

The information given in the output file for this check is:

```nohighlight
uncertTemplSame: {
  <Uncertainty name>: {
    <analysis category>: {
      <process>: {
        "value_d":<value>
        "value_u":<value>
      } 
    }
  }
}
```
Where `value_u` and `value_d` are the values of the 'up' and 'down' normalisation effects.


### Up and Down templates vary the yield in the same direction

Again this check only applies to shape uncertainties - it highlights cases where the 'Up' template and the 'Down' template both have the effect of increasing or decreasing the normalisation of a process.


The information given in the output file for this check is:

```nohighlight
uncertVarySameDirect: {
  <Uncertainty name>: {
    <analysis category>: {
      <process>: {
        "value_d":<value>
        "value_u":<value>
      } 
    }
  }
}
```
Where `value_u` and `value_d` are the values of the 'up' and 'down' normalisation effects.


### Uncertainty probably has no genuine shape effect

In this check, applying only to shape uncertainties, the normalised nominal templates are compared with the normalised templates for the 'up' and 'down' systematic variations. The script calculates
$$ \Sigma_i \frac{2|\text{up}(i) - \text{nominal}(i)|}{|\text{up}(i)| + |\text{nominal}(i)|}$$ and $$ \Sigma_i \frac{2|\text{down}(i) - \text{nominal}(i)|}{|\text{down}(i)| + |\text{nominal}(i)|} $$

where the sums run over all bins in the histograms, and 'nominal', 'up', and 'down' are the central template and up and down varied templates, all normalised.

If both sums are smaller than 0.001, the uncertainty is flagged up as probably not having a genuine shape effect. This means a 0.1% variation in one bin is enough to avoid being reported, but many smaller variations can also sum to be large enough to pass the threshold.
It should be noted that the chosen threshold is somewhat arbitrary: if an uncertainty is flagged up as probably having no genuine shape effect you should take this as a starting point to investigate. 

The information given in the output file for this check is:

```nohighlight
smallShapeEff: {
  <Uncertainty name>: {
    <analysis category>: {
      <process>: {
        "diff_d":<value>
        "diff_u":<value>
      } 
    }
  }
}
```
Where `diff_d` and `diff_u` are the values of the sums described above for the 'down' variation and the 'up' variation.

### Empty process

If a process is listed in the datacard, but the yield is 0, it is flagged up by this check. 

The information given in the output file for this check is:

```nohighlight
emptyProcessShape: {
  <analysis category>: {
    <process1>,
    <process2>,
    <process3>
  }
}
```


### Bins which have signal but no background 

For shape-based analyses, this checks whether there are any bins in the nominal templates which have signal contributions, but no background contributions. 

The information given in the output file for this check is:

```nohighlight
emptyBkgBin: {
  <analysis category>: {
    <bin_nr1>,
    <bin_nr2>,
    <bin_nr3>
  }
}
```


### Small signal process

This reports signal processes which contribute less than the fraction specified by `--reportSigUnder` (default 0.001 = 0.1%) of the total signal in a given category. This produces an alert, not a warning, as it does not hint at a potential problem.
However, in analyses with many signal contributions and with long fitting times, it can be helpful to remove signals from a category in which they do not contribute a significant amount.

The information given in the output file for this check is:

```nohighlight
smallSignalProc: {
  <analysis category>: {
    <process>: {
      "sigrate_tot":<value>
      "procrate":<value>
    } 
  }
}
```

Where `sigrate_tot` is the total signal yield in the analysis category and `procrate` is the yield of signal process `<process>`.

## What to do in case of a warning

These checks are mostly a tool to help you investigate your datacards: a warning does not necessarily mean there is a mistake in your datacard, but you should use it as a starting point to investigate. Empty processes and emtpy shape uncertainties connected to nonempty processes will most likely be unintended. The same holds for cases where the 'up' and 'down' shape templates are identical. If there are bins which contain signal but no background contributions, this should be corrected. See the [FAQ](/part4/usefullinks#faq) for more information on that point.

For other checks it depends on where the check is fired whether there is a problem or not. Some examples:

- An analysis-specific noncloser uncertainty could be larger than 10%. A theoretical uncertainty in the ttbar normalisation probably not.
- In an analysis with a selection that requires the presence of exactly 1 jet, 'up' and 'down' variations in the jet energy uncertainty *could* both change the process normalisation in the same direction. (But they don't have to!)

As always: think about whether you expect a check to yield a warning in case of your analysis, and investigate to make sure.

