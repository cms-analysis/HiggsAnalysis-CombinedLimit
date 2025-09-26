# In-Browser Tools

The following tools can be used in-browser to explore several of the data files used by <span style="font-variant:small-caps;">Combine</span>

## jsROOT (external) <sup>[open](https://root.cern/js/latest/)</sup>

Can display (most) contents of `*.root` files, such as histograms e.g. the [`shapes` files](../part2/settinguptheanalysis.md#binned-shape-analyses).

## Interactive Viewers

* nice display of overly large tables, such that headers & important labels stay on screen
* extensive support for live filtering of table contents
* numeric values are colored, with configurable scale

### Datacard Viewer <sup>[open](view_datacard.html)</sup>

* for `*.txt` datacards
* filter by: bin, process, nuisance-name, -type, and -group
* also shows overview or `rate` values for given bin & process combinations
* also shows: `autoMCStats`, `nuisance edit`, `extArg` & `rateParam`
* can detect & report some layout/formatting issue

### Covariance Viewer <sup>[open](view_cov_json.html)</sup>

* for specific `*.cov.json` files, see following format description
* configurable filtering of row/columns with only small off-diagonal elements
* automatically higlight common parts in nuisance parameter names in different colors for easier reading

#### File Format

It needs to be a valid JSON file, with the following fields:

* `labels`, required
    * type: `Arrray` of N `String`s
    * the nuisance/fit parameter labels
* `cov` or `cor`, at least one of them
    * type: `Array` of N `Array`s of N `Number`s between -1 and 1
    * the covariance/correlation matrix
    * ordering of rows/columns corresponds to that or `labels`
* `qual`, optional
    * type: `Number` or `String`
    * the fit result/quality number

#### Example File

```json
{
  "labels": ["lumi", "BR", "r"],
  "cov": [
    [1   , 0.1, 0.01],
    [0.1 , 1  , 0.1 ],
    [0.01, 0.1, 1   ]
  ]
}
```







