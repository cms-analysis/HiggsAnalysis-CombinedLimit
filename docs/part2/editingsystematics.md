# Editing systematic uncertainties

The systematic uncertainties of an existing datacard are renamed, added, removed or changed with `CombineHarvester`: it reads the datacard, the changes are made through its python or C++ interface, and it writes a new datacard.

!!! info
    A datacard can also carry [`nuisance edit` lines](settinguptheanalysis.md#nuisance-modifiers), which <span style="font-variant:small-caps;">Combine</span> applies each time it reads the datacard. `CombineHarvester` doesn't apply them: it keeps the `nuisance edit` lines of a datacard as they are and writes them at the end of the datacards it writes. The changes on this page are made to the model itself, so the systematic lines of the new datacard have them. The [table at the end of this page](#equivalents-of-nuisance-edit-lines) gives the equivalent of each `nuisance edit` line.

## Reading and writing the datacard

```python
import HiggsAnalysis.CombinedLimit.CombineHarvester.ch as ch

cb = ch.CombineHarvester()
cb.ParseDatacard("datacard.txt", mass="125")

# ... the changes described below ...

cb.WriteDatacard("datacard_new.txt", "datacard_new.root")
```

The `mass` is the value for `$MASS` in the `shapes` lines of the datacard, and is not needed for a datacard without it. `WriteDatacard` writes the shapes of the new datacard to the ROOT file given as its second argument. A counting experiment has no shapes, and is written with `cb.WriteDatacard("datacard_new.txt")`.

The examples below use this datacard:

```nohighlight
imax 2
jmax 1
kmax *
---------------
bin         ch_a ch_b
observation 10   20
---------------
bin      ch_a ch_a ch_b ch_b
process  sig  bkg  sig  bkg
process  0    1    0    1
rate     1    10   2    20
---------------
lumi     lnN  1.02 1.02 1.02 1.02
eff      lnN  1.05 -    1.07 -
norm     lnN  -    1.10 -    1.20
bkg_stat lnN  -    1.04 -    1.03
exp group = lumi eff
```

## Selecting the systematic uncertainties to change

`CombineHarvester` has one entry for each systematic uncertainty of each process in each bin. `cb.cp()` followed by filters selects a part of them:

```python
cb.cp().bin(["ch_b"])                    # the bin ch_b
cb.cp().process(["sig"])                 # the process sig, in all bins
cb.cp().signals()                        # the signal processes; also backgrounds()
cb.cp().syst_name(["eff", "lumi"])       # the systematics with these names
cb.cp().syst_type(["lnN"])               # the systematics of this type
cb.cp().bin(["ch_b"]).process(["sig"])   # filters can be combined
```

The selection shares its entries with `cb`: what is changed in the selection is changed in `cb`. Without `cp()`, a filter removes the other entries from `cb` itself.

## Renaming

To rename a systematic uncertainty everywhere:

```python
cb.RenameSystematic(cb, "lumi", "lumi_13TeV")
```

To rename it for some processes or bins only, which makes it independent of the uncertainty that the other ones keep, select them first:

```python
cb.cp().bin(["ch_b"]).process(["sig"]).RenameSystematic(cb, "eff", "eff_b")
```

The datacard then has an `eff` line for the signal in `ch_a` and an `eff_b` line for the signal in `ch_b`. The first argument of `RenameSystematic` is always the complete `cb`, which holds the nuisance parameters. The new name is in the same groups as the old one, here `exp`.

The histograms or pdfs of a shape uncertainty are written to the new ROOT file with the new name.

To rename all the systematic uncertainties that match a pattern, loop over their names:

```python
import re

for name in cb.syst_name_set():
    if re.fullmatch("eff|norm", name):
        cb.RenameSystematic(cb, name, "CMS_" + name)
```

The name of a `param` line is the name of a variable in the input workspace, and so is the name of a `flatParam`, `rateParam` or `discrete` parameter that comes from a workspace: these are renamed in the workspace.

## Adding

`AddSyst` adds a systematic uncertainty to the selected processes. Its arguments are the complete `cb`, the name, the type and the values:

```python
# The same value for the signal in all bins
cb.cp().process(["sig"]).AddSyst(cb, "theory", "lnN", ch.SystMap()(1.1))

# A value for each bin, the second one asymmetric (down, up)
cb.cp().process(["bkg"]).AddSyst(cb, "fakes", "lnN", ch.SystMap("bin")(["ch_a"], 1.2)(["ch_b"], (0.8, 1.3)))
```

The values of a `SystMap` can depend on any of `"bin"`, `"process"`, `"era"`, `"channel"`, `"analysis"`, `"bin_id"` and `"mass"`: `ch.SystMap("bin", "process")(["ch_a"], ["sig"], 1.1)`. A process that the map has no value for gets no systematic uncertainty. The name can contain `$BIN`, `$PROCESS`, `$ERA`, `$CHANNEL`, `$ANALYSIS`, `$BINID` and `$MASS`: `"fakes_$BIN"` gives the independent uncertainties `fakes_ch_a` and `fakes_ch_b`.

`AddSyst` doesn't look at the systematic uncertainties that a process already has. Change the [value](#changing-the-value-or-the-type) of an existing one instead of adding it again.

## Removing

`FilterSysts` removes the systematic uncertainties for which the function returns `True`:

```python
# From the bin ch_a
cb.FilterSysts(lambda s: s.name() == "norm" and s.bin() == "ch_a")

# Everywhere
cb.FilterSysts(lambda s: s.name() == "bkg_stat")
```

The function can use `s.name()`, `s.type()`, `s.bin()`, `s.process()` and the other properties of the entry.

## Changing the value or the type

`ForEachSyst` calls a function for each selected systematic uncertainty:

```python
# A new value
cb.cp().bin(["ch_a"]).syst_name(["eff"]).ForEachSyst(lambda s: s.set_value_u(1.08))


# An asymmetric value
def set_asymm(s):
    s.set_asymm(True)
    s.set_value_d(0.9)
    s.set_value_u(1.15)


cb.cp().bin(["ch_b"]).syst_name(["norm"]).ForEachSyst(set_asymm)

# The opposite direction
cb.cp().syst_name(["bkg_stat"]).ForEachSyst(lambda s: s.set_value_u(1 / s.value_u()))

# Another type: lnN or lnU, or one of the shape types
cb.cp().syst_name(["norm"]).ForEachSyst(lambda s: s.set_type("lnU"))
```

## Splitting

To split a systematic uncertainty into two independent ones, change its values and add the second one:

```python
cb.cp().syst_name(["norm"]).ForEachSyst(lambda s: s.set_value_u(1.06))
cb.cp().process(["bkg"]).AddSyst(cb, "norm_theory", "lnN", ch.SystMap()(1.08))
```

## Merging

To merge the `lnN` uncertainty `bkg_stat` into `norm`, add its values in quadrature and remove it:

```python
import math

stat = {}
cb.cp().syst_name(["bkg_stat"]).ForEachSyst(lambda s: stat.update({(s.bin(), s.process()): s.value_u()}))


def add_stat(s):
    if (s.bin(), s.process()) in stat:
        s.set_value_u(math.exp(math.hypot(math.log(s.value_u()), math.log(stat[s.bin(), s.process()]))))


cb.cp().syst_name(["norm"]).ForEachSyst(add_stat)
cb.FilterSysts(lambda s: s.name() == "bkg_stat")
```

This is for symmetric values that go in the same direction.

## Freezing

Parameters are frozen when running <span style="font-variant:small-caps;">Combine</span>, with `--freezeParameters name1,name2` or, for the [groups](settinguptheanalysis.md#groups-of-nuisances) of the datacard, with `--freezeNuisanceGroups group`: see the [common command-line options](../part3/runningthetool.md#common-command-line-options).

The groups are changed with `SetGroup` and `RemoveGroup`. The names are regular expressions:

```python
cb.SetGroup("bkg", ["norm", "bkg_.*"])   # adds norm and bkg_stat to the group bkg
cb.RemoveGroup("exp", ["eff"])           # removes eff from the group exp
```

To freeze a parameter in the datacard itself, add a `nuisance edit freeze` line to it:

```python
cb.AddDatacardLineAtEnd("nuisance edit freeze lumi")
```

`WriteDatacard` writes this line by itself for a `rateParam` or `extArg` parameter that was frozen with `cb.GetParameter("name").set_frozen(True)`.

## C++ interface

The C++ interface has the same methods:

```cpp
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/CombineHarvester.h"
#include "HiggsAnalysis/CombinedLimit/interface/CombineHarvester/Systematics.h"

ch::CombineHarvester cb;
cb.ParseDatacard("datacard.txt", "", "", "", 0, "125");

// Rename everywhere, and for the signal in one bin
cb.RenameSystematic(cb, "lumi", "lumi_13TeV");
cb.cp().bin({"ch_b"}).process({"sig"}).RenameSystematic(cb, "eff", "eff_b");
// Add a systematic uncertainty to the signal
cb.cp().process({"sig"}).AddSyst(cb, "theory", "lnN", ch::syst::SystMap<>::init(1.1));
// Remove a systematic uncertainty from one bin
cb.FilterSysts([](ch::Systematic const* s) { return s->name() == "norm" && s->bin() == "ch_a"; });
// Change a value
cb.cp().syst_name({"bkg_stat"}).ForEachSyst([](ch::Systematic* s) { s->set_value_u(1.05); });

cb.WriteDatacard("datacard_new.txt");
```

## Equivalents of `nuisance edit` lines

| Line in the datacard | With `CombineHarvester` |
|---|---|
| `nuisance edit rename old new` | [`RenameSystematic`](#renaming) on `cb` |
| `nuisance edit rename process channel old new` | [`RenameSystematic`](#renaming) on a selection of processes and bins |
| `nuisance edit add process channel name pdf value` | [`AddSyst`](#adding); with the `overwrite` option, [`set_value_u`](#changing-the-value-or-the-type) |
| `nuisance edit drop process channel name` | [`FilterSysts`](#removing) |
| `nuisance edit changepdf name newpdf` | [`set_type`](#changing-the-value-or-the-type) |
| `nuisance edit flip process channel name` | [`set_value_u`](#changing-the-value-or-the-type) with the inverse of the value |
| `nuisance edit split process channel old new1 new2 value1 value2` | [`set_value_u` and `AddSyst`](#splitting) |
| `nuisance edit merge process channel name1 name2` | [the sum in quadrature](#merging) |
| `nuisance edit freeze name` | [`AddDatacardLineAtEnd`](#freezing) with the same line, or `--freezeParameters name` when running <span style="font-variant:small-caps;">Combine</span> |
