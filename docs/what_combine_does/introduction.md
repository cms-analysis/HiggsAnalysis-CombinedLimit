# Introduction And Capabilities

Combine is a tool for making statistical analyses based on a model of expected observations and a dataset.
Example statistical analyses are claiming discovery of a new particle or process, setting limits on the existence of new physics, and measuring cross sections.

The package has no physics-specific knowledge, it is completely agnostic to the interpretation of the analysis being performed, but its usage and development is based around common cases in High Energy Physics.
This documentation is a description of what combine does and how you can use it to run your analyses.

Roughly, combine does three things:

1. Helps you to build a statistical model of expected observations;
2. Runs statistical tests on the model and observed data;
3. Provides tools for validating, inspecting, and understanding the model and the statistical tests.

Combine can be used for analyses in HEP ranging from simple counting experiments, to unfolded measurements, new physics searches, combinations of measurements and EFT fits.

## Model Building 

Combine provides a powerful, human-readable, and lightweight interface for [building likelihood models](../../part2/settinguptheanalysis/#preparing-the-datacard) for both [binned](../../part2/settinguptheanalysis/#binned-shape-analysis) and [unbinned](../../part2/settinguptheanalysis/#unbinned-or-parametric-shape-analysis) data.
The likelihood definition allows the user to define many processes which contribute to the observation, as well as multiple channels which may be fit simultaneously.

Furthermore, combine provides a powerful and intuitive interface for [combining models](../../part2/settinguptheanalysis/#combination-of-multiple-datacards), as it was originally developped for combinations of higgs boson analysis at the CMS experiment.

The interface simplifies many common tasks, while providing many options for customizations.
Common nuisance parameter types are defined for easy use, while user-defined functions can also be provided.
Input histograms defining the model can be provide in root format, or in other tabular formats compatable with pandas.

Custom [physics models](../../part2/physicsmodels/) can be defined in python which determine how the parameters of interest alter the model, and a number of predefined models are provided by default.

A number of tools are also provided for run-time alterations of the model, allowing for straightforward comparisons of alternative models.

## Statistical Tests

Combine can be used for statistical tests in frequentist or bayesian frameworks as well as some hybrid frequentist-bayesian methods.

Combine implements various methods for [commonly used statistical tests](../../part3/commonstatsmethods/) in high energy physics, including for discovery, limit setting, and parameter estimation.
Statistical tests can be customized to use various test statistics and confidence levels, as well as providing different output formats.

A number of asymptotic methods, relying on Wilks' theorem, and valid in appropriate conditions are implemented for fast evaluation.
Generation of pseudo-data from the model can also be performed, and tests are implemented to automatically run over emprical distributions without relying on asymptotic approximations.
Pseudo-data generation and fitting over the pseudo-data can be customized in a number of ways.

## Validation and Inspection

Combine provides tools for [inspecting the model](../../part3/validation/#validating-datacards) for things like potentially problematic input templates.

[Various methods](../../part3/nonstandard/) are provided for inspecting the likelihood function and the performance of the fits.

Methods are provided for comparing pre-fit and postfit results of all values including nuisance parameters, and summaries of the results can produced.

Plotting utilities allow the pre- and post-fit model expectations and their uncertainties to be plotted, as well as plotted summaries of debugging stups such as the nuisance parameter values and likelihood scans.





