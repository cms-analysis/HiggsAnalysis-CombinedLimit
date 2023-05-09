# Introduction And Preliminaries

Combine is a tool for making statistical analysis based on a model of expected observations and a dataset.
Example tatistical analyses are claiming discovery of a new particle or process, setting limits on the existence of new physics, or measuring physics parameters.

The package is has no physics-specific knowledge, it is completely agnostic to the interpretation of the analysis being performed, but its usage and development is based around common cases in High Energy Physics.
This documentation is a description of what combine does and how you can use it to run your analyses.

Roughly, combine does three things:

1. Helps you to build a statistical model of expected observations;
2. Runs statistical tests on the model and observed data;
3. Provides tools for validating, inspecting, and understanding the model and the statistical tests.

Furthermore, combine provides a powerful and intuitive interface for combining models.

