# Statistical Tests and Fitting Concepts

Combine is a likelihood based statistical tool. 
That means that it uses the likelihood function to define statistical tests.


## General Framework

Most statistical analyses involve devising some "**test statistic**", $t$, which is simply any real-valued function of the observed data:

$$ t(\mathrm{data}) \in \mathbb{R} $$

The choice of tests statistic is, however, very important.
Ideally a good test statistic should return different values for likely outcomes as compared to unlikely outcomes.

As it turns out, in many situations, extremely useful test statistics, and sometimes optimal ones, can be constructed from the likelihood function itself.

### Rejecting a model

In order to reject a model, all one has to do is observe data which are sufficiently unlikely given that model.

This is a two step process.

1. Determine the expected distribution of the test statistic, given the model.
2. Check the value of the test statistic for the actual observed data.

If the value of the test statistic is considered sufficiently unlikely given the model, then that model is rejected.

Often, we are not interested in rejecting a model on its own, but actually in rejecting one model in favour of another.
Largely, thats because models can be rejected for many trivial and uninteresting reasons (e.g. the wrong value of the luminosity was used).

In order to be confident that we are doing something interesting in rejecting a model, we like to have an alternative model which explains the data better.
For example, we would like to show that we can reject the standard model, while simultaneously *not* rejecting a model with a hypothetical new particle.

In that case we gain confidence in our hypothetical new particle, and we are more confident that we have not simply mismodelled some uniteresting parameter.

#### Discoveries

In order to claim discovery of a new phenomenon, we would like to find some data which reject the standard model, while failing to reject the standard model + new phenomenon.

#### Limit Setting

In order to set limits on a model, for example on the cross section of a new process, the model rejection procedure can be applied.

In this case, we would like to reject the model where the cross section is larger than a certain value, while failing to reject the model where it the cross section is small or 0.

#### Measurements 

In order to measure some parameter, we would like to find ranges of some physics parameters which are not rejected, while simultaneously rejecting other ranges of those parameters.


## Likelihood related statistics

The likelihood itself often forms a good basis for test building test statistics.

Typically the absolute value of the likelihood itself is not very meaningful.
However, quantities such as the ratio of the likelihood at two different points in parameter space are very informative about the relative merits of those two models.

Furthermore, in some cases analytic formula approximating the distribution of likelihood based test statistics are known.
This is an extremely helpful property, because when the expected distribution of the test statistic is known, only the observed value needs to be calculated and then it can be compared to the expected distribution.

On the other hand, in the general case where the distriubtion of the test statistic is not known analytically, one has to generate pseudodata from the model to generate an expected distribution.
This process is subject to limitations of how many pseudodata samples are produced, and can be computationally very expensive.


### The likelihood Ratio

The likelihood ratio




