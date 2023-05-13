# Statistical Tests 


Combine is a likelihood based statistical tool. 
That means that it uses the likelihood function to define statistical tests.


## General Framework


### Test Statistics

Most statistical analyses involve devising some "**test statistic**", $t$, which is simply any real-valued function of the observed data:

$$ t(\mathrm{data}) \in \mathbb{R} $$

For example, in a simple coin-flipping experiment, the number of heads could be used as the test statistic.

The importance of the test statistic is that it associates with every observation a single value. 
This provides a method of ordering all possible observations.
The ordering allows us to recast the question "how likely was this observation?" in the form of a quantitative question about where the observation falls in the assigned ordering.

The choice of tests statistic is very important as it influences the power of the statistical test.
Ideally a good test statistic should return different values for likely outcomes as compared to unlikely outcomes.

In many situations, extremely useful test statistics, sometimes optimal ones for particular tasks, can be constructed from the likelihood function itself: 

$$t = f(\mathcal{L}(\mathrm{data}))$$


### Statistical tests

Test statistics can be used in order to perform statistical tests, such as parameter estimation, limit setting, or alternative hypothesis testing.

The statistical test is always performed given some model, and allows us to make statements about that model, e.g. that we have found new physics.

The most common statistical test can be characterized by the following steps:

1. Determine the expected distribution of the test statistic, given the model. 
2. Check the value of the test statistic for the actual observed data.
3. Determine the p-value for observing a test statistic at least as extreme as the actual observed value.
4. Draw some conclusion from the p-value.

For example, if the p-value is sufficiently small, the model may be considered rejected;
This would be the case if we were testing a coin for fairness, using the fraction of flips that came up as heads as the test statistic and observed its value to be 0.987.

/// details | **Mathematical details of the general statistical test**

The distribution of the test statistic, $t$ under some model hypothesis $\mathcal{M}$ is:

$$t \stackrel{\mathcal{M}}{\sim} D_{\mathcal{M}}$$

And the observed value of the test statistic is $t_{\mathrm{obs}}$.
The p-value can then be calculated as:

$$p(t_{\mathrm{obs}}) = \int_{\Omega(t_{\mathrm{obs}})} D_{\mathcal{M}}$$

where $\Omega(t_{\mathrm{obs}})$ is a region in the space of the test statistic bounded by $t_{\mathrm{obs}}$.
Its precise definition depends on the test under consideration but it is often either:

$$
[t_{\mathrm{min}}, t_{\mathrm{obs}} ] \mathrm{\ or\ } [t_{\mathrm{obs}}, t_{\mathrm{max}} ] 
$$

Where $t_{\mathrm{min}}$ and  $t_{\mathrm{max}}$ are the lower and upper bounds of the domain of the test statistic.

///

#### Considering Alternative Models in the test

Often, we are not interested in statistical tests performed on a single model in isolation, but in how one model relates to another.

For example, we might find that the data appears very likely given some new physics scenario. 
However, unless that data also appears unlikely given the standard model, this is not a particularly interesting result.

Alternatively if the data appear unlikely given our standard model-based observation model, but we do not have a corresponding new physics scenario for which the data appear likely, it may just be that we've mismodelled some aspect of the detector.

To find new physics, we would like to show that we can reject the standard model, while simultaneously *not* rejecting a model with a hypothetical new particle.
In that case we gain confidence in our hypothetical new particle, and we are more confident that we have not simply mismodelled some other, less interesting aspect of the experiment.

Incorporating alternative models in the statistical testing framework described here can either be done by including them in the definition of the test statistic, or by altering the steps above to include steps based on the alternative model.

## Tests with Likelihood Ratio Test Statistics

The [likelihood function](model_and_likelihood.md#the-likelihood) itself often forms a good basis for building test statistics.

Typically the absolute value of the likelihood itself is not very meaningful as it depends on many fixed aspects we are usually not interested in on their own, like the size of the parameter space and the number of observations.
However, quantities such as the ratio of the likelihood at two different points in parameter space are very informative about the relative merits of those two models.


### The likelihood ratio and likelihood ratio based test statistics

A very useful test statistic is the likelihood ratio of two models: 

$$ \Lambda \equiv \frac{\mathcal{L}_{\mathcal{M}}}{\mathcal{L}_{\mathcal{M}'}} $$ 

For technical and convenience reasons, often the negative logarithm of the likelihood ratio is used: 

$$t \propto -\log(\Lambda) = \log(\mathcal{L}_{\mathcal{M}'}) - \log(\mathcal{L}_{\mathcal{M}})$$

With different proportionality constants being most convenient in different circumstances. 
The negative sign is used by convention since usually the ratios are constructed so that the larger likelihood value must be in the denominator.
This way, $t$ is positive, and larger values of $t$ represent larger differences between the likelihoods of the two models.

#### Sets of test statistics

If the parameters of both likelihoods in the ratio are fixed, either a priori or as a definite function of the data, then that defines a single test statistic.
Often, however, we are interested in testing ["sets" of models](model_and_likelihood.md#sets-of-observation-models), parameterized by some set of values $(\vec{\mu}, \vec{\theta})$.

This is important in limit setting for example, where we perform statistical tests to exclude entire ranges of the parameter space.

In these cases, the likelihood ratio (or a function of it) is not a single test statistic, but a set of test statistics parameterized by the model parameters.
For example, a very useful set of test statistics is:

$$ t_{\vec{x}} \propto -\log(\frac{\mathcal{L}(\vec{x})}{\mathcal{L}(\vec{\hat{x}})}) $$.

Where the likelihood parameters in the bottom are fixed to their maximum likelihood values, but the parameter $\vec{x}$ indexing the test statistic appears in the numerator of the likelihood ratio.

When calculating the p-values for these statistical tests, the p-values are calculated assuming the true model has the parameters used to define the test statistic. 
i.e. as $p_\vec{\mu} \equiv p(t_{\vec{\mu}}(\mathrm{data}) | \mathcal{M}_{\vec{\mu}})$.
In other words, the observed and expected distributions of the test statistics are computed separately at each parameter point $\vec{x}$ being considered.

#### Distributions of likelihood ratio test statistics

Under appropriate conditions, the distribution of $t_\vec{\mu}$ can be approximated analytically, via [Wilks' Theorem](https://en.wikipedia.org/wiki/Wilks'_theorem) or other extensions of that work.
Then, the p-value of the observed test statistic can be calculated from the known form of the expected distribution. 
This is also true for a number of the other test statistics derived from the likelihood ratio.

In the general case, however, the distribution of the test statistic is not known, and it must be estimated.
Typically it is estimated by generating many sets of pseudo-data from the model and using the emprirical distribution of the test statistic.

When the expected distribution is known, computation of the test results is typically much faster, which can be extremely beneficial.
On the other hand, the expected distributions are all aproximations which are sufficiently good under appropriate conditions, and its important to ensure that these approximations work well for the case at hand.


### Parameter Estimation using the likelihood ratio

A common use case for likelihood ratios is estimating the values of some parameters, such as the parameters of interest, $\vec{\mu}$.
The point estimate for the parameters is simply the maximum likelihood estimate, but the likelihood ratio can be used for estimating the uncertainty as a confidence region.

A [confidence region](fitting_concepts.md#frequentist-confidence-regions) for the parameters $\vec{\mu}$ can be defined by using an appropriate test statistic.
Typically, we use the log likelihood ratio:

$$ t_{\vec{\mu}} \propto -\log(\frac{\mathcal{L}(\vec{\mu},\vec{\hat{\theta}}(\vec{\mu}))}{\mathcal{L}(\vec{\hat{\mu}},\vec{\hat{\theta}})}) $$

Where the likelihood in the top is the value of the likelihood at a point $\vec{\mu}$ [profiled over](fitting_concepts.md#profiling) $\vec{\theta}$; and the likelihood on the bottom is at the best fit point.

Then the confience region can be defined as the region where the p-value of the observed test-statistic is less than the confidence level:

$$ \{ \vec{\mu}_{\mathrm{CL}} \} =  \{ \vec{\mu} : p_{\vec{\mu}} \lt \mathrm{CL} \}.$$

This construction will satisfy the frequentist coverage property that the confidence region contains the parameter values used to generate the data in $\mathrm{CL}$ fraction of cases.

In many cases, Wilks' theorem can be used to calculate the p-value and the criteria on $p_{\vec{\mu}}$ can be converted directly into a criterion on $t_{\vec{\mu}}$ itself, $t_{\vec{\mu}} \lt \gamma_{\mathrm{CL}}$.
Where $\gamma_{\mathrm{CL}}$ is a known function of the confidence level which depends on the parameter space being considered.

### Discoveries using the likelihood ratio

A common method for claiming discovery is based on a likelihood ratio test by showing that the new physics model has a "significantly" larger likelihood than the standard model.

This could be done by using the standard log likelihood ratio test statistic:

$$ t_{\mathrm{NP}} = -2\log(\frac{\mathcal{L}(\mu_{\mathrm{NP}} = 0, \vec{\hat{\theta}}(\mu_{\mathrm{NP}} = 0))}{\mathcal{L}(\hat{\mu}_{\mathrm{NP}},\vec{\hat{\theta}})}) $$

Where $\mu_{\mathrm{NP}}$ represents the strength of some new physics quantity, such as the cross section for creation of a new particle.
However, this would also allow for claiming "discovery" in cases where  the best fit value is negative, i.e. $\hat{\mu} \lt 0$, which in particle physics is often an unphysical model, such as a negative cross section.
In order to avoid such a situation, we typically use a modified test statistic:

$$ q_{0} = \begin{cases}
    0 & \hat{\mu} \lt 0 \\ 
    -2\log(\frac{\mathcal{L}(\mathrm{\mu}_{\mathrm{NP}} = 0)}{\mathcal{L}(\hat{\mu}_{\mathrm{NP}})}) & \hat{\mu} \geq 0 
\end{cases}
$$

which excludes the possibility of claiming discovery when the best fit value of $\mu$ is negative.

As with the likelihood ratio test statistic, $t$, defined above, under suitable conditions, analytic expressions for the distribution of $q_0$ are known.

Once the value $q_{0}(\mathrm{data})$ is calculated, it can be compared to the expected distribution of $q_{0}$ under the standard model hypothesis to calculate the p-value.
If the p-value is below some threshold, discovery is often claimed. 
In high-energy physics the standard threshold is $\sim 5\times10^{-7}$.




### Limit Setting using the likelihood ratio

Various test statistics built from likelihood ratios can be used for limit setting, i.e. excluding some parameter values.

The simplest thing one could do to set a limit on a parameter $\mu$ is to find the values of $\mu$ that are outside the confidence regions defined above by using the likelihood ratio test statistic:

$$ t_{\mu} = -2\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})}) $$

However, this could "exclude" $\mu = 0$ or small values of $\mu$ at a typical limit setting confidence level, such as 95%, while still not claiming a discovery.
This is considered undesirable, and often we only want to set upper limits on the value of $\mu$, rather than excluding any possible set of parameters outside our chosen confidence interval.

This can be done using a modified test statistic:

$$ \tilde{t}_{\mu} = -2\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\min(\mu,\hat{\mu}))}) = 
\begin{cases} 
    -2\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})})& \hat{\mu} \lt \mu  \\
    0 &  \mu \leq \hat{\mu}  
\end{cases}
$$

However, this can also have undesirable properties when the best fit value, $\hat{\mu}$, is less than 0.
In that case, we may set limits below 0.
In order to avoid these situations, another modified test statistic can be used:


$$ 
\tilde{q}_{\mu} = \begin{cases} 
    -2\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\mu = 0)})& \hat{\mu} \lt 0  \\
    -2\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})})& 0 \lt \hat{\mu} \lt \mu  \\
    0&  \mu \lt \hat{\mu}  
\end{cases}
$$

In following the [standard test-methodology](#statistical-tests_1) described above, one can then set a limit at a given confidence level $\mathrm{CL}$ by finding the value of $\mu$ for which $p_{\mu} \equiv p(t_{\mu}(\mathrm{data})|\mathcal{M}_{\mu}) = 1 - \mathrm{CL}$. Larger values of $\mu$ will have smaller p-values and are considered excluded at the given confidence level.

However, this procedure is rarely used, in almost every case we use a modified test procedure which uses the $\mathrm{CL}_{s}$ criterion, explained below.

#### The $\mathrm{CL}_{s}$ criterion 

Regardless of which of these test statistics is used, the standard test-methodology has some undesirable properties for limit setting.

Even for an experiment with almost no sensitivity to new physics, 5% of the time the experiment is performed the experimenter will find $p_{\mu} \lt 0.05$ for small values of $\mu$ and set limits on parameter values to which the experiment is not sensitive. 

In order to avoid such situations, the $\mathrm{CL}_{s}$ criterion was developped.
Rather than requiring $p_{\mu} \lt (1-\mathrm{CL})$ to exclude $\mu$, as would be done in the general framework described above, the $\mathrm{CL}_{s}$ criterion requires:

$$ \frac{p_{\mu}}{p_{0}} \lt (1-\mathrm{CL}) $$

Where $p_{\mu}$ is the usual probability of observing the observed value of the test statistic under the signal + background model with signal strength $\mu$ (and $p_0$ for signal strength 0).

Using the $\mathrm{CL}_{s}$ criterion fixes the issue of setting limits much stricter than the experimental sensitivity, because for values of $\vec{\mu}$ to which the experiment is not sensitive $p_{\mu} \approx p_{0}$ and the ratio approaches 1.

Note that this means that a limit set using the $\mathrm{CL}_{s}$ criterion at a given $\mathrm{CL}$ will exclude the true parameter value $\mu$ with a frequency less than the nominal rate of $1-\mathrm{CL}$.
The actual frequency at which it is excluded depends on the sensitivity of the experiment to that parameter value.

### Goodness of fit tests using the likelihood ratio

The likelihood ratio can also be used as a measure of goodness of fit, i.e. describing how well the data match the model for binned data.

A standard likelihood-based measure of the goodness of fit is determined by using the log likelihood ratio with the likelihood in the denominator coming from the **saturated model**.

$$ t_{\mathrm{saturated}} \propto -\log(\frac{\mathcal{L}_{\mathcal{M}}}{\mathcal{L}_{\mathcal{M}_\mathrm{saturated}}}) $$

Here $\mathcal{M}$ is whatever model one is testing the goodness of fit for, and the saturated model is a model for which the prediction matches the observed value in every bin.
Typically, the saturated model would be one in which there are as many free parameters as bins.

This ratio is then providing a comparison between how well the actual data are fit as compared to a hypothetical optimal fit.

Unfortunately, the distribution of $t_{\mathcal{saturated}}$ often is not known a priori and must typically be estimated by generating pseudodata from the model $\mathcal{L}$ and calculating the empirical distribution of the statistic.

Once the distribution is determined a p-value for the statistic can be derived which indicates the probability of observing data with that quality of fit given the model, and therefore serves as a measure of the goodness of fit.

