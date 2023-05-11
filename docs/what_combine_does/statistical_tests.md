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
The ordering allows us to recast the question "how likely was this observation?" in the form of a question about where the observation falls in the assigned ordering.

The choice of tests statistic is very important as it influences the power of the statistical test.
Ideally a good test statistic should return different values for likely outcomes as compared to unlikely outcomes.

In many situations, extremely useful test statistics, sometimes optimal ones, can be constructed from the likelihood function itself: 

$$t = f(\mathcal{L}(\mathrm{data}))$$


### Statistical tests

Test statistics can be used in order to perform statistical tests, such as parameter estimation, limit setting, or alternative hypothesis testing.

The statistical test is always performed given some model, and allows us to make statements about that model, e.g. that we have found new physics.

The most common statistical test can be characterized by the following steps:

1. Determine the expected distribution of the test statistic, given the model. 
2. Check the value of the test statistic for the actual observed data.
3. Determine the p-value for observing a test statistic at least as extreme as the actual observed value.
4. Draw some conclusion from the p-value.

For example, if the p-value is sufficiently small, the model may be considered rejected.

<details>
<summary><b> Show Mathematical details of the general statistical test </b></summary>

#### Mathematical Details of the general statistical test

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

</details>

#### Considering Alternative Models in the test

Often, we are not interested in statistical tests performed on a single model in isolation, but in how one model relates to another.

For example, we might find that the data appears very likely given some new physics scenario. 
However, unless that data also appears unlikely given the standard model, this is not a particularly interesting result.

Alternatively if the data appear unlikely given our standard model-based observation model, but we do not have a corresponding new physics scenario for which the data appear likely, it may just be that we've mismodelled some aspect of the detector.

To find new physics, we would like to show that we can reject the standard model, while simultaneously *not* rejecting a model with a hypothetical new particle.
In that case we gain confidence in our hypothetical new particle, and we are more confident that we have not simply mismodelled some uniteresting parameter.

Incorporating alternative models in the statistical testing framework described here can either be done by including them in the definition of the test statistic, or by altering the steps above to include steps based on the alternative model.

## Likelihood related test statistics

The [likelihood](model_and_likelihood.md#the-likelihood) itself often forms a good basis for building test statistics.

Typically the absolute value of the likelihood itself is not very meaningful.
However, quantities such as the ratio of the likelihood at two different points in parameter space are very informative about the relative merits of those two models.


### The likelihood ratio and likelihood ratio based test statistics

The likelihood ratio of two models: 

$$ \Lambda \equiv \frac{\mathcal{L}_{\mathcal{M}}}{\mathcal{L}_{\mathcal{M}'}} $$ 

is a very useful ingredient in statistics.

The absolute value of the likelihood depends on many things, such as the size of the parameter space being considered and of the possible sets of observations.
However, the likelihood ratio is very informative.
It describes how much more likely the data is under one model than the other.

For technical and convenience reasons, often the logarithm of the likelihood ratio: 

$$\log(\Lambda) = \log(\mathcal{L}_{\mathcal{M}}) - \log(\mathcal{L}_{\mathcal{M}'})$$

or some multiple of it, such as $-2\log(\Lambda)$ is used.

Typically, both models being compared are "composite" models, or ["sets" of models](model_and_likelihood.md#sets-of-observation-models), defined by some set of parameters $(\vec{\mu}, \vec{\theta})$.
In these cases, the likelihood ratio (or a function of it) is not a single test statistic, but a set of test statistics parameterized by the model parameters.
If the parameters of the likelihoods in the ratio are fixed, either a priori or based on the data, then that defines a single test statistic.

Often, the parameters in one of the two likelihood terms is fixed, and the other is not.
In this case a different test statistic is used at every point in the parameter space, i.e. for each set of model parameters. 
For example, one can define the family of test statistics:

$$ t_{\vec{x}} = -2\log(\frac{\mathcal{L}(\vec{x})}{\mathcal{L}(\vec{\hat{x}})}) $$.

Where the likelihood parameters in the bottom are fixed to their maximum likelihood values, but the parameter $\vec{x}$ indexing the test statistic appears in the numerator of the likelihood ratio.

In these cases, when calculating the p-values for statistical tests, the p-value is calculated assuming the true model has the parameters used to define the test statistic. 
i.e. as $p_\vec{\mu} = p(t_{\vec{\mu}}(\mathrm{data}) | \mathcal{M}_{\vec{\mu}})$.

### Parameter Estimation using the likelihood ratio

A common use case for likelihood ratios is estimating the values of some parameters, such as the parameters of interest, $\vec{\mu}$.
The point estimate for the parameters is simply the maximum likelihood estimate, but the likelihood ratio can be used for estimating the uncertainty as a confidence set.

A [confidence set](fitting_concepts.md#uncertainties-in-maxmimum-likelihood-fits-confidence-sets) for the parameters $\vec{\mu}$ can be defined by using an appropriate test statistic.
A good choice is the log likelihood ratio:

$$ t_{\vec{\mu}} \propto -\log(\frac{\mathcal{L}(\vec{\mu},\vec{\hat{\theta}}(\vec{\mu}))}{\mathcal{L}(\vec{\hat{\mu}},\vec{\hat{\theta}})}) $$

Where the likelihood in the top is the value of the likelihood at a point $\vec{\mu}$ [profiled over](fitting_concepts.md#profiling) $\vec{\theta}$; and the likelihood on the bottom is at the best fit point.

Then the confience set is the set where the p-value of the observed test-statistic is less than the confidence level:

$$ \{ \vec{\mu}_{\mathrm{CL}} \} =  \{ \vec{\mu} : p_{\vec{\mu}} \lt \mathrm{CL} \}.$$

Under appropriate conditions, the distribution of $t_\vec{\mu}$ is known analytically, via [Wilks' Theorem](https://en.wikipedia.org/wiki/Wilks'_theorem) or other appropriate extensions.
In these cases, one convert the criteria $p_{\vec{\mu}} \lt \mathrm{CL}$ directly into a criteria on $t_{\vec{\mu}}$ itself, $t_{\vec{\mu}} \lt \gamma_{\mathrm{CL}}$.

In the general case, the distribution $t_{\vec{\mu}}$ may not be known a priori and then it must be estimated from pseudo-experiments.

### Discoveries using the likelihood ratio

A common method for claiming discovery is based on a likelihood ratio test by showing that the new physics model has a "significantly" larger likelihood than the standard model.

This could be done by using the standard log likelihood ratio test statistic:

$$ t_{\mathrm{NP}} \propto -\log(\frac{\mathcal{L}_{\mathrm{SM}}}{\mathcal{L}_{\mathrm{SM+NP}}}) $$

However, this would also allow for claiming "discovery" in cases where $\hat{\mu} \lt 0$, which in particle physics is often an unphysical model, such as a negative cross section.
In order to avoid such a situation, we typically use a modified test statistic:

$$ q_{0} \propto \begin{cases}
    0 & \hat{\mu} \lt 0 \\ 
    -\log(\frac{\mathcal{L}_{\mathrm{SM}}}{\mathcal{L}_{\mathrm{SM+NP}}}) & \hat{\mu} \geq 0 
\end{cases}
$$

As with the likelihood ratio test statistic, $t$, defined above, under suitable conditions, analytic expressions for the distribution of $q_0$ are known.

Once the value $\vec{q}_{0}(\mathrm{data})$ is calculated, it can be compared to the expected distribution of $q_{0}$ under the standard model hypothesis to calculate the p-value.
If the p-value is below some threshold, discovery is often claimed. 
In high-energy physics the standard threshold is $\sim 0.0000005$.




### Limit Setting using the likelihood ratio

Various test statistics built from likelihood ratios can be used for limit setting.

The simplest thing one could do to set a limit on a parameter $\mu$ is to find the values of $\mu$ such that $\mu \not \in \{ \mu_{CL} \}$.
This can be done using the likelihood ratio test statistic as defined above:

$$ t_{\mu} \propto -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})}) $$

However, this could "exclude" $\mu = 0$ or small values of $\mu$ at a typical limit setting confidence level, such as 95%, while still not claiming a discovery.
This is considered undesirable, and often we only want to set upper limits on the value of $\mu$.

This can be done using a modified test statistic:

$$ \tilde{t}_{\mu} \propto -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\min(\mu,\hat{\mu}))}) = 
\begin{cases} 
    -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})})& \hat{\mu} \lt \mu  \\
    0 &  \mu \leq \hat{\mu}  
\end{cases}
$$

However, this can also have undesirable properties when the best fit value, $\hat{\mu}$, is less than 0.
In that case, we may set limits below 0.
In order to avoid these situations, another modified test statistic can be used:


$$ 
\tilde{q}_{\mu} \propto \begin{cases} 
    -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\mu = 0)})& \hat{\mu} \lt 0  \\
    -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})})& 0 \lt \hat{\mu} \lt \mu  \\
    0&  \mu \lt \hat{\mu}  
\end{cases}
$$

In following the [standard test-methodology](#statistical-tests_1) described above, one can then set a limit at a given confidence level $\mathrm{CL}$ by finding the value of $\mu$ for which $p_{\mu} = p(t_{\mu}(\mathrm{data})|\mathcal{M}_{\mu}) = 1 - \mathrm{CL}$. Larger values of $\mu$ are excluded at the given confidence level.

However, this procedure is rarely used, in almost every case we use a modified procedure called the $\mathrm{CL}_{s}$ criterion.

#### The $\mathrm{CL}_{s}$ criterion 

Regardless of which test statistic is used, the standard test-methodology described above would lead to setting limits even in cases where there is no sensitivity.
For example, with the usual 95% confidence level for the limit setting ( requiring $p_{\mu} = p(t_{\mu}(\mathrm{data})|\mathcal{M}_{\mu}) \lt 0.05$) would lead to excluding the true value of the parameter 5% of the time.

This can lead to some very undesirable results in the context of limit setting. 
For an experiment with no sensitivity to new physics $\mu$ is indistinguishable from 0 for that experiment.
Nonetheless, 5% of the time this experiment is performed the experimenter will find $p_{0} \lt 0.05$ and will set a limit of 0, seemingly claiming that there is no new physics.
This, despite the fact that the experiment is entirely insensitive to new physics.

In order to avoid such situations, the $\mathrm{CL}_{s}$ criteria was developped.
Rather than requiring $p_{\mu} \lt (1-\mathrm{CL})$ to exclude $\mu$, as would be done in the general framework described above, the $\mathrm{CL}_{s}$ criterion requires:

$$ \frac{p_{\mu}}{p_{0}} \lt (1-\mathrm{CL}) $$

Where $p_{\mu}$ is the usual probability of observing the observed value of the test statistic under the signal + background model with signal strength $\mu$.

Using the $\mathrm{CL}_{s}$ criteria fixes the issue of setting limits using models with no sensitivty, because in that case $p_{\mu} = p_{0}$ and the criteria can never be satisfied for any $\mathrm{CL} \gt 0$.

Note that this means that a limit set using the $\mathrm{CL}_{s}$ criterion at a given $\mathrm{CL}$ will exclude the true parameter value $\mu$ with a frequency less than the nominal rate of $1-\mathrm{CL}$.
The actual frequency at which it is excluded depends on the sensitivity of the experiment to that parameter value.



### Goodness of fit tests using the likelihood ratio



### Distributions of Test Statistics

In order to perform statistical tests as outlined above, the distribution of the test statistic, given some model $\mathcal{M}$ must be known in order to calculate the p-value.

In some cases analytic formula approximating the distribution of likelihood based test statistics are known.
This is an extremely helpful property, because when the expected distribution of the test statistic is known, only the observed value needs to be calculated and then it can be compared to the expected distribution.

On the other hand, in the general case where the distriubtion of the test statistic is not known a priori, one has to generate pseudodata from the model to estimate the expected distribution.
This process is subject to limitations of how many pseudodata samples are produced, and can be computationally very expensive.
