# Statistical Tests 


Combine is a likelihood based statistical tool. 
That means that it uses the likelihood function to define statistical tests.

Combine provides a number of customization options for each test; as always it is up to the user to chose an appropriate test and options.


## General Framework


### Statistical tests

Combine implements a number of different customizable [statistical tests](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#section.40.3). 
These tests can be used for purposes such as determining the significance of some new physics model over the standard model, setting limits, estimating parameters, and checking goodness of fit.

These tests are all performed on a given model (null hypothesis), and often require additional specification of an alternative model.
The statistical test then typically requires defining some "**test statistic**", $t$, which is simply any real-valued function of the observed data:

$$ t(\mathrm{data}) \in \mathbb{R} $$

For example, in a simple coin-flipping experiment, the number of heads could be used as the test statistic.

The distribution of the test statistic should be estimated under the null hypothesis (and the alternative hypothesis, if applicable).
Then the value of the test statistic on the actual observed data, $t^{\mathrm{obs}}$ is compared with its expected value under the relevant hypotheses.

<div class="p value explainer", id="pvalexplainer">
<!--- this div exists as a convenient link target to the p-value explainer below. It is placed slightly above the explainer, because otherwise the explainer is covered by the page header-->
</div>

This comparison, which depends on the test in question, defines the results of the test, which may be simple binary results (e.g. this model point is rejected at a given confidence level), or continuous (e.g. defining the degree to which the data are considered surprising, given the model).
Often, as either a final result or as an intermediate step, the p-value of the observed test statistic under a given hypothesis is calculated.



/// details | **How p-values are calculated**


The distribution of the test statistic, $t$ under some model hypothesis $\mathcal{M}$ is:

$$t \stackrel{\mathcal{M}}{\sim} D_{\mathcal{M}}$$

And the observed value of the test statistic is $t_{\mathrm{obs}}$.
The p-value of the observed result gives the probability of having observed a test statistic at least as extreme as the actual observation. For example, this may be:

$$p = \int_{t_{\mathrm{min}}}^{t_\mathrm{obs}} D_{\mathcal{M}} \mathrm{d}t$$

In some cases, the bounds of the integral may be modified, such as $( t_{\mathrm{obs}}, t_{\mathrm{max}} )$ or $(-t_{\mathrm{obs}}, t_{\mathrm{obs}} )$, depending on the details of the test being performed.
And specifically, for the distribution in question, whether an observed value in the right tail, left tail, or either tail of the distribution is considered as unexpected.

The p-values using the left-tail and right tail are related to each other via $p_{\mathrm{left}} = 1 - p_{\mathrm{right}}$.

///

#### Test Statistics

The test statistic can be any real valued function of the data.
While in principle, many valid test statistics can be used, the choice of tests statistic is very important as it influences the power of the statistical test.

By associating a single real value with every observation, the test statistic allows us to recast the question "how likely was this observation?" in the form of a quantitative question about where the value of the test statistic.
Ideally a good test statistic should return different values for likely outcomes as compared to unlikely outcomes and the expected distributions under the null and alternate hypotheses should be well-separated.

In many situations, extremely useful test statistics, sometimes optimal ones for particular tasks, can be constructed from the likelihood function itself: 

$$ t(\mathrm{data}) = f(\mathcal{L}) $$

Even for a given statistical tests, several likelihood-based test-statistics may be suitable, and for some tests combine implements multiple test-statistics from which the user can chose.

## Tests with Likelihood Ratio Test Statistics

The [likelihood function](../../what_combine_does/model_and_likelihood/#the-likelihood) itself often forms a good basis for building test statistics.

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
Often, however, we are interested in testing ["sets" of models](../../what_combine_does/model_and_likelihood/#sets-of-observation-models), parameterized by some set of values $(\vec{\mu}, \vec{\nu})$.

This is important in limit setting for example, where we perform statistical tests to exclude entire ranges of the parameter space.

In these cases, the likelihood ratio (or a function of it) is not a single test statistic, but a set of test statistics parameterized by the model parameters.
For example, a very useful set of test statistics is:

$$ t_{\vec{x}} \propto -\log(\frac{\mathcal{L}(\vec{x})}{\mathcal{L}(\vec{\hat{x}})}) $$.

Where the likelihood parameters in the bottom are fixed to their maximum likelihood values, but the parameter $\vec{x}$ indexing the test statistic appears in the numerator of the likelihood ratio.

When calculating the p-values for these statistical tests, the p-values are calculated at each point in parameter space taking those parammeter values as the null hypothesis,
i.e. as $p_\vec{\mu} \equiv p(t_{\vec{\mu}}(\mathrm{data}) ; \mathcal{M}_{\vec{\mu}})$.
In other words, the observed and expected distributions of the test statistics are computed separately at each parameter point $\vec{x}$ being considered.

#### Expected distributions of likelihood ratio test statistics

Under [appropriate conditions](https://arxiv.org/abs/1911.10237), the distribution of $t_\vec{\mu}$ can be approximated analytically, via [Wilks' Theorem](https://en.wikipedia.org/wiki/Wilks'_theorem) or other extensions of that work.
Then, the p-value of the observed test statistic can be calculated from the known form of the expected distribution. 
This is also true for a number of the other test statistics derived from the likelihood ratio, where [asymptotic approximations have been derived](https://arxiv.org/abs/1007.1727).

Combine provides asymptotic methods, for [limit setting](../../part3/commonstatsmethods/#asymptotic-frequentist-limits) and [significance tests](../../part3/commonstatsmethods/#asymptotic-significances), which make used of these approximations for fast calculations.

In the general case, however, the distribution of the test statistic is not known, and it must be estimated.
Typically it is estimated by generating many sets of pseudo-data from the model and using the emprirical distribution of the test statistic.

Combine provides methods for [limit setting](../../part3/commonstatsmethods/#computing-limits-with-toys) and [significance tests](../../part3/commonstatsmethods/#computing-significances-with-toys) which use pseudodata generation to estimate the expected test-statistic distributions, and therefore don't depend on the asymptotic approximation.
Methods are also provided for [generating pseudodata](../../part3/runningthetool/#toy-data-generation) without running a particular test, which can be saved and used for estimating expected distributions.

### Parameter Estimation using the likelihood ratio

A common use case for likelihood ratios is estimating the values of some parameters, such as the parameters of interest, $\vec{\mu}$.
The point estimate for the parameters is simply the maximum likelihood estimate, but the likelihood ratio can be used for estimating the uncertainty as a confidence region.

A [confidence region](../../what_combine_does/fitting_concepts/#frequentist-confidence-regions) for the parameters $\vec{\mu}$ can be defined by using an appropriate test statistic.
Typically, we use the log likelihood ratio:

$$ t_{\vec{\mu}} \propto -\log(\frac{\mathcal{L}(\vec{\mu},\vec{\hat{\nu}}(\vec{\mu}))}{\mathcal{L}(\vec{\hat{\mu}},\vec{\hat{\nu}})}) $$

Where the likelihood in the top is the value of the likelihood at a point $\vec{\mu}$ [profiled over](../../what_combine_does/fitting_concepts/#profiling) $\vec{\nu}$; and the likelihood on the bottom is at the best fit point.

Then the confidence region can be defined as the region where the p-value of the observed test-statistic is less than the confidence level:

$$ \{ \vec{\mu}_{\mathrm{CL}} \} =  \{ \vec{\mu} : p_{\vec{\mu}} \lt \mathrm{CL} \}.$$

This construction will satisfy the frequentist coverage property that the confidence region contains the parameter values used to generate the data in $\mathrm{CL}$ fraction of cases.

In many cases, Wilks' theorem can be used to calculate the p-value and the criteria on $p_{\vec{\mu}}$ can be converted directly into a criterion on $t_{\vec{\mu}}$ itself, $t_{\vec{\mu}} \lt \gamma_{\mathrm{CL}}$.
Where $\gamma_{\mathrm{CL}}$ is a known function of the confidence level which depends on the parameter space being considered.

### Discoveries using the likelihood ratio

A common method for claiming discovery is based on a likelihood ratio test by showing that the new physics model has a "significantly" larger likelihood than the standard model.

This could be done by using the standard log likelihood ratio test statistic:

$$ t_{\mathrm{NP}} = -2\log(\frac{\mathcal{L}(\mu_{\mathrm{NP}} = 0, \vec{\hat{\nu}}(\mu_{\mathrm{NP}} = 0))}{\mathcal{L}(\hat{\mu}_{\mathrm{NP}},\vec{\hat{\nu}})}) $$

Where $\mu_{\mathrm{NP}}$ represents the strength of some new physics quantity, such as the cross section for creation of a new particle.
However, this would also allow for claiming "discovery" in cases where  the best fit value is negative, i.e. $\hat{\mu} \lt 0$, which in particle physics is often an unphysical model, such as a negative cross section.
In order to avoid such a situation, we typically use a modified test statistic:

$$ q_{0} = \begin{cases}
    0 & \hat{\mu} \lt 0 \\ 
    -2\log(\frac{\mathcal{L}(\mathrm{\mu}_{\mathrm{NP}} = 0)}{\mathcal{L}(\hat{\mu}_{\mathrm{NP}})}) & \hat{\mu} \geq 0 
\end{cases}
$$

which excludes the possibility of claiming discovery when the best fit value of $\mu$ is negative.

As with the likelihood ratio test statistic, $t$, defined above, under suitable conditions, [analytic expressions for the distribution of $q_0$](https://ar5iv.labs.arxiv.org/html/1007.1727#S3.SS5) are known.

Once the value $q_{0}(\mathrm{data})$ is calculated, it can be compared to the expected distribution of $q_{0}$ under the standard model hypothesis to calculate the p-value.
If the p-value is below some threshold, discovery is often claimed. 
In high-energy physics the standard threshold is $\sim 5\times10^{-7}$.




### Limit Setting using the likelihood ratio

Various test statistics built from likelihood ratios can be used for limit setting, i.e. excluding some parameter values.

One could do to set limits on a parameter $\mu$ by finding the values of $\mu$ that are outside the confidence regions defined above by using the likelihood ratio test statistic:

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

Which also has a [known distribution](https://ar5iv.labs.arxiv.org/html/1007.1727#S3.SS7) under appropriate conditions, or can be estimated from pseudo-experiments. One can then set a limit at a given confidence level, $\mathrm{CL}$, by finding the value of $\mu$ for which $p_{\mu} \equiv p(t_{\mu}(\mathrm{data});\mathcal{M}_{\mu}) = 1 - \mathrm{CL}$. Larger values of $\mu$ will have smaller p-values and are considered excluded at the given confidence level.

However, this procedure is rarely used, in almost every case we use a modified test procedure which uses the $\mathrm{CL}_{s}$ criterion, explained below.

#### The $\mathrm{CL}_{s}$ criterion 

Regardless of which of these test statistics is used, the standard test-methodology has some undesirable properties for limit setting.

Even for an experiment with almost no sensitivity to new physics, 5% of the time the experiment is performed the experimenter will find $p_{\mu} \lt 0.05$ for small values of $\mu$ and set limits on parameter values to which the experiment is not sensitive. 

In order to avoid such situations the $\mathrm{CL}_{s}$ criterion was developped, as explained in these [two](https://cdsweb.cern.ch/record/451614) [papers](https://arxiv.org/abs/hep-ex/9902006).
Rather than requiring $p_{\mu} \lt (1-\mathrm{CL})$ to exclude $\mu$, as would be done in the general framework described above, the $\mathrm{CL}_{s}$ criterion requires:



$$ \frac{p_{\mu}}{1-p_{b}} \lt (1-\mathrm{CL}) $$

Where $p_{\mu}$ is the usual probability of observing the observed value of the test statistic under the signal + background model with signal strength $\mu$, and $p_{b}$ is the p-value for the background-only hypothesis, with the p-value defined [using the opposite tail](#pvalexplainer) from the definition of $p_{\mu}$.

Using the $\mathrm{CL}_{s}$ criterion fixes the issue of setting limits much stricter than the experimental sensitivity, because for values of $\mu$ to which the experiment is not sensitive the distribution of the test statistic under the signal hypothesis is nearly the same as under the background hypothesis. Therefore, given the use of opposite tails in the p-value definition, $p_{\mu} \approx 1-p_{b}$, and the ratio approaches 1.

Note that this means that a limit set using the $\mathrm{CL}_{s}$ criterion at a given $\mathrm{CL}$ will exclude the true parameter value $\mu$ with a frequency less than the nominal rate of $1-\mathrm{CL}$.
The actual frequency at which it is excluded depends on the sensitivity of the experiment to that parameter value.

### Goodness of fit tests using the likelihood ratio

The likelihood ratio can also be used [as a measure of goodness of fit](https://www.physics.ucla.edu/~cousins/stats/cousins_saturated.pdf), i.e. describing how well the data match the model for binned data.

A standard likelihood-based measure of the goodness of fit is determined by using the log likelihood ratio with the likelihood in the denominator coming from the **saturated model**.

$$ t_{\mathrm{saturated}} \propto -\log(\frac{\mathcal{L}_{\mathcal{M}}}{\mathcal{L}_{\mathcal{M}_\mathrm{saturated}}}) $$

Here $\mathcal{M}$ is whatever model one is testing the goodness of fit for, and the saturated model is a model for which the prediction matches the observed value in every bin.
Typically, the saturated model would be one in which there are as many free parameters as bins.

This ratio is then providing a comparison between how well the actual data are fit as compared to a hypothetical optimal fit.

Unfortunately, the distribution of $t_{\mathcal{saturated}}$ often is not known a priori and must typically be estimated by generating pseudodata from the model $\mathcal{L}$ and calculating the empirical distribution of the statistic.

Once the distribution is determined a p-value for the statistic can be derived which indicates the probability of observing data with that quality of fit given the model, and therefore serves as a measure of the goodness of fit.

### Channel Compatibility test using the likelihood ratio

When performing an anlysis across many different channels (for example, different Higgs decay modes), it is often interesting to check the level of compatibility of the various channels.

Combine implements a [channel compatibility test](../../part3/commonstatsmethods/#channel-compatibility), by considering the a model, $\mathcal{M}_{\mathrm{c-independent}}$, in which the signal is independent in every channel.
As a test statistic, this test uses the likelihood ratio between the best fit value of the nominal model and the model with independent signal strength for each channel:

$$ t = -\log(\frac{\mathcal{L}_{\mathcal{M}}(\vec{\hat{\mu}},\vec{\hat{\nu}})}{\mathcal{L}_{\mathcal{M}_{\mathrm{c-indep}}}(\vec{\hat{\mu}}_{c1}, \vec{\hat{\mu}}_{c2}, ..., \vec{\hat{\nu}})}) $$

The distribution of the test statistic is not known a priori, and needs to be calculated by generating pseudo-data samples.

## Other Statistical Tests

While combine is a likelihood based statistical framework, it does not require that all statistical tests use the likelihood ratio.

### Other Goodness of Fit Tests

As well as the saturated goodness of fit test, defined above, combine implements [Kolmogorov-Smirnov](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) and [Anderson-Darling](https://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test) goodness of fit tests.

For the Kolomogorov-Smirnov (KS) test, the test statistic is the maximum absolute difference between the cumulative distribution function between the data and the model:

$$ D = \max_{x} | F_{\mathcal{M}}(x) - F_{\mathrm{data}}(x) | $$

Where $F(x)$ is the Cumulative Distribution Function (i.e. cumulative sum) of the model or data at point $\vec{x}$.

For the Anderson-Darling (AD) test, the test statistic is the integral of the square of difference between the two cumulative distribution functions:

$$ A^2 = \int_{x_{\mathrm{min}}}^{x_{\mathrm{max}}} \frac{ (F_{\mathcal{M}}(x) - F_{\mathrm{data}}(x))^2}{ F_\mathcal{M}(x) (1 - F_{\mathcal{M}}(x)) } \mathrm{d}F_\mathcal{M}(x) $$

Notably, both the Anderson-Darling and Kolmogorov-Smirnov test rely on the cumulative distribution.
Because the ordering of different channels of a model is not well defined, the tests themselves are not well defined over multiple channels.

