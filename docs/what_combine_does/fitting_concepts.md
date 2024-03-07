# Likelihood based fitting

"Fitting" simply means estimating some parameters of a model (or really a [set of models](../../what_combine_does/model_and_likelihood/#sets-of-observation-models)) based on data.
Likelihood-based fitting does this through the [likelihood function](../../what_combine_does/model_and_likelihood/#the-likelihood).

In frequentist frameworks, this typically means doing [maximum likelihood estimation](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.2).
In bayesian frameworks, usually [posterior distributions](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.6) of the parameters are calculated from the likelihood.

## Fitting Frameworks

Likelihood fits typically either follow a frequentist framework of maximum likelihood estimation, or a bayesian framework of updating estimates to find posterior distributions given the data.

### Maximum Likelihood fits

A [maximum likelihood fit](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.2) means finding the values of the model parameters $(\vec{\mu}, \vec{\nu})$ which maximize the likelihood, $\mathcal{L}(\vec{\mu},\vec{\nu};\mathrm{data})$
The values which maximize the likelihood, are the parameter estimates, denoted with a "hat" ($\hat{}$):

$$(\vec{\hat{\mu}}, \vec{\hat{\nu}}) \equiv \underset{\vec{\mu},\vec{\nu}}{\operatorname{argmax}} \mathcal{L}(\vec{\mu}, \vec{\nu};\mathrm{data})$$ 

These values provide **point estimates** for the parameter values.

Because the likelihood is equal to the probability of observing the data given the model, the maximum likelihood estimate finds the parameter values for which the data is most probable.

### Bayesian Posterior Calculation

In a bayesian framework, the likelihood represents the probability of observing the data given the model and some prior probability distribution over the model parameters.

Beliefs about the values of the parameters are updated based on the data to provide a [posterior distributions](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.6)

$$ p(\vec{\nu};\mathrm{data}) = \frac{ p(\mathrm{data};\vec{\nu}) p(\vec{\nu}) }{ p(\mathrm{data}) } = \frac{ \mathcal{L}_{\mathrm{data}}(\vec{\nu};\mathrm{data}) \mathcal{L}_{\mathrm{constraint}}(\vec{\nu}) }{ \int_{\vec{\nu'}} \mathcal{L}_{\mathrm{data}}(\vec{\nu'};\mathrm{data}) \mathcal{L}_{\mathrm{constraint}}(\vec{\nu'}) }$$ 

The posterior distribution p$(\vec{\nu};\mathrm{data})$ defines the updated belief about the parameters $\vec{\nu}$.

## Methods for considering subsets of models

Often, one is interested in some particular aspect of a model. 
This may be for example information related to the parameters of interest, but not the nuisance parameters.
In this case, one needs a method for specifying precisely what is meant by a model considering only those parameters of interest.

There are several methods for considering sub models which each have their own interpretations and use cases.

### Conditioning 

Conditional Sub-models can be made by simply restricting the values of some parameters. 
The conditional likelihood of the parameters $\vec{\mu}$ conditioned on particular values of the parameters $\vec{\nu}$ is:

$$ \mathcal{L}(\vec{\mu},\vec{\nu}) \xrightarrow{\mathrm{conditioned\ on\ } \vec{\nu} = \vec{\nu}_0} \mathcal{L}(\vec{\mu}) = \mathcal{L}(\vec{\mu},\vec{\nu}_0) $$

### Profiling

The profiled likelihood $\mathcal{L}(\vec{\mu})$ is defined from the full likelihood, $\mathcal{L}(\vec{\mu},\vec{\nu})$, such that for every point $\vec{\mu}$ it is equal to the full likelihood at $\vec{\mu}$ maximized over $\vec{\nu}$.

$$ \mathcal{L}(\vec{\mu},\vec{\nu}) \xrightarrow{\mathrm{profiling\ } \vec{\nu}} \mathcal{L}({\vec{\mu}}) = \max_{\vec{\nu}} \mathcal{L}(\vec{\mu},\vec{\nu})$$

In some sense, the profiled likelihood is the best estimate of the likelihood at every point $\vec{\mu}$, it is sometimes also denoted $\mathcal{L}(\vec{\mu},\vec{\hat{\nu}}(\vec{\mu}))$.

### Marginalization 

Marginalization is a procedure for producing a probability distribution $p(\vec{\mu};\mathrm{data})$ for a set of parameters $\vec{\mu}$, which are only a subset of the parameters in the full distribution $p(\vec{\mu},\vec{\nu};\mathrm{data})$.
The marginal probability density $p(\vec{\mu})$ is defined such that for every point $\vec{\mu}$ it is equal to the probability at $\vec{\mu}$ integrated over $\vec{\nu}$.

$$ p(\vec{\mu},\vec{\nu}) \xrightarrow{\mathrm{marginalizing\ } \vec{\nu}} p({\vec{\mu}}) = \int_{\vec{\nu}} p(\vec{\mu},\vec{\nu})$$

The marginalized probability $p(\vec{\mu})$ is the probability for the parameter values $\vec{\mu}$ taking into account all possible values of $\vec{\nu}$.

Marginalized likelihoods can also be defined, by their relationship to the probability distributions:

$$ \mathcal{L}(\vec{\mu};\mathrm{data}) = p(\mathrm{data};\vec{\mu}) $$


## Parameter Uncertainties 

Parameter uncertainties describe regions of parameter values which are considered reasonable parameter values, rather than single estimates.
These can be defined either in terms of frequentist **confidence regions** or bayesian **credibility regions**.

In both cases the region is defined by a confidence or credibility level $CL$, which quantifies the meaning of the region.
For frequentist confidence regions, the confidence level $CL$ describes how often the confidence region will contain the true parameter values if the model is a sufficiently accurate approximation of the truth.
For bayesian credibility regions, the credibility level $CL$ describes the bayesian probability that the true parameter value is in that region for under the given model.


The confidence or credibility regions are described by a set of points $\{ \vec{\mu} \}_{\mathrm{CL}}$ which meet some criteria.
In most situations of interest, the credibility region or confidence region for a single parameter, $\mu$, is effectively described by an interval:

$$ \{ \mu \}_{\mathrm{CL}} = [ \mu^{-}_{\mathrm{CL}}, \mu^{+}_{\mathrm{CL}} ] $$

Typically indicated as:

$$ \mu = X^{+\mathrm{up}}_{-\mathrm{down}} $$


### Frequentist Confidence Regions

Frequentist confidence regions are random variables of the observed data. 
These are very often the construction used to define the uncertainties reported on a parameter.

If the same experiment is repeated multiple times, different data will be osbserved each time and a different confidence set $\{ \vec{\mu}\}_{\mathrm{CL}}^{i}$ will be found for each experiment.
If the data are generated by the model with some set of values $\vec{\mu}_{\mathrm{gen}}$, then the fraction of the regions $\{ \vec{\mu}\}_{\mathrm{CL}}^{i}$ which contain the values $\vec{\mu}_{\mathrm{gen}}$ will be equal to the confidence level ${\mathrm{CL}}$. 
The fraction of intervals which contain the generating parameter value is referred to as the "coverage".

From first principles, the intervals can be constructed using the [Neyman construction](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsubsection.40.4.2.1).

In practice, the likelihood can be used to construct confidence regions for a set of parameters $\vec{\mu}$ by using the profile likelikhood ratio:

$$ \Lambda \equiv \frac{\mathcal{L}(\vec{\mu},\vec{\hat{\mu}}(\vec{\mu}))}{\mathcal{L}(\vec{\hat{\mu}},\vec{\hat{\mu}})} $$

i.e. the ratio of the profile likelihood at point $\vec{\mu}$ to the maxmimum likelihood. For technical reasons, the negative logarithm of this quantity is typically used in practice.

Each point $\vec{\mu}$ can be tested to see if it is in the confidence region, by checking the value of the likelihood ratio at that point and comparing it to the expected distribution if that point were the true generating value of the data.

$$ \{ \vec{\mu} \}_{\mathrm{CL}} = \{ \vec{\mu} : -\log(\Lambda) \lt  \gamma_{\mathrm{CL}}(\vec{\mu}) \} $$

The cutoff value $\gamma_{\mathrm{CL}}$ must be chosen to match this desired coverage of the confidence set. 

Under some conditions, the value of $\gamma_{\mathrm{CL}}$ is known analytically for any desired confidence level, and is independent of $\vec{\mu}$, which greatly simplifies estimating confidence regions.

/// details | **Constructing Frequentist Confidence Regions in Practice**

When a single fit is performed by some numerical minimization program and parameter values are reported along with some uncertainty values, they are usually reported as frequentist intervals.
The [MINUIT minimizer](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf) which evaluates likelihood functions has two methods for [estimating parameter uncertainties](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf#section.2.5).

These two methods are the most commonly used methods for estimating confidence regions in a fit; they are the [**minos** method](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf#subsection.2.5.3), and the [**hessian** method](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf#subsection.2.5.2).
In both cases, Wilk's theorem is assumed to hold at all points in parameter space, such the $\gamma_{\mathrm{CL}}$ is independent of $\vec{\mu}$.

When $\gamma_{\mathrm{CL}}$ is independent of $\vec{\mu}$ the problem simplifies to finding the boundaries where $-\log(\Lambda) = \gamma_{\mathrm{CL}}$.
This boundary point is referred to as the "crossing", i.e. where $-\log(\Lambda)$ crosses the threshold value.

#### The Minos method for estimating confidence regions

In the minos method, once the best fit point $\vec{\hat{\mu}}$ is determined, the confidence region for any parameter $\mu_i$ can be found by moving away from its best fit value $\hat{\mu}_i$.
At each value of $\mu_i$, the other parameters are profiled, and $-\log{\Lambda}$ is calculated.

Following this procedure, $\mu_i$ is searched for the boundary of the confidence regions, where $-\log{\Lambda} = \gamma_{\mathrm{CL}}$.

The search is performed in both directions, away from the best fit value of the parameter and the two crossings are taken as the borders of the confidence region.

This procedure has to be followed sepately for each parameter $\mu_i$ for which a confidence interval is calculated.

#### The Hessian method for estimating confidence regions

The Hessian method relies on the second derivatives (i.e. the [hessian](https://en.wikipedia.org/wiki/Hessian_matrix)) of the likelihood at the best fit point.

By assuming that the shape of the likelihood function is well described by its second-order approximation, the values at which $-\log(\Lambda) = \gamma_{\mathrm{CL}}$ can be calculated analytically without the need for a seach

$$ \mu_i^{\mathrm{crossing}} - \hat{\mu} \propto (\frac{\partial^2{\mathcal{L}(\vec{\hat{\mu}})}}{\partial\mu_i^2})^{-1} $$

By computing and then inverting the full hessian matrix, all individual confidence regions and the full covariance matrix are determined.
By construction, this method always reports symmetric confidence intervals, as it assumes that the likelihood is well described by a second order expansion.

///

### Bayesian Credibility Regions

Often the full posterior probability distribution is summarized in terms of some **credible region** which contains some specified portion of the posterior probability of the parameter.

$$ \{ \vec{\mu} \}_{\mathrm{CL}} =  \{ \vec{\mu} : \vec{\mu} \in \Omega, \int_{\Omega} p(\vec{\mu};\mathrm{data}) = \mathrm{CL}  \}$$

The credible region represents a region in which the bayesian probability of the parameter being in that region is equal to the chosen Credibility Level.

