# Likelihood based fitting

"Fitting" simply means estimating some parameters of a model (or really a [set of models](../../what_combine_does/model_and_likelihood/#sets-of-observation-models)) based on data.
Likelihood-based fitting does this through the [likelihood function](../../what_combine_does/model_and_likelihood/#the-likelihood).

In frequentist frameworks, this typically means doing [maximum likelihood estimation](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.2).
In bayesian frameworks, usually [posterior distributions](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.6) of the parameters are calculated from the likelihood.

## Fitting Frameworks

Likelihood fits typically either follow a frequentist framework of maximum likelihood estimation, or a bayesian framework of updating estimates to find posterior distributions given the data.

### Maximum Likelihood fits

A [maximum likelihood fit](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.2) means finding the values of the model parameters $(\vec{\mu}, \vec{\theta})$ which maximize the likelihood, $\mathcal{L}(\vec{\mu},\vec{\theta}|\mathrm{data})$
The values which maximize the likelihood, are the parameter estimates, denoted with a "hat" ($\hat{}$):

$$(\vec{\hat{\mu}}, \vec{\hat{\theta}}) \equiv \underset{\vec{\mu},\vec{\theta}}{\operatorname{argmax}} \mathcal{L}(\vec{\mu}, \vec{\theta}|\mathrm{data})$$ 

These values provide **point estimates** for the parameter values.

Because the likelihood is equal to the probability of observing the data given the model, the maximum likelihood estimate finds the parameter values for which the data is most probable.

### Bayesian Posterior Calculation

In a bayesian framework, the likelihood represents the probability of observing the data given the model and some prior probability distribution over the model parameters.

Beliefs about the values of the parameters are updated based on the data to provide a [posterior distributions](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsection.40.2.6)

$$ p(\vec{\theta}|\mathrm{data}) = \frac{ p(\mathrm{data}|\vec{\theta}) p(\vec{\theta}) }{ p(\mathrm{data}) } = \frac{ \mathcal{L}_{\mathrm{data}}(\vec{\theta}|\mathrm{data}) \mathcal{L}_{\mathrm{constraint}}(\vec{\theta}) }{ \int_{\vec{\theta'}} \mathcal{L}_{\mathrm{data}}(\vec{\theta'}|\mathrm{data}) \mathcal{L}_{\mathrm{constraint}}(\vec{\theta'}) }$$ 

The posterior distribution p$(\vec{\theta}|\mathrm{data})$ defines the updated belief about the parameters $\vec{\theta}$.

## Methods for considering subsets of models

Often, one is interested in some particular aspect of a model. 
This may be for example information related to the parameters of interest, but not the nuisance parameters.
In this case, one needs a method for specifying precisely what is meant by a model considering only those parameters of interest.

There are several methods for considering sub models which each have their own interpretations and use cases.

### Conditioning 

Conditional Sub-models can be made by simply restricting the values of some parameters. 
The conditional likelihood of the parameters $\vec{x}$ conditioned on particular values of the parameters $\vec{y}$ is:

$$ \mathcal{L}(\vec{x},\vec{y}) \xrightarrow{\mathrm{conditioned\ on\ } \vec{y} = \vec{y}_0} \mathcal{L}(\vec{x}) = \mathcal{L}(\vec{x},\vec{y}_0) $$

### Profiling

The profiled likelihood $\mathcal{L}(\vec{x})$ is defined from the full likelihood, $\mathcal{L}(\vec{x},\vec{y})$, such that for every point $\vec{x}$ it is equal to the full likelihood at $\vec{x}$ maximized over $\vec{y}$.

$$ \mathcal{L}(\vec{x},\vec{y}) \xrightarrow{\mathrm{profiling\ } \vec{y}} \mathcal{L}({\vec{x}}) = \max_{\vec{y}} \mathcal{L}(\vec{x},\vec{y})$$

In some sense, the profiled likelihood is the best estimate of the likelihood at every point $\vec{x}$, it is sometimes also denoted $\mathcal{L}(\vec{x},\vec{\hat{y}}(\vec{x}))$.

### Marginalization 

Marginalization is a procedure for producing a probability distribution $p(\vec{x}|\mathrm{data})$ for a set of parameters $\vec{x}$, which are only a subset of the parameters in the full distribution $p(\vec{x},\vec{y}|\mathrm{data})$.
The marginal probability density $p(\vec{x})$ is defined such that for every point $\vec{x}$ it is equal to the probability at $\vec{x}$ integrated over $\vec{y}$.

$$ p(\vec{x},\vec{y}) \xrightarrow{\mathrm{marginalizing\ } \vec{y}} p({\vec{x}}) = \int_{\vec{y}} p(\vec{x},\vec{y})$$

The marginalized probability $p(\vec{x})$ is the probability for the parameter values $\vec{x}$ taking into account all possible values of $\vec{y}$.

Marginalized likelihoods can also be defined, by their relationship to the probability distributions:

$$ \mathcal{L}(\vec{x}|\mathrm{data}) = p(\mathrm{data}|\vec{x}) $$


## Parameter Uncertainties 

Parameter uncertainties describe regions of parameter values which are considered reasonable parameter values, rather than single estimates.
These can be defined either in terms of frequentist **confidence regions** or bayesian **credibility regions**.

In both cases the region is defined by a confidence or credibility level $CL$, which quantifies the meaning of the region.
For frequentist confidence regions, the confidence level $CL$ describes how often the confidence region will contain the true parameter values.
For bayesian credibility regions, the credibility level $CL$ describes the bayesian probability that the true parameter value is in that region.


The confidence or credibility regions are described by a set of points $\{ \vec{\theta} \}_{\mathrm{CL}}$ which meet some criteria.
In most situations of interest, the credibility region or confidence region for a single parameter, $\theta$, is effectively described by an interval:

$$ \{ \theta \}_{\mathrm{CL}} = [ \theta^{-}_{\mathrm{CL}}, \theta^{+}_{\mathrm{CL}} ] $$

Typically indicated as:

$$ \theta = X^{+\mathrm{up}}_{-\mathrm{down}} $$


### Frequentist Confidence Regions

Frequentist confidence regions are random variables of the observed data. 
These are very often the construction used to define the uncertainties reported on a parameter.

If the same experiment is repeated multiple times, different data will be osbserved each time and a different confidence set $\{ \vec{\theta}\}_{\mathrm{CL}}^{i}$ will be found for each experiment.
The confidence regions, If the data are generated by the model with some set of values $\vec{\theta}_{\mathrm{gen}}$, then the fraction of the intervals $\{ \vec{\theta}\}_{\mathrm{CL}}^{i}$ which contain the values $\vec{\theta}_{\mathrm{gen}}$ will be equal to the confidence level ${\mathrm{CL}}$. 
The fraction of intervals which contain the generating parameter value is referred to as the "coverage".

From first principles, the intervals can be constructed using the [Neyman construction](https://pdg.lbl.gov/2022/web/viewer.html?file=../reviews/rpp2022-rev-statistics.pdf#subsubsection.40.4.2.1).

In practice, the likelihood can be used to construct confidence regions for a set of parameters $\vec{\mu}$ by using the profile likelikhood ratio:

$$ \Lambda \equiv \frac{\mathcal{L}(\vec{\mu},\vec{\hat{\theta}}(\vec{\mu}))}{\mathcal{L}(\vec{\hat{\mu}},\vec{\hat{\theta}})} $$

i.e. the ratio of the profile likelihood at point $\vec{\mu}$ to the maxmimum likelihood. For technical reasons, the negative logarithm of this quantity is typically used in practice.

Each point $\vec{\theta}$ can be tested to see if it is in the confidence region, by checking the value of the likelihood ratio at that point and comparing it to the expected distribution if that point were the true generating value of the data.

$$ \{ \vec{\theta} \}_{\mathrm{CL}} = \{ \vec{\theta} : -\log(\Lambda) \lt  \gamma_{\mathrm{CL}}(\vec{\theta}) \} $$

The cutoff value $\gamma_{\mathrm{CL}}$ must be chosen to match this desired coverage of the confidence set. 

Under some conditions, the value of $\gamma_{\mathrm{CL}}$ is known analytically for any desired confidence level, and is independent of $\vec{\theta}$, which greatly simplifies estimating confidence regions.

/// details | **Constructing Frequentist Confidence Regions in Practice**

When a single fit is performed by some numerical minimization program and parameter values are reported along with some uncertainty values, they are usually reported as frequentist intervals.
The [MINUIT minimizer](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf) which evaluates likelihood functions has two methods for [estimating parameter uncertainties](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf#section.2.5).

These two methods are the most commonly used methods for estimating confidence regions in a fit; they are the [**minos** method](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf#subsection.2.5.3), and the [**hessian** method](https://root.cern/root/htmldoc/guides/minuit2/Minuit2.pdf#subsection.2.5.2).
In both cases, Wilk's theorem is assumed to hold at all points in parameter space, such the $\gamma_{\mathrm{CL}}$ is independent of $\vec{\theta}$.

When $\gamma_{\mathrm{CL}}$ is independent of $\vec{\theta}$ the problem simplifies to finding the boundaries where $-\log(\Lambda) = \gamma_{\mathrm{CL}}$.
This boundary point is referred to as the "crossing", i.e. where $-\log(\Lambda)$ crosses the threshold value.

#### The Minos method for estimating confidence regions

In the minos method, once the best fit point $\vec{\hat{\theta}}$ is determined, the confidence region for any parameter $\theta_i$ can be found by moving away from its best fit value $\hat{\theta}_i$.
At each value of $\theta_i$, the other parameters are profiled, and $-\log{\Lambda}$ is calculated.

Following this procedure, $\theta_i$ is searched for the boundary of the confidence regions, where $-\log{\Lambda} = \gamma_{\mathrm{CL}}$.

The search is performed in both directions, away from the best fit value of the parameter and the two crossings are taken as the borders of the confidence region.

This procedure has to be followed sepately for each parameter $\theta_i$ for which a confidence interval is calculated.

#### The Hessian method for estimating confidence regions

The Hessian method relies on the second derivatives (i.e. the [hessian](https://en.wikipedia.org/wiki/Hessian_matrix)) of the likelihood at the best fit point.

By assuming that the shape of the likelihood function is well described by its second-order approximation, the values at which $-\log(\Lambda) = \gamma_{\mathrm{CL}}$ can be calculated analytically without the need for a seach

$$ \theta_i^{\mathrm{crossing}} - \hat{\theta} \propto (\frac{\partial^2{\mathcal{L}(\vec{\hat{\theta}})}}{\partial\theta_i^2})^{-1} $$

By computing and then inverting the full hessian matrix, all individual confidence regions and the full covariance matrix are determined.
By construction, this method always reports symmetric confidence intervals, as it assumes that the likelihood is well described by a second order expansion.

///

### Bayesian Credibility Regions

Often the full posterior probability distribution is summarized in terms of some **credible region** which contains some specified portion of the posterior probability of the parameter.

$$ \{ \vec{\theta} \}_{\mathrm{CL}} =  \{ \vec{\theta} : \vec{\theta} \in \Omega, \int_{\Omega} p(\vec{\theta}|\mathrm{data}) = \mathrm{CL}  \}$$

The credible region represents a region in which the bayesian probability of the parameter being in that region is equal to the chosen Credibility Level.

