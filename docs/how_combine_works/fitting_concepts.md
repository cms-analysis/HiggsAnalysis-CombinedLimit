# Likelihood based fitting

"Fitting" simply means estimating some parameters of a model (or really a [set of models](how_combine_works/model_and_likelihood.md#sets-of-observation-models)) based on data.
Likelihood-based fitting does this through the [likelihood function](how_combine_works/model_and_likelihood.md#the-likelihood).

In frequentist frameworks, this typically means doing maximum likelihood estimation.
In bayesian frameworks, usually [posterior distributios](https://en.wikipedia.org/wiki/Posterior_probability) of the parameters are calculated from the likelihood.

## Fitting Frameworks

Likelihood fits typically either follow a frequentist framework of maximum likelihood estimation, or a bayesian framework of updating estimates to find posterior distributions given the data.

### Maximum Likelihood fits

A [maximum likelihood fit](https://en.wikipedia.org/wiki/Maximum_likelihood_estimation) means finding the values of the model parameters $(\vec{\mu}, \vec{\theta})$ which maximize the likelihood, $\mathcal{L}(\vec{\mu},\vec{\theta})$
The values which maximize the likelihood, 

$$(\vec{\hat{\mu}}, \vec{\hat{\theta}}) \equiv \mathrm{argmax\ } \mathcal{L}(\vec{\mu}, \vec{\theta})$$ 

are taken to be the best estimates of the parameter values.
These values provide **point estimates** for the parameter values.

Because the likelihood is equal to the probability of observing the data given the model, the maximum likelihood estimate finds the parameter values for which the data is most probable.

### Bayesian Posterior Calculation

In a bayesian framework, the likelihood represents the probability of observing the data given the model and some prior probability distribution over the model parameters.

Beliefs about the values of the parameters are updated based on the data to provide a posterior distribution:

$$ p(\vec{\theta}|\mathrm{data}) = \frac{ p(\mathrm{data}|\vec{\theta}) p(\vec{\theta}) }{ p(\mathrm{data}) } = \frac{ \mathcal{L}_{\mathrm{data}}(\vec{\theta}|\mathrm{data}) \mathcal{L}_{\mathrm{constraint}}(\vec{\theta}) }{ \int_{\vec{\theta'}} \mathcal{L}_{\mathrm{data}}(\vec{\theta'}|\mathrm{data}) \mathcal{L}_{\mathrm{constraint}}(\vec{\theta'}) }$$ 

The posterior distribution p$(\vec{\theta}|\mathrm{data})$ defines the updated belief about the parameters $\vec{\theta}$.

## Methods for Making Sub-Models

Often, one is interest in some particular aspect of a model. 
This may be for example information related to the parameters of interest, but not the nuisance parameters.
In this case, one needs a method for specifying precisely what is meant by a model considering only those parameters of interest.

There are several methods for considering sub models which each have their own interpretations and use cases.

### Conditioning 

Conditional Sub-models can be made by simply restricting the values of some parameters. 
The conditional likelihood of the parameters $\vec{x}$ conditioned on particular values of the parameters $\vec{y}$ is:

$$ \mathcal{L}(\vec{x},\vec{y}) \xrightarrow{\mathrm{conditioned\ on\ } \vec{y} = \vec{y}_0} \mathcal{L}(\vec{x}) = \mathcal{L}(\vec{x},\vec{y}_0) $$

### Profiling

Profiling is a procedure for producing a likelihood $\mathcal{L}(\vec{x})$ for a set of parameters $\vec{x}$, which are only a subset of the parameters in the full likelihood $\mathcal{L}(\vec{x},\vec{y})$.
The profiled likelihood $\mathcal{L}(\vec{x})$ is defined such that for every point $\vec{x}$ it is equal to the full likelihood at $\vec{x}$ maximized over $\vec{y}$.

$$ \mathcal{L}(\vec{x},\vec{y}) \xrightarrow{\mathrm{profiling\ } \vec{y}} \mathcal{L}({\vec{x}}) = \max_{\vec{y}} \mathcal{L}(\vec{x},\vec{y})$$

In some sense, the profiled likelihood is the best estimate of the likelihood at every point $\vec{x}$, it is sometimes also denoted $\mathcal{L}(\vec{x},\vec{\hat{y}}(\vec{x}))$.

### Marginalization 

Marginalization is a procedure for producing a probability distribution $p(\vec{x}|\mathrm{data})$ for a set of parameters $\vec{x}$, which are only a subset of the parameters in the full distribution $p(\vec{x},\vec{y}|\mathrm{data})$.
The marginal probability density $p(\vec{x})$ is defined such that for every point $\vec{x}$ it is equal to the probability at $\vec{x}$ integrated over $\vec{y}$.

$$ p(\vec{x},\vec{y}) \xrightarrow{\mathrm{marginalizing\ } \vec{y}} p({\vec{x}}) = \int_{\vec{y}} p(\vec{x},\vec{y})$$

The marginalized probability $p(\vec{x})$ is the probability for the parameter values $\vec{x}$ taking into account all possible values of the $\vec{y}$.


## Fit Uncertainties 

Fit uncertainties describe regions of parameter values which are considered reasonable parameter values, rather than single estimates.
These can be defined either in terms of frequentist **confidence regions** or bayesian **credibility regions**.

In both cases the region is defined by a confidence or credibility level $CL$, which quantifies the meaning of the region.
For frequentist confidence regions, the confidence level $CL$ describes how often the confidence region will contain the true parameter values.
For bayesian credibility regions, the credibility level $CL$ describes the bayesian probability that the true parameter value is in that region.

### Frequentist Confidence Regions

The likelihood can also be used to estimate uncertainties on the parameters. 
By determining how much the likelihood changes when the parameter values are changed, we can find the parameter values for which the data are more likely than some predetermined cutoff value, $\gamma$.

$$ \{ \vec{\theta}_{\mathrm{CL}} \} = \{ \vec{\theta} : \mathcal{L}(\vec{\theta}|\mathrm{data}) \gt  \gamma_{\mathrm{CL}} \} $$

If the same experiment is repeated multiple times, different data will be osbserved each time and a different confidence set $\{ \vec{\theta}_{\mathrm{CL}} \}^{i}$ will be found for each experiment.
If the data are generated by the model with some set of values $\vec{\theta}_{\mathrm{gen}}$, the set $\{ \vec{\theta}_{\mathrm{CL}} \}^{i}$ will contain the values $\vec{\theta}_{\mathrm{gen}}$ with some fixed frequency.
The value $\mathrm{CL}$ represents the frequency with which the confidence sets contain the generating parameters, and the cutoff value $\gamma_{\mathrm{CL}}$ must be chosen to match this desired frequency.


### Bayesian Credibility Regions

Often the full posteriror probability distribution is summarized in terms of some **credible region** which contains some specified portion of the probability.

$$ \{ \vec{\theta}_{\mathrm{CL}} \} =  \{ \vec{\theta} \in \Omega : \int_{\Omega} p(\vec{\theta}|\mathrm{data}) = \mathrm{CL}  \}$$

where the label $\Omega$ is used to indicate that there are an infinite number of such regions which can be chosen.
The credible region represents a region in which the bayesian probability of the parameter being in that region is equal to the chosen Credibility Level.

### Describing Regions with Intervals

In most situations of interest, the credibility region or confidince region for a single parameter, $\theta$, is effectively described by an interval:

$$ \{ \theta_{\mathrm{CL}} \} = [ \theta^{-}_{\mathrm{CL}}, \theta^{+}_{\mathrm{CL}} ] $$

Often this is indicated as:

$$ \theta = X^{+\mathrm{up}}_{-\mathrm{down}} $$

