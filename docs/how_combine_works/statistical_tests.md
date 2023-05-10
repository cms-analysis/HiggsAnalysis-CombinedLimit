# Statistical Tests


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

This is a three step process.

1. Determine the expected distribution of the test statistic, given the model. 
2. Check the value of the test statistic for the actual observed data.
3. Determine the p-value for observing a test statistic at least as extreme as the actual observed value .

If the p-value is sufficiently small, the model is considered rejected.

<details>
<summary><b> Show Mathematical details of the general statistical test </b></summary>

#### Mathematical Details of the general statistical test

The distribution of the test statistic, $t$ under some model hypothesis $\mathcal{M}$ is:

$$t \stackrel{\mathcal{M}}{\sim} D_{\mathcal{M}}$$

And the observed value of the test statistic is $t_{\mathrm{obs}}$.
The p-value can the be calculated as:

$$p(t_{\mathrm{obs}}) = \int_{\Omega(t_{\mathrm{obs}})} D_{\mathcal{M}}$$

where $\Omega(t_{\mathrm{obs}})$ is a region in the space of the test statistic bounded by $t_{\mathrm{obs}$.
Its precise definition depends on the test under consideration but in general it is either:

$$
[t_{\mathrm{min}}, t_{\mathrm{obs}} ] \mathrm{\ or\ } [t_{\mathrm{obs}}, t_{\mathrm{max}} ] 
$$

Where $t_{\mathrm{min}}$ and  $t_{\mathrm{max}}$ are the low and upper bounds of the domain of the test statistic.

</details>

#### Alternative Models

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


## Likelihood related test statistics

The [likelihood](model_and_likelihood.md#the-likelihood) itself often forms a good basis for test building test statistics.

Typically the absolute value of the likelihood itself is not very meaningful.
However, quantities such as the ratio of the likelihood at two different points in parameter space are very informative about the relative merits of those two models.

Furthermore, in some cases analytic formula approximating the distribution of likelihood based test statistics are known.
This is an extremely helpful property, because when the expected distribution of the test statistic is known, only the observed value needs to be calculated and then it can be compared to the expected distribution.

On the other hand, in the general case where the distriubtion of the test statistic is not known analytically, one has to generate pseudodata from the model to generate an expected distribution.
This process is subject to limitations of how many pseudodata samples are produced, and can be computationally very expensive.


### The likelihood ratio and likelihood ratio based test statistics

The likelihood ratio of two models: 

$$ \Lambda \equiv \frac{\mathcal{L}_{\mathcal{M}}}{\mathcal{L}_{\mathcal{M}'}} $$ 

is a very useful ingredient in statistics.

The absolute value of the likelihood depends on many things, such as the size of the parameter space being considered and of the possible sets of observations.
However, if the models $\mathcal{M}$ and $\mathcal{M}'$ have the same or appropriately similar parameter spaces, and describe the same dataset, their ratio is very informative.
It describes how much more likely the data is under one model than the other.

For technical and convenience reasons, often the logarithm of the likelihood ratio: 

$$\log(\Lambda) = \log(\mathcal{L}_{\mathcal{M}}) - \log(\mathcal{L}_{\mathcal{M}'})$$

or some multiple of it, such as $-2\log(\Lambda)$ is used.

Typically, both models under are "composite" models, or ["sets" of models](model_and_likelihood.md#sets-of-observation-models), defined by some set of parameters $(\vec{\mu}, \vec{\theta})$.
In these cases, we are typically interested in the likelihood ratio for the models in each set which best fit the data:

$$ \mathcal{L}_{\mathcal{M}}(\vec{\hat{\mu}},\vec{\hat{\theta}})\mathrm{,\ and\ }\mathcal{L}_{\mathcal{M}'}(\vec{\hat{\mu}},\vec{\hat{\theta}})$$

### Parameter Estimation using the likelihood ratio

A common use case for likelihood ratios is estimating the values of the parameters of interest, $\vec{\mu}$.
The point estimate for the parameters is simply the maximum likelihood estimate, but the likelihood ratio can be used for estimating the uncertainty as a confidence set.

A [confidence set](fitting_concepts.md#uncertainties-in-maxmimum-likelihood-fits-confidence-sets) for the parameter $\vec{\mu}$ can be defined by comparing models with $\vec{\mu}$ floating, to one where they are fixed to their best fit values $\vec{\mu} = \vec{\hat{\mu}}$.

The confidence set is determined as the region for which $-\log(\Lambda) \gt c$, where $c$ is some constant.

$$ \{ \vec{\mu}_{CL} \} =  \{ \vec{\mu} : -\log(\frac{\mathcal{L}(\vec{\mu},\vec{\hat{\theta}})}{\mathcal{L}(\vec{\hat{\mu}},\vec{\hat{\theta}})}) \gt c_\mathrm{CL} \}.$$

In other words, the confidence set is defined by chosing as the test statistic, $t$, the log likelihood ratio between the [profiled likelihood](fitting_concepts.md#profiling) $\mathcal{L}(\vec{\mu})$ and the likelihood at the best fit point.

In order to determine the correct value of $c_{\mathrm{CL}}$ for a chosen confidence level, one must know the expected distribution of $t$.

Under appropriate conditions, the distribution of $t$ is known analytically, via [Wilks' Theorem](https://en.wikipedia.org/wiki/Wilks'_theorem).
In the general case, the

### Discoveries using the likelihood ratio

A common method for claiming discovery is based on a likelihood ratio test.
Essentially, by showing that the best parameter estimate for some new physics parameter excludes the standard model value of that parameter.

Typicaly, in order to claim discovery, a strong criterion is used for the appropriate confidence level.
In particle physics, we usually require that $\mu_{SM} \not \in \{ \mu_{0.9999995} \}$.

Practically speaking, this can be done by using a test statistic:

$$ t \propto -\log(\frac{\mathcal{L}_{\mathrm{SM}}}{\mathcal{L}_{\mathrm{SM+NP}}}) $$

However, this also allows for claiming "discovery" in cases where $\vec{\mu} \lt 0$, which in particle physics is often an unphysical model.
In order to avoid such a situation, we typically use a modified test statistic:

$$ q_{0} \propto \begin{cases}
    0 & \hat{\mu} \lt 0 \\ 
    -\log(\frac{\mathcal{L}_{\mathrm{SM}}}{\mathcal{L}_{\mathrm{SM+NP}}}) & \hat{\mu} \geq 0 
\end{cases}
$$


### Limit Setting using the likelihood ratio

Various test statistics built from likelihood ratios can be used for limit setting.

The simplest thing to do to set a limit on a parameter $\mu$ is to find the values of $\mu$ such that $\mu \not \in \{ \mu_{CL} \}$.

This can be done using the likelihood ratio as defined above. 

$$ t \propto -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})}) $$

However, these excluded values $\mu$ may include values closer to the standard model than the best fit value of $\mu$,  $\mu_{\mathrm{SM}} \lt \mu \lt \hat{\mu}$.
Generally speaking, we prefer not to exclude these values.

Therefore, we typically use a modified test statistic:

$$ t' \propto -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\min(\mu,\hat{\mu}))}) = 
\begin{cases} 
    -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})})& \hat{\mu} \lt \mu  \\
    -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\mu)})&  \mu \lt \hat{\mu}  
\end{cases}
$$

However, this can also have undesirable properties when the best fit value $\hat{\mu}$ is less than 0.
In that case, we may set limits below 0.
In order to avoid these situations, another modified test statistic can account for this:


$$ 
\tilde{q}_{\mu} \propto \begin{cases} 
    0& \hat{\mu} \lt 0  \\
    -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\hat{\mu})})& 0 \lt \hat{\mu} \lt \mu  \\
    -\log(\frac{\mathcal{L}(\mu)}{\mathcal{L}(\mu)})& 0 \lt \mu \lt \hat{\mu}  
\end{cases}
$$






### Goodness of fit tests using the likelihood ratio




