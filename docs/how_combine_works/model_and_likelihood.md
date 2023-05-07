# Models and Likelihoods

## The Model

The Model $\mathcal{M}(\vec{mu},\vec{theta})$ defines the expected observations given the input parameters of interest $\vec{\mu}$ and nuisance parameters $\vec{\theta}$.

Combine is designed for essentially counting experiments, where the number of events with particular features are counted.  
Three basic types of model can be provided to combine:

1. A simple counting experiment, in which only the expected number of observed events is specified and how this number depends on a set of parameters.
2. A template analysis, in which observed events are divided into various regions (e.g. histogram bins defining regions of invariant mass of two photons in an event), and the expected observation for each bin is defined
3. A parametric analysis, in which continuous functions define the probability density of observing events with a given set of characterstics (e.g. where the highest energy photon has a given pseudorapidity and transvere momentum).

These last two types of analysis are often collectively referred to as "shape-based" analyses. 

Furthermore, models can be composed of multiple channels and processes.

**Channels**: The model is a set of models for multiple channels. $\mathcal{M} = \{ \mathcal{m}_{c1}, \mathcal{m}_{c2}, .... \mathcal{m}_{cN}\}a$. 
    The expected observations which define the model is then the union of the sets of expected observations in each individual channel.

**Processes**: The model consists of a sum over different processes, each processes having its own model defining it. $\mathcal{m}_c =  \Sum_p \mathcal{m}_{c,p}$
    the expected observations is then the sum of the expected observations for each of the processes.

Combining full models together is then possible by combining their channels together, assuming that the channels are each independent.

Generally we are not interested in one specific model with perfectly known parameters. We are typically interested in a set of models, parameterized by a set of real numbers representing possible versions of the model.

**Model Parameters**: A number of parameters can modify the model expectation. These include the parameters of interest ( $\vec{\mu}$, those being measured) as well as nuisance parameters ($\vec{\theta}$), which may not be of interest but still affect the model expectation.

The full model therefore defines the expected observations over all the channels, given all the processes and model parameters:

$$\mathcal{M} = \{ \mathcal{m}_{c1}(\vec{\mu},\vec{\theta}), \mathcal{m}_{c2}(\vec{\mu},\vec{theta}), ... }$$

Combine provides tools and interfaces for defining the model as arbitrary functions of the input parameters.
In practice, however, there are a number of most commonly used functional forms which define how the expected events depend on the model parameters.
These are discussed in detail in the context of the full likelihood below.

## The Likelihood 

Combine builds a likelihood model through the `text2workspace.py` script, given the model which is defined by the `Datacard`.
The `Datacard` defines the observables being considered, the observations, the processes which contribute, the uncertainties, and various relationships between each of these.
The `text2workspace.py` script takes this datacard and any other input files it references (for example root files which may contain histograms) and builds a `RooWorkspace` which defines the likelihood model.

The likelihood takes the general form:

$ \mathcal{L} =  \mathcal{L}_\mathrm{data} \cdot \mathcal{L}_\mathrm{constraint}$

Where $ \mathcal{L}_\mathrm{data} $ is the likelihood of observing the data given your model, and $\mathcal{L}_\mathrm{constraint}$ represent some external constraints.
The data term determine how probable any set of data is given the model for particular values of the model parameters. 
The constraint term does not depend on the observed data, but rather encodes other constraints which may be constraints from previous measurements (such as Jet Energy Scales) or prior beliefs about the value some parameter in the model should have. i.e. they encode the a priori likelihood of the given parameters.

In general each is, in turn, composed of many sublikelihoods. 
For example in a binned shape analysis every bin has its own term representing the probability of observing the actual observed number of events given the expected number of events.
Similarly, for every systematic uncertainty there will be some constraint term associated with the uncertainty.

This model is extremely general. Without yet specifying any function forms for these likelihoods, this model can represent anything at all.
However, as with the model itself, there are typical forms that the likelihoods take which will cover most use cases, and for which combine is primarily designed.

Combine is designed for constructing data likelihoods which refer to observation of some count of events, though they can be binned or unbinned.
The constraint terms also often take on common forms, such as that of a gaussian distribution.

### Binned Likelihood Models

For a binned likelihood, the probability of observing a certain number of counts, given a model takes on a simple form. For each bin:

$ \mathcal{L}_\mathrm{bin}(\vec{\mu},\vec{\theta};\mathrm{data}) = \mathrm{Poiss}(\mathrm{obs}; n_\mathrm{exp})$.

i.e. it is a poisson distribution with the mean given by the expected number of events in that bin. 
The full data likelihood is simply the product of each of the bins' likelihoods:

$ \mathcal{L}_\mathrm{data} = \Prod_\mathrm{bins} \mathcal{L}_\mathrm{bin}.

This is the underlying likelihood model used for every binned analysis. 
The freedom in the analysis comes in how $n_\mathrm{exp}$ depends on the model parameters, and the constraints that are place on those parameters.

### Unbinned Likelihood Models

For unbinned likelihood models, a likelihood can be given to each data point. It is just equal to probability density function at that point, $\vec{x}$.

$\mathcal{L}_\mathrm{data} = \Prod_{i} \mathrm{pdf}(\vec{x}_i | \vec{\mu}, \vec{\theta} )$

### Constraint Likelihoods

In principle, the constraint terms in the likelihood can be arbitrary. 
They simply encode the probability of model parameters taking on a certain value.
This might represent the probability that a cross section has a certain value, or that the true luminosity is shifted by some amount as compared to the nominal value.
In practice, most constraint terms are gaussian, of the form $\mathcal{N}(0,1)$ where the model parameter, often denoted $\theta$, has been defined such that 0 is its nominal value and $\pm1$ represent its value shifted to the edges of the 68\% confidence range; i.e. these two values define the bounds of the region in which you expect to find the parameter 68\% of the time.

It is not *necessary* that these constraint terms take on such a gaussian form. They can be arbitrary functions.
However, it is very convenient to use gaussian constraints for a number of reasons\*. 

Furthermore, the form of the constraint terms, on their own is not in and of itself meaningful; what is meaningful is the relationship between the constraint term and how the model expectation is altered by the parameter.
Any co-ordinate transformation of the parameter values can be absorbed into the definition of the parameter; i.e. the likelihood model is independent of how it is parameterized.
A reparameterization would change the mathematical form of the constraint term, but would also simultaneously change how the model depends on the parameter in such a way that the total likelihood is unchanged.
e.g. if you define  $\theta = \sigma(tt)$ or $\theta = \sigma(tt) - \sigma_0$ you will change the form of the constraint term, but the you will not change the overall likelihood.


### How the model is constructed from the inputs

Most of the magic and variability of the different likelihoods is encoded in what I have refered to as the "model"; i.e. the relationship between $n_\mathrm{exp}$ and the model parameters.



**Modifiers**: A number of parameters can modify the model expectation. These include the parameters of interest (those being measured) as well as nuisance parameters, which may not be of interest but still affect the model expectation.

There is a lot of flexibility in defining custom modifiers both for parameters of interest and nuisance parameters, such that they may alter the model in almost arbitrary ways.
However, there are typical modifiers which cover most use cases.

**Parameter of interest**: By default, the parameter of interest simply scales the overall expectation from any signal processes. $m_{c,\mathrm{sig}(\mu) = \mu * m_{c,\mathrm{sig}}^0

**Multiplicative Modifiers**: Generally speaking, there are many modifiers which modify the model in a multiplicative fashion. $m_{c,p}(\theta) = f(\theta) m_{c,p}^0$
    note that this is still quite general, as $f$ can take on arbitrary function forms.
    However, the most common forms are:
        - $f = \theta$
        - $f = \kappa^{\theta} $

**Additive Modifiers**: We also often use modifiers which modify the model in an additive way. $m_{c,p}(\theta) = m_{c,p}^0 + f(\theta)$
    where again, $f$ can take o arbitrary functional forms, leaving this quite general.
    The most common forms are:
        

The resolution order of these parameters is important, since $ f\cdot m + g \neq f\cdot(m+g)$. 




