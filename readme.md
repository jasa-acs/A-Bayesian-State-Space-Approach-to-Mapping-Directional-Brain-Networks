# A Bayesian State-Space Approach to Mapping Directional Brain Networks

# Author Contributions Checklist Form

## Data

### Abstract (Mandatory)

The purposes of this paper are to present a new statistical approach for analysis of intracranial
electroencephalographic (iEEG) data and to use our approach to uncover the normal and
abnormal directional brain networks of epileptic patients over the course of seizure development.
iEEG data are high-dimensional multivariate time series recordings of many small regions'
neuronal activities at a high temporal resolution (millisecond scale) and spatial resolution (about
10 mm in diameter) and with a strong signal-to-noise ratio, in contrast to popular functional
magnetic resonance imaging (fMRI) with a low temporal resolution and scalp EEG with a low
spatial resolution. As such, iEEG data provide valuable information about directional brain
networks.

### Availability (Mandatory)

iEEG data used in the paper were collected from an epileptic patient who received epilepsy
diagnosis and treatment at the University of Virginia (UVA) Hospital. The authors analyzed the
de-identified iEEG data of the patient. The de-identified data from the UVA hospital are stored
on a password protected server, which is managed by UVA Health System IT. The authors are
not allowed to share the data with a third party.

The UVA hospital owns all the iEEG data. The third party can request the data directly from the
UVA Health System. For questions regarding the data, please contact the co-author Mark Quigg
via email at MSQ6G@hscmail.mcc.virginia.edu .

## Code

### Abstract (Mandatory)

We build a state-space multivariate autoregression (SSMAR) for iEEG data to model the
underlying directional brain network system. We identify connected brain regions (i.e., mapping
the brain network) through estimating the SSMAR parameters that denote directional
connectivity. In contrast to most existing network models that were developed mainly for
observed network edges, we develop a Bayesian framework to estimate the proposed high-
dimensional model, infer directional connections, and identify clusters for the unobserved
network edges. The computing codes are to implement the proposed Bayesian method.

### Description (Mandatory)

We developed a MATLAB software package for implementing our proposed method.
The MATLAB used in our work has a Total Academic Headcount License purchased by
the University of Virginia.

https://github.com/StatDeptZhang/A-Bayesian-State-Space-Approach-to-Mapping-
Directional-Brain-Networks

We used two MATLAB toolboxes, including Curve Fitting Toolbox and Statistics and
Machine Learning Toolbox.

## Instructions for Use

### Reproducibility (Mandatory)

Figures 2-3 are to be reproduced.

