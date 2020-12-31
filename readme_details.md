# A Bayesian State-Space Approach to Mapping Directional Brain Networks

## Description

We develop these functions to implement our proposed Bayesian method with a SBM-motivated prior (BSBM) in our paper *A Bayesian State-Space Approach to Mapping Directional Brain Networks*, submitted
to *Journal of the American Statistical Association*. We also include a control panel (BSBM_ContPnl.m) and the simulated datasets we used in our paper to help readers understand the usages of the functions, reproduce the results of the simulation study, and analyze their data.

The codes are also avaialbe from https://github.com/StatDeptZhang/A-Bayesian-State-Space-Approach-to-Mapping-Directional-Brain-Networks 

## Installation

Please download all the functions and add them to search path before using any of them. We also provide a function addContainingDirAndSubDir.m in our toolbox to do so. You may want to put our toolbox in your current working directory and execute addContainingDirAndSubDir() before using any other functions.

The Matlab codes are written under Matlab version 9.5.0.944444 (R2018b). The R codes are only for generating figures and are written under R version 3.5.1 (2018-07-02). The following R packages are needed:

- R.matlab: [Henrik Bengtsson (2016). R.matlab: Read and Write MAT Files and Call MATLAB from Within R.](https://CRAN.R-project.org/package=R.matlab);
- gplots: [Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas
  Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz and Bill Venables (2019). gplots:
  Various R Programming Tools for Plotting Data.](https://CRAN.R-project.org/package=gplots);
- ggplot2: [Hadley Wickham (2016). ggplot2: Elegant Graphics for Data Analysis.](https://ggplot2.tidyverse.org);
- network: ‪[Carter T. Butts (2015). network: Classes for Relational Data.](http://CRAN.R-project.org/package=network).

## Usage

Detailed guidelines on using these functions can be found in BSBM\_manual.pdf.

- BSBM_ContPnl.m: The control panel of BSBM. It shows how to analyze the simulated data step by step and serves as an example of using our toolbox. It also includes the calculation of the clustering probabilities and the network edge probabilities. The analyses based on the clustering probabilities and the network edge probabilities, including generating figures, are described in Section 2.4, Section 3, and Section 4 of our paper in detail. The running time of 1000 MCMC iterations for data at 50 channels and 1000 time points is about 40 minutes;

- BSBM_EM.m: The primary function of the EM algorithm, which is used to set starting values and hyperparameter for the MCMC algorithm. Other functions in the subfolder *EM* are the updating steps in the EM algorithm; 
- BSBM_MCMC.m: The primary function of the MCMC algorithm, which is developed for the posterior inference. Other functions in the subfolder *MCMC* are the simulation steps in the MCMC algorithm.

- BSBM\_Network\_Example.R: An example of generating network plots from the outputs of BSBM_ConPnl.m. This file shows how we generated Figure 2(c) in our paper. The values of the thresholds in this file are calculated following the procedure described in Section 2.4. The functions called in BSBM\_Network\_Example.R are written in a separate file BSBM\_Network\_Source.R for conciseness. 


## Data
iEEG data analyzed in the paper were collected from an epileptic patient who received epilepsy diagnosis and treatment at the University of Virginia (UVA) Hospital. The authors analyzed the de-identified iEEG data of the patient. The de-identified data from the UVA hospital are stored on a password protected server, which is managed by UVA Health System IT. The authors are not allowed to share the data with a third party. 

The UVA hospital owns all the iEEG data. The third party can request the data directly from the UVA Health System. For questions regarding the data, please contact the co-author Mark Quigg via email at MSQ6G@hscmail.mcc.virginia.edu .

This toolbox contains the simulated datasets we used in our paper.

- Simulation1.mat: The simulated data presented in Section 3.1, with 1000 time points;
- Simulation2_1.mat: The simulated data presented in Section 3.2, with 2714 time points;
- Simulation2_2.mat: The simulated data presented in Section 3.2, with 1000 time points.

Each MAT file contains two variables:

- Y: The simulated data, which is a real d\*T matrix. d is the number of channels and T is the number of time points. For example, in Simulation1.mat, Y is a 50\*1000 matrix;
- A\_true: The true network structure, which is a binary d\*d matrix. "A\_true(i,j) = 1" indicates that there is a directional connection from channel j to channel i. For example, in Simulation1.mat, A\_true is a 50\*50 matrix.

The models used to generate these simulated datasets are described in Section 3.1 and Section 3.2 of our paper in detail.


## Reference
Huazhang Li, Yaotian Wang, Guofen Yan, Yinge Sun, Seiji Tanabe, Chang-Chia Liu, Mark S. Quigg, Tingting Zhang (2020).  A Bayesian State-Space Approach to Mapping
Directional Brain Networks, *Journal of the American Statistical Association*. Under Review. 
