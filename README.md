# Likelihood-based-inference
This repository contains several Julia codes associated with "Likelihood-based inference, identifiability, and prediction using count data from lattice-based random walk models" by Yihan Liu, David J. Warne, and Matthew J. Simpson. The preprint is available at https://arxiv.org/abs/2406.16296.

# MotivatingSimulations
 This file contains two sub-files: SinglePopulation and TwoSubpopulations.
## SinglePopulation from MotivatingSimulations
 DiscreteModel.jl is responsible for generating all data needed to create Figure 2(a), and Figure2a.jl is responsible for creating Figure 2(a) using the data.
## TwoSubpopulations from MotivatingSimulations
 DiscreteModel.jl is responsible for generating all data needed to create Figure 2(b), and Figure2b.jl is responsible for creating Figure 2(b) using the data.
# Case1
 This file contains three sub-files: Data, AdditiveGaussian, and Multinomial.
## Data from Case1
 DiscreteModel.jl is responsible for generating all data needed to create Figure 3(a)-(d) and generate result for both the additive Gaussian measurement error model and the multinomial measurement error model. Figure3a-d.jl is responsible for creating Figure 3(a)-(d) using the data.
## AdditiveGaussian from Case1
 AdditiveGaussian.jl is used to generate the results for the additive Gaussian measurement error model for Case 1 using the data generated in DiscreteModel.jl.
## Multinomial from Case1
 Multinomial.jl is used to generate the results for the multinomial measurement error model for Case 1 using the data generated in DiscreteModel.jl.
# Case2
 This file contains three sub-files: Data, AdditiveGaussian, and Multinomial.
## Data from Case2
 DiscreteModel.jl is responsible for generating all data needed to create Figure 3(e)-(j) and generate result for both the additive Gaussian measurement error model and the multinomial measurement error model. Figure3e-j.jl is responsible for creating Figure 3(e)-(j) using the data.
## AdditiveGaussian from Case2
 AdditiveGaussian.jl is used to generate the results for the additive Gaussian measurement error model for Case 2 using the data generated in DiscreteModel.jl.
## Multinomial from Case2
 Multinomial.jl is used to generate the results for the multinomial measurement error model for Case 2 using the data generated in DiscreteModel.jl.
