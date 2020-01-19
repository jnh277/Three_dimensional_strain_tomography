# Three dimensional strain tomography
Demonstration of full-field triaxial strain tomography from neutron transmission strain images. Strain images were RADEN energy-resolved-neutron-imaging instrument at the Japan Proton Accelerator Research Complex (J-PARC) in Japan (Jan 2019). This repository contains companion code to the paper available at https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.3.113803 and https://arxiv.org/abs/1906.08506, and has matlab example scripts to reconstruct and plot the comparisons given in the paper.

The strain image images have been preprocessed and is included in 'strain_image_data_set.m'. This involved Bragg-edge fitting, computing the relative strain values, and determining the LRT geometry (exit and entry locations of each ray through the sample). This preprocessed data is included in the repository. The original strain images data is very large but can be made available on request.

Three example matlab scripts are included:
- 'example_plane1.m': Reconstructs the strain field for cutting plane 1 and produces Figures 3 and 5 from the paper.
- 'example_plane2.m': Reconstructs the strain field for cutting plane 2 and produces Figures 4 and 6 from the paper.
- 'example_plane3.m': Reconstructs the strain field for cutting plane 3 and produces Figure 7 from the paper.


Note: The reconstruction process is quite memory intensive and may not run well on some personal computers. Additionally, hyperparameters for the Gaussian process covariance functions have been prechosen by performing an optimisation to maximise the Marginal Log Likelihood. This process takes several days and so is not included in the example scripts. 
