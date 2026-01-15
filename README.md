[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

# Virtual Neuromodulation Toolbox
The Virtual Neuromodulation Toolbox for MATLAB.

## Introduction
Whole-brain (BOLD signal) data-driven model assumes no specific brain tissue structure, one voxel of 4 × 4 × 4 mm gray matter, and 25445 voxels based on the Allen human brain atlas. 
All voxels are fully connected, and vector auto-regression (VAR) surrogate model is used to learn group data and generate group surrogate data.

<div align="center">
<img src="data/img/fig1.jpg" width="70%">
</div>


Then, we extended output of VAR surrogate model for virtual neuromodulation.
<div align="center">
<img src="data/img/fig2.jpg" width="50%">
</div>


where x<sub>i</sub>(t) is output of voxel i, y<sub>i</sub>(t) is VAR surrogate output of voxel i, z<sub>i</sub>(t),u<sub>i</sub>(t) is modulation term of voxel i, and c<sub>i</sub>,σ<sub>X</sub>∈R. 
σ<sub>X</sub> is calculated as the standard deviation from the entire voxel time-series.
z<sub>i</sub>(t) is constructed by convolution of the canonical Hemodynamic Response Function (HRF) and the Box-car task design.
u<sub>i</sub>(t) performs direct adjustment to the output of VAR surrogate.
In other words, if there is a virtual neuromodulation stimulus, it adjusts to prioritize neuromodulation stimulus, such as u<sub>i</sub>(t)=1-0.5∙z<sub>i</sub>(t).
Using the above-mentioned virtual neuromodulation, BOLD signal addition, i.e., DBS treatment, can be virtually performed for specific voxels. 


<b>Command line tools</b>

| name | description |
|:---|:---|
| mtess | Calculate and plot MTESS for a group of multivariate time-series data. |
| gsdgm | Generate a whole-brain data-driven model based on the group surrogate model (VAR surrogate).|
| vneumod | Generate virtual neuromodulation time-series surrogate data based on the whole-brain data-driven model.|

## Requirements: Software
* MATLAB R2019b or later
* Deep Learning Toolbox ver12.1 or later
* Fuzzy Logic Toolbox ver2.6 or later
* Econometrics Toolbox ver5.3 or later
* Parallel Computing Toolbox ver7.1 or later
* [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn)

Please download the [VARDNN Toolbox](https://github.com/takuto-okuno-riken/vardnn) and "Add Path" in the MATLAB before using Virtual Neuromodulation toolbox.


## Command Line Tools Demos
<b>Demo 1</b><br>


## Command Line Tools Reference
<b>mtess command</b><br>


## Citing Virtual Neuromodulation Toolbox
If you find Virtual Neuromodulation Toolbox useful in your research, please cite it as follows: 

Takuto Okuno, Alexander Woodward, Yuji Takahashi, Hideyuki Okano, Junichi Hata (20XX)
["A digital brain study: Is deep brain stimulation activating associative circuit?"](https://www.google.com/), ,


