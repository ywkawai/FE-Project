# FE-project 

## What is FE-project?
This project provides a library and some sample programs for some discontinuous Galerkin methods. 
Futhermore, some dynamical cores of atmospheric and oceanic model with DGM will be provided. 

- Example of simulation results by nonhydrostatic atmospheric models with nodal DGM

![A simulation of density current](https://github.com/ywkawai/FE-project/wiki/gallery/atm_nonhydro2d/density_current/density_current.gif)

![A simulation of idealized baroclinic instability](https://github.com/ywkawai/FE-project/wiki/gallery/atm_nonhydro3d/baroc_wave/BarocWave.gif)

 For more animations, please see 'FE-project gallery' channel on YouTube (url: https://www.youtube.com/channel/UCO17OQtKHwkkQwmHD9y9mQg/featured) or Wiki of FE-project on github. 


## Models with FE library
### Shallow water model
- Global shallow water model using cubed-sphere mesh

### Nonhydrostatic atmospheric model
- Simple 2D model with only dynamical process
- Regional/Global model with dynamical process and some physical processes

## Simple samples for intrduction to DGM. 
### 1D problems 
  - Linear advection equation
    - various profiles of advected quantity
  - linear advection-diffusion equation

### 2D problems 
  - Linear advection equation in a rectangle domain
    - various profiles of advected quantity and flow
  - Linear advection equation in a cubed sphere domain

### 3D problems 
  - Linear advection equation in a cubic domain
    - various profiles of advected quantity
  - Linear advection equation in a cubed sphere domain
  - Euler equation in a cubic domain
    - Test the propagation of sound waves with HEVI temporal methods

 
## Build FE library
Please see ``INSTALL.md''.  


