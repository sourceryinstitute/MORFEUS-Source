---
title: High-level FV API documentation
---


<details><summary><b>Table of Contents</b></summary>

[TOC]

</details>

Overview
--------

Here an overview of the Morfeus finite volume (FV) framework is given.
Included is a brief description of how to setup a new solver using the framework.
Key objects/classes, and procedures are outlined.
An overview of reading input files, reading grid files,
initializing the library and data structures and, setting up PDEs, integrating the solution forward in time
and writing outputs is given.

Solver and Problem Setup
------------------------

@todo
This needs more text added here.

Timestepping and Integration
----------------------------

@todo
This needs more text added here.
Pressumably something better than forward Euler is being used...

\begin{equation}
u_i^{n+1} = u_i^n + \Delta_t f(t,u_i^n)
\end{equation}

High Level Classes and Objects
------------------------------

Below is a list of the important high-level objects and classes, with a brief discussion and
links to detailed API documentation and their implementations.

### mesh


### matptr


### bc_poly


### scalar_source


### scalar_field


### vector_field


### scalar_pde


### vector_pde


### vector


### iterating


### output


High Level Procedures and Methods
---------------------------------

### start_psblas


### create_bc


### create_material


### create_scalar_field


### create_vector_field


### create_pde


### create_source


### create_iterating


### create_output


### compute_channel_parameters


### pde_ddt


### pde_laplacian


### pde_source


### solve_pde


### symmetric_


### adjust_strain


### compute_stress


### pde_div


### apply_flow_rule


### adjust_displacements



### coolant_enthalpy_rise


### corrosion_thickness


### update_field


### stop_iterating


### write_vtk_morfeus


### wrfreq


### write_scalar_pde


### free_pde


### free_field


### free_bc


### free_mesh


### stop_psblas
