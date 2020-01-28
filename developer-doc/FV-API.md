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
A complete list of all FD and FV classes and types can be found on the [types list page].

[types list page]: ../lists/types.html

### [[mesh(type)]]


### [[matptr(type)]]


### [[bc_poly(type)]]


### [[scalar_source(type)]]


### [[scalar_field(type)]]


### [[vector_field(type)]]


### [[scalar_pde(type)]]


### [[vector_pde(type)]]


### [[vector(type)]]


### [[iterating(type)]]


### [[output(type)]]


High Level Procedures and Methods
---------------------------------


All Morfeus FV and FD procedures are listed on the [procedures list page], but below is a curated list of those that correspond to high-level
operations in FV. These should help you write your own programs and kernels using Morfeus FV.

[procedures list page]: ../lists/procedures.html

### [[start_psblas(proc)]]


### [[create_bc(proc)]]


### [[create_material(proc)]]


### [[create_scalar_field(proc)]]


### [[create_vector_field(proc)]]


### [[create_pde(proc)]]


### [[create_source(proc)]]


### [[create_iterating(proc)]]


### [[create_output(proc)]]


### [[compute_channel_parameters(proc)]]


### [[pde_ddt(proc)]]


### [[pde_laplacian(proc)]]


### [[pde_source(proc)]]


### [[solve_pde(proc)]]


<!-- These are all internal procedures, a bunch can/should be moved to Morfeus FV

### [[symmetric_(proc)]]


### [[adjust_strain(proc)]]


### [[compute_stress(proc)]]



### [[pde_div(proc)]]


### [[apply_flow_rule(proc)]]


### [[adjust_displacements(proc)]]



### [[coolant_enthalpy_rise(proc)]]


### [[corrosion_thickness(proc)]]

-->

### [[update_field(proc)]]


### [[stop_iterating(proc)]]


### [[write_vtk_morfeus(proc)]]


### [[wrfreq(proc)]]


### [[write_scalar_pde(proc)]]


### [[free_pde(proc)]]


### [[free_field(proc)]]


### [[free_bc(proc)]]


### [[free_mesh(proc)]]


### [[stop_psblas(proc)]]
