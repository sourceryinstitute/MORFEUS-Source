---
title: High-level FD API documentation
---


<details><summary><b>Table of Contents</b></summary>

[TOC]

</details>

@warning The block-structure parallel FD solver is still a work in progress

Overview
--------

Here an overview of the Morfeus finite difference (FD) framework is given.
Included is a brief description of how to setup a new solver using the framework.
Key objects/classes, and procedures are outlined.
An overview of reading input files, reading grid files,
initializing the library and data structures and, setting up PDEs, integrating the solution forward in time
and writing outputs is given.

Solver and Problem Setup
------------------------

@todo
This needs more text added here.

Time-stepping and Integration
----------------------------

@todo
This needs more text added here.
Presumably something better than forward Euler is being used...

\begin{equation}
u_i^{n+1} = u_i^n + \Delta_t f(t,u_i^n)
\end{equation}

High Level Classes and Objects
------------------------------

Below is a list of the important high-level objects and classes, with a brief discussion and
links to detailed API documentation and their implementations.
A complete list of all FD and FV classes and types can be found on the [types list page].

[types list page]: ../lists/types.html

@todo
This needs updating/finishing


### [[object(type)]]

`[[object(type)]]` is an abstract type to ensure descendents declare the required deferred bindings to

#### Methods

* `[[object(type):mark_as_defined]]` : Set the object as defined
* `[[object(type):user_defined]]` : Query whether the object is defined

### [[grid(type)]]

`[[grid(type)]]` Is the Morfeus universal base type for all grids, even for FV grids. It extends `[[object(type)]]`.

#### Methods

* `[[grid(type):set_units]]` : Set the physical units used by the grid
* `[[grid(type):get_units]]` : Get the physical units used by the grid

### [[structured_grid(type)]]

This is the main structured grid class used by Morfeus FD. It extends the `[[grid(type)]]` class.
Methods defined for this class can get information about the global grid as well as this images particular block.

#### Methods

@todo
This section needs to be filled in/expanded!

* `[[structured_grid(type):clone]]` :
* `[[structured_grid(type):space_dimension]]` :
* `[[structured_grid(type):set_metadata]]` :

### [[cartesian_grid(type)]]

@todo
add description and how to use, etc.

#### Methods

@todo
discussion of important methods provided and how to use them, etc.

### [[curvilinear_grid(type)]]

@todo
add description and how to use, etc.

#### Methods

@todo
discussion of important methods provided and how to use them, etc.

### [[co_object(type)]]


### [[surfaces(type)]]


### [[material_t(type)]]


### [[thickness_t(type)]]


### [[problem_discretization(type)]]


### [[package(type)]]


### [[material_t(type)]]


### [[flux_planes(type)]]


### [[plate_3D(type)]]


### [[geometry(type)]]


### [[subdomain_t(type)]]


### etc., etc.


High Level Procedures and Methods
---------------------------------

All Morfeus FV and FD procedures are listed on the [procedures list page], but below is a curated list of those that correspond to high-level
operations in FD. These should help you write your own programs and kernels using Morfeus FD.

[procedures list page]: ../lists/procedures.html

### [[assert(proc)]]

`[[assert(proc)]]` Is an assertion utility used in design by contract (DBC) for enforcing pre-conditions, post-conditions and invariants, as well in testing.


@todo
Add other stand-alone procedures used to instantiate objects of for other purposes where a TBP method is not used.
