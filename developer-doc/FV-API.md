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

![Sphere geometry and mesh for finite volume solver](https://github.com/sourceryinstitute/MORFEUS-Source/raw/fv-main-page/developer-doc/media/sphere.png)

<!-- 'fv-main-page' should be 'master' -->


Solver Description 
------------------

The morfeus solver consists of two parts:

  - A common library of routines for reading input files, creating the grid, discretization of the equations, solving the equations using PSBLAS, and plotting of the results. 
  
  - A problem-specific solver built using the routines from the morfeus finite volume library.


Timestepping and Integration
----------------------------

Morfeus is a finite volume solver for transport equations. The solver uses  explicit finite difference scheme in time and cell-based finite volume scheme to solve the transport equations. Morfeus is built on top on PSBLAS using an object-oriented paradigm in Fortran. For a more complete description of the solver implementation, refer to [Toninel thesis].

[Toninel thesis]: http://people.uniroma2.it/salvatore.filippone/nemo/toninel_phd.pdf

Input Files 
------------

Morfeus requires two input files. The first is a geometry file in GAMBIT neutral file format or in EXODUS II format; the second is a file called fast.json that contains the problem description, i.e. the description of materials, boundary conditions, solver parameters such as convergence criteria, and output parameters. The `fast.json` file should be present in the same folder as the mesh-file and the solver. 



```
{
  "MORFEUS_FV": {
    "MESH": {
      "mesh-dir": "./",
      "mesh-file": "coaxial_cable.e",
      "scale": 1,
      "gps-renumbering": false,
      "partitioning-scheme": "Block",
      "write-matrix-pattern-diagnostics": false
    },
    "MATERIALS": {
      "copper-core": 100,
      "dielectric-insulator": 4,
      "copper-sheath": 100,
      "plastic-sheath": 402,
      "tape": 100
    },
    "PDES": {
      "Energy": {
        "convergence-method": "BICGSTAB",
        "preconditioning-method": "BJAC",
        "solve-epsilon": 1e-08,
        "max-solve-iterations": 100,
        "write-matrix-system-diagnostics": false
      },
      "Concentration": {
        "convergence-method": "BICGSTAB",
        "preconditioning-method": "BJAC",
        "solve-epsilon": 1e-08,
        "max-solve-iterations": 100,
        "write-matrix-system-diagnostics": false
      },
      "Displacement": {
        "convergence-method": "BICGSTAB",
        "preconditioning-method": "BJAC",
        "solve-epsilon": 1e-06,
        "max-solve-iterations": 100,
        "write-matrix-system-diagnostics": false
      }
    },
    "iterations": {
      "time": {
        "max-steps": 200,
        "delta-t": 1
      },
      "big-solver": {
        "tolerance": 1e-08,
        "max-steps": 100
      },
      "output" : {
        "max-steps": 5
      }
    },
    "output": {
      "format": "vtk",
      "base-path": "out_nemo"
    },
    "Source-terms": {
      "temperature": {
        "sc": {
          "value": 10,
          "units": "[W/m^3]"
        },
        "sp": {
          "value": 0,
          "units": "[W/m^3 K]"
        }
      },
      "concentration": {
        "sc": {
          "value": 10,
          "units": "[W/m^3]"
        },
        "sp": {
          "value": 0,
          "units": "[W/m^3 K]"
        }
      }
    },
    "BCS": [
      {
        "type": "wall",
        "description": "copper-core 0",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "dielectric-insulator 0",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "copper-sheath 0",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "plastic-sheath 0",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "tape 0",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "copper-core 1",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "dielectric-insulator 1",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "copper-sheath 1",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "plastic-sheath 1",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "tape 1",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "copper-core surface",
        "temperature": {
          "id" : 1,
          "value" : 550.0
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      },
      {
        "type": "wall",
        "description": "tape surface",
        "temperature": {
          "id" : 2
        },
        "stress": {
          "id" : 1
        },
        "velocity": {
          "id" : 1
        }
      }
    ]
  }
}

```


High Level Objects
-------------------
The example code listing shows the input file for solving the transport equations inside a coaxial cable. The input file is written using the `json` format and is read using the `json-fortran` library. The file contains a top-level object called `MORFEUS` which contains child objects that define various problem parameters. The child objects can occur in any order within the top-level  order without affecting the functionality of the solver.


### Mesh object
The `MESH` object in the input file describes the name and directory of the file containing the mesh geometry. The coaxial cable geometry consists of several concentric materials. All computations in morfeus are done in MKS units. Any scaling of the geometry is also specified in this object. Mesh renumbering, partitioning scheme, and mesh diagnostics output for debugging can also be specified in this object.
```
"MESH": {
  "mesh-dir": "./",
  "mesh-file": "coaxial_cable.e",
  "scale": 1,
  "gps-renumbering": false,
  partitioning-scheme": "Block",
  "write-matrix-pattern-diagnostics": false
}
```

### Materials object
The `MATERIALS` object of the input file describes the materials used in the problem. Each material corresponds to a group of cells in the geometry file, and the materials should be listed in the same order as the cell groups in the geometry file.
```
"MATERIALS": {
  "copper-core": 100,
  "dielectric-insulator": 4,
  "copper-sheath": 100,
  "plastic-sheath": 402,
  "tape": 100
}
```

### PDES object
The `PDES` objects describes the solver parameters for the different PDEs. It contains child objects for each of the  energy, concentration and displacement PDEs that need to be solved. The `Energy`  object describes the parameters used for the heat transfer equation, the `Concentration`  object describes the parameters used for the diffusion equations solver, and the `Displacement` object describes the parameters for the solid mechanics equations solver. The main parameters defined in each of these objects are the numerical solver method for the matrix problem, the matrix preconditioner used, convergence criteria, and maximum number of interactions.  For illustration, only the `Energy` PDE is shown below.
```
"Energy": {
  "convergence-method": "BICGSTAB",
  "preconditioning-method": "BJAC",
  "solve-epsilon": 1e-08,
  "max-solve-iterations": 100,
  "write-matrix-system-diagnostics": false
}
```

### Iterations object
The `iterations` object contains three child objects which define the solver iterations and convergence parameters.

* Time object
The `time` object describes the number of time steps to be executed and the $\delta t$ in seconds.

* Big-solver object
The `big-solver` object describes the convergence tolerance for the iterative solver, and the maximum number of steps to be executed before moving to the next time-step.

* Output object
The `output` object describes how frequently (per how many time steps) the results should be output in the chosen format.

```
"iterations": {
  "time": {
    "max-steps": 200,
    "delta-t": 1
  },
  "big-solver": {
    "tolerance": 1e-08,
    "max-steps": 100
  },
  "output" : {
    "max-steps": 5
  }
}
```

### Output object
The `output` object describes the output filename and format.
```
"output": {
  "format": "vtk",
  "base-path": "out_nemo"
}
```

### Source terms object
The `Source-terms` objects contains two child objects that contain the source terms for the `Energy` PDE and the `Concentration` PDE.
```
"Source-terms": {
  "temperature": {
    "sc": {
      "value": 10,
      "units": "[W/m^3]"
    },
    "sp": {
      "value": 0,
      "units": "[W/m^3 K]"
    }
  },
  "concentration": {
    "sc": {
      "value": 10,
      "units": "[]"
    },
    "sp": {
      "value": 0,
      "units": "[1/K]"
    }
  }
}
```


### Boundary conditions object
The `BCS` object contains a list of child objects for each boundary surface. The number of child objects should match the boundary surfaces in the geometry file and should be in the same order as the surfaces in the geometry file.
```
"BCS": [
  {
    "type": "wall",
    "description": "copper-core 0",
    "temperature": {
      "id" : 1,
      "value" : 500
    },
    "stress": {
      "id" : 1
    },
    "velocity": {
      "id" : 1
    }
  },
  {
    "type": "wall",
    "description": "dielectric-insulator 0",
    "temperature": {
      "id" : 2
    },
    "stress": {
      "id" : 1
    },
    "velocity": {
      "id" : 1
    }
  }
]
```


Table: Temperature boundary conditions

| BC ID <br>| Description       | Value 1          | Value 2          
|-------|-------------------|------------------|------------------:
| 1     | Fixed temperature | Temperature in K |                  
| 2     | Adiabatic BC      |                  |                  
| 3     | Fixed flux        | flux in W/m2     |                  
| 4     | Convection        | Coeff in W/m2K   | Temperature in K 


Table: Stress boundary conditions

| BC ID <br>  | Description   | Value   
|------------ |-------------------  |----------------:
| 1   | Stress free   |   
| 2   | Prescribed Stress   | Stress in N/m2  


Table: Velocity boundary conditions

| BC ID <br>| Description <br>| Value 
|-------|-------------|-------:
| 1 | No slip | 
| 2 | Free slip | 
| 3 | Sliding | 
| 4 | Moving | Velocity in m/s 
| 5 | Free sliding | 


High Level Classes
------------------

Below is a list of the important high-level objects and classes, with a brief discussion and
links to detailed API documentation and their implementations.
A complete list of all FD and FV classes and types can be found on the [types list page].

[types list page]: ../lists/types.html

### [[object(type)]]

`[[object(type)]]` is an abstract type to ensure descendents declare the required deferred bindings to

#### Methods

* `[[object(type):mark_as_defined]]` : Set the object as defined
* `[[object(type):user_defined]]` : Query whether the object is defined

### [[grid(type)]]

`[[grid(type)]]` ss the Morfeus universal base type for all grids, for FD and FV grids. It extends `[[object(type)]]`.

#### Methods

* `[[grid(type):set_units]]` : Set the physical units used by the grid
* `[[grid(type):get_units]]` : Get the physical units used by the grid


### [[field(type)]]

`[[field(type)]]` is the Morfeus base type for all fields. It extends `[[grid(type)]]`.

#### Methods

* `[[field(type):create_field]]` : Class constructor
* `[[field(type):free_field]]` : Class destructor
* `[[field(type):on_faces_]]` :
* `[[field(type):mat_]]` :
* `[[field(type):bc_]]` :
* `[[field(type):fld_size]]` :
* `[[field(type):get_material]]` :
* `[[field(type):msh_]]` :
* `[[field(type):dim_]]` :
* `[[field(type):name_]]` :
* `[[field(type):set_field_dim]]` :
* `[[field(type):set_field_on_faces]]` :
* `[[field(type):check_field_operands]]` :
* `[[field(type):nemo_sizeof]]` :
* `[[field(type):get_mesh]]` :
* `[[field(type):check_mesh_consistency]]` :


### [[material(type)]]

`[[material(type)]]` is a class to describe material state and specify state equations. 

#### Methods

* `[[material(type):name_]]` : 
* `[[material(type):mat_id]]` :
* `[[material(type):nemo_sizeof]]` :

### [[matptr(type)]]

`[[matptr(type)]]` is a class with a single `[[material(type)]]` type object

### [[mesh(type)]]

`[[mesh(type)]]` is a class to define and manipulate data describing the discretization of space into connected finite-volume cells and surfaces

#### Methods

* `[[mesh(type):create_mesh]]` : 
* `[[mesh(type):free_mesh]]` : 
* `[[mesh(type):check_mesh_unused_el]]` : 
* `[[mesh(type):nemo_sizeof]]` : 

### [[bc_poly(type)]]

`[[bc_poly(type)]]` Implement a runtime polymorphism pattern [1] to emulate an abstract parent of the `[[bc_math(type)]]`/`[[bc_wall(type)]]` classes. 

[1] Akin, E. (2003) Object-oriented programming via Fortran 90/95. Cambridge University Press.

#### Methods

* `[[bc_poly(type):get_abc]]` : 
* `[[bc_poly(type):set_bc]]` : 
* `[[bc_poly(type):surface_motion]]` : 
* `[[bc_poly(type):vertex_motion]]` : 
* `[[bc_poly(type):get_displacement]]` : 
* `[[bc_poly(type):get_velocity]]` : 
* `[[bc_poly(type):nemo_sizeof]]` : 

### [[scalar_source(type)]]

`[[scalar_source(type)]]`

### [[scalar_field(type)]]

`[[scalar_field(type)]]`

### [[scalar_pde(type)]]

`[[scalar_pde(type)]]`

#### Methods

* `[[scalar_pde(type):create_pde]]` :
* `[[scalar_pde(type):free_pde]]` :
* `[[scalar_pde(type):write_scalar_pde]]` :
* `[[scalar_pde(type):geins_pde]]` :
* `[[scalar_pde(type):nemo_sizeof]]` :
* `[[scalar_pde(type):reinit_pde]]` :
* `[[scalar_pde(type):asb_pde_]]` :
* `[[scalar_pde(type):solve_pde]]` :

### [[vector(type)]]

`[[vector(type)]]`

### [[vector_field(type)]]

`[[vector_field(type)]]`

### [[vector_pde(type)]]

`[[vector_pde(type)]]`

#### Methods

* `[[vector_pde(type):create_pde]]` :
* `[[vector_pde(type):free_pde]]` :
* `[[vector_pde(type):write_vector_pde]]` :
* `[[vector_pde(type):geins_pde]]` :
* `[[vector_pde(type):nemo_sizeof]]` :
* `[[vector_pde(type):reinit_pde]]` :
* `[[vector_pde(type):asb_pde_]]` :
* `[[vector_pde(type):solve_pde]]` :

### [[iterating(type)]]

`[[iterating(type)]]`


High Level Procedures and Methods
---------------------------------

All Morfeus FV and FD procedures are listed on the [procedures list page], but below is a curated list of those that correspond to high-level
operations in FV. These should help you write your own programs and kernels using Morfeus FV.

[procedures list page]: ../lists/procedures.html

### [[assert(proc)]]

`[[assert(proc)]]` Is an assertion utility used in design by contract (DBC) for enforcing pre-conditions, post-conditions and invariants, as well in testing.


@todo
Add other stand-alone procedures used to instantiate objects of for other purposes where a TBP method is not used.
