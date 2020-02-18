---
project: Morfeus Framework
summary: A parallel finite volume and finite difference PDE toolkit
src_dir: src/FD
         src/FV/src
exclude_dir: src/FD/tests
             src/FV/src/unit-tests
             API-doc
output_dir: API-doc
media_dir: developer-doc/media
extra_filetypes: txt #
page_dir: developer-doc
preprocess: true
macro: FORD
preprocessor: gfortran-8 -E
display: public
         protected
source: true
graph: true
md_extensions: markdown.extensions.toc
coloured_edges: true
sort: permission-alpha
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            iso_c_binding:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html#ISO_005fC_005fBINDING
lower: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/sourceryinstitute/MORFEUS-Source
project_download: https://github.com/sourceryinstitute/MORFEUS-Source/releases
license: bsd
author: Guide Star Engineering, LLC
email: damian@gsellc.com
author_description: An engineering services company specializing in hardware & software engineering design, R&D, testing, and systems integration.
author_pic: media/GSElogo.png
favicon: developer-doc/media/gsellcfavicons.png
github: https://github.com/sourceryinstitute/MORFEUS-Source
website: https://gsellc.com
---

[_____ Comments _______]:#
[source: display source code corresponding to item being documented]:#
[graph: generate call graphs, module dependency graphs, derive type composition/inheritance graphs ]:#
[sort: different sorting schemes for the modules or procedures or programs or derived types (alpha = alphabetical see wiki).]:#
[extra_mods: documentation for intrinsic modules]:#

[This document is a FORD project file, formatted with Pythonic Markdown                                      ]:#
[See https://github.com/Fortran-FOSS-programmers/ford/wiki/Project-File-Options for more info on writing FORD project files]:#

[TOC]

@warning While the finite volume (FV) capability of this package is relatively mature, the finite difference (FD) capability is still a work in progress.

Morfeus Developer Documentation
===============================

Welcome to the Morefeus Framework developer documentation.
This online documentation is automatically generated from inline comments and static analysis
using the [FORD] tool.
Please keep reading to learn how to best use this documentation to become more familiar with
the Morfeus Framework and how to get started building an application or modifying one that uses
this framework.

[FORD]: https://github.com/Fortran-FOSS-Programmers/ford#readme


Organization
------------

The [FORD] tool is used to document Modern Fortran source code.
At first browsing the documentation can be overwhelming, but understanding the structure of the documentation and
consulting the getting started and high-level API documentation should provide an appropriate orientation.

### Top Navigation Bar

The chief means of navigating through the source code is by using the black navigation bar at the top of this landing page.
In addition to the search box on the right side of the top navbar, the following links are available:

* [Morfeus Framework]:  
  This link takes the viewer back to this developer documentation homepage.
* [Developer Documentation]:  
  Visit the high-level developer documentation which outlines:
  * how to get started
  * procedures and classes most relavent to creating or modifying a PDE solver
  * instructions for building the Morfeus Framework library
  * a high-level overview of the CMake based build system
* [Source Files]:  
  This landing page enumerates the source files associated with the project,
  includes a graph depicting their interdependencies and links to their dedicated pages.
  You probably do not want to start here as much of this information is redundant with
  the [Modules] landing page.
* [Modules]:  
  The landing page enumerating and describing all the modules and submodules in the project.
  More usefull than the files page, as it groups submodules with their parent modules and
  provides some description for what each one does.
* [Procedures]:  
  Enumerates all procedures including generic overloaded interfaces, operators, module procedures, functions and subroutines.
  The table includes a link to which module the procedure is defined in, and a description of what type of procedure it is.
  `Interface` denotes an overloaded procedure or operator, or a module procedure
  with the interface defined in a module and the implementation defined in a submodule.
  This is a decent place to start if you are looking for a particular procedure,
  and where to find the specific implementation.
* [Abstract Interfaces]:  
  Enumerates abstract procedures defined, with links to the modules in which they are defined.
  Abstract procedures only serve as prototypes for other code to implement, or are used to
  specify deferred type bound procedures.
* [Derived Types]:  
  A list enumerating the derived types or classes defined within Morfeus.
  The list includes links to the type definitions, a link to any type they extend and a link to the module in which they are defined.
  This is one of the better places to start when trying to understand a particular object, its methods (type bound procedures),
  where it is defined, and whether there are any types that extend it or types which it extends.

[Morfeus Framework]: https://sourceryinstitute.github.io/MORFEUS-Source/index.html
[Developer Documentation]: https://sourceryinstitute.github.io/MORFEUS-Source/page/index.html
[Source Files]: https://sourceryinstitute.github.io/MORFEUS-Source/lists/files.html
[Modules]: https://sourceryinstitute.github.io/MORFEUS-Source/lists/modules.html
[Procedures]: https://sourceryinstitute.github.io/MORFEUS-Source/lists/procedures.html
[Abstract Interfaces]: https://sourceryinstitute.github.io/MORFEUS-Source/lists/absint.html


### [Developer Documentation]

Navigating to the [Developer Documentation] page using the top navbar provides information, written in prose using Markdown,
to help orient developers to the Morfeus Framework.
The two most important pages to visit in this section are the [Getting Started] page and the [High-level FV API and Usage] page.

The [Getting Started] page contains information about:

* Prerequisites
  * software
  * compilers
  * OSes
* Configuring with [CMake]
* Building (compiling) the library
* Running tests

In addition, projects utilizing the Morfeus Framework are encouraged to provide instructions for obtaining,
configuring, and compiling the project in question in a `README` file distributed with the source code.
These files are typically rendered to HTML and displayed prominently on version control websites like
[Github], [Gitlab] or [BitBucket] to name a few.

The [High-level FV API and Usage] page provides:
* an overview and references for the numerical algorithms used
* input files and entries therein to control the mesh, governing PDEs, boundary conditions and input/output
* a list of high-level classes used in constructing a solver with links to the detailed [FORD] documentation for each
  * sub-lists of the most importand methods each class provides and links to their detailed documentation
* a list of other high-level procedures that are not type bound procedures (methods)


[Getting Started]: https://sourceryinstitute.github.io/MORFEUS-Source/page/getting-started.html
[High-level FV API and Usage]: https://sourceryinstitute.github.io/MORFEUS-Source/page/FV-API.html
[CMake]: https://cmake.org/cmake/help/latest/
[Github]: https://github.com
[Gitlab]: https://gitlab.com
[BitBucket]: https://bitbucket.org


Getting Help
------------

If you encounter a problem, have a suggestion, or want to ask a question,
we encourage you to post an issue in [this projects Github repository] by
[opening a new issue]. Every effort will be made to respond to your inquiry in a timely fashion.

[this projects Github repository]: https://github.com/sourceryinstitute/MORFEUS-Source
[opening a new issue]: https://github.com/sourceryinstitute/MORFEUS-Source/issues/new
