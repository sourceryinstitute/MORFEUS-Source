---
project: morfeus-fd
summary: A parallel finite difference solver
src_dir: grid
src_dir: tests
src_dir: utilities
output_dir: doc
preprocess: true
display: public
        protected
        private
source: true
graph: true
sort: alpha
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/sourceryinstitute/morfeus
project_download: https://github.com/sourceryinstitute/morfeus/releases
license: bsd
author: Sourcery Institute
email: damian@sourceryinstitute.org
author_description: A California public-benefit nonprofit corporation engaged in research, education, and consulting in computational science. engineering, and mathematics
author_pic: http://www.sourceryinstitute.org/uploads/4/9/9/6/49967347/sourcery-logo-rgb-hi-rez-1.png
github: https://github.com/sourceryinstitute
website: http://www.sourceryinstitute.org
---

[source: display source code corresponding to item being documented]:#
[graph: generate call graphs, module dependency graphs, derive type composition/inheritance graphs ]:#
[sort: different sorting schemes for the modules or procedures or programs or derived types (alpha = alphabetical see wiki).]:#
[extra_mods: documentation for intrinsic modules]:#

[This document is a FORD project file, formatted with Pythonic Markdown                                      ]:#
[See https://github.com/cmacmackin/ford/wiki/Project-File-Options for more info on writing FORD project files]:#

--------------------

Compilers
---------
@warning
This archive makes extensive use of Fortran 2018.   We recommend compiling wiht with the latest available version of any of the following recommended compilers:

* GNU Compiler Collection (GCC)
* Intel
* Cray
