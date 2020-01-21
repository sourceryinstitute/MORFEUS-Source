---
project: Morfeus-FD
summary: A parallel finite difference solver for block-structured grids
macro: FORD
src_dir: src/FD
exclude_dir: src/FD/tests
             API-doc
output_dir: API-doc
media_dir: documentation/media
page_dir: documentation
project_url: https://sourceryinstitute.github.io/MORFEUS-Source
preprocess: true
display: public
         protected
         private
source: true
graph: true
md_extensions: markdown.extensions.toc
coloured_edges: true
sort: alpha
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
            iso_c_binding:https://gcc.gnu.org/onlinedocs/gcc-4.9.4/gfortran/ISO_005fC_005fBINDING.html
lower: true
print_creation_date: true
creation_date: %Y-%m-%d %H:%M %z
project_github: https://github.com/sourceryinstitute/MORFEUS-Source
project_download: https://github.com/sourceryinstitute/MORFEUS-Source/releases
license: bsd
author: Guide Star Engineering, LLC
email: damian@gsellc.com
author_description: An engineering services company specializing in hardware & software engineering design, R&D, testing, and systems integration.
author_pic: https://gsellc.com/wp-content/uploads/2019/07/logo21.png
github: https://github.com/sourceryinstitute/MORFEUS-Source
website: http://www.gsellc.com
---

[_____ Comments _______]:#
[source: display source code corresponding to item being documented]:#
[graph: generate call graphs, module dependency graphs, derive type composition/inheritance graphs ]:#
[sort: different sorting schemes for the modules or procedures or programs or derived types (alpha = alphabetical see wiki).]:#
[extra_mods: documentation for intrinsic modules]:#

[This document is a FORD project file, formatted with Pythonic Markdown                                      ]:#
[See https://github.com/Fortran-FOSS-programmers/ford/wiki/Project-File-Options for more info on writing FORD project files]:#

[TOC]

{!src/FD/README.md!}
