---
project: Morfeus-FD
summary: A parallel finite difference solver for block-structured grids
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
         private
source: true
graph: true
md_extensions: markdown.extensions.toc
coloured_edges: true
sort: alpha
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

{!src/FD/README.md!}
