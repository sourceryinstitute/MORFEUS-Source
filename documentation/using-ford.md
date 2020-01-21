---
title: Documenting New and Modified Code
---

# FORD the Fortran Documentation generator

## Getting started

[FORD] has is [quite well documented]. Documentation is written in
markdown syntax using the [python markdown] parser and
generator. There are some important and often subtle differences from
Github flavored markdown ([GFM]). In addition, other markdown source
can be [included] in any FORD markdown with the custom FORD syntax
`\{!include_file.md!\}` which gets injected before the document is
processed. (This allows code to be included this way as well as [meta
data].)

[FORD] parses and understands modern Fortran and was written in a way
to try to make [writing documentation] easy. The project options are
mostly set using the [project file], although some options may be
overridden on the [command line]. The usual [FORD] documentation
comment string is `!!` with documentation following immediately after
the entity being documented. This documentation sentinel can be
changed in the [project file].

## Gotchas

### Meta Data

The trickiest thing to deal with is probably [meta data]. Meta data
must *NOT* have an proceeding FORD documentation strings for a given
entity even if they are empty. For example:

```fortran

SUBROUTINE bad_sub1()
!!
!! author: me (my company)
!! date: 1/1/3091
!! regular markdown
...
END SUBROUTINE

SUBROUTINE bad_sub2()
!! This is my sub
!! author: J. Doe
!! date: May 1908
!!
!! more normal markdown
...
END SUBROUTINE

SUBROUTINE good1()
!! Author: me (my organization)
!! Date: June 2018
!!
!! good1() frobnicate the widget <!--This becomes the brief/summary-->
!!
!! More main body descriptive text.
...
END SUBROUTINE

SUBROUTINE good2()
!! summary: good2 provides not the widget you want, but the widget you need
!!     Meta data may be continued on the next IFF there are 4 spaces or more of indentation
!!     The summary: meta data is basically equivalent to Doxygen's !>@brief but can just be specified
!!     as normal FORD markdown after the meta data, as shown in [good1].
!! author: me
!!
!! more info, this won't be included in the summary/brief text
...
END SUBROUTINE
```

As can be seen above, the summary (`!>@brief` in Doxygen parlance)
text may be specified either using the `!! summary: ...` meta data
entry __OR__ as the first paragraph of the regular FORD markdown
documentation text.

### Comments (in body text and meta data)

The second trickiest thing might be comments that you do not wish to
render with the documentation (but will be visible if someone reads
the HTML source or the original code source file). Comments CANNOT BE
INCLUDED IN THE META DATA SECTION. In the normal markdown section HTML
comments, e.g., `<!-- your comment goes here -->` or the dummy link
declaration comments:

```
!! [Comment text can go here, linking to an empty anchor]:#
```

<!-- here is an HTML style comment that is only visible when editing this wiki page's source -->
[Here is a dummy anchor/link markdown style comment, only visible when editing this wiki page's source]:#

Finally, a minor annoyance is that there is no convenient way to
document many variables declared on the same line. The best solution
is usually to split their declarations into multiple individual
declarations.

### FORD errors and bugs

[FORD] is not always maintained as actively as one would hope, but
@cmacmackin has expressed interest in passing the torch, possibly to
the @Fortran-FOSS-Programmers group. As such, when you encounter a bug
it is important to try to reproduce it easily and provide a Python
backtrace. To get python backtraces from [FORD] the `--debug` command
line flag may be passed when generating documentation.

### Slow output generation

A few steps may be taken to speed up running [FORD].

 - Install the `lxml` python package; this helps speed up generating search data
 - Skip generating search data:
   - Pass `-DALWAYS_SKIP_SEARCH_GEN=YES` to CMake when configuring MORFEUS
   - Pass `--no-search` to [FORD] when running manually
 - Edit the project file meta data field for graph generation: `graph: true` --> `graph: false`
 - Don't run [FORD] at all!
   - Pass `-DSKIP_DOC_GEN=YES` to CMake during configuration
   - Don't have `ford` on your `PATH`. (Don't install it, or move it, or edit your `PATH`.)

## Example code

For an example of the suggested FORD documentation style, please see
the [style guide] page in this wiki.

[FORD]: https://github.com/Fortran-FOSS-Programmers/ford
[quite well documented]: https://github.com/Fortran-FOSS-Programmers/ford/wiki
[python markdown]: https://python-markdown.github.io/#features
[GFM]: https://github.github.com/gfm/
[meta data]: https://github.com/Fortran-FOSS-Programmers/ford/wiki/Documentation-Meta-Data
[project file]: https://github.com/Fortran-FOSS-Programmers/ford/wiki/Project-File-Options
[command line]: https://github.com/Fortran-FOSS-Programmers/ford/wiki/Command-Line-Options
[included]: https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation#include-capabilities
[writing documentation]: https://github.com/Fortran-FOSS-Programmers/ford/wiki/Writing-Documentation
[style guide]: ./style-guide
