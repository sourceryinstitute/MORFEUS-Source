---
title: Developer Documentation
---

<details><summary><b>Table of Contents</b></summary>

[TOC]

</details>

About
-----

This page and directory contain links to various high-level [Morfeus] developer documentation pages.
This content is written in Pythonic Markdown and rendered to HTML by [FORD] and served at [https://sourceryinstitute.github.io/MORFEUS-Source].
For more details on writing documentation with [FORD], please see the [FORD wiki].

Contents
--------

* __[Getting Started]__: Information required for building [Morfeus] from source
    * __[Prerequisites]__: Required software libraries and tooling
    * __[Configuring the Library]__: Running [CMake] to configure the [Morfeus] build, find compilers and optional and required prerequisites
    * __[Building the Library]__: Compiling [Morfeus] from source
    * __[Running Tests]__: How to run unit and integration tests with [CTest]
* __[Documenting New and Modified Code]__: How to write [FORD] documentation, and what documentation is required in contributed code
* __[Style Guide]__: Style requirements, many of which are enforced with [EditorConfig] or [findent].
* __[Example Usage]__: Example code and discussion on how to use [Morfeus] in your application


[Morfeus]: https://github.com/sourceryinstitute/MORFEUS-Source#readme
[https://sourceryinstitute.github.io/MORFEUS-Source]: https://sourceryinstitute.github.io/MORFEUS-Source
[FORD]: https://github.com/Fortran-FOSS-Programmers/ford#readme
[FORD wiki]: https://github.com/Fortran-FOSS-Programmers/ford/wiki
[EditorConfig]: https://editorconfig.org/
[findent]: https://www.ratrabbit.nl/ratrabbit/content/findent/introduction
[CMake]: https://cmake.org/
[CTest]: https://cmake.org/cmake/help/latest/manual/ctest.1.html
[Getting Started]: ./getting-started.html
[Prerequisites]: ./getting-started.html#prerequisites
[Configuring the Library]: ./getting-started.html#configuring-the-library
[Building the Library]: ./getting-started.html#building-the-library
[Running Tests]: ./getting-started.html#running-tests
[Documenting New and Modified Code]: ./using-ford.html
[Style Guide]: ./style-guide.html
[Example Usage]: ./examples/
