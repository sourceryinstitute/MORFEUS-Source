<div align="center">

# MORFEUS

**M**ulti-physics **O**bject-oriented **R**econfigurable **F**luid **E**nvironment for **U**nified **S**imulations

</div>

<details><summary><b>Table of Contents</b></summary>
<p>

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [MORFEUS](#morfeus)
    - [What is MORFEUS?](#what-is-morfeus)
        - [Capabilities](#capabilities)
        - [Compatibility](#compatibility)
    - [Repository Mirror and Source](#repository-mirror-and-source)
        - [Mirroring and Updating Mechanism](#mirroring-and-updating-mechanism)
        - [Development Workflow](#development-workflow)
    - [Developer Contributing Checklist](#developer-contributing-checklist)

<!-- markdown-toc end -->

</p>
</details>

## What is MORFEUS?

### Capabilities

MORFEUS is an open source framework for the solution of partial differential equations
(PDEs) written in modern Fortran. It is object oriented, and attempts to provide usefull
abstractions for the solution of PDEs that are easy to use and performant.

MORFEUS consists of two main solution approaches:

  - Finite Volume ([FV]) which can handle very complex geometris and grids
  - Finite Difference ([FD]) which is less mature and under active development, but provides
    some base objects and functionality consumed by [FV]. The [FD] solver operates on
	block-structured grids.

### Compatibility

| Compiler  | Linux | macOS | Windows (64bit) |
|----------:|:-----:|:-----:|:---------------:|
| GCC 8     |  ✔    |   ✔  |   ✘             |
| Intel 18  |  ✔    |   ✘  |   ✔             |


## Repository Mirror and Source

:warning: This repository is automatically mirrored. Ensure that all
development happens from [the source repository][SI source repository]!
The [Sourcery Institute Source Repository][SI source repository] is
located at:

https://github.com/sourceryinstitute/MORFEUS-Source


__YOU SHOULD NOT COMMIT TO OR OTHERWISE MODIFY THE MIRROR
REPOSITORY!!!__ You should not have write access to the repository to
start with, so it should be hard to push a commit there
accidentally. In addition, any wiki edits made to the Sourcery
Institute source repository will automatically be forwarded to the
mirrored wiki.

### Mirroring and Updating Mechanism

The mirror is updated automatically whenever commits are pushed to
the SI source repository using [GitHub Actions]. It may take a few
minutes for the mirror to update. If it takes any longer than that,
then please check the [status of the GitHub Actions], and if there's a
problem fix it, or notify [@zbeekman].


### Development Workflow

When developing [MORFEUS][SI Source Repository] as it is consumed
by a parent project using a mirror, ensure that you set the submodule remote to
point to this repository, *NOT* a mirror repository:

``` bash
git remote -v # check first that your remote points to the mirror or
              # source
git remote set-url origin git@github.com:sourceryinstitute/MORFEUS-Source.git
```

This will ensure that you are pulling the same commits and branches
that the mirror repositories will be using, while your pushed branches
& commits will go to the [SI source repository].

The update cycle of the repositories looks something like this:

```
                <mirror repo>
                         __
                        /\
                          \
                           \
                            \ (GitHub Actions mirror)
                             \
                              \
                               \
                                \
<local-repo> <------------------> <SI source>
            (git push/pull origin)
```

## Developer Contributing Checklist

Here is a developer checklist to remind current developers and new developers
of steps they must take to ensure they are correctly contributing to MORFEUS.

  - [ ] Ensure that the `origin` remote of a larger package they are working on points to
        [this repository][SI source repository] and not a mirror. (See [above][set origin].)
  - [ ] Ensure that all new files added have the appropriate copyright statement for any
        contract work.
  - [ ] Ensure that you have [setup your editor] to respect the projects [style guidelines],
        which are partially codified via [EditorConfig].
  - [ ] Ensure that you are [documenting your source code] with appropriate [FORD] comments.
  - [ ] Ensure that you are following the [style guidelines].
  - [ ] Double check that you are not adding sensitive or proprietary information.
  - [ ] Double check that you did not overwrite a [git submodule] with its contents.


[FV]: https://github.com/sourceryinstitute/MORFEUS-Source/tree/master/src/FV
[FD]: https://github.com/sourceryinstitute/MORFEUS-Source/tree/master/src/FD
[@zbeekman]: https://github.com/zbeekman
[GitHub Actions]: https://github.com/features/actions
[status of the GitHub Actions]: https://github.com/sourceryinstitute/MORFEUS-Source/actions
[SI source repository]: https://github.com/sourceryinstitute/MORFEUS-Source
[set origin]: https://github.com/sourceryinstitute/MORFEUS-Source/#development-workflow
[setup your editor]: https://editorconfig.org/#download
[style guidelines]: https://github.com/sourceryinstitute/MORFEUS-Source/wiki/Code-Style-Guidelines
[EditorConfig]: https://editorconfig.org
[documenting your source code]: https://github.com/sourceryinstitute/MORFEUS-Source/wiki/Documenting-your-code-with-FORD
[FORD]: https://github.com/Fortran-FOSS-Programmers/ford
[git submodule]: https://git-scm.com/book/en/v2/Git-Tools-Submodules
