<div align="center">

# MORFEUS

**M**ulti-physics **O**bject-oriented ...
</div>

<!-- toc -->

- [Repository Mirror and Source](#repository-mirror-and-source)
  * [NRC Mirror](#nrc-mirror)
  * [Sourcery Institute Source Repository](#sourcery-institute-source-repository)
  * [Mirroring and Updating Mechanism](#mirroring-and-updating-mechanism)
- [Workflow and Roadmap](#workflow-and-roadmap)
  * [Workflow](#workflow)
  * [Roadmap](#roadmap)

<!-- tocstop -->

## Repository Mirror and Source

### [NRC Mirror]
This repository is automatically mirrored to
https://github.com/nrc-fuels/MORFEUS-mirror so that NRC staff always
have access to the source code while the repository is
private. Eventually the respository will be made open source, and this
will not be an issue, although the mirror will likely continue to
exist to ensure continuity of access.

__YOU SHOULD NOT COMMIT TO OR OTHERWISE MODIFY THE NRC MIRROR
REPOSITORY!!!__ You should not have write access to the repository to
start with, so it should be hard to push a commit there
accidentally. In addition, any wiki edits made to the Sourcery
Institute source repository will automatically be forwarded to the NRC
mirror wiki.

### [Sourcery Institute Source Repository][SI source repository]
Sourcery Institute and GSE staff working on the project have access to
the (for now) private MORFEUS source repsitory at
https://github.com/sourceryinstitute/MORFEUS-Source. The plan is to
eventually open source this, once appropriate checks have been made
and the source is deemed ready for public scrutiny.

__NRC Staff *DO NOT* have access to the source repository__ (until it
is open-sourced). Therefore, in the short term, they will not be able
to see any issues, pull requests, etc. on that repository. __However,
*NRC STAFF ALONG WITH THE GENERAL PUBLIC WILL LIKELY BE ABLE TO SEE
THIS CONTENT AT A LATER DATE!*__ Any conversation in issues and PRs
should be considered public, but temporarily inaccessible.

### Mirroring and Updating Mechanism
The NRC mirror is updated automatically whenever commits are pushed to
the SI source repository using [GitHub Actions]. It may take a few
minutes for the mirror to update. If it takes any longer than that,
then please check the [status of the GitHub Actions], and if there's a
problem fix it, or notify @zbeekman.


## Workflow and Roadmap

### Workflow
Until MORFEUS is made open source, the private NRC repositories will be pulling
in MORFEUS as a submodule from the [NRC mirror]. This means that, in
general, you should pull from the [NRC mirror], but when developing,
you will need to push to the [SI source repository]. To ensure that
you have git setup correctly, first clone the repository, (or navigate
to the git submodule directory, `submodules/morfeus`) and run the
following commands:

``` bash
git remote -v # check first that your remote points to the mirror or
              # source
git remote set-url --push origin git@github.com:sourceryinstitute/MORFEUS-Source.git
git remote set-url --pull origin git@github.com:nrc-fuels/MORFEUS-mirror.git
```

This will ensure that you are pulling the same commits and branches
that nrc-fuels repositories will be using, while your pushed branches
& commits will go to the [SI source repository]. Once MORFEUS is open
sourced, then everything can be set to point to the [SI source
repository].

The update cycle of the repositories looks something like this:

```
                <NRC mirror>
                   /     __
                  /     /\
(git pull origin)/        \
                /          \
               /            \ (GitHub Actions mirror)
              /              \
             /                \
            /                  \
          \/_                   \
<local-repo> ------------------> <SI source>
              (git push orign)
```

### Roadmap
Intagrating MORFEUS with the existing code will happen in a two phase
process before MORFEUS is open sourced.

1. Directories `src/morfeus_{fd,fv,fe}` will be deleted in the other
   code, and CMake will be adjusted so that
   `submodules/morfues/src/{FD,FV,FE}` will be used as the existing
   source.
1. Next, MORFEUS will be given its own build system, and submodules
   for the necessary TPLs, and it will be built via a CMake Super
   Build process alongside other TPLs using its native build
   system. The other code(s) will consume morfeus as an External
   Project.

After this work has been completed, so that MORFEUS can stand on its
own, and checks for sensitive info are complete, and the authorization
to open source the software has been confirmed, then MORFEUS will be
ready to be open sourced, and consumed by the general public.

[GitHub Actions]: https://github.com/features/actions
[status of the GitHub Actions]: https://github.com/sourceryinstitute/MORFEUS-Source/actions
[NRC Mirror]: https://github.com/nrc-fuels/MORFEUS-mirror
[SI source repository]: https://github.com/sourceryinstitute/MORFEUS-Source
