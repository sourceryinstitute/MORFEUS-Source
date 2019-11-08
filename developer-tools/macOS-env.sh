#!/usr/bin/env bash
#
#
#
# Source this script on macOS to setup environment for compiling MORFEUS and its prerequisites
#


# Unlink GCC > 8 if linked
brew unlink gcc || true

# Unlink OpenMPI if linked
brew unlink openmpi || true

# Link MPICH if not linked
brew link mpich || true

# Setup all the compilers
# Find the GCC@8 opt directory
_gcc8_loc="$(brew --prefix)/opt/gcc@8/bin"

# Set compilers
export FC="${_gcc8_loc}/gfortran-8"
export CC="${_gcc8_loc}/gcc-8"
export CXX="${_gcc8_loc}/g++-8"

# Set compilers to use for MPICH
export MPICH_FC="${FC}"
export MPICH_CC="${CC}"
export MPICH_CXX="${CXX}"

ccache_opt="$(brew --prefix)/opt/ccache/libexec"

if [[ -x "${ccache_opt}/gcc-8" && -x "${ccache_opt}/g++-8" ]] ; then
  export CC="${ccache_opt}/gcc-8"
  export CXX="${ccache_opt}/g++-8"
fi
