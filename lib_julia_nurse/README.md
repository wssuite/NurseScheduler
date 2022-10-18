Install
------------------
### Parameters
We use the package CxxWrap.jl (and its C++ interface libcxxwrap_julia) to write this wrapper. This last library requires C++17 to be build.

### Julia
CxxWrap.jl
through the Julia package manager
### C++
libcxxwrap_julia as advised on the git repo.

Build
----------------------
In the folder julia-wrap, create a CMakeDefinitionsLists.txt file with the following variables:
- JlCxx_DIR to libcxxwrap-julia-build.
- NURSE_DIR to pNurseScheduler.
- NURSE_LIB to pNurseScheduler/lib or wherever the library file is.