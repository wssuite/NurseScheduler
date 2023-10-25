# INSTALLATION GUIDE

The nurse rotation scheduler requires two external libraries. Before compiling the code, COIN BCP must be installed and 
the corresponding environment variables must be set or pkg-config may be used. The submodule coin might be used to install BCP.
This document provides a manual to install the libraries and compile the code.

You can run the INSTALL.sh script for performing this first step.

## 1. COIN-OR BCP

  ### 1.a. Use the INSTALL.sh script
  To compile the source from the repository in the submodule. It's a modified BCP version.
  Please refer to the README.md of the submodule once initialized with: 
  ``git submodule update --init --recursive``

  ### 1.b. Download the latest versions of BCP url 
  Our published results have been found with BCP-1.4: http://www.coin-or.org/download/source/Bcp. Then, extract the sources in a directory of your choice /path/to/coinor/directory/.

  Or, if svn is available, download directly the latest stable versions using the following commands in the terminal:
  ````bash
  cd  /path/to/coinor/directory/
  svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
  ````
  #### Build BCP (https://projects.coin-or.org/BuildTools/wiki/downloadUnix):
  ````bash
  cd /path/to/coinor/directory/Bcp-1.4
  mkdir build
  cd build
  ../configure --enable-shared=yes --enable-static=yes --with-cpx=no
  make
  make test
  make install
  ````
  We add the option --with-cpx=no because some incompatibilities have been noticed between CPLEX 12.7+ and BCP 1.4. The option can be withdrawn if no recent version of CPLEX is installed.

  On Windows, BCP can be installed using the gcc distribution packaged in MSYS2 (https://www.msys2.org/).

  ##### For developers: if you need to debug the code and enter BCP.
  Build BCP once again with the following commands:
  ````bash
  cd /path/to/coinor/directory/Bcp-1.4
  mkdir debug
  cd debug
  ../configure --enable-debug --enable-shared=yes --enable-static=yes --with-cpx=no
  make
  make test
  make install
  ````

## 2. Boost Graph library (OPTIONAL)

1. Download Boost Graph Library (BGL) at http://sourceforge.net/projects/boost/files/latest/download?source=files.
   At least, boost 1.58 is needed.

2. The part of the BGL we use is in the header only library, so there is no need to build anything. Uncompress the archive and move the folder to any directory /path/to/boost/directory/boost_X_XX (where X_XX stands for the version of boost you downloaded).

## 3. Build the nurse scheduling code:
The code can be built with cmake. However, the paths to BCP need to be provided if you do not use pkg_config.
You can either set the variable when running cmake: ``cmake -DBCPDIROPT=/path/to/coinor/directory/Bcp-1.4/build``.

Or, you can specify the values by creating at the root the file CMakeDefinitionsLists.txt to provide your own definitions (you need to define only one of the last variables for BCP):
  ````bash
  set(BCPDIROPT /path/to/coinor/directory/Bcp-1.4/build)
  set(BCPDIRDBG /path/to/coinor/directory/Bcp-1.4/debug)
  ````
You can also use the flag MBCP with the variable USE_MBCP ``set(USE_MBCP True)`` 
if you installed the version in the submodule. There is a small enhancement that allows you to pause and restart BCP 
when solving different problems (especially useful when dividing the problems into several components).

Also, if boost is installed in a classical path, cmake will be able to find the path automatically. 
You need also to define the variable USE_BOOST to compile with it, and you may provide the path to the boost library.
````bash
  set(USE_BOOST True)
  set(BOOST_DIR /path/to/boost/directory/boost_X_XX)
````

#### Finally, go to the root directory of the nurse scheduling project and type:
````bash
mkdir build && cd build && cmake ..
make
````
You can use ``cmake -DCMAKE_BUILD_TYPE=Debug ..`` to compile in debug.

## 4. Mac M1 processors
If you are using a M1 processor, you will have some issues.
Unfortunately, we are not able to compile BCP on an arm architecture. 
Therefore, BCP must be compiled on a x86_64 architecture: 
you will just need to set the right compiler by default in your environment (check ``gcc -v``).
However, cmake will use by default the system architecture, and thus the x86_64.cmake 
configuration file must be given to cmake with the command "-DCMAKE_TOOLCHAIN_FILE=x86_64.cmake".

## 5. Optional flags

### External solvers
If you have compiled BCP with Cplex or Gurobi, you need to provide the right paths to CMake.
Currently, Gurobi is handled and tested. Add with the right values:
```bash
set(GUROBI_INCLUDE_DIR /Library/gurobi901/mac64/include)
set(GUROBI_LIB gurobi)
```

### ASAN
When developing and debugging using an address sanitizer may be very handy, add this line to activate asan:
```bash
SET(ASAN True)
```

### Control
If you add the flag CTR (``SET(CTR True)``), the code will solve the pricing problem 
with different approaches if enable to ensure the found solutions are the ight ones.
