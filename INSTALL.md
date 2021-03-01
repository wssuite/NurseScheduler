# INSTALLATION GUIDE

The nurse rotation scheduler requires two external libraries. Before compiling the code, BOOST and COIN BCP  must be installed and the corresponding environment variables must be set.
This document provides a manual to install the libraries and compile the code.

You can run the INSTALL.sh script for performing the first two steps (give as first argument the directory where you wish to install the software, otherwise the current directory will be used).

#### 1. Boost Graph library

  1. Download Boost Graph Library (BGL) at http://sourceforge.net/projects/boost/files/latest/download?source=files.
  At least, boost 1.58 is needed.

  2. The part of the BGL we use is in the header only library, so there is no need to build anything. Uncompress the archive and move the folder to any directory /path/to/boost/directory/boost_X_XX (where X_XX stands for the version of boost you downloaded).


#### 2. COIN-OR BCP

  #####- Download the latest versions of BCP url 
  Our published results have been found with BCP-1.4: http://www.coin-or.org/download/source/Bcp. Then, extract the sources in a directory of your choice /path/to/coinor/directory/.

  Or, if svn is available, download directly the latest stable versions using the following commands in the terminal:
  ````bash
  cd  /path/to/coinor/directory/
  svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
  ````
  #####- Build BCP (https://projects.coin-or.org/BuildTools/wiki/downloadUnix):
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

  #####- For developpers: if you need to debug the code and enter BCP.
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

#### 3. Build the nurse scheduling code:
The code can be built either with a traditional Makefile or with cmake. In both cases, the paths to boost and BCP need to be provided.

  - ##### From the Makefile:
  First, add the following environment variables (you need to define only one of the last variables for BCP):
  ````bash
  export BOOST_DIR="/path/to/boost/directory/boost_X_XX"
  export BCPDIROPT="/path/to/coinor/directory/Bcp-1.4/build"
  export BCPDIRDBG="/path/to/coinor/directory/Bcp-1.4/debug"
  ````
  Then, go to the root directory of the nurse scheduling project and type:
  ````bash
  make
  ````

  - ##### From cmake:
  First, create at the root the file CMakeDefinitionsLists.txt to provide your own paths (you need to define only one of the last variables for BCP). Also, if boost is installed in a classical path, cmake will be able to find the path automatically. 
  ````bash
  set(BOOST_DIR /path/to/boost/directory/boost_X_XX)
  set(BCPDIROPT /path/to/coinor/directory/Bcp-1.4/build)
  set(BCPDIRDBG /path/to/coinor/directory/Bcp-1.4/debug)
  ````
  
  Then, go to the root directory of the nurse scheduling project and type:
  ````bash
  mkdir build && cd build && cmake ..
  make
  ````
