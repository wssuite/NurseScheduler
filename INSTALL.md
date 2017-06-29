++------------------------------------------------------++
++------------------------------------------------------++
||                                                      ||
||                                                      ||
||                                                      ||
||    INSTALL GUIDE FOR THE STATIC NURSE SCHEDULER      ||
||                                                      ||
++------------------------------------------------------++
++------------------------------------------------------++

The nurse rotation scheduler requires two external libraries. Before compiling the code, BOOST and COIN BCP  must be installed and the corresponding environment variables must be set. 
This document provides a manual to install the libraries and compile the code.

I. Boost Graph library
  |
  | 1.Download Boost Graph Library (BGL) at 
  |
  |   http://sourceforge.net/projects/boost/files/latest/download?source=files
  |
  |
  | 2. The part of the BGL we use is in the header only library, so there is no need to build anything. 
  |    Uncompress the archive and move the folder to any directory /path/to/boost/directory/boost_X_XX (where X_XX stands for the version of boost you downloaded)
  |
  | 3. Add the following environment variable 
  |
       BOOST_DIR="/path/to/boost/directory/boost_X_XX"
       export BOOST_DIR


II. COIN-OR BCP
  |
  | 1a. Download the latest versions of BCP url (our published results have been found with BCP-1.4):
  |
       http://www.coin-or.org/download/source/Bcp

  |  And extract the sources in a directory of your choice /path/to/coinor/directory/

  | 1b. Or, if svn is available, download directly the latest stable versions using the following commands in the terminal:
      cd  /path/to/coinor/directory/
      svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
  |
  | 2. Build BCP (https://projects.coin-or.org/BuildTools/wiki/downloadUnix):

  cd /path/to/coinor/directory/Bcp-1.4
  mkdir build
  cd build
  ../configure --enable-shared=no --enable-static=yes --with-cpx=no
  make
  make test
  make install

  | We add the option --with-cpx=no because some incompatibilities have been noticed between CPLEX 12.7+ and BCP 1.4. The option can be withdrawn if no recent version of CPLEX is installed. 

  |3. For developpers: If you need to debug the code and enter BCP 
  |
  | Build BCP once again with the following commands:
  cd /path/to/coinor/directory/Bcp-1.4
  mkdir debug
  cd debug
  ../configure --enable-debug --enable-shared=no --enable-static=yes --with-cpx=no
  make
  make test
  make install
  |
  | 4. Modify the environment
       export BCPDIROPT=/path/to/coinor/directory/Bcp-1.4/build
       export BCPDIRDBG=/path/to/coinor/directory/Bcp-1.4/debug
  |
  +--------------------------------------------------------------------------------------


III. Build the static nurse scheduling code: 
  | Go to the root directory of the nurse scheduling project and type make 
