++------------------------------------------------------++
++------------------------------------------------------++
||                                                      ||
||  README :                                            ||
||                                                      ||
||       INSTALLATION DE BOOST POUR ROSTER DES NURSES   ||
||                                                      ||
++------------------------------------------------------++
++------------------------------------------------------++

The nurse rotation scheduler requires three external libraries. Before compiling the code, BOOST, COIN BCP and COIN CBC must be installed and the corresponding environment variables must be set. 
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
  |
  +--------------------------------------------------------------------------------------


II. Modifier l'environnement
  |
  | 1. Ajouter les lignes suivantes au .bashrc (ou .profile, c'est selon) :
  |
       BOOST_DIR="/mon/adresse/boost_1_57_0"
       export BOOST_DIR
  |
  +--------------------------------------------------------------------------------------


III. Modifier votre environnement de développement (ça c'est un peu le problème de chacun, il suffit de lui ajouter le chemin que vous avez choisi...)


IV. Me dire si ça merde.
++------------------------------------------------------++
++------------------------------------------------------++
||                                                      ||
||  README :                                            ||
||                                                      ||
||       INSTALLATION DES OUTILS DE COIN_OR             ||
||             POUR ROSTER DES NURSES                   ||
||                                                      ||
++------------------------------------------------------++
++------------------------------------------------------++



I. Build the source codes of COIN-OR BCP
  |
  | 1a. Download the latest versions of BCP url (our published results have been found with BCP-1.4):
  |
       http://www.coin-or.org/download/source/Bcp

  |  And extract the sources in a directory of your choice /path/to/coinor/directory/

  | 1b. Or, if svn is available, download directly the latest stable versions using the following commands in the terminal:
      cd  /path/to/coinor/directory/
      svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
  |
  | 2. Compiler BCP (https://projects.coin-or.org/BuildTools/wiki/downloadUnix):

  cd /path/to/coinor/directory/Bcp-1.4
  mkdir build
  cd build
  ../configure --enable-shared=no --enable-static=yes --with-cpx=no
  make
  make test
  make install

  | We add the option --with-cpx=no because some incompatibilities have been noticed between CPLEX 12.7+ and BCP 1.4. The option can be withdrawn if no recent version of CPLEX is installed. 

  II. Si besoin de DEBUGER
  |
  | 1. Recompiler CBC et BCP avec la commande :
  cd /le/repertoire/dinstallation/des/librairies/Bcp-1.4
  mkdir debug
  cd debug
  ../configure --enable-debug --enable-shared=no --enable-static=yes
  make
  make test
  make install

  cd /le/repertoire/dinstallation/des/librairies/Cbc-2.9
  mkdir build
  cd build
  ../configure --enable-debug
  make
  make test
  make install
  |
  III. Pour utiliser CPLEX avec BCP
  |
  | 1. Dans le .bashrc ou .profile sur Mac, ajouter les chemins suivants :
  |
  |  export CPLEX_HOME="/repertoire/dinstallation/dilog/CPLEX_Studio126/cplex"
  |  export CPXINCDIR="${CPLEX_HOME}/include/ilcplex"
  |  export CPXLIBDIR="${CPLEX_HOME}/lib/x86-64_linux/static_pic"
  |  export CPXLIB="-L${CPLEX_HOME}/lib/x86-64_linux/static_pic -lcplex -lm -lpthread -DIL_STD"
  |  (attention à remplacer x86-64_linux par le bon dossier selon l'OS d'installation, c'est x86-64_osx sur mac)

  | 2. Relancer toute la configuration en ajoutant l'option --disable-cplex-libcheck à la configuration

  | 3. Compiler/Installer
  |
  IV. Pour utiliser GUROBI avec BCP (cumulable avec CPLEX pour faire des tests)
  |
  | 1. Télécharger Gurobi sur le site : http://www.gurobi.com/downloads/download-center)
  |    Pour Linux, il s'agit de gurobi6.5.0_linux64.tar.gz
  |
  | 2. Installer Gurobi sur votre ordi (dépend de l'OS).
  |
  |    Pour Linux, les étapes d'installation sont détaillées ici :
  |    https://www.gurobi.com/documentation/6.5/quickstart_linux/software_installation_guid.html#section:Installation
  |
  |    En résumé (sur linux... pour mac je ne sais pas), il s'agit de
  |
  |      (i)  décompresser le fichier (typiquement, dans le répertoire /opt)
  |           tar xvfz gurobi6.5.0_linux64.tar.gz
  |
  | 3. Dans le .bashrc ou .profile sur Mac, ajouter les chemins suivants :
  |
  | Sur Linux :
  |  export GUROBI_HOME="/opt/gurobi650/linux64"
  |  export PATH="${PATH}:${GUROBI_HOME}/bin"
  |  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
  |  export GRBINCDIR="${GUROBI_HOME}/include"
  |  export GRBLIB="${GUROBI_HOME}/lib/libgurobi.so.6.5.0"
  |
  | Sur Mac :
  |  export GUROBI_HOME="/Library/gurobi650/mac64"
  |  export PATH="${PATH}:${GUROBI_HOME}/bin"
  |  export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
  |  export GRBINCDIR="${GUROBI_HOME}/include"
  |  export GRBLIB="${GUROBI_HOME}/lib/libgurobi65.so"

  | 4. Relancer toute la configuration

  | 5. Compiler/Installer
  +--------------------------------------------------------------------------------------


V. Modifier l'environnement
  |
  | 1. Ajouter les lignes suivantes au .bashrc (ou .profile, c'est selon) :
  |
       export BCPDIROPT=/le/repertoire/dinstallation/des/librairies/Bcp-1.4/build
       export CBCDIROPT=/le/repertoire/dinstallation/des/librairies/Cbc-2.9/build

  | 2. Si on veut passer en mode debug, simplement changer ces  variables d'environnement :

       export BCPDIRDBG=/le/repertoire/dinstallation/des/librairies/Bcp-1.4/debug
       export CBCDIRDBG=/le/repertoire/dinstallation/des/librairies/Cbc-2.9/debug
  |
  +--------------------------------------------------------------------------------------


IV. Le fichier Makefile dans /le/repertoire/de/Git/RosterDesNurses a été écrit à partir des exemples de COIN pour fonctionner avec ces réglages. Donc juste faire make pour compiler le code.
