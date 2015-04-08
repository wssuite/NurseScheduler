++------------------------------------------------------++
++------------------------------------------------------++
||                                                      ||
||  README :                                            ||
||                                                      ||
||       INSTALLATION DES OUTILS DE COIN_OR             ||
||             POUR ROSTER DES NURSES   	              ||
||                                                      ||
++------------------------------------------------------++
++------------------------------------------------------++



I. Compiler les sources
  |
  | 1a.Télécharger BCP et CBC aux adresses suivantes :
  |
       http://www.coin-or.org/download/source/Bcp
       http://www.coin-or.org/download/source/Cbc

  |  Puis extraire dans le répertoire d'installation.

  | 1b. Ou, mieux, utiliser svn pour télécharger directement les dernières versions stables : 

      cd /le/repertoire/dinstallation/des/librairies
      svn co https://projects.coin-or.org/svn/Cbc/stable/2.9 Cbc-2.9
      svn co https://projects.coin-or.org/svn/Bcp/stable/1.4 Bcp-1.4
  |
  | 3. Compiler CBC et BCP (https://projects.coin-or.org/BuildTools/wiki/downloadUnix):
	cd /le/repertoire/dinstallation/des/librairies/Bcp-1.4
	mkdir build
	cd build
	../configure
	make
	make test
	make install

  cd /le/repertoire/dinstallation/des/librairies/Cbc-2.9
  mkdir build
  cd build
  ../configure
  make
  make test
  make install

  II. Si besoin de DEBUGER
  |
  | 1. Recompiler CBC et BCP avec la commande :
  cd /le/repertoire/dinstallation/des/librairies/Bcp-1.4
  mkdir debug
  cd debug
  ../configure --enable-debug
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
  +--------------------------------------------------------------------------------------


III. Modifier l'environnement
  |
  | 1. Ajouter les lignes suivantes au .bashrc (ou .profile, c'est selon) :
  |
       export BCPDIR=/le/repertoire/dinstallation/des/librairies/Bcp-1.4/build
       export CBCDIR=/le/repertoire/dinstallation/des/librairies/Cbc-2.9/build

  | 2. Si on veut passer en mode debug, simplement changer ces  variables d'environnement :

       export BCPDIR=/le/repertoire/dinstallation/des/librairies/Bcp-1.4/debug
       export CBCDIR=/le/repertoire/dinstallation/des/librairies/Cbc-2.9/debug
  |
  +--------------------------------------------------------------------------------------


IV. Le fichier Makefile dans /le/repertoire/de/Git/RosterDesNurses a été écrit à partir des exemples de COIN pour fonctionner avec ces réglages. Donc juste faire make pour compiler le code.

  |
  +--------------------------------------------------------------------------------------
