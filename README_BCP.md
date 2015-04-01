++------------------------------------------------------++
++------------------------------------------------------++
||                                                      ||
||  README :                                            ||
||                                                      ||
||       INSTALLATION DE BCP POUR ROSTER DES NURSES   	||
||                                                      ||
++------------------------------------------------------++
++------------------------------------------------------++



I. SCIP source
  |
  | 1.Télécharger BCP à l'adresse suivante:
  |
       http://www.coin-or.org/download/source/Bcp/
  |
  |
  | 2. L'extraire dans votre répertoire favori (adresse : /home/local/)
  |
  |
  | 3. Compiler BCP (https://projects.coin-or.org/BuildTools/wiki/downloadUnix):
	cd /home/local/Bcp-1.4.0/
	mkdir build
	cd build
	../configure
	make
	make test
	make install
  |
  +--------------------------------------------------------------------------------------


II. Modifier l'environnement
  |
  | 1. Ajouter les lignes suivantes au .bashrc (ou .profile, c'est selon) :
  |
       BCPDIR=/local/Bcp-1.4.0/build
       export BCPDIR
  |
  +--------------------------------------------------------------------------------------


III. Utiliser le Makefile donné par BCP et modifié pour notre projet (il est à la racine de Git)

IV. Si besoin de DEBUGER
  |
  | 1. Recompiler BCP avec la commande :
	make clean
          ../configure --enable-debug
	make
	make install
  |
  +--------------------------------------------------------------------------------------
