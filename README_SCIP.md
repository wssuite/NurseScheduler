++------------------------------------------------------++
++------------------------------------------------------++
||                                                      ||
||  README :                                            ||
||                                                      ||
||       INSTALLATION DE SCIP POUR ROSTER DES NURSES   ||
||                                                      ||
++------------------------------------------------------++
++------------------------------------------------------++



I. SCIP source
  |
  | 1.Télécharger SCIP Optimization Suite à l'adresse suivante:
  |
       http://scip.zib.de/download.php?fname=scipoptsuite-3.1.1.tgz
  |
  |
  | 2. L'extraire dans votre répertoire favori (adresse : /home/local/)
  |

  |
  |
  | 3. Compiler SCIP:
	cd /home/local/scipoptsuite-3.1.1/
	make
  | Si make ne fonctionne pas, essayer:
	make LEGACY=true SPX_LEGACY=true

  |
  |
  | 4. Faites les tests:
	make test
  |
  +--------------------------------------------------------------------------------------


II. Modifier l'environnement
  |
  | 1. Ajouter les lignes suivantes au .bashrc (ou .profile, c'est selon) :
  |
       SCIPDIR=/local/scipoptsuite-3.1.1/scip-3.1.1
       export SCIPDIR
  |
  +--------------------------------------------------------------------------------------


III. Utiliser le Makefile donné par SCIP et modifié pour notre projet (il est à la racine de Git)

IV. Si besoin de DEBUGER
  |
  | 1. Recompiler SCIP avec la commande :
  |       make OPT=dbg
  |
  | 2. Compiler le projet avec l'option OPT=dbg (ex: make OPT=dbg)
  |
  +--------------------------------------------------------------------------------------

V. Me dire si ça merde.




