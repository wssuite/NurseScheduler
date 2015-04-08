#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic Licence.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for C++ vehicle routing example using SCIP for branch-and-price
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Andreas Bley
#@author  Marc Pfetsch


#-----------------------------------------------------------------------------
# paths
#-----------------------------------------------------------------------------
#Should be define in bashrc or profile
#SCIPDIR         =       ../..

#-----------------------------------------------------------------------------
# add user flags
#-----------------------------------------------------------------------------
INCLUDESFLAGS	=	-I$(BOOST_DIR) -I$(CBCDIR)/include/coin

USRFLAGS	=	
USROFLAGS	=
USRCFLAGS	=
USRCXXFLAGS	=	-g -O0 -w -fPIC -fexceptions -DNDEBUG -DIL_STD  $(INCLUDESFLAGS)
USRLDFLAGS	= 	-g -O0
OS = $(shell uname -s)
ifeq ($(OS),Linux)
	USRLDFLAGS += -lrt
endif

USRARFLAGS	=
USRDFLAGS	=
LIBS		=



#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
include $(SCIPDIR)/make/make.project

#-----------------------------------------------------------------------------
# include project Makefile from BCP
#-----------------------------------------------------------------------------
include make.coin

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

MAINNAME	=	roster
MAINOBJ		=	main.o main_test.o MyTools.o Demand.o Nurse.o Scenario.o ReadWrite.o Roster.o MasterProblem.o SubProblem.o Solver.o Greedy.o RotationPricer.o ScipModeler.o BcpModeler.o CbcModeler.o
MAINSRC		=	$(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))
MAINDEP		=	$(SRCDIR)/depend.cppmain.$(OPT)


MAIN		=	$(MAINNAME).$(BASE).$(LPS)$(EXEEXTENSION)
MAINFILE	=	$(BINDIR)/$(MAIN)
MAINSHORTLINK	=	$(BINDIR)/$(MAINNAME)
MAINOBJFILES	=	$(addprefix $(OBJDIR)/,$(MAINOBJ))

# <<<<<<< HEAD
# =======
# INCLUDESFLAGS = -I$(BOOST_DIR)
# CXXFLAGS	+= -g -O0 -m64 -w -fPIC -fexceptions  -DIL_STD -std=c++98 $(INCLUDESFLAGS)

# >>>>>>> Fin-du-greedy-et-ouputs
#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------

ifeq ($(VERBOSE),false)
.SILENT:	$(MAINFILE) $(MAINOBJFILES) $(MAINSHORTLINK)
endif

.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

.PHONY: lint
lint:		$(MAINSRC)
		-rm -f lint.out
		$(SHELL) -ec 'for i in $^; \
			do \
			echo $$i; \
			$(LINT) $(SCIPDIR)/lint/scip.lnt +os\(lint.out\) -u -zero \
# =======
# 			$(LINT) -I$(SCIPDIR) lint/main-gcc.lnt +os\(lint.out\) -u -zero \
# >>>>>>> Fin-du-greedy-et-ouputs
			$(FLAGS) -UNDEBUG -UWITH_READLINE -UROUNDING_FE $$i; \
			done'

.PHONY: scip
scip:
		@$(MAKE) -C $(SCIPDIR) libs $^

.PHONY: doc
doc:
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scip.css);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/scippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/pictures/miniscippy.png);
		@-(cd doc && ln -fs ../$(SCIPDIR)/doc/scipfooter.html footer.html);
		cd doc; $(DOXY) $(MAINNAME).dxy

$(MAINSHORTLINK):	$(MAINFILE)
		@rm -f $@
		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/*.o
		-rmdir $(OBJDIR)
endif
		rm -f $(BINDIR)/$(MAINNAME)

.PHONY: test
test:           $(MAINFILE)
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/evalcheck_cluster.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/check.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/getlastprob.awk);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_set.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_logfiles.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/configuration_tmpfile_setup_scip.sh);
		@-(cd check && ln -fs ../$(SCIPDIR)/check/run.sh);
		cd check; \
		$(SHELL) ./check.sh $(TEST) $(MAINFILE) $(SETTINGS) $(notdir $(MAINFILE)) $(TIME) $(NODES) $(MEM) $(THREADS) $(FEASTOL) $(DISPFREQ) $(CONTINUE) $(LOCK) "example" $(LPS) $(VALGRIND) $(CLIENTTMPDIR) $(OPTCOMMAND);

.PHONY: depend
depend:		$(SCIPDIR)
		$(SHELL) -ec '$(DCXX) $(FLAGS) $(DFLAGS) $(MAINSRC) \
		| sed '\''s|^\([0-9A-Za-z\_]\{1,\}\)\.o *: *$(SRCDIR)/\([0-9A-Za-z\_]*\).cpp|$$\(OBJDIR\)/\2.o: $(SRCDIR)/\2.cpp|g'\'' \
		>$(MAINDEP)'

-include	$(MAINDEP)

$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
				@echo 		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
		$(LIBS)	\
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@
		$(LINKCXX) $(MAINOBJFILES) \
		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
                $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
		$(LIBS)	\
                $(OFLAGS) $(LPSLDFLAGS) \
		$(LDFLAGS) $(LINKCXX_o)$@
# =======
# 		$(CXXFLAGS) \
# 		$(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) \
#                 $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX) \
#                 $(OFLAGS) $(LPSLDFLAGS) \
# 		$(LDFLAGS) $(LINKCXX_o)$@ 
# >>>>>>> Fin-du-greedy-et-ouputs

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $@"
		@echo $(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
