#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program 
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#@file    Makefile
#@brief   Makefile for C++ nurse rostering branch-and-price
#@author  Antoine Legrain
#@author  Jeremy Omer
#@author  Samuel Rosat


#-----------------------------------------------------------------------------
# OPTIONS
#-----------------------------------------------------------------------------
USE_SCIP = FALSE
USE_BCP = TRUE
USE_CBC = FALSE
DEBUG  = FALSE

#-----------------------------------------------------------------------------
# default flags
#-----------------------------------------------------------------------------
INCLUDESFLAGS =
FLAGS =  
CFLAGS   =
CXXFLAGS=
LDFLAGS  =  
LINKCXX =
LIBS  =

#-----------------------------------------------------------------------------
# include default project Makefile from SCIP
#-----------------------------------------------------------------------------
ifeq ($(USE_SCIP), TRUE)
   ifeq ($(DEBUG), TRUE)
      OPT=dbg
   endif
   include $(SCIPDIR)/make/make.project
endif

#-----------------------------------------------------------------------------
# include project Makefile from COIN
# define BCPDIR and CBCDIR (if not defined)
#-----------------------------------------------------------------------------
ifeq ($(USE_BCP), TRUE)
   ifeq ($(DEBUG), TRUE)
      BCPDIR = $(BCPDIRDBG)
    else
      BCPDIR = $(BCPDIROPT)
    endif
   include make.coin
endif

ifeq ($(USE_BCP), TRUE)
   ifeq ($(DEBUG), TRUE)
      BCPDIR = $(BCPDIRDBG)
    else
      BCPDIR = $(BCPDIROPT)
    endif
   include make.coin
endif

#-----------------------------------------------------------------------------
# add user flags
#-----------------------------------------------------------------------------
INCLUDESFLAGS  += -I$(BOOST_DIR)
CXXFLAGS    += -w -fPIC -fexceptions -std=c++11  -DNDEBUG -DIL_STD  $(INCLUDESFLAGS)
ifeq ($(DEBUG), TRUE)
   CXXFLAGS += -g -O0
   LDFLAGS  += -g -O0
else
   CXXFLAGS += -O3
endif    
OS = $(shell uname -s)
ifeq ($(OS),Linux)
   LDFLAGS += -lrt
endif

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------
BINDIR      =  bin
SRCDIR      =  src
OBJDIR      =  obj

MAINNAME    = roster
MAINOBJ     = main.o # OptimalSolver.o
MAINOBJ	   += main_test.o MyTools.o Demand.o Nurse.o Scenario.o ReadWrite.o DemandGenerator.o Roster.o MasterProblem.o SubProblem.o Solver.o Greedy.o StochasticSolver.o RotationPricer.o

ifeq ($(USE_SCIP), TRUE)
   MAINOBJ  += ScipModeler.o 
endif
ifeq ($(USE_BCP), TRUE)
   MAINOBJ  += BcpModeler.o
endif
ifeq ($(USE_CBC), TRUE)
   MAINOBJ  += CbcModeler.o
endif
MAINSRC     =  $(addprefix $(SRCDIR)/,$(MAINOBJ:.o=.cpp))

MAIN     	=  $(MAINNAME)
MAINFILE 	=  $(BINDIR)/$(MAIN)
MAINOBJFILES	=  $(addprefix $(OBJDIR)/,$(MAINOBJ))

#-----------------------------------------------------------------------------
# Default compiler parameters
#-----------------------------------------------------------------------------
CXX      =  g++
CXX_c    =  -c # the trailing space is important
CXX_o    =  -o # the trailing space is important
LINKCXX     =  g++
LINKCXX_L   =  -L
LINKCXX_l   =  -l
LINKCXX_o   =  -o # the trailing space is important

#-----------------------------------------------------------------------------
# SCIP Flags
#-----------------------------------------------------------------------------
ifeq ($(USE_SCIP), TRUE)
   LIBS +=  $(LINKCXX_L)$(SCIPDIR)/lib $(LINKCXX_l)$(SCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(OBJSCIPLIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(LPILIB)$(LINKLIBSUFFIX) $(LINKCXX_l)$(NLPILIB)$(LINKLIBSUFFIX)
endif

#-----------------------------------------------------------------------------
# Rules
#-----------------------------------------------------------------------------
.PHONY: all
all:            $(SCIPDIR) $(MAINFILE) $(MAINSHORTLINK)

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
ifneq ($(BINDIR),)
		-rm -f $(BINDIR)/$(MAINNAME)
#		-rmdir $(BINDIR)
endif


$(MAINFILE):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(MAINOBJFILES)
		@echo "-> linking $@"
				@echo 		$(LINKCXX) $(MAINOBJFILES) $(LIBS) $(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCXX_o)$@
		$(LINKCXX) $(MAINOBJFILES) $(LIBS) $(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCXX_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $(OBJDIR) $@"
		@echo $(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
