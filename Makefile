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
USE_BCP = TRUE
USE_BCP_GUROBI = FALSE
USE_BCP_CPLEX = FALSE
USE_CBC = FALSE
DEBUG  = FALSE
MEMORY_PROFILE = FALSE

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
# include project Makefile from COIN
# define BCPDIR and CBCDIR (if not defined)
#-----------------------------------------------------------------------------
ifeq ($(USE_BCP), TRUE)
   ifeq ($(DEBUG), TRUE)
      BCPDIR = $(BCPDIRDBG)
   else
      BCPDIR = $(BCPDIROPT)
   endif
#   ifeq ($(USE_BCP_GUROBI), TRUE)
#      BCPDIR = $(BCPGRBDIR)
#   endif
   include make.coin
endif

ifeq ($(USE_CBC), TRUE)
   ifeq ($(DEBUG), TRUE)
      CBCDIR = $(CBCDIRDBG)
    else
      CBCDIR = $(CBCDIROPT)
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
ifeq ($(MEMORY_PROFILE), TRUE)
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------
BINDIR      =  bin
SRCDIR      =  src
OBJDIR      =  obj

EXEC_DET = staticscheduler
OBJ_DET = DeterministicMain.o DeterministicMain_test.o
COMMONOBJ  = InputPaths.o GlobalStats.o InitializeSolver.o MyTools.o Demand.o Nurse.o Scenario.o ReadWrite.o Roster.o MasterProblem.o SubProblem.o Solver.o  RotationPricer.o TreeManager.o 
#Greedy.o

ifeq ($(USE_SCIP), TRUE)
	CXXFLAGS += -DUSE_SCIP
endif
ifeq ($(USE_BCP), TRUE)
   COMMONOBJ  += BcpModeler.o
endif
ifeq ($(USE_CBC), TRUE)
	 CXXFLAGS += -DUSE_CBC
endif
ifeq ($(USE_BCP_GUROBI), TRUE)
	 CXXFLAGS += -DUSE_GUROBI
endif
ifeq ($(USE_BCP_CPLEX), TRUE)
	 CXXFLAGS += -DUSE_CPLEX
endif

OBJ_DET += $(COMMONOBJ) DeterministicSolver.o
OBJFILES_DET = $(addprefix $(OBJDIR)/,$(OBJ_DET))

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
all: $(SCIPDIR) $(EXEC_STO) $(EXEC_DET) # $(MAINSHORTLINK)

#$(MAINSHORTLINK):	$(MAINFILE)
#		@rm -f $@
#		cd $(dir $@) && ln -s $(notdir $(MAINFILE)) $(notdir $@)

$(OBJDIR):
		@-mkdir -p $(OBJDIR)

$(BINDIR):
		@-mkdir -p $(BINDIR)

.PHONY: clean
clean:		$(OBJDIR)
ifneq ($(OBJDIR),)
		-rm -f $(OBJDIR)/*
endif
ifneq ($(BINDIR),)
		-rm -f $(BINDIR)/*
endif

$(EXEC_DET):	$(BINDIR) $(OBJDIR) $(SCIPLIBFILE) $(LPILIBFILE) $(NLPILIBFILE) $(OBJFILES_DET)
		@echo "-> linking $@"
		@echo 		$(LINKCXX) $(OBJFILES_DET) $(LIBS) $(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCXX_o)$@
		$(LINKCXX) $(OBJFILES_DET) $(LIBS) $(OFLAGS) $(LPSLDFLAGS) $(LDFLAGS) $(LINKCXX_o)$@
		-mv $(EXEC_DET) $(BINDIR)

$(OBJDIR)/%.o:	$(SRCDIR)/%.c
		@echo "-> compiling $@"
		$(CC) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CFLAGS) -c $< $(CC_o)$@

$(OBJDIR)/%.o:	$(SRCDIR)/%.cpp
		@echo "-> compiling $(OBJDIR) $@"
		$(CXX) $(FLAGS) $(OFLAGS) $(BINOFLAGS) $(CXXFLAGS) -c $< $(CXX_o)$@

#---- EOF --------------------------------------------------------------------
