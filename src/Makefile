#!/bin/bash

SYSTEM = x86-64_linux
LIBFORMAT = static_pic

INCLUDESFLAGS = -I$(BOOST_DIR)

CXX		= g++
#CXXFLAGS	= -m64 -w -fPIC -fexceptions -DNDEBUG -DIL_STD $(INCLUDESFLAGS)
CXXFLAGS	= -g -m64 -w -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++11 $(INCLUDESFLAGS)
OBJECTS		= MyTools.o Demand.o Nurse.o Scenario.o ReadWrite.o Roster.o SubProblem.o Solver.o Greedy.o main.o
TARGET		= ../bin/roster

$(TARGET):	$(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(CXXFLAGS) $(OTHERFLAGS)

clean:
	rm -f $(OBJECTS) $(TARGET)
