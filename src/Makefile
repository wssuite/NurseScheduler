#!/bin/bash

SYSTEM = x86-64_linux
LIBFORMAT = static_pic

#INCLUDESFLAGS = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -I$(LEMONDIR)/include -I$(CLIQUERDIR) -I$(OTHERDIR)

CXX		= g++
#CXXFLAGS	= -m64 -w -fPIC -fexceptions -DNDEBUG -DIL_STD $(INCLUDESFLAGS)
CXXFLAGS	= -g -m64 -w -fPIC -fexceptions -DNDEBUG -DIL_STD #$(INCLUDESFLAGS)
OBJECTS		= MyTools.o Nurse.o Scenario.o ReadWrite.o main.o
TARGET		= roster

$(TARGET):	$(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(CXXFLAGS) #$(OTHERFLAGS)

clean:
	rm -f $(OBJECTS) $(TARGET)
