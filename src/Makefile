#!/bin/bash

SYSTEM = x86-64_linux
LIBFORMAT = static_pic

#OTHERDIR = /home/samuel/Dropbox/PhD/dev/includes

#OTHERFLAGS = -L/home/samuel/Dropbox/PhD/dev/lib

#INCLUDESFLAGS = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -I$(LEMONDIR)/include -I$(CLIQUERDIR) -I$(OTHERDIR)

CXX		= g++
#CXXFLAGS	= -m64 -w -fPIC -fexceptions -DNDEBUG -DIL_STD $(INCLUDESFLAGS)
CXXFLAGS	= -m64 -w -fPIC -fexceptions -DNDEBUG -DIL_STD #$(INCLUDESFLAGS)
OBJECTS		= ReadWrite.o main.o
TARGET		= roster

$(TARGET):	$(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(CXXFLAGS) #$(OTHERFLAGS)

clean:
	rm -f $(OBJECTS) $(TARGET)
