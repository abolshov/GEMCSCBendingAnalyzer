#!/bin/bash
g++ -g -o GEM_fitter.o GEM_fitter.cpp $(root-config --cflags --glibs --ldflags) -lMinuit
./GEM_fitter.o
