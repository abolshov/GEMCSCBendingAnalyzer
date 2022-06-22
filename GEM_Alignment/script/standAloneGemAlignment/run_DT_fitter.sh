g++ -g -o DT_fitter.o DT_fitter.cpp $(root-config --cflags --glibs --ldflags) -lMinuit
./DT_fitter.o
