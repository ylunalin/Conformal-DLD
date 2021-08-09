# This is a sample config file for MacPorts installed GNU complier
# To use copy into ../ and delete the suffix .mac
# Replace the compiler commands and the flags with ones on your computing system
# Compiler and compilation flags
cc=gcc-mp-9
cxx=g++-mp-9 -fopenmp
cflags=-Wall -march=native -ansi -pedantic -O3 -fopenmp -std=c++11 -fPIC -lstdc++fs -Wno-variadic-macros -g
