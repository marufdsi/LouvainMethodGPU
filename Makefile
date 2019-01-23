#!/bin/bash

CUDAVERSION=8.0

CC=/usr/local/cuda-$(CUDAVERSION)/bin/nvcc 
CPP = g++
CFLAGS= -I/usr/local/cuda-$(CUDAVERSION)/include -O3 -std=c++11


DFLAGS= -D RUNONGPU
CUDAFLAGS= -arch sm_35 

DEPS = communityGPU.h  graphGPU.h  graphHOST.h openaddressing.h

OBJ = binWiseGaussSeidel.o communityGPU.o preprocessing.o  aggregateCommunity.o coreutility.o independentKernels.o gatherInformation.o graphHOST.o graphGPU.o main.o assignGraph.o computeModularity.o computeTime.o


LIBS= -L/usr/local/cuda-$(CUDAVERSION)/lib64 -lcudart 


EXEC=run_CU_community

all:$(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LIBS) 

%.o: %.cu $(DEPS)
	$(CC) -o $@ -c $< $(CFLAGS) $(DFLAGS) $(CUDAFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) -o $@ -c $< $(CFLAGS) 


clean:
	rm -f *.o *~ $(EXEC)

