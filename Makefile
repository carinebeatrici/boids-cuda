CUDACC = /usr/local/cuda/bin/nvcc
CC = gcc
SRC = src
INC = include
OBJ = obj
CFLAGS = -fbounds-check -Wall -Wextra -I$(INC) -lm -g
CUDAFLAGS =  -I$(INC) -I/usr/local/cuda/include -lm -lcuda -G
EXEC = gpu-todos


XOPT="-arch compute_30 -code sm_30"
if [ "$EXE" = "" ]; then EXE="GPU-BOX";fi

all:  
	$(CC) -c $(SRC)/alocation.c -o $(OBJ)/alocation.o $(CFLAGS)
	$(CC) -c $(SRC)/free.c -o $(OBJ)/free.o $(CFLAGS)
	$(CUDACC) -c $(SRC)/initialize.cu -o $(OBJ)/initialize.o $(CUDAFLAGS) 
	$(CUDACC) -c $(SRC)/main.cu -o $(OBJ)/main.o $(CUDAFLAGS)
	$(CUDACC)  $(OBJ)/main.o  $(OBJ)/free.o $(OBJ)/alocation.o $(OBJ)/initialize.o -o $(EXEC) $(CUDAFLAGS)
	
clean:
	rm src/*~	
	rm $(INC)/*~
	rm $(OBJ)/*

