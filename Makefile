CC = g++ -O3 -Wno-deprecated

CFLAGS = -Wall -fpic -g

INCLUDE = -I

SRC = ./src
LIB = ./lib
BIN = ./bin
OBJ = ./obj
STREE = $(LIB)/stree
FRACT = $(LIB)/fractional

all: folder multigene.o

folder:
	mkdir -p $(OBJ) $(BIN)

teste.o: dinkelbach.o thermalign.o
	$(CC) $(CFLAGS) -o $(BIN)/teste $(SRC)/teste.cpp $(OBJ)/*.o $(INCLUDE) $(FRACT)

multigene.so:
	g++ -shared -o ./gui/Flask/$@ $(OBJ)/*.o

multigene.o: dinkelbach.o thermalign.o stree_strmat.o stree_ukkonen.o 
	$(CC) $(CFLAGS) -o $(OBJ)/$@ -c $(SRC)/multigene.cpp $(OBJ)/*.o $(INCLUDE) $(STREE) $(INCLUDE) $(FRACT) 

dinkelbach.o: galign.o nnparams.o
	$(CC) $(CFLAGS) -o $(OBJ)/$@ -c $(FRACT)/dinkelbach.cpp

thermalign.o: nnparams.o
	$(CC) $(CFLAGS) -o $(OBJ)/$@ -c $(FRACT)/thermalign.cpp

galign.o: nnparams.o
	$(CC) $(CFLAGS) -o $(OBJ)/$@ -c $(FRACT)/galign.cpp

nnparams.o:
	$(CC) $(CFLAGS) -o $(OBJ)/$@ -c $(FRACT)/nnparams.cpp

stree_strmat.o:
	$(CC) $(CFLAGS) -o $(OBJ)/$@ -c $(STREE)/stree_strmat.c 

stree_ukkonen.o:
	$(CC) $(CFLAGS) -o $(OBJ)/$@ -c $(STREE)/stree_ukkonen.c 

clean:
	rm $(BIN)/* $(OBJ)/*.o
