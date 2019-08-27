CC = g++ -O3 -Wno-deprecated
CFLAGS = -Wall -g
INCLUDE = -I

SRC = ./src
LIB = ./lib
BIN = ./bin
OBJ = ./obj
STREE = ${LIB}/stree
FRACT = ${LIB}/fractional

$(shell mkdir -p $(OBJ))
$(shell mkdir -p $(BIN))

all: multigene.o

multigene.o: dinkelbach.o thermalign.o stree_strmat.o stree_ukkonen.o 
	${CC} ${CFLAGS} -o ${BIN}/multigene ${SRC}/multigene.cpp ${OBJ}/*.o ${INCLUDE} ${STREE} ${INCLUDE} ${FRACT} 

teste.o: dinkelbach.o thermalign.o
	${CC} ${CFLAGS} -o ${BIN}/teste ${SRC}/teste.cpp ${OBJ}/*.o ${INCLUDE} ${FRACT}

stree_strmat.o:
	${CC} ${CFLAGS} -o ${OBJ}/stree_strmat.o -c ${STREE}/stree_strmat.c 

stree_ukkonen.o:
	${CC} ${CFLAGS} -o ${OBJ}/stree_ukkonen.o -c ${STREE}/stree_ukkonen.c 

dinkelbach.o: galign.o nnparams.o
	${CC} ${CFLAGS} -o ${OBJ}/dinkelbach.o -c ${FRACT}/dinkelbach.cpp

thermalign.o: nnparams.o
	${CC} ${CFLAGS} -o ${OBJ}/thermalign.o -c ${FRACT}/thermalign.cpp

galign.o: nnparams.o
	${CC} ${CFLAGS} -o ${OBJ}/galign.o -c ${FRACT}/galign.cpp

nnparams.o:
	${CC} ${CFLAGS} -o ${OBJ}/nnparams.o -c ${FRACT}/nnparams.cpp

clean:
	rm ${BIN}/* ${OBJ}/*.o
