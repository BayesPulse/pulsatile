################################################################################
# Makefile for decon
################################################################################

SRC  = deconvolution_main.c format_data.c linklistv3.c hash.c\
	     randgen.c birthdeath_strauss.c cholesky.c normal.c
OBJ  = deconvolution_main.o format_data.o linklistv3.o hash.o\
   	   randgen.o birthdeath_strauss.o cholesky.o normal.o
DEPS = deconvolution_main.h format_data.h linklistv3.h hash.h\
	     randgen.h birthdeath_strauss.h cholesky.h normal.h
LIBS = -lm

# allows using git version in C code
# see http://goo.gl/PHJ9Zk
# removed options: --dirty 
GIT_VERSION := $(shell git describe --abbrev=4 --always --tags)


# Debugging version
# -pg for gprof
#  -g for gdb debugging
#  -Wall for display all warnings
#CFLAGS = -g -I. -DVERSION=\"$(GIT_VERSION)\"
#decon: ${OBJ}
#	gcc  -o decon ${OBJ} ${LIBS} ${CFLAGS} -lprofiler
#
#${OBJ} : ${SRC}
#	gcc -Wall -c  ${SRC} ${CFLAGS} -lprofiler


# Production version
CFLAGS = -O3 -I. -DVERSION=\"$(GIT_VERSION)\"
decon: ${OBJ}
	gcc -o decon ${OBJ} ${LIBS} ${CFLAGS}

${OBJ} : ${SRC}
	gcc -c  ${SRC} ${CFLAGS}


## 'make clean' will remove obj files
.PHONY: clean

clean:
	rm ${OBJ}
