CCLINUX = g++ -O3 -std=c++11 -Wall
CCMAC = g++ -fast -std=c++11 -Wall
# Compiler 
CC = ${CCLINUX}
# Remove file command
RM=@rm -f
# Library path flags
LIBS=
# comment the following line and uncomment the next one if your gcc version does not include tr1
XFLAGS = 
#XFLAGS = -Ithird-party -Ithird-party/boost/tr1
INC=-I ../include
OBJDIR = .
SRCDIR = .
SOURCES = `ls *.cpp`
# Object files
OBJS=$(subst .cpp,.o,$(SOURCES))
EXEDIR = ../bin
# Executables
EXE = filter viterbi likelihood 


all: ${EXE}

%.o: %.cpp
	@echo "building $@ ..."
	${CC} ${XFLAGS} -c -o $@ $<

viterbi: 	viterbi.o utils.o inference.o 
	@echo 'linking viterbi2...'
	${CC} -o ${EXEDIR}/viterbi $^ ${LIBS}

viterbi2: 	inference.o utils.o viterbi2.o
	@echo 'linking viterbi2...'
	${CC} -o ${EXEDIR}/viterbi2 $^ ${LIBS}

enum_mpe: 	inference.o utils.o enum_mpe.o
	@echo 'linking enumerative...'
	${CC} -o ${EXEDIR}/enum_mpe S^ ${LIBS}

likelihood:	inference.o utils.o likelihood.o
	@echo 'linking likelihood...'
	${CC} -o ${EXEDIR}/likelihood $^ ${LIBS}

filter:		inference.o utils.o filtering.o
	@echo 'linking filter...'
	${CC} -o ${EXEDIR}/filter $^ ${LIBS}

filter2:	inference.o utils.o filtering2.o
	@echo 'linking filter2...'
	${CC} -o ${EXEDIR}/filter2 $^ ${LIBS}


# CLEAN

.PHONY: clean

clean:
	@echo "removing objects..."
	${RM} ${OBJS}

dist-clean: clean
	@echo "removing executables..."
	${RM} ${EXE}
