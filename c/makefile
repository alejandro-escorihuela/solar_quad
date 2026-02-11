# 06-06-2024
# alex

OBJ = vector.so solar.so kpert.so evol.so # Objectes
ICC = /opt/intel/oneapi/compiler/2024.0/bin/icx  # Compilador C intel
# Per al compilador d'intel primer:
# $ export LD_LIBRARY_PATH="/opt/intel/oneapi/compiler/latest/lib:/usr/local/lib"
GCC = gcc # Compilador C gnu
COM = $(GCC)
LIB = -I/usr/lib/gcc/x86_64-linux-gnu/12/include/ -lm -lquadmath # Llibreries (-l*, -L*, -I*)
AVS = -W -Wall -Wextra -Wno-format # -pedantic i -ansi no funcionen bé amb -lquadmath
OPT = -march=native -O2 -ffast-math -pipe -std=c17 -shared -fPIC # Optimitzacio
OPC = $(OPT) $(AVS) # Opcions

all: $(OBJ)

vector.so: vector.c vector.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@ $(LIB)

solar.so: solar.c solar.h vector.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@ $(LIB)

kpert.so: kpert.c kpert.h vector.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@ $(LIB)

evol.so: evol.c evol.h vector.h solar.h kpert.h
	@echo Compilant $<
	$(COM) $(OPC) $< -o $@ $(LIB)

clean:
	@echo Netejant...
	rm -rf $(OBJ) *_d.so *~
