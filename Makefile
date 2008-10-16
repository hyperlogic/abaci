# build mathtest for Mac OS 10.5.2 (leopard)

GCC = g++
SYMBOLS = true
OPT = true

# symbol flag
ifeq ($(SYMBOLS),true)
  SYMBOL_CFLAGS = -g
else
  SYMBOL_CFLAGS =
endif

# opt flag
ifeq ($(OPT),true)
  OPT_CFLAGS = -O3
else
  OPT_CFLAGS =
endif

# for reference, this is how you generate disassembly in intel syntax (as opposed to at&t style)
# -mintel-syntax -S

CONFIG_CFLAGS = $(OPT_CFLAGS) $(SYMBOL_CFLAGS)

# compile flags
CFLAGS = $(CONFIG_CFLAGS) -Wall -DDARWIN

# linker flags
LFLAGS = -lstdc++ -lm -framework CoreServices

OBJ = math3d.o mathtest.o

# TODO: fix dependencies
DEPENDS = math3d.h

mathtest : $(OBJ)
	$(GCC) $(OBJ) -o mathtest $(LFLAGS)

mathtest.o: mathtest.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

math3d.o: math3d.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

clean:
	rm $(OBJ) mathtest
