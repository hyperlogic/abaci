# build mathtest for Mac OS 10.5.2 (leopard)

GCC = g++
SYMBOLS = false
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
CFLAGS = $(CONFIG_CFLAGS) -Wall -DDARWIN -msse3

# linker flags
LFLAGS = -lstdc++ -lm -framework CoreServices

OBJ = math3d.o mathunittest.o unittest.o

# TODO: fix dependencies
DEPENDS = math3d.h math3dinlines.cpp unittest.h

mathunittest : $(OBJ)
	$(GCC) $(OBJ) -o mathunittest $(LFLAGS)

mathunittest.o: mathunittest.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

unittest.o: unittest.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

math3d.o: math3d.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

clean:
	rm $(OBJ) mathunittest
