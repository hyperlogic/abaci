# build abacitest for Mac OS 10.5.2 (leopard)

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

CONFIG_CFLAGS = $(OPT_CFLAGS) $(SYMBOL_CFLAGS)

# compile flags
CFLAGS = $(CONFIG_CFLAGS) -Wall -DDARWIN -DCMD_LINE -msse3 -I../src/

# linker flags
LFLAGS = -lstdc++ -lm -framework CoreServices

OBJ = abacitest.o unittest.o

# TODO: fix dependencies
DEPENDS = ../src/abaci.h ../src/abaciinlines.cpp unittest.h

abacitest : $(OBJ)
	$(GCC) $(OBJ) -o abacitest $(LFLAGS)

abacitest.o: abacitest.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

unittest.o: unittest.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

clean:
	rm $(OBJ) abacitest
