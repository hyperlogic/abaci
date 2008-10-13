# build mathtest for Mac OS 10.5.2 (leopard)

GCC = g++
DEBUG = false

# config specific flags
ifeq ($(DEBUG),true)
  CONFIG_CFLAGS = -g
else
  CONFIG_CFLAGS = -O3
endif

# Mac OS X
# uses SDL, SDL_image, OpenGL, FreeType2, Chicken & Bullet
CFLAGS = $(CONFIG_CFLAGS) -Wall -DDARWIN

# er? do i have to do a -Wl, to pass options to the linker?
LFLAGS = -lstdc++ -lm

OBJ = math3d.o mathtest.o

# TODO: fix dependencies
DEPENDS = math3d.h

mathtest : $(OBJ)
	$(GCC) $(OBJ) -o mathtest $(LFLAGS)

math3d.o: math3d.cpp $(DEPENDS)
	$(GCC) $(CFLAGS) -c $<

clean:
	rm $(OBJ) mathtest
