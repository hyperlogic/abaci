CC = gcc

CFLAGS = -Wall -I ../ -I /opt/local/include
LFLAGS = -shared -L/opt/local/lib -lstdc++ -undefined dynamic_lookup
OBJ = abaci_ext.o vec2_ext.o
DEPENDS = abaci_ext.h vec2_ext.h ../abaci.h ../abaciinlines.cpp ../abaci.cpp

# link
abaci.so: $(OBJ)
	$(CC) $(OBJ) -o abaci.so $(LFLAGS)

# compile
abaci_ext.o: abaci_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

vec2_ext.o: vec2_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm $(OBJ) abaci.so