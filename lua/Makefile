CC = gcc

CFLAGS = -Wall -I ../ -I /opt/local/include
LFLAGS = -shared -L/opt/local/lib -lstdc++ -undefined dynamic_lookup
OBJ = abaci_ext.o vec2_ext.o vec3_ext.o vec4_ext.o quat_ext.o complex_ext.o
DEPENDS = abaci_ext.h vec2_ext.h vec3_ext.h vec4_ext.h quat_ext.h ../abaci.h ../abaciinlines.cpp ../abaci.cpp

# link
abaci.so: $(OBJ)
	$(CC) $(OBJ) -o abaci.so $(LFLAGS)

# compile
abaci_ext.o: abaci_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

vec2_ext.o: vec2_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

vec3_ext.o: vec3_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

vec4_ext.o: vec4_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

quat_ext.o: quat_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

complex_ext.o: complex_ext.cpp $(DEPENDS)
	$(CC) $(CFLAGS) -c $<

clean:
	rm $(OBJ) abaci.so