#ifndef ABACI_EXT_H
#define ABACI_EXT_H

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

extern "C" int luaopen_abaci(lua_State* L);


// helper macros

#define check_vec2(L, n) (Vector2f*)luaL_checkudata(L, n, "abaci.vec2")
#define new_vec2(L, name)                                             \
	Vector2f* name = (Vector2f*)lua_newuserdata(L, sizeof(Vector2f)); \
	luaL_getmetatable(L, "abaci.vec2");                               \
	lua_setmetatable(L, -2)

#define check_vec3(L, n) (Vector3f*)luaL_checkudata(L, n, "abaci.vec3")
#define new_vec3(L, name)                                             \
	Vector3f* name = (Vector3f*)lua_newuserdata(L, sizeof(Vector3f)); \
	luaL_getmetatable(L, "abaci.vec3");                               \
	lua_setmetatable(L, -2)

#define check_vec4(L, n) (Vector4f*)luaL_checkudata(L, n, "abaci.vec4")
#define new_vec4(L, name)                                             \
	Vector4f* name = (Vector4f*)lua_newuserdata(L, sizeof(Vector4f)); \
	luaL_getmetatable(L, "abaci.vec4");                               \
	lua_setmetatable(L, -2)

#define check_quat(L, n) (Quatf*)luaL_checkudata(L, n, "abaci.quat")
#define new_quat(L, name)                                             \
	Quatf* name = (Quatf*)lua_newuserdata(L, sizeof(Quatf));          \
	luaL_getmetatable(L, "abaci.quat");                               \
	lua_setmetatable(L, -2)

#define check_complex(L, n) (Complexf*)luaL_checkudata(L, n, "abaci.complex")
#define new_complex(L, name)                                          \
	Complexf* name = (Complexf*)lua_newuserdata(L, sizeof(Complexf)); \
	luaL_getmetatable(L, "abaci.complex");                            \
	lua_setmetatable(L, -2)

#endif
