// lua extentions

#include "abaci.h"

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

#define check_vec2(L) (Vector2f*)luaL_checkudata(L, 1, "abaci.vec2")

static int vec2_new(lua_State* L)
{
	double x = luaL_checknumber(L, 1);
	double y = luaL_checknumber(L, 2);
	Vector2f* v = (Vector2f*)lua_newuserdata(L, sizeof(Vector2f));
	v->Set((float)x, (float)y);

	luaL_getmetatable(L, "abaci.vec2");
	lua_setmetatable(L, -2);

	return 1;
}

static int vec2_len(lua_State* L)
{
	Vector2f* v = check_vec2(L);
	lua_pushnumber(L, (double)v->Len());
	return 1;
}

static int vec2_getx(lua_State* L)
{
	Vector2f* v = check_vec2(L);
	lua_pushnumber(L, (double)v->x);
	return 1;
}

static int vec2_gety(lua_State* L)
{
	Vector2f* v = check_vec2(L);
	lua_pushnumber(L, (double)v->y);
	return 1;
}

static const luaL_Reg vec2_lib_f [] = {
	{"new", vec2_new},
	{NULL, NULL}
};

static const luaL_Reg vec2_lib_m [] = {
	{"len", vec2_len},
	{"getx", vec2_getx},
	{"gety", vec2_gety},
	{NULL, NULL}
};

int luaopen_abaci(lua_State* L)
{
	luaL_newmetatable(L, "abaci.vec2");

	lua_pushvalue(L, -1);
	lua_setfield(L, -2, "__index");
	luaL_register(L, NULL, vec2_lib_m);
	
	luaL_register(L, "vec2", vec2_lib_f);
	return 1;
}
