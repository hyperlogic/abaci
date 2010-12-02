// lua extentions

#include "abaci.h"

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

static int deg_to_rad(lua_State* L)
{
	double d = luaL_checknumber(L, 1);
	lua_pushnumber(L, DegToRad(d));
	return 1;
}

static int rad_to_deg(lua_State* L)
{
	double r = luaL_checknumber(L, 1);
	lua_pushnumber(L, RadToDeg(r));
	return 1;
}

static int clamp(lua_State* L)
{
	double v = luaL_checknumber(L, 1);
	double min = luaL_checknumber(L, 2);
	double max = luaL_checknumber(L, 3);
	lua_pushnumber(L, Clamp(v, min, max));
	return 1;
}

static int limit_pi(lua_State* L)
{
	double t = luaL_checknumber(L, 1);
	lua_pushnumber(L, LimitPi(t));
	return 1;
}

static int mod_two_pi(lua_State* L)
{
	double t = luaL_checknumber(L, 1);
	lua_pushnumber(L, Mod2Pi(t));
	return 1;
}

static int fuzzy_equals(lua_State* L)
{
	double a = luaL_checknumber(L, 1);
	double b = luaL_checknumber(L, 1);
	lua_pushboolean(L, FuzzyEqual(a, b));
	return 1;
}

static int lerp(lua_State* L)
{
	double a = luaL_checknumber(L, 1);
	double b = luaL_checknumber(L, 2);
	double t = luaL_checknumber(L, 3);
	lua_pushnumber(L, Lerp(a, b, t));
	return 1;
}

static int rand_int(lua_State* L)
{
	int min = luaL_checkinteger(L, 1);
	int max = luaL_checkinteger(L, 2);
	lua_pushinteger(L, RandomInt(min, max));
	return 1;
}

static int rand(lua_State* L)
{
	double min = luaL_checknumber(L, 1);
	double max = luaL_checknumber(L, 2);
	lua_pushnumber(L, RandomScalar(min, max));
	return 1;
}

static const luaL_Reg abaci_lib_f [] = {
	{"deg_to_rad", deg_to_rad},
	{"rad_to_deg", rad_to_deg},
	{"clamp", clamp},
	{"limit_pi", limit_pi},
	{"mod_two_pi", mod_two_pi},
	{"fuzzy_equals", fuzzy_equals},
	{"lerp", lerp},
	{"rand_int", rand_int},
	{"rand", rand},
	{NULL, NULL}
};

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

extern "C" int luaopen_abaci(lua_State* L)
{
//	luaL_newmetatable(L, "abaci.abaci");
	luaL_register(L, "abaci", abaci_lib_f);

	luaL_newmetatable(L, "abaci.vec2");

	lua_pushvalue(L, -1);
	lua_setfield(L, -2, "__index");
	luaL_register(L, NULL, vec2_lib_m);

	
	luaL_register(L, "vec2", vec2_lib_f);

	return 1;
}
