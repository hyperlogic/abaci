// lua extentions

#include "abaci.h"
#include <stdio.h>
#include <string.h>
#include "vec2_ext.h"


//
// abaci functions
//

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
	double b = luaL_checknumber(L, 2);
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


/* stack notes: (i always forget this)
   -1 is top, -2 is just under that.
   1 is bottom of stack, 2 is the one above that.

   3 | c | -1
   2 | b | -2
   1 | a | -3
     -----
*/

extern "C" int luaopen_abaci(lua_State* L)
{
    // register abaci module table, and push it.
	luaL_register(L, "abaci", abaci_lib_f);

    init_vec2(L);

	return 1;
}
