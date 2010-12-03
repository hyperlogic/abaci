// lua extentions

#include "abaci.h"
#include <stdio.h>
#include <string.h>

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

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

//
// vec2 functions
//

#define check_vec2(L, n) (Vector2f*)luaL_checkudata(L, n, "abaci.vec2")
#define new_vec2(L, name)                                             \
	Vector2f* name = (Vector2f*)lua_newuserdata(L, sizeof(Vector2f)); \
	luaL_getmetatable(L, "abaci.vec2");                               \
	lua_setmetatable(L, -2)

static int vec2_new(lua_State* L)
{
	lua_Number x = luaL_checknumber(L, 1);
	lua_Number y = luaL_checknumber(L, 2);
    new_vec2(L, result);
	result->Set((float)x, (float)y);
	return 1;
}

static int vec2_random_unit(lua_State* L)
{
    new_vec2(L, result);
    *result = Vector2f::RandomUnitVector();
    return 1;
}

static int vec2_lerp(lua_State* L)
{
    Vector2f* a = check_vec2(L, 1);
    Vector2f* b = check_vec2(L, 2);
    lua_Number t = luaL_checknumber(L, 3);

    new_vec2(L, result);
    *result = Lerp(*a, *b, (float)t);
    return 1;
}

static const luaL_Reg vec2_class_funcs [] = {
	{"new", vec2_new},
    {"random_unit", vec2_random_unit},
    {"lerp", vec2_lerp},
	{NULL, NULL}
};

//
// vec2 methods
//

static int vec2_len(lua_State* L)
{
    Vector2f* self = check_vec2(L, 1);
    lua_pushnumber(L, self->Len());
    return 1;
}

static int vec2_unit(lua_State* L)
{
    Vector2f* self = check_vec2(L, 1);
    new_vec2(L, result);
    *result = self->Unit();
    return 1;
}

static int vec2_min_len(lua_State* L)
{
    Vector2f* self = check_vec2(L, 1);
    lua_Number len = luaL_checknumber(L, 2);

    new_vec2(L, result);
    *result = self->MinLen(len);
    return 1;
}

static const luaL_Reg vec2_method_funcs [] = {
    {"len", vec2_len},
    {"unit", vec2_unit},
	{"min_len", vec2_min_len},
	{NULL, NULL}
};

//
// vec2 meta methods
//

static int vec2_index(lua_State* L)
{
    Vector2f* v = check_vec2(L, 1);
    luaL_checkany(L, 2);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i >= 0 && i < 2)
        {
            lua_pushnumber(L, (*v)[i]);
            return 1;
        }

        // error
        lua_pushstring(L, "vec2: num index out of range, must be 0 or 1");
        lua_error(L);
    }
    else if (key_type == LUA_TSTRING)
    {
        size_t size;
        const char* s = lua_tolstring(L, 2, &size);
        if (size == 1)
        {
            switch (s[0])
            {
            case 'x':
                lua_pushnumber(L, v->x);
                return 1;
            case 'y':
                lua_pushnumber(L, v->y);
                return 1;
            }
        }
        else
        {
            // search for string key in vec2_method_funcs array.
            int i = 0;
            while (vec2_method_funcs[i].name)
            {
                if (!strcmp(s, vec2_method_funcs[i].name))
                {
                    lua_pushcfunction(L, vec2_method_funcs[i].func);
                    return 1;
                }
                i++;
            }
        }

        // error
        lua_pushstring(L, "vec2: unknown string key, must be x or y");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "vec2: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int vec2_newindex(lua_State* L)
{
    Vector2f* v = check_vec2(L, 1);
    luaL_checkany(L, 2);
    lua_Number value = luaL_checknumber(L, 3);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i >= 0 && i < 2)
        {
            (*v)[i] = value;
            lua_pushnumber(L, (*v)[i]);
            return 0;
        }

        // error
        lua_pushstring(L, "vec2: num index out of range, must be 0 or 1");
        lua_error(L);
    }
    else if (key_type == LUA_TSTRING)
    {
        size_t size;
        const char* s = lua_tolstring(L, 2, &size);
        if (size == 1)
        {
            switch (s[0])
            {
            case 'x':
                v->x = value;
                return 0;
            case 'y':
                v->y = value;
                return 0;
            }
        }

        // error
        lua_pushstring(L, "vec2: unknown string key, must be x or y");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "vec2: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int vec2_tostring(lua_State* L)
{
    Vector2f* v = check_vec2(L, 1);
    char temp[64];
    snprintf(temp, 64, "vec2(%.5f, %.5f)", v->x, v->y);
    lua_pushstring(L, temp);
    return 1;
}

static int vec2_lenop(lua_State* L)
{
    lua_pushinteger(L, 2);
    return 1;
}

#define BINARY_VEC_OP_FUNC(name, op)            \
    static int vec2_##name(lua_State* L)        \
    {                                           \
        Vector2f* a = check_vec2(L, 1);         \
        Vector2f* b = check_vec2(L, 2);         \
        new_vec2(L, result);                    \
        *result = *a op *b;                     \
        return 1;                               \
    }

BINARY_VEC_OP_FUNC(add, +)
BINARY_VEC_OP_FUNC(sub, -)
BINARY_VEC_OP_FUNC(mul, *)
BINARY_VEC_OP_FUNC(div, /)

static int vec2_dot(lua_State* L)
{
    Vector2f* a = check_vec2(L, 1);
    Vector2f* b = check_vec2(L, 2);
    lua_pushnumber(L, Dot(*a, *b));
    return 1;
}

static int vec2_unm(lua_State* L)
{
    Vector2f* a = check_vec2(L, 1);
    new_vec2(L, result);
    *result = -*a;
    return 1;
}

static const luaL_Reg vec2_meta_funcs [] = {
    {"__index", vec2_index},
    {"__newindex", vec2_newindex},
	{"__len", vec2_lenop},
    {"__tostring", vec2_tostring},
    {"__add", vec2_add},
    {"__sub", vec2_sub},
    {"__mul", vec2_mul},
    {"__div", vec2_div},
    {"__pow", vec2_dot},  // ^ as dot product
    {"__unm", vec2_unm},  // unary minus
	{NULL, NULL}
};

/* stack notes: (i always forget this)
   -1 is top, -2 is just under that.
   1 is bottom of stack, 2 is the one above that.

   3 | c | -1
   2 | b | -2
   1 | a | -3
     -----

   Conundrums  (lean toward :)
   * should length of a vector be v.len or v:len()
   * should unit vector be v.unit or v:unit()
   * should min_len be or v.min_len(10), v:min_len(10)

*/

extern "C" int luaopen_abaci(lua_State* L)
{
    // register abaci module table, and push it.
	luaL_register(L, "abaci", abaci_lib_f);

    // registers all of vec2_class_funcs functions in the abaci.vec2 table.
    lua_newtable(L);
    luaL_register(L, NULL, vec2_class_funcs);
    lua_setfield(L, -2, "vec2");

    // metatable for use with vec2 userdata.
	luaL_newmetatable(L, "abaci.vec2");

    /*
    // vec2_mt.__index = vec2_index
    lua_pushcfunction(L, vec2_index);
	lua_setfield(L, -2, "__index");

    // vec2_mt.__newindex = vec2_newindex
    lua_pushcfunction(L, vec2_newindex);
	lua_setfield(L, -2, "__newindex");

    // vec2_mt.__tostring = vec2_tostring
    lua_pushcfunction(L, vec2_tostring);
	lua_setfield(L, -2, "__tostring");
    */

    // registers all vec2_meta_funcs functions in vec2_mt
	luaL_register(L, NULL, vec2_meta_funcs);
    lua_pop(L, 1);

	return 1;
}
