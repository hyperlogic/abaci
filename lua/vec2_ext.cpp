#include "vec2_ext.h"
#include "abaci.h"
#include <string.h>

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
        if (i > 0 && i <= 2)
        {
            lua_pushnumber(L, (*v)[i-1]);
            return 1;
        }

        // error
        lua_pushstring(L, "vec2: num index out of range, must be 1 or 2");
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
        if (i > 0 && i <= 2)
        {
            (*v)[i] = value;
            return 0;
        }

        // error
        lua_pushstring(L, "vec2: num index out of range, must be 1 or 2");
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
    const int kLen = 2;
    lua_pushinteger(L, kLen);
    return 1;
}

// TODO: support numbers! i.e. vec * num
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

// assumes the "abaci" table is the top of the stack.
void init_vec2(lua_State* L)
{
    // registers all of vec2_class_funcs functions in the abaci.vec2 table.
    lua_newtable(L);
    luaL_register(L, NULL, vec2_class_funcs);
    lua_setfield(L, -2, "vec2");

    // metatable for use with vec2 userdata.
	luaL_newmetatable(L, "abaci.vec2");

    // registers all vec2_meta_funcs functions in vec2_mt
	luaL_register(L, NULL, vec2_meta_funcs);
    lua_pop(L, 1);
}
