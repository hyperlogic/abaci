#include "complex_ext.h"
#include "abaci.h"
#include <string.h>

//
// complex functions
//

static int complex_new(lua_State* L)
{
    luaL_checkany(L, 1);

    int key_type = lua_type(L, 1);
    if (key_type == LUA_TNUMBER)
    {
        lua_Number r = luaL_checknumber(L, 1);
        lua_Number i = luaL_checknumber(L, 2);
        new_complex(L, result);
        result->Set((float)r, (float)i);
        return 1;
    }
    else if (key_type == LUA_TUSERDATA)
    {
        Vector2f* v = check_vec2(L, 1);
        new_complex(L, result);
        result->Set((float)v->x, (float)v->y);
        return 1;
    }

    // error
    lua_pushstring(L, "complex: new expected vec2 or 2 numbers");
    lua_error(L);
    return 0;
}

static int complex_exp(lua_State* L)
{
    Complexf* z = check_complex(L, 1);
    new_complex(L, result);
    *result = Exp(*z);
    return 1;
}

static int complex_expi(lua_State* L)
{
    lua_Number x = luaL_checknumber(L, 1);
    new_complex(L, result);
    *result = ExpI((float)x);
    return 1;
}

static int complex_log(lua_State* L)
{
    Complexf* z = check_complex(L, 1);
    new_complex(L, result);
    *result = Log(*z);
    return 1;
}

static const luaL_Reg complex_class_funcs [] = {
	{"new", complex_new},
    {"exp", complex_exp},
    {"expi", complex_expi},
    {"log", complex_log},
	{NULL, NULL}
};

//
// complex methods
//

static int complex_len(lua_State* L)
{
    Complexf* self = check_complex(L, 1);
    lua_pushnumber(L, self->Len());
    return 1;
}

static int complex_unit(lua_State* L)
{
    Complexf* self = check_complex(L, 1);
    new_complex(L, result);
    *result = self->Unit();
    return 1;
}

static int complex_min_len(lua_State* L)
{
    Complexf* self = check_complex(L, 1);
    lua_Number len = luaL_checknumber(L, 2);

    new_complex(L, result);
    *result = self->MinLen(len);
    return 1;
}

static int complex_conj(lua_State* L)
{
    Complexf* self = check_complex(L, 1);
    new_complex(L, result);
    *result = ~(*self);
    return 1;
}

static const luaL_Reg complex_method_funcs [] = {
    {"len", complex_len},
    {"unit", complex_unit},
	{"min_len", complex_min_len},
    {"conj", complex_conj},
	{NULL, NULL}
};

//
// complex meta methods
//

static int complex_index(lua_State* L)
{
    Complexf* z = check_complex(L, 1);
    luaL_checkany(L, 2);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i == 1)
        {
            lua_pushnumber(L, z->r);
            return 1;
        }
        else if (i == 2)
        {
            lua_pushnumber(L, z->i);
            return 1;
        }

        // error
        lua_pushstring(L, "complex: num index out of range, must be 1 or 2");
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
            case 'r':
                lua_pushnumber(L, z->r);
                return 1;
            case 'i':
                lua_pushnumber(L, z->i);
                return 1;
            }
        }
        else
        {
            // search for string key in complex_method_funcs array.
            int i = 0;
            while (complex_method_funcs[i].name)
            {
                if (!strcmp(s, complex_method_funcs[i].name))
                {
                    lua_pushcfunction(L, complex_method_funcs[i].func);
                    return 1;
                }
                i++;
            }
        }

        // error
        lua_pushstring(L, "complex: unknown string key");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "complex: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int complex_newindex(lua_State* L)
{
    Complexf* z = check_complex(L, 1);
    luaL_checkany(L, 2);
    lua_Number value = luaL_checknumber(L, 3);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i == 0)
        {
            z->r = value;
            return 0;
        }
        else if (i == 1)
        {
            z->i = value;
            return 0;
        }

        // error
        lua_pushstring(L, "complex: num index out of range, must be 1 or 2");
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
            case 'r':
                z->r = value;
                return 0;
            case 'i':
                z->i = value;
                return 0;
            }
        }

        // error
        lua_pushstring(L, "complex: unknown string key, expecting r or i");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "complex: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int complex_tostring(lua_State* L)
{
    Complexf* z = check_complex(L, 1);
    char temp[64];
    snprintf(temp, 64, "complex(%.5f, %.5f)", z->r, z->i);
    lua_pushstring(L, temp);
    return 1;
}

static int complex_lenop(lua_State* L)
{
    const int kLen = 2;
    lua_pushinteger(L, kLen);
    return 1;
}

// TODO: support numbers! i.e. vec * num
#define BINARY_VEC_OP_FUNC(name, op)            \
    static int complex_##name(lua_State* L)     \
    {                                           \
        Complexf* a = check_complex(L, 1);      \
        Complexf* b = check_complex(L, 2);      \
        new_complex(L, result);                 \
        *result = *a op *b;                     \
        return 1;                               \
    }

BINARY_VEC_OP_FUNC(add, +)
BINARY_VEC_OP_FUNC(sub, -)
BINARY_VEC_OP_FUNC(mul, *)
BINARY_VEC_OP_FUNC(div, /)

static int complex_dot(lua_State* L)
{
    Complexf* a = check_complex(L, 1);
    Complexf* b = check_complex(L, 2);
    lua_pushnumber(L, Dot(*a, *b));
    return 1;
}

static int complex_unm(lua_State* L)
{
    Complexf* a = check_complex(L, 1);
    new_complex(L, result);
    *result = -*a;
    return 1;
}

static const luaL_Reg complex_meta_funcs [] = {
    {"__index", complex_index},
    {"__newindex", complex_newindex},
	{"__len", complex_lenop},
    {"__tostring", complex_tostring},
    {"__add", complex_add},
    {"__sub", complex_sub},
    {"__mul", complex_mul},
    {"__div", complex_div},
    {"__pow", complex_dot},  // ^ as dot product
    {"__unm", complex_unm},  // unary minus
	{NULL, NULL}
};

// assumes the "abaci" table is the top of the stack.
void init_complex(lua_State* L)
{
    // registers all of complex_class_funcs functions in the abaci.complex table.
    lua_newtable(L);
    luaL_register(L, NULL, complex_class_funcs);
    lua_setfield(L, -2, "complex");

    // metatable for use with complex userdata.
	luaL_newmetatable(L, "abaci.complex");

    // registers all complex_meta_funcs functions in complex_mt
	luaL_register(L, NULL, complex_meta_funcs);
    lua_pop(L, 1);
}
