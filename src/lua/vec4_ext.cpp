#include "vec4_ext.h"
#include "abaci.h"
#include <string.h>

//
// vec4 functions
//

static int vec4_new(lua_State* L)
{
    int argCount = lua_gettop(L);
    if (argCount == 1)
    {
        luaL_checktype(L, 1, LUA_TTABLE);
        lua_rawgeti(L, 1, 1);
        lua_rawgeti(L, 1, 2);
        lua_rawgeti(L, 1, 3);
        lua_rawgeti(L, 1, 4);
        lua_Number x = luaL_checknumber(L, 2);
        lua_Number y = luaL_checknumber(L, 3);
        lua_Number z = luaL_checknumber(L, 4);
        lua_Number w = luaL_checknumber(L, 5);
        new_vec4(L, result);
        result->Set((float)x, (float)y, (float)z, (float)w);
        return 1;
    }
    else if (argCount == 4)
    {
        lua_Number x = luaL_checknumber(L, 1);
        lua_Number y = luaL_checknumber(L, 2);
        lua_Number z = luaL_checknumber(L, 3);
        lua_Number w = luaL_checknumber(L, 4);
        new_vec4(L, result);
        result->Set((float)x, (float)y, (float)z, (float)w);
        return 1;
    }
    else
    {
        lua_pushstring(L, "vec4.new: expected four numbers or a single table");
        lua_error(L);
        return 0;
    }
}

static int vec4_lerp(lua_State* L)
{
    Vector4f* a = check_vec4(L, 1);
    Vector4f* b = check_vec4(L, 2);
    lua_Number t = luaL_checknumber(L, 3);

    new_vec4(L, result);
    *result = Lerp(*a, *b, (float)t);
    return 1;
}

static const luaL_Reg vec4_class_funcs [] = {
    {"new", vec4_new},
    {"lerp", vec4_lerp},
    {NULL, NULL}
};

//
// vec4 methods
//

static int vec4_len(lua_State* L)
{
    Vector4f* self = check_vec4(L, 1);
    lua_pushnumber(L, self->Len());
    return 1;
}

static int vec4_unit(lua_State* L)
{
    Vector4f* self = check_vec4(L, 1);
    new_vec4(L, result);
    *result = self->Unit();
    return 1;
}

static int vec4_min_len(lua_State* L)
{
    Vector4f* self = check_vec4(L, 1);
    lua_Number len = luaL_checknumber(L, 2);

    new_vec4(L, result);
    *result = self->MinLen(len);
    return 1;
}

static const luaL_Reg vec4_method_funcs [] = {
    {"len", vec4_len},
    {"unit", vec4_unit},
    {"min_len", vec4_min_len},
    {NULL, NULL}
};

//
// vec4 meta methods
//

static int vec4_index(lua_State* L)
{
    Vector4f* v = check_vec4(L, 1);
    luaL_checkany(L, 2);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i > 0 && i <= 4)
        {
            lua_pushnumber(L, (*v)[i-1]);
            return 1;
        }

        // error
        lua_pushstring(L, "vec4: num index out of range, must be 1, 2, 3 or 4");
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
            case 'z':
                lua_pushnumber(L, v->z);
                return 1;
            case 'w':
                lua_pushnumber(L, v->w);
                return 1;
            }
        }
        else
        {
            // search for string key in vec4_method_funcs array.
            int i = 0;
            while (vec4_method_funcs[i].name)
            {
                if (!strcmp(s, vec4_method_funcs[i].name))
                {
                    lua_pushcfunction(L, vec4_method_funcs[i].func);
                    return 1;
                }
                i++;
            }
        }

        // error
        lua_pushstring(L, "vec4: unknown string key");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "vec4: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int vec4_newindex(lua_State* L)
{
    Vector4f* v = check_vec4(L, 1);
    luaL_checkany(L, 2);
    lua_Number value = luaL_checknumber(L, 3);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i > 0 && i <= 4)
        {
            (*v)[i] = value;
            return 0;
        }

        // error
        lua_pushstring(L, "vec4: num index out of range, must be 1, 2, 3 or 4");
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
            case 'z':
                v->z = value;
                return 0;
            case 'w':
                v->w = value;
                return 0;
            }
        }

        // error
        lua_pushstring(L, "vec4: unknown string key, must be x, y, z or w");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "vec4: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int vec4_tostring(lua_State* L)
{
    Vector4f* v = check_vec4(L, 1);
    char temp[64];
    snprintf(temp, 64, "vec4(%.5f, %.5f, %.5f, %.5f)", v->x, v->y, v->z, v->w);
    lua_pushstring(L, temp);
    return 1;
}

static int vec4_lenop(lua_State* L)
{
    const int kLen = 4;
    lua_pushinteger(L, kLen);
    return 1;
}

BINARY_VEC_OP_FUNC(add, Vector4f, vec4, operator+)
BINARY_VEC_OP_FUNC(sub, Vector4f, vec4, operator-)
BINARY_VEC_OP_FUNC2(mul, Vector4f, vec4, operator*)
BINARY_VEC_OP_FUNC2(div, Vector4f, vec4, operator/)

static int vec4_dot(lua_State* L)
{
    Vector4f* a = check_vec4(L, 1);
    Vector4f* b = check_vec4(L, 2);
    lua_pushnumber(L, Dot(*a, *b));
    return 1;
}

static int vec4_unm(lua_State* L)
{
    Vector4f* a = check_vec4(L, 1);
    new_vec4(L, result);
    *result = -*a;
    return 1;
}

static const luaL_Reg vec4_meta_funcs [] = {
    {"__index", vec4_index},
    {"__newindex", vec4_newindex},
    {"__len", vec4_lenop},
    {"__tostring", vec4_tostring},
    {"__add", vec4_add},
    {"__sub", vec4_sub},
    {"__mul", vec4_mul},
    {"__div", vec4_div},
    {"__pow", vec4_dot},  // ^ as dot product
    {"__unm", vec4_unm},  // unary minus
    {NULL, NULL}
};

// assumes the "abaci" table is the top of the stack.
void init_vec4(lua_State* L)
{
    // registers all of vec4_class_funcs functions in the abaci.vec4 table.
    lua_newtable(L);
    luaL_register(L, NULL, vec4_class_funcs);
    lua_setfield(L, -2, "vec4");

    // metatable for use with vec4 userdata.
    luaL_newmetatable(L, "abaci.vec4");

    // registers all vec4_meta_funcs functions in vec4_mt
    luaL_register(L, NULL, vec4_meta_funcs);
    lua_pop(L, 1);
}
