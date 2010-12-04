#include "vec3_ext.h"
#include "abaci.h"
#include <string.h>

//
// vec3 functions
//

static int vec3_new(lua_State* L)
{
    int argCount = lua_gettop(L);
    if (argCount == 1)
    {
        luaL_checktype(L, 1, LUA_TTABLE);
        lua_rawgeti(L, 1, 1);
        lua_rawgeti(L, 1, 2);
        lua_rawgeti(L, 1, 3);
        lua_Number x = luaL_checknumber(L, 2);
        lua_Number y = luaL_checknumber(L, 3);
        lua_Number z = luaL_checknumber(L, 4);
        new_vec3(L, result);
        result->Set((float)x, (float)y, (float)z);
        return 1;
    }
    else if (argCount == 3)
    {
        lua_Number x = luaL_checknumber(L, 1);
        lua_Number y = luaL_checknumber(L, 2);
        lua_Number z = luaL_checknumber(L, 3);
        new_vec3(L, result);
        result->Set((float)x, (float)y, (float)z);
        return 1;
    }
    else
    {
        lua_pushstring(L, "vec3.new: expected three numbers or a single table");
        lua_error(L);
        return 0;
    }
}

static int vec3_random_unit(lua_State* L)
{
    new_vec3(L, result);
    *result = Vector3f::RandomUnitVector();
    return 1;
}

static int vec3_lerp(lua_State* L)
{
    Vector3f* a = check_vec3(L, 1);
    Vector3f* b = check_vec3(L, 2);
    lua_Number t = luaL_checknumber(L, 3);

    new_vec3(L, result);
    *result = Lerp(*a, *b, (float)t);
    return 1;
}

static const luaL_Reg vec3_class_funcs [] = {
	{"new", vec3_new},
    {"random_unit", vec3_random_unit},
    {"lerp", vec3_lerp},
	{NULL, NULL}
};

//
// vec3 methods
//

static int vec3_len(lua_State* L)
{
    Vector3f* self = check_vec3(L, 1);
    lua_pushnumber(L, self->Len());
    return 1;
}

static int vec3_unit(lua_State* L)
{
    Vector3f* self = check_vec3(L, 1);
    new_vec3(L, result);
    *result = self->Unit();
    return 1;
}

static int vec3_min_len(lua_State* L)
{
    Vector3f* self = check_vec3(L, 1);
    lua_Number len = luaL_checknumber(L, 2);

    new_vec3(L, result);
    *result = self->MinLen(len);
    return 1;
}

static const luaL_Reg vec3_method_funcs [] = {
    {"len", vec3_len},
    {"unit", vec3_unit},
	{"min_len", vec3_min_len},
	{NULL, NULL}
};

//
// vec3 meta methods
//

static int vec3_index(lua_State* L)
{
    Vector3f* v = check_vec3(L, 1);
    luaL_checkany(L, 2);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i > 0 && i <= 3)
        {
            lua_pushnumber(L, (*v)[i-1]);
            return 1;
        }

        // error
        lua_pushstring(L, "vec3: num index out of range, must be 1, 2 or 3");
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
            }
        }
        else
        {
            // search for string key in vec3_method_funcs array.
            int i = 0;
            while (vec3_method_funcs[i].name)
            {
                if (!strcmp(s, vec3_method_funcs[i].name))
                {
                    lua_pushcfunction(L, vec3_method_funcs[i].func);
                    return 1;
                }
                i++;
            }
        }

        // error
        lua_pushstring(L, "vec3: unknown string key");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "vec3: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int vec3_newindex(lua_State* L)
{
    Vector3f* v = check_vec3(L, 1);
    luaL_checkany(L, 2);
    lua_Number value = luaL_checknumber(L, 3);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i > 0 && i <= 3)
        {
            (*v)[i] = value;
            return 0;
        }

        // error
        lua_pushstring(L, "vec3: num index out of range, must be 0, 1 or 2");
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
            }
        }

        // error
        lua_pushstring(L, "vec3: unknown string key, expecting x, y or z");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "vec3: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int vec3_tostring(lua_State* L)
{
    Vector3f* v = check_vec3(L, 1);
    char temp[64];
    snprintf(temp, 64, "vec3(%.5f, %.5f, %.5f)", v->x, v->y, v->z);
    lua_pushstring(L, temp);
    return 1;
}

static int vec3_lenop(lua_State* L)
{
    const int kLen = 3;
    lua_pushinteger(L, kLen);
    return 1;
}

BINARY_VEC_OP_FUNC(add, Vector3f, vec3, operator+)
BINARY_VEC_OP_FUNC(sub, Vector3f, vec3, operator-)
BINARY_VEC_OP_FUNC2(mul, Vector3f, vec3, operator*)
BINARY_VEC_OP_FUNC2(div, Vector3f, vec3, operator/)
BINARY_VEC_OP_FUNC(cross, Vector3f, vec3, Cross)

static int vec3_dot(lua_State* L)
{
    Vector3f* a = check_vec3(L, 1);
    Vector3f* b = check_vec3(L, 2);
    lua_pushnumber(L, Dot(*a, *b));
    return 1;
}

static int vec3_unm(lua_State* L)
{
    Vector3f* a = check_vec3(L, 1);
    new_vec3(L, result);
    *result = -*a;
    return 1;
}

static const luaL_Reg vec3_meta_funcs [] = {
    {"__index", vec3_index},
    {"__newindex", vec3_newindex},
	{"__len", vec3_lenop},
    {"__tostring", vec3_tostring},
    {"__add", vec3_add},
    {"__sub", vec3_sub},
    {"__mul", vec3_mul},
    {"__div", vec3_div},
    {"__pow", vec3_dot},    // ^ as dot product
    {"__unm", vec3_unm},    // unary minus
    {"__mod", vec3_cross},  // % as cross product
	{NULL, NULL}
};

// assumes the "abaci" table is the top of the stack.
void init_vec3(lua_State* L)
{
    // registers all of vec3_class_funcs functions in the abaci.vec3 table.
    lua_newtable(L);
    luaL_register(L, NULL, vec3_class_funcs);
    lua_setfield(L, -2, "vec3");

    // metatable for use with vec3 userdata.
	luaL_newmetatable(L, "abaci.vec3");

    // registers all vec3_meta_funcs functions in vec3_mt
	luaL_register(L, NULL, vec3_meta_funcs);
    lua_pop(L, 1);
}
