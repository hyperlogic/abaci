#include "quat_ext.h"
#include "abaci.h"
#include <string.h>

//
// vec4 functions
//

static int quat_new(lua_State* L)
{
    lua_Number i = luaL_checknumber(L, 1);
    lua_Number j = luaL_checknumber(L, 2);
    lua_Number k = luaL_checknumber(L, 3);
    lua_Number r = luaL_checknumber(L, 4);
    new_quat(L, result);
    result->Set((float)i, (float)j, (float)k, (float)r);
    return 1;
}

static int quat_identity(lua_State* L)
{
    new_quat(L, result);
    *result = Quatf::Identity();
    return 1;
}

static int quat_lerp(lua_State* L)
{
    Quatf* a = check_quat(L, 1);
    Quatf* b = check_quat(L, 2);
    lua_Number t = luaL_checknumber(L, 3);

    new_quat(L, result);
    *result = Lerp(*a, *b, (float)t);
    return 1;
}

static int quat_axis_angle(lua_State* L)
{
    luaL_checkany(L, 1);
    luaL_checkany(L, 2);

    int key_type = lua_type(L, 1);
    if (key_type == LUA_TNUMBER)
    {
        lua_Number axis_x = luaL_checknumber(L, 1);
        lua_Number axis_y = luaL_checknumber(L, 2);
        lua_Number axis_z = luaL_checknumber(L, 3);
        lua_Number theta = luaL_checknumber(L, 4);
        new_quat(L, result);
        *result = Quatf::AxisAngle(Vector3f(axis_x, axis_y, axis_z), theta);
        return 1;
    }
    else if (key_type == LUA_TUSERDATA)
    {
        Vector3f* axis = check_vec3(L, 1);
        lua_Number theta = luaL_checknumber(L, 2);
        new_quat(L, result);
        *result = Quatf::AxisAngle(*axis, theta);
        return 1;
    }
    else
    {
        // error
        lua_pushstring(L, "quat: expected 4 numbers, or 1 vec3 and 1 number");
        lua_error(L);
        return 0;
    }
}

static const luaL_Reg quat_class_funcs [] = {
    {"new", quat_new},
    {"identity", quat_identity},
    {"lerp", quat_lerp},
    {"axis_angle", quat_axis_angle},
    {NULL, NULL}
};

//
// quat methods
//

static int quat_len(lua_State* L)
{
    Quatf* self = check_quat(L, 1);
    lua_pushnumber(L, self->Len());
    return 1;
}

static int quat_unit(lua_State* L)
{
    Quatf* self = check_quat(L, 1);
    new_quat(L, result);
    *result = self->Unit();
    return 1;
}

static int quat_conj(lua_State* L)
{
    Quatf* self = check_quat(L, 1);
    new_quat(L, result);
    *result = ~(*self);
    return 1;
}

static int quat_rotate(lua_State* L)
{
    Quatf* self = check_quat(L, 1);
    Vector3f* x = check_vec3(L, 2);
    new_vec3(L, result);
    *result = self->Rotate(*x);
    return 1;
}

static int quat_exp(lua_State* L)
{
    Quatf* x = check_quat(L, 1);
    new_quat(L, result);
    *result = Exp(*x);
    return 1;
}

static int quat_log(lua_State* L)
{
    Quatf* x = check_quat(L, 1);
    new_quat(L, result);
    *result = Log(*x);
    return 1;
}

static const luaL_Reg quat_method_funcs [] = {
    {"len", quat_len},
    {"unit", quat_unit},
    {"conj", quat_conj},
    {"rotate", quat_rotate},
    {"exp", quat_exp},
    {"log", quat_log},
    {NULL, NULL}
};

//
// quat meta methods
//

static int quat_index(lua_State* L)
{
    Quatf* v = check_quat(L, 1);
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
        lua_pushstring(L, "quat: num index out of range, must be 1, 2, 3 or 4");
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
            case 'i':
                lua_pushnumber(L, v->i);
                return 1;
            case 'j':
                lua_pushnumber(L, v->j);
                return 1;
            case 'k':
                lua_pushnumber(L, v->k);
                return 1;
            case 'r':
                lua_pushnumber(L, v->r);
                return 1;
            }
        }
        else
        {
            // search for string key in quat_method_funcs array.
            int i = 0;
            while (quat_method_funcs[i].name)
            {
                if (!strcmp(s, quat_method_funcs[i].name))
                {
                    lua_pushcfunction(L, quat_method_funcs[i].func);
                    return 1;
                }
                i++;
            }
        }

        // error
        lua_pushstring(L, "quat: unknown string key");
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

static int quat_newindex(lua_State* L)
{
    Quatf* v = check_quat(L, 1);
    luaL_checkany(L, 2);
    lua_Number value = luaL_checknumber(L, 3);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Integer i = lua_tointeger(L, 2);
        if (i > 0 && i <= 4)
        {
            (*v)[i] = value;
            lua_pushnumber(L, (*v)[i-1]);
            return 0;
        }

        // error
        lua_pushstring(L, "quat: num index out of range, must be 1, 2, 3 or 4");
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
            case 'i':
                v->i = value;
                return 0;
            case 'j':
                v->j = value;
                return 0;
            case 'k':
                v->k = value;
                return 0;
            case 'r':
                v->r = value;
                return 0;
            }
        }

        // error
        lua_pushstring(L, "quat: unknown string key, must be i, j, k or r");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "quat: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int quat_tostring(lua_State* L)
{
    Quatf* v = check_quat(L, 1);
    char temp[64];
    snprintf(temp, 64, "quat(%.5f, %.5f, %.5f, %.5f)", v->i, v->j, v->k, v->r);
    lua_pushstring(L, temp);
    return 1;
}

static int quat_lenop(lua_State* L)
{
    const int kLen = 4;
    lua_pushinteger(L, kLen);
    return 1;
}

BINARY_VEC_OP_FUNC(add, Quatf, quat, operator+)
BINARY_VEC_OP_FUNC(sub, Quatf, quat, operator-)
BINARY_VEC_OP_FUNC(mul, Quatf, quat, operator*)

static int quat_dot(lua_State* L)
{
    Quatf* a = check_quat(L, 1);
    Quatf* b = check_quat(L, 2);
    lua_pushnumber(L, Dot(*a, *b));
    return 1;
}

static int quat_unm(lua_State* L)
{
    Quatf* a = check_quat(L, 1);
    new_quat(L, result);
    *result = -*a;
    return 1;
}

static const luaL_Reg quat_meta_funcs [] = {
    {"__index", quat_index},
    {"__newindex", quat_newindex},
    {"__len", quat_lenop},
    {"__tostring", quat_tostring},
    {"__add", quat_add},
    {"__sub", quat_sub},
    {"__mul", quat_mul},
    {"__pow", quat_dot},  // ^ as dot product
    {"__unm", quat_unm},  // unary minus
    {NULL, NULL}
};

// assumes the "abaci" table is the top of the stack.
void init_quat(lua_State* L)
{
    // registers all of quat_class_funcs functions in the abaci.quat table.
    lua_newtable(L);
    luaL_register(L, NULL, quat_class_funcs);
    lua_setfield(L, -2, "quat");

    // metatable for use with vec4 userdata.
    luaL_newmetatable(L, "abaci.quat");

    // registers all quat_meta_funcs functions in quat_mt
    luaL_register(L, NULL, quat_meta_funcs);
    lua_pop(L, 1);
}
