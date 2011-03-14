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

#define check_matrix(L, n) (Matrixf*)luaL_checkudata(L, n, "abaci.matrix")
#define new_matrix(L, name)                                           \
    Matrixf* name = (Matrixf*)lua_newuserdata(L, sizeof(Matrixf));    \
    luaL_getmetatable(L, "abaci.matrix");                             \
    lua_setmetatable(L, -2)


// this version only supports vec * vec
#define BINARY_VEC_OP_FUNC(func, c_type, lua_type, op)              \
    static int lua_type##_##func(lua_State* L)                      \
    {                                                               \
        c_type * a = check_##lua_type(L, 1);                        \
        c_type * b = check_##lua_type(L, 2);                        \
        new_##lua_type(L, result);                                  \
        *result = op(*a, *b);                                       \
        return 1;                                                   \
    }

// this version supports vec * num and num * vec
#define BINARY_VEC_OP_FUNC2(func, c_type, lua_typename, op)          \
    static int lua_typename##_##func(lua_State* L)                   \
    {                                                                \
        luaL_checkany(L, 1);                                         \
        luaL_checkany(L, 2);                                         \
        int a_type = lua_type(L, 1);                                 \
        int b_type = lua_type(L, 2);                                 \
        if (a_type == LUA_TUSERDATA && b_type == LUA_TUSERDATA)      \
        {                                                            \
            c_type * a = check_##lua_typename(L, 1);                 \
            c_type * b = check_##lua_typename(L, 2);                 \
            new_##lua_typename(L, result);                           \
            *result = op(*a, *b);                                    \
            return 1;                                                \
        }                                                            \
        else if (a_type == LUA_TNUMBER && b_type == LUA_TUSERDATA)   \
        {                                                            \
            lua_Number a = luaL_checknumber(L, 1);                   \
            c_type * b = check_##lua_typename(L, 2);                 \
            new_##lua_typename(L, result);                           \
            *result = op((float)a, *b);                              \
            return 1;                                                \
        }                                                            \
        else if (a_type == LUA_TUSERDATA && b_type == LUA_TNUMBER)   \
        {                                                            \
            c_type * a = check_##lua_typename(L, 1);                 \
            lua_Number b = luaL_checknumber(L, 2);                   \
            new_##lua_typename(L, result);                           \
            *result = op(*a, (float)b);                              \
            return 1;                                                \
        }                                                            \
        else                                                         \
        {                                                            \
            lua_pushstring(L, " lua_type : bad args");               \
            lua_error(L);                                            \
            return 0;                                                \
        }                                                            \
    }


#endif
