#include "matrix_ext.h"
#include "abaci.h"
#include <string.h>

//
// matrix functions
//

static int matrix_ident(lua_State* L)
{
    new_matrix(L, result);
	*result = Matrixf::Identity();
	return 1;
}

static int matrix_rows(lua_State* L)
{
    Vector4f* row0 = check_vec4(L, 1);
    Vector4f* row1 = check_vec4(L, 2);
    Vector4f* row2 = check_vec4(L, 3);
    Vector4f* row3 = check_vec4(L, 4);
    new_matrix(L, result);
	*result = Matrixf::Rows(*row0, *row1, *row2, *row3);
	return 1;
}

static int matrix_axes(lua_State* L)
{
    int argCount = lua_gettop(L);
    Vector3f* xAxis = check_vec3(L, 1);
    Vector3f* yAxis = check_vec3(L, 2);
    Vector3f* zAxis = check_vec3(L, 3);

    // trans is optional
    Vector3f zero(0, 0, 0);
    Vector3f* trans = &zero;
    if (argCount == 4)
        trans = check_vec3(L, 4);
    
    new_matrix(L, result);
    *result = Matrixf::Axes(*xAxis, *yAxis, *zAxis, *trans);
    return 1;
}

static int matrix_trans(lua_State* L)
{
    Vector3f* trans = check_vec3(L, 1);
    new_matrix(L, result);
	*result = Matrixf::Trans(*trans);
	return 1;
}

static int matrix_quat(lua_State* L)
{
    Quatf* quat = check_quat(L, 1);

    // optional trans
    Vector3f zero(0,0,0);
    Vector3f* trans = &zero;
    int argCount = lua_gettop(L);
    if (argCount == 2)
        trans = check_vec3(L, 2);

    new_matrix(L, result);
	*result = Matrixf::QuatTrans(*quat, *trans);
	return 1;
}

static int matrix_axis_angle(lua_State* L)
{
    Vector3f* axis = check_vec3(L, 1);
    lua_Number theta = luaL_checknumber(L, 2);

    // optional trans
    Vector3f zero(0,0,0);
    Vector3f* trans = &zero;
    int argCount = lua_gettop(L);
    if (argCount == 3)
        trans = check_vec3(L, 3);

    new_matrix(L, result);
	*result = Matrixf::QuatTrans(Quatf::AxisAngle(*axis, theta), *trans);
	return 1;
}

static int matrix_scale_quat(lua_State* L)
{
    Vector3f one(1,1,1);
    Vector3f* scale = &one;
    luaL_checkany(L, 1);
    int key_type = lua_type(L, 1);
    if (key_type == LUA_TNUMBER)
    {
        lua_Number factor = luaL_checknumber(L, 1);
        scale->Set(factor, factor, factor);
    }
    else if (key_type == LUA_TUSERDATA)
    {
        scale = check_vec3(L, 1);
    }
    else
    {
        lua_pushstring(L, "matrix: first arg should be a num or a vec3");
        lua_error(L);
        return 0;
    }
    
    Quatf* quat = check_quat(L, 2);

    // optional trans
    Vector3f zero(0,0,0);
    Vector3f* trans = &zero;
    int argCount = lua_gettop(L);
    if (argCount == 3)
        trans = check_vec3(L, 3);

    new_matrix(L, result);
	*result = Matrixf::ScaleQuatTrans(*scale, *quat, *trans);
	return 1;
}

static int matrix_frustum(lua_State* L)
{
    lua_Number fovy = luaL_checknumber(L, 1);
    lua_Number aspect = luaL_checknumber(L, 2);
    lua_Number nearVal = luaL_checknumber(L, 3);
    lua_Number farVal = luaL_checknumber(L, 4);
    
    new_matrix(L, result);
	*result = Matrixf::Frustum(fovy, aspect, nearVal, farVal);
	return 1;
}

static int matrix_ortho(lua_State* L)
{
    lua_Number left = luaL_checknumber(L, 1);
    lua_Number right = luaL_checknumber(L, 2);
    lua_Number bottom = luaL_checknumber(L, 3);
    lua_Number top = luaL_checknumber(L, 4);
    lua_Number nearVal = luaL_checknumber(L, 5);
    lua_Number farVal = luaL_checknumber(L, 6);
    
    new_matrix(L, result);
	*result = Matrixf::Ortho(left, right, bottom, top, nearVal, farVal);
	return 1;
}

static int matrix_look_at(lua_State* L)
{
    Vector3f* eye = check_vec3(L, 1);
    Vector3f* target = check_vec3(L, 2);
    Vector3f* up = check_vec3(L, 3);
    new_matrix(L, result);
	*result = Matrixf::LookAt(*eye, *target, *up);
	return 1;
}

static const luaL_Reg matrix_class_funcs [] = {
	{"new", matrix_ident},
    {"identity", matrix_ident},
    {"rows", matrix_rows},
    {"axes", matrix_axes},
    {"trans", matrix_trans},
    {"quat", matrix_quat},
    {"axis_angle", matrix_axis_angle},
    {"scale_quat", matrix_scale_quat},
    {"frustum", matrix_frustum},
    {"ortho", matrix_ortho},
    {"look_at", matrix_look_at},
	{NULL, NULL}
};

//
// matrix methods
//

#define GETTER(name, c_name, l_type)            \
    static int matrix_##name (lua_State* L)     \
    {                                           \
        Matrixf* self = check_matrix(L, 1);     \
        new_##l_type(L, result);                \
        *result = self->c_name();               \
        return 1;                               \
    }

GETTER(get_xaxis, GetXAxis, vec3);
GETTER(get_yaxis, GetYAxis, vec3);
GETTER(get_zaxis, GetZAxis, vec3);
GETTER(get_trans, GetTrans, vec3);
GETTER(get_quat, GetQuat, quat);

static int matrix_get_row(lua_State* L)
{
    Matrixf* self = check_matrix(L, 1);
    int i = luaL_checkinteger(L, 2);
    if (i > 0 && i <= 4)
    {
        new_vec4(L, result);
        *result = self->GetRow(i-1);
        return 1;
    }

    // error
    lua_pushstring(L, "matrix: expected index between 1 and 4");
    lua_error(L);
    return 0;
}

GETTER(transpose, Transpose, matrix);
GETTER(inverse, FullInverse, matrix);
GETTER(ortho_inverse, OrthoInverse, matrix);

#define VEC3_SETTER(name, c_name)               \
    static int matrix_##name (lua_State* L)     \
    {                                           \
        Matrixf* self = check_matrix(L, 1);     \
        Vector3f* arg = check_vec3(L, 2);       \
        self->c_name(*arg);                     \
        return 0;                               \
    }

VEC3_SETTER(set_xaxis, SetXAxis)
VEC3_SETTER(set_yaxis, SetYAxis)
VEC3_SETTER(set_zaxis, SetZAxis)
VEC3_SETTER(set_trans, SetTrans)

static int matrix_set_scale(lua_State* L)
{
    Matrixf* self = check_matrix(L, 1);
    luaL_checkany(L, 2);
    int key_type = lua_type(L, 2);
    if (key_type == LUA_TNUMBER)
    {
        lua_Number scale = luaL_checknumber(L, 2);
        self->SetScale(scale);
        return 0;
    }
    else if (key_type == LUA_TUSERDATA)
    {
        Vector3f* scale = check_vec3(L, 2);
        self->SetScale(*scale);
        return 0;
    }

    // error
    lua_pushstring(L, "matrix: expected num or vec3");
    lua_error(L);
    return 0;
}

static int matrix_mul3x3(lua_State* L)
{
    Matrixf* self = check_matrix(L, 1);
    Vector3f* arg = check_vec3(L, 2);
    new_vec3(L, result);
    *result = self->Mul3x3(*arg);
    return 1;
}

static int matrix_mul3x4(lua_State* L)
{
    Matrixf* self = check_matrix(L, 1);
    Vector3f* arg = check_vec3(L, 2);
    new_vec3(L, result);
    *result = self->Mul3x4(*arg);
    return 1;
}

static int matrix_mul4x4(lua_State* L)
{
    Matrixf* self = check_matrix(L, 1);
    Vector4f* arg = check_vec4(L, 2);
    new_vec4(L, result);
    *result = self->Mul4x4(*arg);
    return 1;
}

static const luaL_Reg matrix_method_funcs [] = {
    {"get_xaxis", matrix_get_xaxis},
    {"get_yaxis", matrix_get_yaxis},
    {"get_zaxis", matrix_get_zaxis},
    {"get_trans", matrix_get_trans},
    {"get_quat", matrix_get_quat},
    {"get_row", matrix_get_row},
    {"transpose", matrix_transpose},
    {"inverse", matrix_inverse},
    {"ortho_inverse", matrix_ortho_inverse},
    {"set_scale", matrix_set_scale},
    {"set_xaxis", matrix_set_xaxis},
    {"set_yaxis", matrix_set_yaxis},
    {"set_zaxis", matrix_set_zaxis},
    {"set_trans", matrix_set_trans},
    {"mul3x3", matrix_mul3x3},
    {"mul3x4", matrix_mul3x4},
    {"mul4x4", matrix_mul4x4},
	{NULL, NULL}
};

//
// matrix meta methods
//

static int matrix_index(lua_State* L)
{
    //Matrixf* m = check_matrix(L, 1);
    luaL_checkany(L, 2);

    int key_type = lua_type(L, 2);
    if (key_type == LUA_TSTRING)
    {
        size_t size;
        const char* s = lua_tolstring(L, 2, &size);

        // search for string key in matrix_method_funcs array.
        int i = 0;
        while (matrix_method_funcs[i].name)
        {
            if (!strcmp(s, matrix_method_funcs[i].name))
            {
                lua_pushcfunction(L, matrix_method_funcs[i].func);
                return 1;
            }
            i++;
        }

        // error
        lua_pushstring(L, "matrix: unknown string key");
        lua_error(L);
    }
    else
    {
        // error
        lua_pushstring(L, "matrix: unsupported key");
        lua_error(L);
    }

    return 0;
}

static int matrix_newindex(lua_State* L)
{
    // error
    lua_pushstring(L, "matrix: read only attributes");
    lua_error(L);

    return 0;
}

static int matrix_tostring(lua_State* L)
{
    Matrixf* m = check_matrix(L, 1);
    char temp[512];
    Vector4f row0 = m->GetRow(0);
    Vector4f row1 = m->GetRow(1);
    Vector4f row2 = m->GetRow(2);
    Vector4f row3 = m->GetRow(3);
    snprintf(temp, 512, "matrix((%.3f, %.3f, %.3f, %.3f), (%.3f, %.3f, %.3f, %.3f), (%.3f, %.3f, %.3f, %.3f), (%.3f, %.3f, %.3f, %.3f))",
             row0.x, row0.y, row0.z, row0.w, 
             row1.x, row1.y, row1.z, row1.w, 
             row2.x, row2.y, row2.z, row2.w, 
             row3.x, row3.y, row3.z, row3.w);
    lua_pushstring(L, temp);
    return 1;
}

BINARY_VEC_OP_FUNC(add, Matrixf, matrix, operator+)
BINARY_VEC_OP_FUNC(sub, Matrixf, matrix, operator-)
BINARY_VEC_OP_FUNC(mul, Matrixf, matrix, operator*)

static const luaL_Reg matrix_meta_funcs [] = {
    {"__index", matrix_index},
    {"__newindex", matrix_newindex},
    {"__tostring", matrix_tostring},
    {"__add", matrix_add},
    {"__sub", matrix_sub},
    {"__mul", matrix_mul},
	{NULL, NULL}
};

// assumes the "abaci" table is the top of the stack.
void init_matrix(lua_State* L)
{
    // registers all of matrix_class_funcs functions in the abaci.matrix table.
    lua_newtable(L);
    luaL_register(L, NULL, matrix_class_funcs);
    lua_setfield(L, -2, "matrix");

    // metatable for use with matrix userdata.
	luaL_newmetatable(L, "abaci.matrix");

    // registers all matrix_meta_funcs functions in matrix_mt
	luaL_register(L, NULL, matrix_meta_funcs);
    lua_pop(L, 1);
}
