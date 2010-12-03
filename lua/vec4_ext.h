#ifndef VEC4_EXT_H
#define VEC4_EXT_H

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

void init_vec4(lua_State* L);

#endif
