#ifndef COMPLEX_EXT_H
#define COMPLEX_EXT_H

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

void init_complex(lua_State* L);

#endif
