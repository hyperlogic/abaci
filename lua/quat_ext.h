#ifndef QUAT_EXT_H
#define QUAT_EXT_H

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

void init_quat(lua_State* L);

#endif
