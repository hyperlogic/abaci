#ifndef ABACI_EXT_H
#define ABACI_EXT_H

extern "C" {
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}

int luaopen_abaci(lua_State* L);

#endif
