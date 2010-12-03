require "abaci"
require "math"
local vec2 = abaci.vec2
local vec3 = abaci.vec3
local vec4 = abaci.vec4
local eq = abaci.fuzzy_equals

print("abaci test")

--
-- abaci
--

function abaci_test()

    -- test fuzzy_equals it's used by everything else below
    assert(eq(1.00001, 1.0))
    assert(not eq(1.001, 1.0))

    -- deg_to_rad
    assert(eq(math.pi, abaci.deg_to_rad(180)))

    -- rad_to_deg
    assert(eq(180, abaci.rad_to_deg(math.pi)))

    -- clamp
    min = 1.0
    max = 3.0
    assert(eq(2.0, abaci.clamp(2.0, min, max)))
    assert(eq(min, abaci.clamp(0.0, min, max)))
    assert(eq(max, abaci.clamp(9.0, min, max)))

    -- limit_pi
    delta = 0.1
    assert(eq(-math.pi + delta, abaci.limit_pi(math.pi + delta)))
    assert(eq(math.pi - delta, abaci.limit_pi(-math.pi - delta)))

    -- mod_two_pi
    assert(eq(delta, abaci.mod_two_pi(2 * math.pi + delta)))
    assert(eq(2 * math.pi - delta, abaci.mod_two_pi(-delta)))

    -- lerp
    assert(eq(10.0, abaci.lerp(0, 100, 0.1)))
    assert(eq(90.0, abaci.lerp(100, 0, 0.1)))

    -- rand_int
    min = 10
    max = 100
    for i=1,10000 do
        r = abaci.rand_int(min, max)
        assert(r >= min and r <= max)
    end

    -- rand
    min = 10
    max = 100
    for i=1,10000 do
        r = abaci.rand(min, max)
        assert(r >= min and r <= max)
    end
end

--
-- vec2
--
function vec2_test()
    -- new
    a = vec2.new(1, 2)

    -- __index
    assert(a.x == 1)
    assert(a.y == 2)

    -- __newindex
    a.x = 0.1
    a.y = 0.2
    assert(eq(0.1, a.x))
    assert(eq(0.2, a.y))

    -- len
    b = vec2.new(0.3, 0.5)
    assert(eq(math.sqrt(b.x * b.x + b.y * b.y), b:len()))
    assert(eq(math.sqrt(a.x * a.x + a.y * a.y), a:len()))

    -- random_unit
    b = vec2.random_unit()
    assert(eq(1, math.sqrt(b.x * b.x + b.y * b.y)))

    -- lerp
    l = vec2.lerp(a, b, 0.3)
    assert(eq(abaci.lerp(a.x, b.x, 0.3), l.x))
    assert(eq(abaci.lerp(a.y, b.y, 0.3), l.y))

    -- min_len
    c = b:min_len(0.5)
    d = b:min_len(100.0)
    assert(eq(c:len(), 0.5))
    assert(eq(d:len(), b:len()))

    --print("a = ", a)
    --print("b = ", b)

    -- add
    sum = a + b
    assert(eq(sum.x, a.x + b.x))
    assert(eq(sum.y, a.y + b.y))

    -- sub
    dif = a - b
    assert(eq(dif.x, a.x - b.x))
    assert(eq(dif.y, a.y - b.y))

    -- mul
    prod = a * b
    assert(eq(prod.x, a.x * b.x))
    assert(eq(prod.y, a.y * b.y))

    -- div
    quo = a / b
    assert(eq(quo.x, a.x / b.x))
    assert(eq(quo.y, a.y / b.y))

    -- dot
    dot = a ^ b
    assert(eq(dot, a.x * b.x + a.y * b.y))

    -- unm
    neg = -b
    assert(eq(-b.x, neg.x))
    assert(eq(-b.y, neg.y))

    -- len
    assert(eq(#a, 2))
end

--
-- vec3
--
function vec3_test()

    -- new
    a = vec3.new(1, 2, 3)

    -- __index
    assert(a.x == 1)
    assert(a.y == 2)
    assert(a.z == 3)

    -- __newindex
    a.x = 0.1
    a.y = 0.2
    a.z = 0.3
    assert(eq(0.1, a.x))
    assert(eq(0.2, a.y))
    assert(eq(0.3, a.z))

    -- len
    b = vec3.new(0.3, 0.5, 0.7)
    assert(eq(math.sqrt(b.x * b.x + b.y * b.y + b.z * b.z), b:len()))
    assert(eq(math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z), a:len()))

    -- random_unit
    b = vec3.random_unit()
    assert(eq(1, math.sqrt(b.x * b.x + b.y * b.y + b.z * b.z)))

    -- lerp
    l = vec3.lerp(a, b, 0.3)
    assert(eq(abaci.lerp(a.x, b.x, 0.3), l.x))
    assert(eq(abaci.lerp(a.y, b.y, 0.3), l.y))
    assert(eq(abaci.lerp(a.z, b.z, 0.3), l.z))

    -- min_len
    c = b:min_len(0.5)
    d = b:min_len(100.0)
    assert(eq(c:len(), 0.5))
    assert(eq(d:len(), b:len()))

    --print("a = ", a)
    --print("b = ", b)

    -- add
    sum = a + b
    assert(eq(sum.x, a.x + b.x))
    assert(eq(sum.y, a.y + b.y))
    assert(eq(sum.z, a.z + b.z))

    -- sub
    dif = a - b
    assert(eq(dif.x, a.x - b.x))
    assert(eq(dif.y, a.y - b.y))
    assert(eq(dif.z, a.z - b.z))

    -- mul
    prod = a * b
    assert(eq(prod.x, a.x * b.x))
    assert(eq(prod.y, a.y * b.y))
    assert(eq(prod.z, a.z * b.z))

    -- div
    quo = a / b
    assert(eq(quo.x, a.x / b.x))
    assert(eq(quo.y, a.y / b.y))
    assert(eq(quo.z, a.z / b.z))

    -- dot
    dot = a ^ b
    assert(eq(dot, a.x * b.x + a.y * b.y + a.z * b.z))

    -- unm
    neg = -b
    assert(eq(-b.x, neg.x))
    assert(eq(-b.y, neg.y))
    assert(eq(-b.z, neg.z))

    -- cross
    cross = a % b
    -- a × b = (a2b3 − a3b2) i + (a3b1 − a1b3) j + (a1b2 − a2b1)
    assert(eq(a[2]*b[3] - a[3]*b[2], cross.x))
    assert(eq(a[3]*b[1] - a[1]*b[3], cross.y))
    assert(eq(a[1]*b[2] - a[2]*b[1], cross.z))

    -- len
    assert(eq(#a, 3))
end

--
-- vec4
--
function vec4_test()
    -- new
    a = vec4.new(1, 2, 3, 4)

    -- __index
    assert(a.x == 1)
    assert(a.y == 2)
    assert(a.z == 3)
    assert(a.w == 4)

    -- __newindex
    a.x = 0.1
    a.y = 0.2
    a.z = 0.3
    a.w = 0.4
    assert(eq(0.1, a.x))
    assert(eq(0.2, a.y))
    assert(eq(0.3, a.z))
    assert(eq(0.4, a.w))

    -- len
    b = vec4.new(0.3, 0.5, 0.7, 0.9)
    assert(eq(math.sqrt(b.x * b.x + b.y * b.y + b.z * b.z + b.w * b.w), b:len()))
    assert(eq(math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w), a:len()))

    -- lerp
    l = vec4.lerp(a, b, 0.3)
    assert(eq(abaci.lerp(a.x, b.x, 0.3), l.x))
    assert(eq(abaci.lerp(a.y, b.y, 0.3), l.y))
    assert(eq(abaci.lerp(a.z, b.z, 0.3), l.z))
    assert(eq(abaci.lerp(a.w, b.w, 0.3), l.w))

    -- min_len
    c = b:min_len(0.5)
    d = b:min_len(100.0)
    assert(eq(c:len(), 0.5))
    assert(eq(d:len(), b:len()))

    --print("a = ", a)
    --print("b = ", b)

    -- add
    sum = a + b
    assert(eq(sum.x, a.x + b.x))
    assert(eq(sum.y, a.y + b.y))
    assert(eq(sum.z, a.z + b.z))
    assert(eq(sum.w, a.w + b.w))

    -- sub
    dif = a - b
    assert(eq(dif.x, a.x - b.x))
    assert(eq(dif.y, a.y - b.y))
    assert(eq(dif.z, a.z - b.z))
    assert(eq(dif.w, a.w - b.w))

    -- mul
    prod = a * b
    assert(eq(prod.x, a.x * b.x))
    assert(eq(prod.y, a.y * b.y))
    assert(eq(prod.z, a.z * b.z))
    assert(eq(prod.w, a.w * b.w))

    -- div
    quo = a / b
    assert(eq(quo.x, a.x / b.x))
    assert(eq(quo.y, a.y / b.y))
    assert(eq(quo.z, a.z / b.z))
    assert(eq(quo.w, a.w / b.w))

    -- dot
    dot = a ^ b
    assert(eq(dot, a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w))

    -- unm
    neg = -b
    assert(eq(-b.x, neg.x))
    assert(eq(-b.y, neg.y))
    assert(eq(-b.z, neg.z))
    assert(eq(-b.w, neg.w))

    -- len
    assert(eq(#a, 4))
end

abaci_test()
vec2_test()
vec3_test()
vec4_test()

print("Success!")
