require "abaci"
require "math"
local vec2 = abaci.vec2
local vec3 = abaci.vec3
local vec4 = abaci.vec4
local quat = abaci.quat
local complex = abaci.complex
local eq = abaci.fuzzy_equal

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

    -- test mul by scalar
    prod = a * 2
    assert(eq(prod.x, a.x * 2))
    assert(eq(prod.y, a.y * 2))
    prod = 2 * a
    assert(eq(prod.x, a.x * 2))
    assert(eq(prod.y, a.y * 2))

    -- div
    quo = a / b
    assert(eq(quo.x, a.x / b.x))
    assert(eq(quo.y, a.y / b.y))

    -- test div by scalar
    prod = a / 2
    assert(eq(prod.x, a.x / 2))
    assert(eq(prod.y, a.y / 2))
    prod = 2 / a
    assert(eq(prod.x, 2 / a.x))
    assert(eq(prod.y, 2 / a.y))

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

--
-- quat
--
function quat_test()
    -- new
    a = quat.new(1, 2, 3, 4)

    -- __index
    assert(a.i == 1)
    assert(a.j == 2)
    assert(a.k == 3)
    assert(a.r == 4)

    -- __newindex
    a.i = 0.1
    a.j = 0.2
    a.k = 0.3
    a.r = 0.4
    assert(eq(0.1, a.i))
    assert(eq(0.2, a.j))
    assert(eq(0.3, a.k))
    assert(eq(0.4, a.r))

    -- len
    b = quat.new(0.3, 0.5, 0.7, 0.9)
    assert(eq(math.sqrt(b.i * b.i + b.j * b.j + b.k * b.k + b.r * b.r), b:len()))
    assert(eq(math.sqrt(a.i * a.i + a.j * a.j + a.k * a.k + a.r * a.r), a:len()))

    -- lerp
    l = quat.lerp(a, b, 0.3)
    assert(eq(abaci.lerp(a.i, b.i, 0.3), l.i))
    assert(eq(abaci.lerp(a.j, b.j, 0.3), l.j))
    assert(eq(abaci.lerp(a.k, b.k, 0.3), l.k))
    assert(eq(abaci.lerp(a.r, b.r, 0.3), l.r))

    --print("a = ", a)
    --print("b = ", b)

    -- add
    sum = a + b
    assert(eq(sum.i, a.i + b.i))
    assert(eq(sum.j, a.j + b.j))
    assert(eq(sum.k, a.k + b.k))
    assert(eq(sum.r, a.r + b.r))

    -- sub
    dif = a - b
    assert(eq(dif.i, a.i - b.i))
    assert(eq(dif.j, a.j - b.j))
    assert(eq(dif.k, a.k - b.k))
    assert(eq(dif.r, a.r - b.r))

    -- mul
    prod = a * b
	assert(eq(prod.i, a.i * b.r + a.j * b.k - a.k * b.j + a.r * b.i))
    assert(eq(prod.j, -a.i * b.k + a.j * b.r + a.k * b.i + a.r * b.j))
    assert(eq(prod.k, a.i * b.j - a.j * b.i + a.k * b.r + a.r * b.k))
    assert(eq(prod.r, -a.i * b.i - a.j * b.j - a.k * b.k + a.r * b.r))

    -- dot
    dot = a ^ b
    assert(eq(dot, a.i * b.i + a.j * b.j + a.k * b.k + a.r * b.r))

    -- unm
    neg = -b
    assert(eq(-b.i, neg.i))
    assert(eq(-b.j, neg.j))
    assert(eq(-b.k, neg.k))
    assert(eq(-b.r, neg.r))

    -- len
    assert(eq(#a, 4))

    -- conj
    a_conj = a:conj()
    assert(-a.i, a_conj.i)
    assert(-a.j, a_conj.j)
    assert(-a.k, a_conj.k)
    assert(a.r, a_conj.r)

    -- axis angle
    rot_x_90 = quat.axis_angle(vec3.new(1, 0, 0), math.pi/2)
    rot_y_90 = quat.axis_angle(0, 1, 0, math.pi/2)
    rot_z_90 = quat.axis_angle(vec3.new(0, 0, 1), math.pi/2)

    assert(eq(rot_x_90.i, math.sin(math.pi/4)))
    assert(eq(rot_x_90.j, 0))
    assert(eq(rot_x_90.k, 0))
    assert(eq(rot_x_90.r, math.cos(math.pi/4)))

    assert(eq(rot_y_90.i, 0))
    assert(eq(rot_y_90.j, math.sin(math.pi/4)))
    assert(eq(rot_y_90.k, 0))
    assert(eq(rot_y_90.r, math.cos(math.pi/4)))

    assert(eq(rot_z_90.i, 0))
    assert(eq(rot_z_90.j, 0))
    assert(eq(rot_z_90.k, math.sin(math.pi/4)))
    assert(eq(rot_y_90.r, math.cos(math.pi/4)))

    -- log
    q = quat.log(rot_x_90)
    assert(eq(q.i, math.pi/2))
    assert(eq(q.j, 0))
    assert(eq(q.k, 0))
    assert(eq(q.r, 0))
    
    -- exp
    q = quat.exp(q)
    assert(eq(q.i, rot_x_90.i))
    assert(eq(q.j, rot_x_90.j))
    assert(eq(q.k, rot_x_90.k))
    assert(eq(q.r, rot_x_90.r))

    -- rotate
    x_axis = vec3.new(1, 0, 0)
    x_prime = rot_y_90:rotate(x_axis)
    assert(eq(0, x_prime.x))
    assert(eq(0, x_prime.y))
    assert(eq(-1, x_prime.z))
end

--
-- complex
--
function complex_test()
    -- new
    a = complex.new(1, 2)
    b = complex.new(vec2.new(1, 2))

    -- __index
    assert(a.r == 1)
    assert(a.i == 2)
    assert(b.r == 1)
    assert(b.i == 2)

    -- __newindex
    a.r = 0.1
    a.i = 0.2
    assert(eq(0.1, a.r))
    assert(eq(0.2, a.i))

    -- len
    b = complex.new(0.3, 0.5)
    assert(eq(math.sqrt(b.r * b.r + b.i * b.i), b:len()))
    assert(eq(math.sqrt(a.r * a.r + a.i * a.i), a:len()))

    -- min_len
    c = b:min_len(0.5)
    d = b:min_len(100.0)
    assert(eq(c:len(), 0.5))
    assert(eq(d:len(), b:len()))

    -- add
    sum = a + b
    assert(eq(sum.r, a.r + b.r))
    assert(eq(sum.i, a.i + b.i))

    -- sub
    dif = a - b
    assert(eq(dif.r, a.r - b.r))
    assert(eq(dif.i, a.i - b.i))

    -- mul
    prod = a * b

    -- (a + bi) * (c + di) = (ac - bd) + (bc + ad)i
    aa, bb = a.r, a.i
    cc, dd = b.r, b.i

    assert(eq(prod.r, aa*cc - bb*dd))
    assert(eq(prod.i, bb*cc + aa*dd))

    -- div
    quo = a / b
    -- (a + bi) / (c + di) = ((ac + bd) / (c^2 + d^2)) + ((bc - ad) / (c^2 + d^2))i
    aa, bb = a.r, a.i
    cc, dd = b.r, b.i
    denom = cc * cc + dd * dd

    assert(eq(quo.r, (aa*cc + bb*dd)/denom))
    assert(eq(quo.i, (bb*cc - aa*dd)/denom))

    -- dot
    dot = a ^ b
    assert(eq(dot, a.r * b.r + a.i * b.i))

    -- unm
    neg = -b
    assert(eq(-b.r, neg.r))
    assert(eq(-b.i, neg.i))

    -- conj
    conj = b:conj()
    assert(eq(b.r, conj.r))
    assert(eq(-b.i, conj.i))

    -- len
    assert(eq(#a, 2))

    -- exp  todo
    -- log  todo
    -- expi
    expi = complex.expi(math.pi/4)
    assert(eq(1, expi:len()))
    assert(eq(1/math.sqrt(2), expi.r))
    assert(eq(1/math.sqrt(2), expi.i))
end

abaci_test()
vec2_test()
vec3_test()
vec4_test()
quat_test()
complex_test()

print("Success!")
