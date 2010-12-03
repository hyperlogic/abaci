require "abaci"
require "math"
local vec2 = abaci.vec2
local eq = abaci.fuzzy_equals

print("abaci test")

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


print("Success!")
