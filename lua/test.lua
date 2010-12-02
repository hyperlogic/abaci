require "abaci"

print("hello world")

pi = abaci.deg_to_rad(180)
print(string.format("pi = %.5f", pi))

a = vec2.new(1, 2)

print(string.format("a = (%.5f, %.5f), len = %.5f", a:getx(), a:gety(), a:len()));
