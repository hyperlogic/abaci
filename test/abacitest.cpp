#include "abaci.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "unittest.h"

#ifdef ABACI_NAMESPACE
using namespace ABACI_NAMESPACE;
#endif

// static test containers
std::vector<float> s_floatVec;
std::vector<Vector2f> s_vector2fVec;
std::vector<Vector3f> s_vector3fVec;
std::vector<Vector4f> s_vector4fVec;
std::vector<Quatf> s_quatfVec;
std::vector<Matrixf> s_matrixfVec;

bool FuzzyFloatTest(float rhs, float lhs, float epsilon=0.0001f)
{
    // NaN's are equal, and so are Infs
    if ((isnan(rhs) && isnan(lhs)) ||
        (isinf(rhs) && isinf(lhs)))
        return true;
    else
        return FuzzyEqual(rhs, lhs, epsilon);
}

bool FloatTest(float rhs, float lhs)
{
    // For testing purposes, two Nans are "equal" even if they have different bits.
    if (isnan(rhs) && isnan(lhs))
        return true;

    return rhs == lhs;
}

// 0..1
float RandomNormalizedFloat()
{
    return (float)rand() / RAND_MAX;
}

// -range to range
float RandomFloat(float range = 100.0f)
{
    return (((float)rand() / (RAND_MAX / 2)) - 1.0f) * range;
}

void InitTestData()
{
    srand((unsigned int)666);

    // add random values
    const int RANDOM_DATA_SIZE = 128;
    for (int i = 0; i < RANDOM_DATA_SIZE; ++i)
    {
        s_floatVec.push_back(RandomFloat());
        s_vector2fVec.push_back(Vector2f(RandomFloat(), RandomFloat()));
        s_vector3fVec.push_back(Vector3f(RandomFloat(), RandomFloat(), RandomFloat()));
        s_vector4fVec.push_back(Vector4f(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat()));

        Vector4f row0(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
        Vector4f row1(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
        Vector4f row2(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
        Vector4f row3(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
        s_matrixfVec.push_back(Matrixf::Rows(row0, row1, row2, row3));
    }

    // zero
    s_floatVec.push_back(0);
    s_vector2fVec.push_back(Vector2f(0, 0));
    s_vector3fVec.push_back(Vector3f(0, 0, 0));
    s_vector4fVec.push_back(Vector4f(0, 0, 0, 0));
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(0, 0, 0, 0), Vector4f(0, 0, 0, 0), Vector4f(0, 0, 0, 0), Vector4f(0, 0, 0, 0)));

    // ident-matrix
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(1, 0, 0, 0), Vector4f(0, 1, 0, 0), Vector4f(0, 0, 1, 0), Vector4f(0, 0, 0, 1)));

    // FLT_MAX
    s_floatVec.push_back(FLT_MAX);
    s_vector2fVec.push_back(Vector2f(FLT_MAX, FLT_MAX));
    s_vector3fVec.push_back(Vector3f(FLT_MAX, FLT_MAX, FLT_MAX));
    s_vector4fVec.push_back(Vector4f(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX));
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX),
                                         Vector4f(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX),
                                         Vector4f(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX),
                                         Vector4f(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX)));

    s_floatVec.push_back(-FLT_MAX);
    s_vector2fVec.push_back(Vector2f(-FLT_MAX, -FLT_MAX));
    s_vector3fVec.push_back(Vector3f(-FLT_MAX, -FLT_MAX, -FLT_MAX));
    s_vector4fVec.push_back(Vector4f(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX));
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
                                         Vector4f(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
                                         Vector4f(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX),
                                         Vector4f(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)));

    // FLT_MIN
    s_floatVec.push_back(FLT_MIN);
    s_vector2fVec.push_back(Vector2f(FLT_MIN, FLT_MIN));
    s_vector3fVec.push_back(Vector3f(FLT_MIN, FLT_MIN, FLT_MIN));
    s_vector4fVec.push_back(Vector4f(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN));
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN),
                                         Vector4f(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN),
                                         Vector4f(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN),
                                         Vector4f(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN)));

    s_floatVec.push_back(-FLT_MIN);
    s_vector2fVec.push_back(Vector2f(-FLT_MIN, -FLT_MIN));
    s_vector3fVec.push_back(Vector3f(-FLT_MIN, -FLT_MIN, -FLT_MIN));
    s_vector4fVec.push_back(Vector4f(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN));
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN),
                                         Vector4f(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN),
                                         Vector4f(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN),
                                         Vector4f(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN)));

    // inf
    float inf = FLT_MAX * FLT_MAX;
    s_floatVec.push_back(inf);
    s_vector2fVec.push_back(Vector2f(inf, inf));
    s_vector3fVec.push_back(Vector3f(inf, inf, inf));
    s_vector4fVec.push_back(Vector4f(inf, inf, inf, inf));
    /*
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(inf, inf, inf, inf),
                                         Vector4f(inf, inf, inf, inf),
                                         Vector4f(inf, inf, inf, inf),
                                         Vector4f(inf, inf, inf, inf)));
    */

    // -inf
    inf = -FLT_MAX * FLT_MAX;
    s_floatVec.push_back(inf);
    s_vector2fVec.push_back(Vector2f(inf, inf));
    s_vector3fVec.push_back(Vector3f(inf, inf, inf));
    s_vector4fVec.push_back(Vector4f(inf, inf, inf, inf));
    /*
    s_matrixfVec.push_back(Matrixf::Rows(Vector4f(inf, inf, inf, inf),
                                         Vector4f(inf, inf, inf, inf),
                                         Vector4f(inf, inf, inf, inf),
                                         Vector4f(inf, inf, inf, inf)));
    */

    s_matrixfVec.push_back(Matrixf::Rows(
                               Vector4f(0.79450,         0.00000,         0.00000,        -0.00000),
                               Vector4f(0.00000,         1.19175,         0.00000,         0.00000),
                               Vector4f(0.00000,         0.00000,        -1.00200,         0.80180),
                               Vector4f(0.00000,         0.00000,        -1.00000,         1.00000)));

    // Generate the quats from the vec4s
    for(unsigned int i = 0; i < s_vector4fVec.size(); ++i)
    {
        Vector4f v = s_vector4fVec[i];
        v = v.Unit();
        s_quatfVec.push_back(Quatf(v.x, v.y, v.z, v.w));
    }

    // TODO: nans, q-nans & denormalized...
}

//
// Random tests
//

class RandomIntTest : public TestCase
{
public:
    RandomIntTest() : TestCase("RandomInt") {}
    ~RandomIntTest() {}
    bool Test() const
    {
        const int RangeMin = 0;
        const int RangeMax = 9;

        int bins[10];
        for (int i = 0; i < 10; ++i)
            bins[i] = 0;

        for (int i = 0; i < 10000; ++i)
        {
            int r = RandomInt(RangeMin, RangeMax);
            bins[r]++;
            if (r < RangeMin || r > RangeMax)
                return false;
        }

        // check to see if numbers are more or less uniformly distrubuted
        for (int i = 0; i < 10; ++i)
        {
            if (bins[i] < 950)
                return false;
        }
        return true;
    }
};

class RandomFloatTest : public TestCase
{
public:
    RandomFloatTest() : TestCase("RandomFloat") {}
    ~RandomFloatTest() {}
    bool Test() const
    {
        const float RangeMin = 100;
        const float RangeMax = 300;
        for (int i = 0; i < 1024; ++i)
        {
            float r = RandomFloat(RangeMin, RangeMax);
            if (r < RangeMin || r > RangeMax)
                return false;
        }
        return true;
    }
};

//
// Float tests
//

template <class UnaryOp>
class FloatUnaryOpTest : public TestCase
{
public:
    FloatUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
    ~FloatUnaryOpTest() {}

    bool Test() const
    {
        UnaryOp op;
        for (unsigned int i = 0; i < s_floatVec.size(); ++i)
        {
            float a = s_floatVec[i];
            if (!op(a))
            {
                printf("a = %.5f", a);
                return false;
            }
        }
        return true;
    }
};

class FloatDegToRad
{
public:
    static const char* GetName() { return "DegToRad"; }
    bool operator() (float a)
    {
        float r = a * (PI / 180.0f);
        bool retval = FuzzyFloatTest(DegToRad(a), r);
        if (!retval)
            printf("r = %.5f\n", r);
        return retval;
    }
};

class FloatRadToDeg
{
public:
    static const char* GetName() { return "RadToDeg"; }
    bool operator() (float a)
    {
        float r = a * (180.0f / PI);
        return FuzzyFloatTest(RadToDeg(a), r);
    }
};

class FloatLimitPi
{
public:
    static const char* GetName() { return "LimitPi"; }
    bool operator() (float a)
    {
        float r = a;
        if ((a > PI) || (a < -PI))
            r = (a + PI) - 2*PI * floor((a + PI)/(2*PI)) - PI;
        bool result = FuzzyFloatTest(LimitPi(a), r);
        if (!result)
            printf("a = %.5f, r = %.5f\n", a, r);

        return result;
    }
};

class LimitPiTest : public TestCase
{
public:
    LimitPiTest() : TestCase("LimitPiTest") {}
    ~LimitPiTest() {}

    bool Test() const
    {
        float delta = 0.1f;
        float r = -PI + delta;
        float t = LimitPi(PI + delta);
        bool pass1 = FuzzyFloatTest(r, t);
        if (!pass1)
            printf("r = %.5f, t = %.5f\n", r, t);

        r = PI - delta;
        t = LimitPi(-PI - delta);
        bool pass2 = FuzzyFloatTest(r, t);
        if (!pass2)
            printf("r = %.5f, t = %.5f\n", r, t);

        return pass1 && pass2;
    }
};

class FloatMod2Pi
{
public:
    static const char* GetName() { return "Mod2Pi"; }
    bool operator() (float a)
    {
        // give flt_max a pass.
        if (a == FLT_MAX || a == -FLT_MAX)
            return true;

        float r;
        if (a < 0)
            r = (2.0f * PI) + fmod(a, 2.0f * PI);
        else
            r = fmod(a, 2.0f * PI);
        bool result = FuzzyFloatTest(Mod2Pi(a), r);
        if (!result)
            printf("a = %.5f, r = %.5f, Mod2Pi(a) = %.5f\n", a, r, Mod2Pi(a));
        return result;
    }
};

template <class BinaryOp>
class FloatBinaryOpTest : public TestCase
{
public:
    FloatBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
    ~FloatBinaryOpTest() {}

    bool Test() const
    {
        BinaryOp op;
        for (unsigned int i = 0; i < s_vector2fVec.size(); ++i)
        {
            for (unsigned int j = 0; j < s_vector2fVec.size(); ++j)
            {
                float a = s_floatVec[i];
                float b = s_floatVec[j];
                if (!op(a, b))
                {
                    printf("a = %.5f b = %.5f", a, b);
                    return false;
                }
            }
        }
        return true;
    }
};

class FloatFuzzyEqual
{
public:
    static const char* GetName() { return "FuzzyEqual"; }
    bool operator() (float a, float b)
    {
        const float epsilon = 0.001f;
        bool r = fabs(a - b) <= epsilon;
        return !(r ^ FuzzyEqual(a, b, epsilon));
    }
};

class FloatLerp
{
public:
    static const char* GetName() { return "Lerp"; }
    bool operator() (float a, float b)
    {
        for (float t = 0.0f; t <= 1.0f; t += 0.01f)
        {
            float r = a + ((b - a) * t);
            if (!FuzzyFloatTest(r, Lerp(a, b, t)))
                return false;
        }
        return true;
    }
};

template <class TernaryOp>
class FloatTernaryOpTest : public TestCase
{
public:
    FloatTernaryOpTest() : TestCase(TernaryOp::GetName()) {}
    ~FloatTernaryOpTest() {}

    bool Test() const
    {
        TernaryOp op;
        for (unsigned int i = 0; i < s_vector2fVec.size(); ++i)
        {
            for (unsigned int j = 0; j < s_vector2fVec.size(); ++j)
            {
                for (unsigned int k = 0; k < s_vector2fVec.size(); ++k)
                {
                    float a = s_floatVec[i];
                    float b = s_floatVec[j];
                    float c = s_floatVec[k];
                    if (!op(a, b, c))
                    {
                        printf("a = %.5f b = %.5f, c = %.5f", a, b, c);
                        return false;
                    }
                }
            }
        }
        return true;
    }
};

class FloatClamp
{
public:
    static const char* GetName() { return "Clamp"; }
    bool operator() (float a, float min, float max)
    {
        float r = a;

        // don't check degenerate values. I'm assuming the user does the right thing.
        if (min > max)
            return true;

        if (a < min)
            r = min;
        else if (a > max)
            r = max;

        return FloatTest(Clamp(a, min, max), r);
    }
};

//
// Vector2f tests
//

template <class UnaryOp>
class Vector2fUnaryOpTest : public TestCase
{
public:
    Vector2fUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
    ~Vector2fUnaryOpTest() {}

    bool Test() const
    {
        UnaryOp op;
        for (unsigned int i = 0; i < s_vector2fVec.size(); ++i)
        {
            Vector2f a = s_vector2fVec[i];
            if (!op(a))
            {
                printf("a = (%.5f, %.5f)", a.x, a.y);
                return false;
            }
        }
        return true;
    }
};

class Vector2fNegation
{
public:
    static const char* GetName() { return "Negation"; }
    bool operator() (const Vector2f& a)
    {
        float ax = a.x, ay = a.y;
        float rx = -ax, ry = -ay;
        Vector2f r = -a;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y);
    }
};

class Vector2fLength
{
public:
    static const char* GetName() { return "Length"; }
    bool operator() (const Vector2f& a)
    {
        float ax = a.x, ay = a.y;
        float r = sqrt((ax * ax) + (ay * ay));
        float len = a.Len();
        return FloatTest(r, len);
    }
};

class Vector2fSet
{
public:
    static const char* GetName() { return "Set"; }
    bool operator() (const Vector2f& a)
    {
        Vector2f v(0.0f, 0.0f);
        v.Set(a);
        return FloatTest(v.x, a.x) && FloatTest(v.y, a.y);
    }
};


class Vector2fUnitVec
{
public:
    static const char* GetName() { return "UnitVec"; }
    bool operator() (const Vector2f& a)
    {
        float ax = a.x, ay = a.y;
        float len = sqrt((ax * ax) + (ay * ay));
        float rx = ax / len, ry = ay / len;
        Vector2f r = a.Unit();
        return FloatTest(rx, r.x) && FloatTest(ry, r.y);
    }
};

class Vector2fScalarMultiplication
{
public:
    static const char* GetName() { return "Scalar Multiplication"; }
    bool operator() (const Vector2f& a)
    {
        float ax = a.x, ay = a.y;
        for (unsigned int i = 0; i < s_floatVec.size(); ++i)
        {
            // test (vector * scalar) and (scalar * vector)
            float s = s_floatVec[i];
            float r1x = ax * s, r1y = ay * s;
            float r2x = s * ax, r2y = s * ay;
            Vector2f r1 = a * s;
            Vector2f r2 = s * a;
            if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r1.x, r2.x) && FloatTest(r1.y, r2.y)))
            {
                printf("s = %.5f", s);
                return false;
            }
        }
        return true;
    }
};

class Vector2fScalarDivision
{
public:
    static const char* GetName() { return "Scalar Division"; }
    bool operator() (const Vector2f& a)
    {
        float ax = a.x, ay = a.y;
        for (unsigned int i = 0; i < s_floatVec.size(); ++i)
        {
            // test (vector / scalar) and (scalar / vector)
            float s = s_floatVec[i];
            float r1x = ax / s, r1y = ay / s;
            float r2x = s / ax, r2y = s / ay;
            Vector2f r1 = a / s;
            Vector2f r2 = s / a;

            if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y)))
            {
                printf("s = %.5f\n", s);
                return false;
            }
        }
        return true;
    }
};

template <class BinaryOp>
class Vector2fBinaryOpTest : public TestCase
{
public:
    Vector2fBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
    ~Vector2fBinaryOpTest() {}

    bool Test() const
    {
        BinaryOp op;
        for (unsigned int i = 0; i < s_vector2fVec.size(); ++i)
        {
            for (unsigned int j = 0; j < s_vector2fVec.size(); ++j)
            {
                Vector2f a = s_vector2fVec[i];
                Vector2f b = s_vector2fVec[j];
                if (!op(a, b))
                {
                    printf("a = (%.5f, %.5f), b = (%.5f, %.5f)", a.x, a.y, b.x, b.y);
                    return false;
                }
            }
        }
        return true;
    }
};

class Vector2fAddition
{
public:
    static const char* GetName() { return "Addition"; }
    bool operator()(const Vector2f& a, const Vector2f& b) const
    {
        float ax = a.x, ay = a.y;
        float bx = b.x, by = b.y;
        float rx = ax + bx, ry = ay + by;
        Vector2f r = a + b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y);
    }
};

class Vector2fSubtraction
{
public:
    static const char* GetName() { return "Subtraction"; }
    bool operator()(const Vector2f& a, const Vector2f& b) const
    {
        float ax = a.x, ay = a.y;
        float bx = b.x, by = b.y;
        float rx = ax - bx, ry = ay - by;
        Vector2f r = a - b;

        return FuzzyFloatTest(rx, r.x) && FuzzyFloatTest(ry, r.y);
    }
};

class Vector2fDotProduct
{
public:
    static const char* GetName() { return "Dot Product"; }
    bool operator()(const Vector2f& a, const Vector2f& b) const
    {
        float ax = a.x, ay = a.y;
        float bx = b.x, by = b.y;
        float r = (ax * bx) + (ay * by);
        float dot = Dot(a, b);
        return FloatTest(dot, r);
    }
};

class Vector2fCompMul
{
public:
    static const char* GetName() { return "Comp Mul"; }
    bool operator()(const Vector2f& a, const Vector2f& b) const
    {
        float ax = a.x, ay = a.y;
        float bx = b.x, by = b.y;
        float rx = (ax * bx), ry = (ay * by);
        Vector2f r = a * b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y);
    }
};

class Vector2fCompDiv
{
public:
    static const char* GetName() { return "Comp Div"; }
    bool operator()(const Vector2f& a, const Vector2f& b) const
    {
        float ax = a.x, ay = a.y;
        float bx = b.x, by = b.y;
        float rx = (ax / bx), ry = (ay / by);
        Vector2f r = a / b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y);
    }
};

class Vector2fSetZeroTest : public TestCase
{
public:
    Vector2fSetZeroTest() : TestCase("SetZero") {}
    ~Vector2fSetZeroTest() {}

    bool Test() const
    {
        Vector2f v(1,1);
        v.SetZero();
        return (v.x == 0.0f) && (v.y == 0.0f);
    }
};


//
// Vector3f tests
//

template <class UnaryOp>
class Vector3fUnaryOpTest : public TestCase
{
public:
    Vector3fUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
    ~Vector3fUnaryOpTest() {}

    bool Test() const
    {
        UnaryOp op;
        for (unsigned int i = 0; i < s_vector3fVec.size(); ++i)
        {
            Vector3f a = s_vector3fVec[i];
            if (!op(a))
            {
                printf("a = (%.5f, %.5f, %.5f)", a.x, a.y, a.z);
                return false;
            }
        }
        return true;
    }
};

class Vector3fNegation
{
public:
    static const char* GetName() { return "Negation"; }
    bool operator() (const Vector3f& a)
    {
        float ax = a.x, ay = a.y, az = a.z;
        float rx = -ax, ry = -ay, rz = -az;
        Vector3f r = -a;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
    }
};

class Vector3fLength
{
public:
    static const char* GetName() { return "Length"; }
    bool operator() (const Vector3f& a)
    {
        float ax = a.x, ay = a.y, az = a.z;
        float r = sqrt((ax * ax) + (ay * ay) + (az * az));
        float len = a.Len();
        return FloatTest(r, len);
    }
};

class Vector3fUnitVec
{
public:
    static const char* GetName() { return "UnitVec"; }
    bool operator() (const Vector3f& a)
    {
        float ax = a.x, ay = a.y, az = a.z;
        float len = sqrt((ax * ax) + (ay * ay) + (az * az));
        float rx = ax / len, ry = ay / len, rz = az / len;
        Vector3f r = a.Unit();
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz ,r.z);
    }
};

class Vector3fScalarMultiplication
{
public:
    static const char* GetName() { return "Scalar Multiplication"; }
    bool operator() (const Vector3f& a)
    {
        float ax = a.x, ay = a.y, az = a.z;
        for (unsigned int i = 0; i < s_floatVec.size(); ++i)
        {
            // test (vector * scalar) and (scalar * vector)
            float s = s_floatVec[i];
            float r1x = ax * s, r1y = ay * s, r1z = az * s;
            float r2x = s * ax, r2y = s * ay, r2z = s * az;
            Vector3f r1 = a * s;
            Vector3f r2 = s * a;
            if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) &&
                  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) &&
                  FloatTest(r1.x, r2.x) && FloatTest(r1.y, r2.y) && FloatTest(r1.z, r2.z)))
            {
                printf("s = %.5f", s);
                return false;
            }
        }
        return true;
    }
};

class Vector3fScalarDivision
{
public:
    static const char* GetName() { return "Scalar Division"; }
    bool operator() (const Vector3f& a)
    {
        float ax = a.x, ay = a.y, az = a.z;
        for (unsigned int i = 0; i < s_floatVec.size(); ++i)
        {
            // test (vector / scalar) and (scalar / vector)
            float s = s_floatVec[i];
            float r1x = ax / s, r1y = ay / s, r1z = az / s;
            float r2x = s / ax, r2y = s / ay, r2z = s / az;
            Vector3f r1 = a / s;
            Vector3f r2 = s / a;

            if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) &&
                  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z)))
            {
                printf("s = %.5f\n", s);
                return false;
            }
        }
        return true;
    }
};

template <class BinaryOp>
class Vector3fBinaryOpTest : public TestCase
{
public:
    Vector3fBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
    ~Vector3fBinaryOpTest() {}

    bool Test() const
    {
        BinaryOp op;
        for (unsigned int i = 0; i < s_vector3fVec.size(); ++i)
        {
            for (unsigned int j = 0; j < s_vector3fVec.size(); ++j)
            {
                Vector3f a = s_vector3fVec[i];
                Vector3f b = s_vector3fVec[j];
                if (!op(a, b))
                {
                    printf("a = (%.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f)", a.x, a.y, a.z, b.x, b.y, b.z);
                    return false;
                }
            }
        }
        return true;
    }
};

class Vector3fAddition
{
public:
    static const char* GetName() { return "Addition"; }
    bool operator()(const Vector3f& a, const Vector3f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z;
        float bx = b.x, by = b.y, bz = b.z;
        float rx = ax + bx, ry = ay + by, rz = az + bz;
        Vector3f r = a + b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
    }
};

class Vector3fSubtraction
{
public:
    static const char* GetName() { return "Subtraction"; }
    bool operator()(const Vector3f& a, const Vector3f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z;
        float bx = b.x, by = b.y, bz = b.z;
        float rx = ax - bx, ry = ay - by, rz = az - bz;
        Vector3f r = a - b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
    }
};

class Vector3fDotProduct
{
public:
    static const char* GetName() { return "Dot Product"; }
    bool operator()(const Vector3f& a, const Vector3f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z;
        float bx = b.x, by = b.y, bz = b.z;
        float r = (ax * bx) + (ay * by) + (az * bz);
        float dot = Dot(a, b);
        return FloatTest(dot, r);
    }
};

class Vector3fCompMul
{
public:
    static const char* GetName() { return "Comp Mul"; }
    bool operator()(const Vector3f& a, const Vector3f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z;
        float bx = b.x, by = b.y, bz = b.z;
        float rx = (ax * bx), ry = (ay * by), rz = (az * bz);
        Vector3f r = a * b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
    }
};

class Vector3fCompDiv
{
public:
    static const char* GetName() { return "Comp Div"; }
    bool operator()(const Vector3f& a, const Vector3f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z;
        float bx = b.x, by = b.y, bz = b.z;
        float rx = (ax / bx), ry = (ay / by), rz = (az / bz);
        Vector3f r = a / b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
    }
};

class Vector3fCrossProduct
{
public:
    static const char* GetName() { return "Cross Product"; }
    bool operator()(const Vector3f& a, const Vector3f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z;
        float bx = b.x, by = b.y, bz = b.z;

        float rx = (ay * bz) - (az * by);
        float ry = (az * bx) - (ax * bz);
        float rz = (ax * by) - (ay * bx);

        Vector3f r = Cross(a, b);
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
    }
};

//
// Vector4f tests
//

template <class UnaryOp>
class Vector4fUnaryOpTest : public TestCase
{
public:
    Vector4fUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
    ~Vector4fUnaryOpTest() {}

    bool Test() const
    {
        UnaryOp op;
        for (unsigned int i = 0; i < s_vector4fVec.size(); ++i)
        {
            Vector4f a = s_vector4fVec[i];
            if (!op(a))
            {
                printf("a = (%.5f, %.5f, %.5f, %.5f)", a.x, a.y, a.z, a.w);
                return false;
            }
        }
        return true;
    }
};

class Vector4fNegation
{
public:
    static const char* GetName() { return "Negation"; }
    bool operator() (const Vector4f& a)
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float rx = -ax, ry = -ay, rz = -az, rw = -aw;
        Vector4f r = -a;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
    }
};

class Vector4fLength
{
public:
    static const char* GetName() { return "Length"; }
    bool operator() (const Vector4f& a)
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float r = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
        float len = a.Len();
        return FloatTest(r, len);
    }
};

class Vector4fUnitVec
{
public:
    static const char* GetName() { return "UnitVec"; }
    bool operator() (const Vector4f& a)
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float len = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
        float rx = ax / len, ry = ay / len, rz = az / len, rw = aw / len;
        Vector4f r = a.Unit();
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
    }
};

class Vector4fScalarMultiplication
{
public:
    static const char* GetName() { return "Scalar Multiplication"; }
    bool operator() (const Vector4f& a)
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        for (unsigned int i = 0; i < s_floatVec.size(); ++i)
        {
            // test (vector * scalar) and (scalar * vector)
            float s = s_floatVec[i];
            float r1x = ax * s, r1y = ay * s, r1z = az * s, r1w = aw * s;
            float r2x = s * ax, r2y = s * ay, r2z = s * az, r2w = s * aw;
            Vector4f r1 = a * s;
            Vector4f r2 = s * a;
            if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && FloatTest(r1w, r1.w) &&
                  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) && FloatTest(r2w, r2.w) &&
                  FloatTest(r1.x, r2.x) && FloatTest(r1.y, r2.y) && FloatTest(r1.z, r2.z) && FloatTest(r1.w, r2.w)))
            {
                printf("s = %.5f", s);
                return false;
            }
        }
        return true;
    }
};

class Vector4fScalarDivision
{
public:
    static const char* GetName() { return "Scalar Division"; }
    bool operator() (const Vector4f& a)
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        for (unsigned int i = 0; i < s_floatVec.size(); ++i)
        {
            // test (vector / scalar) and (scalar / vector)
            float s = s_floatVec[i];
            float r1x = ax / s, r1y = ay / s, r1z = az / s, r1w = aw / s;
            float r2x = s / ax, r2y = s / ay, r2z = s / az, r2w = s / aw;
            Vector4f r1 = a / s;
            Vector4f r2 = s / a;

            if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && FloatTest(r1w, r1.w) &&
                  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) && FloatTest(r2w, r2.w)))
            {
                printf("s = %.5f\n", s);
                return false;
            }
        }
        return true;
    }
};

template <class BinaryOp>
class Vector4fBinaryOpTest : public TestCase
{
public:
    Vector4fBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
    ~Vector4fBinaryOpTest() {}

    bool Test() const
    {
        BinaryOp op;
        for (unsigned int i = 0; i < s_vector4fVec.size(); ++i)
        {
            for (unsigned int j = 0; j < s_vector4fVec.size(); ++j)
            {
                Vector4f a = s_vector4fVec[i];
                Vector4f b = s_vector4fVec[j];
                if (!op(a, b))
                {
                    printf("a = (%.5f, %.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f, %.5f)", a.x, a.y, a.z, a.w, b.x, b.y, b.z, b.w);
                    return false;
                }
            }
        }
        return true;
    }
};

class Vector4fAddition
{
public:
    static const char* GetName() { return "Addition"; }
    bool operator()(const Vector4f& a, const Vector4f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float bx = b.x, by = b.y, bz = b.z, bw = b.w;
        float rx = ax + bx, ry = ay + by, rz = az + bz, rw = aw + bw;
        Vector4f r = a + b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
    }
};

class Vector4fSubtraction
{
public:
    static const char* GetName() { return "Subtraction"; }
    bool operator()(const Vector4f& a, const Vector4f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float bx = b.x, by = b.y, bz = b.z, bw = b.w;
        float rx = ax - bx, ry = ay - by, rz = az - bz, rw = aw - bw;
        Vector4f r = a - b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
    }
};

class Vector4fDotProduct
{
public:
    static const char* GetName() { return "Dot Product"; }
    bool operator()(const Vector4f& a, const Vector4f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float bx = b.x, by = b.y, bz = b.z, bw = b.w;
        float r = (ax * bx) + (ay * by) + (az * bz) + (aw * bw);
        float dot = Dot(a, b);
        return FloatTest(dot, r);
    }
};

class Vector4fCompMul
{
public:
    static const char* GetName() { return "Comp Mul"; }
    bool operator()(const Vector4f& a, const Vector4f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float bx = b.x, by = b.y, bz = b.z, bw = b.w;
        float rx = (ax * bx), ry = (ay * by), rz = (az * bz), rw = (aw * bw);
        Vector4f r = a * b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
    }
};

class Vector4fCompDiv
{
public:
    static const char* GetName() { return "Comp Div"; }
    bool operator()(const Vector4f& a, const Vector4f& b) const
    {
        float ax = a.x, ay = a.y, az = a.z, aw = a.w;
        float bx = b.x, by = b.y, bz = b.z, bw = b.w;
        float rx = (ax / bx), ry = (ay / by), rz = (az / bz), rw = (aw / bw);
        Vector4f r = a / b;
        return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
    }
};

//
// Quatf tests
//

template <class UnaryOp>
class QuatfUnaryOpTest : public TestCase
{
public:
    QuatfUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
    ~QuatfUnaryOpTest() {}

    bool Test() const
    {
        UnaryOp op;
        for (unsigned int i = 0; i < s_quatfVec.size(); ++i)
        {
            Quatf a = s_quatfVec[i];
            if (!op(a))
            {
                printf("a = (%.5f, %.5f, %.5f, %.5f)", a.i, a.j, a.k, a.r);
                return false;
            }
        }
        return true;
    }
};

void QuatfMul(float* rx, float* ry, float* rz, float* rw, float ax, float ay, float az, float aw, float bx, float by, float bz, float bw)
{
    *rx =  ax * bw + ay * bz - az * by + aw * bx;
    *ry = -ax * bz + ay * bw + az * bx + aw * by;
    *rz =  ax * by - ay * bx + az * bw + aw * bz;
    *rw = -ax * bx - ay * by - az * bz + aw * bw;
}

class QuatfRotate
{
public:
    static const char* GetName() { return "Rotate"; }
    bool operator() (const Quatf& a)
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float rx, ry, rz, rw;

        for (unsigned int i = 0; i < s_vector3fVec.size(); ++i)
        {
            Vector3f v = s_vector3fVec[i];
            float bx = v.x, by = v.y, bz = v.z, bw = 0.0f;

            QuatfMul(&rx, &ry, &rz, &rw, ax, ay, az, aw, bx, by, bz, bw);
            QuatfMul(&rx, &ry, &rz, &rw, rx, ry, rz, rw, -ax, -ay, -az, aw);

            Vector3f r = a.Rotate(v);
            if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z)))
            {
                printf("v = (%.5f, %.5f, %.5f)", v.x, v.y, v.z);
                return false;
            }
        }
        return true;
    }
};

class QuatfConjugate
{
public:
    static const char* GetName() { return "Conjugate"; }
    bool operator() (const Quatf& a)
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float rx = -ax, ry = -ay, rz = -az, rw = aw;
        Quatf r = ~a;
        return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
    }
};

class QuatfNegation
{
public:
    static const char* GetName() { return "Negation"; }
    bool operator() (const Quatf& a)
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float rx = -ax, ry = -ay, rz = -az, rw = -aw;
        Quatf r = -a;
        return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
    }
};

class QuatfLength
{
public:
    static const char* GetName() { return "Length"; }
    bool operator() (const Quatf& a)
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float r = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
        float len = a.Len();
        return FloatTest(r, len);
    }
};

class QuatfUnitVec
{
public:
    static const char* GetName() { return "UnitVec"; }
    bool operator() (const Quatf& a)
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float len = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
        float rx = ax / len, ry = ay / len, rz = az / len, rw = aw / len;
        Quatf r = a.Unit();
        return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
    }
};

void QuatfExp(float* rx, float* ry, float* rz, float* rw, float qx, float qy, float qz, float qw)
{
    float angle = sqrt((qx * qx) + (qy * qy) + (qz * qz));
    float sin_a = sin(angle / 2.0f);
    // prevent div by zero
    if (angle > 0.0001f)
    {
        *rx = (qx / angle) * sin_a;
        *ry = (qy / angle) * sin_a;
        *rz = (qz / angle) * sin_a;
    }
    else
        *rx = *ry = *rz = 0.0f;

    *rw = cos(angle / 2.0f);
}

class QuatfExponential
{
public:
    static const char* GetName() { return "Quatf Exponential"; }
    bool operator() (const Quatf& a)
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float rx, ry, rz, rw;
        QuatfExp(&rx, &ry, &rz, &rw, ax, ay, az, aw);
        Quatf r = Exp(a);
        bool rval = FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
        if (!rval)
        {
            printf("rx = %.5f, ry = %.5f, rz = %.5f, rw = %.5f\n", rx, ry, rz, rw);
            printf("r = %.5f, %.5f, %.5f, %.5f\n", r.i, r.j, r.k, r.r);
        }
        return rval;
    }
};

void QuatfLog(float* rx, float* ry, float* rz, float* rw, float qx, float qy, float qz, float qw)
{
    float cos_a = qw;
    if (cos_a > 1.0f) cos_a = 1.0f;
    if (cos_a < -1.0f) cos_a = -1.0f;

    float sin_a = (float)sqrt(1.0f - cos_a * cos_a);

    if (fabs(sin_a) < 0.0005f)
        sin_a = 1.0f;
    else
        sin_a = 1.f/sin_a;

    float angle = 2.0f * (float)acos(cos_a);

    *rx = qx * sin_a * angle;
    *ry = qy * sin_a * angle;
    *rz = qz * sin_a * angle;
    *rw = 0.0f;
}

class QuatfLogarithm
{
public:
    static const char* GetName() { return "Quatf Logarithm"; }
    bool operator() (const Quatf& a)
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float rx, ry, rz, rw;
        QuatfLog(&rx, &ry, &rz, &rw, ax, ay, az, aw);
        Quatf r = Log(a);
        bool retval = FuzzyFloatTest(rx, r.i) && FuzzyFloatTest(ry, r.j) && FuzzyFloatTest(rz, r.k) && FuzzyFloatTest(rw, r.r);
        if (!retval)
            printf("rx = %.5f, ry = %.5f, rz = %.5f, rw = %.5f\n", rx, ry, rz, rw);
        return retval;
    }
};

template <class BinaryOp>
class QuatfBinaryOpTest : public TestCase
{
public:
    QuatfBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
    ~QuatfBinaryOpTest() {}

    bool Test() const
    {
        BinaryOp op;
        for (unsigned int i = 0; i < s_quatfVec.size(); ++i)
        {
            for (unsigned int j = 0; j < s_quatfVec.size(); ++j)
            {
                Quatf a = s_quatfVec[i];
                Quatf b = s_quatfVec[j];
                if (!op(a, b))
                {
                    printf("a = (%.5f, %.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f, %.5f)", a.i, a.j, a.k, a.r, b.i, b.j, b.k, b.r);
                    return false;
                }
            }
        }
        return true;
    }
};

class QuatfAddition
{
public:
    static const char* GetName() { return "Addition"; }
    bool operator()(const Quatf& a, const Quatf& b) const
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float bx = b.i, by = b.j, bz = b.k, bw = b.r;
        float rx = ax + bx, ry = ay + by, rz = az + bz, rw = aw + bw;
        Quatf r = a + b;
        return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
    }
};

class QuatfSubtraction
{
public:
    static const char* GetName() { return "Subtraction"; }
    bool operator()(const Quatf& a, const Quatf& b) const
    {
        float ax = a.i, ay = a.j, az = a.k, aw = a.r;
        float bx = b.i, by = b.j, bz = b.k, bw = b.r;
        float rx = ax - bx, ry = ay - by, rz = az - bz, rw = aw - bw;
        Quatf r = a - b;
        return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
    }
};

//
// Matrixf tests
//


template <class UnaryOp>
class MatrixfUnaryOpTest : public TestCase
{
public:
    MatrixfUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
    ~MatrixfUnaryOpTest() {}

    bool Test() const
    {
        UnaryOp op;
        for (unsigned int i = 0; i < s_matrixfVec.size(); ++i)
        {
            Matrixf m = s_matrixfVec[i];
            if (!op(m))
            {
                printf("m = \n");
                PrintMatrix(m);
                return false;
            }
        }
        return true;
    }
};

void MatrixfToFloatVec(float *fv, const Matrixf& m)
{
    for(int r = 0; r < 4; ++r)
    {
        for(int c = 0; c < 4; ++c)
        {
            fv[4 * r + c] = m.Elem(r,c);
        }
    }
}

class MatrixfTransform3x3
{
public:
    static const char* GetName() { return "Transform3x3"; }
    bool operator() (const Matrixf& a)
    {
        float fv[16];
        MatrixfToFloatVec(fv, a);

        for (unsigned int i = 0; i < s_vector3fVec.size(); ++i)
        {
            Vector3f v = s_vector3fVec[i];
            float bx = v.x, by = v.y, bz = v.z;

            float rx = fv[0] * bx + fv[1] * by + fv[2] * bz;
            float ry = fv[4] * bx + fv[5] * by + fv[6] * bz;
            float rz = fv[8] * bx + fv[9] * by + fv[10] * bz;

            Vector3f r = a.Mul3x3(v);
            if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z)))
            {
                printf("v = (%.5f, %.5f, %.5f)", v.x, v.y, v.z);
                return false;
            }
        }
        return true;
    }
};

class MatrixfTransform3x4
{
public:
    static const char* GetName() { return "Transform3x4"; }
    bool operator() (const Matrixf& a)
    {
        float fv[16];
        MatrixfToFloatVec(fv, a);

        for (unsigned int i = 0; i < s_vector3fVec.size(); ++i)
        {
            Vector3f v = s_vector3fVec[i];
            float bx = v.x, by = v.y, bz = v.z;

            float rx = fv[0] * bx + fv[1] * by + fv[2] * bz + fv[3];
            float ry = fv[4] * bx + fv[5] * by + fv[6] * bz + fv[7];
            float rz = fv[8] * bx + fv[9] * by + fv[10] * bz + fv[11];

            Vector3f r = a.Mul3x4(v);
            if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z)))
            {
                printf("v = (%.5f, %.5f, %.5f)", v.x, v.y, v.z);
                return false;
            }
        }
        return true;
    }
};

class MatrixfTransform4x4
{
public:
    static const char* GetName() { return "Transform4x4"; }
    bool operator() (const Matrixf& a)
    {
        float fv[16];
        MatrixfToFloatVec(fv, a);

        for (unsigned int i = 0; i < s_vector4fVec.size(); ++i)
        {
            Vector4f v = s_vector4fVec[i];
            float bx = v.x, by = v.y, bz = v.z, bw = v.w;

            float rx = fv[0] * bx + fv[1] * by + fv[2] * bz + fv[3] * bw;
            float ry = fv[4] * bx + fv[5] * by + fv[6] * bz + fv[7] * bw;
            float rz = fv[8] * bx + fv[9] * by + fv[10] * bz + fv[11] * bw;
            float rw = fv[12] * bx + fv[13] * by + fv[14] * bz + fv[15] * bw;

            Vector4f r = a.Mul4x4(v);
            if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w)))
            {
                printf("v = (%.5f, %.5f, %.5f, %.5f)", v.x, v.y, v.z, v.w);
                return false;
            }
        }
        return true;
    }
};

template <class BinaryOp>
class MatrixfBinaryOpTest : public TestCase
{
public:
    MatrixfBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
    ~MatrixfBinaryOpTest() {}

    bool Test() const
    {
        BinaryOp op;
        for (unsigned int i = 0; i < s_matrixfVec.size(); ++i)
        {
            for (unsigned int j = 0; j < s_matrixfVec.size(); ++j)
            {
                Matrixf a = s_matrixfVec[i];
                Matrixf b = s_matrixfVec[j];
                if (!op(a, b))
                {
                    printf("a = \n");
                    PrintMatrix(a);

                    printf("b = \n");
                    PrintMatrix(b);
                    return false;
                }
            }
        }
        return true;
    }
};



class MatrixfMultiplication
{
public:
    static const char* GetName() { return "Multiplication"; }
    bool operator()(const Matrixf& a, const Matrixf& b) const
    {
        float afv[16], bfv[16];
        MatrixfToFloatVec(afv, a);
        MatrixfToFloatVec(bfv, b);

        float rfv[16];
        static int ri[4][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}};
        static int ci[4][4] = {{0, 4, 8, 12}, {1, 5, 9, 13}, {2, 6, 10, 14}, {3, 7, 11, 15}};
        for (int r = 0; r < 4; ++r)
        {
            for (int c = 0; c < 4; ++c)
            {
                float terms[4];
                for (int i = 0; i < 4; ++i)
                    terms[i] = afv[ri[r][i]] * bfv[ci[c][i]];
                rfv[r * 4 + c] = terms[0] + terms[1] + terms[2] + terms[3];
            }
        }

        Matrixf r = a * b;

        return FloatTest(rfv[0], r.Elem(0,0)) && FloatTest(rfv[1], r.Elem(0,1)) && FloatTest(rfv[2], r.Elem(0,2)) && FloatTest(rfv[3], r.Elem(0,3)) &&
            FloatTest(rfv[4], r.Elem(1,0)) && FloatTest(rfv[5], r.Elem(1,1)) && FloatTest(rfv[6], r.Elem(1,2)) && FloatTest(rfv[7], r.Elem(1,3)) &&
            FloatTest(rfv[8], r.Elem(2,0)) && FloatTest(rfv[9], r.Elem(2,1)) && FloatTest(rfv[10], r.Elem(2,2)) && FloatTest(rfv[11], r.Elem(2,3)) &&
            FloatTest(rfv[12], r.Elem(3,0)) && FloatTest(rfv[13], r.Elem(3,1)) && FloatTest(rfv[14], r.Elem(3,2)) && FloatTest(rfv[15], r.Elem(3,3));
    }
};


/*
 * Compute inverse of 4x4 transformation matrix.
 * Code contributed by Jacques Leroy jle@star.be
 * Return GL_TRUE for success, GL_FALSE for failure (singular matrix)
 */
static bool
invert_matrix(const float * m, float * out)
{
/* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { float *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(m,r,c) (m)[(c)*4+(r)]

   float wtmp[4][8];
   float m0, m1, m2, m3, s;
   float *r0, *r1, *r2, *r3;

   r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

   r0[0] = MAT(m, 0, 0), r0[1] = MAT(m, 0, 1),
      r0[2] = MAT(m, 0, 2), r0[3] = MAT(m, 0, 3),
      r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0,
      r1[0] = MAT(m, 1, 0), r1[1] = MAT(m, 1, 1),
      r1[2] = MAT(m, 1, 2), r1[3] = MAT(m, 1, 3),
      r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0,
      r2[0] = MAT(m, 2, 0), r2[1] = MAT(m, 2, 1),
      r2[2] = MAT(m, 2, 2), r2[3] = MAT(m, 2, 3),
      r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0,
      r3[0] = MAT(m, 3, 0), r3[1] = MAT(m, 3, 1),
      r3[2] = MAT(m, 3, 2), r3[3] = MAT(m, 3, 3),
      r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;

   /* choose pivot - or die */
   if (fabs(r3[0]) > fabs(r2[0]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[0]) > fabs(r1[0]))
      SWAP_ROWS(r2, r1);
   if (fabs(r1[0]) > fabs(r0[0]))
      SWAP_ROWS(r1, r0);
   if (0.0 == r0[0])
      return false;

   /* eliminate first variable     */
   m1 = r1[0] / r0[0];
   m2 = r2[0] / r0[0];
   m3 = r3[0] / r0[0];
   s = r0[1];
   r1[1] -= m1 * s;
   r2[1] -= m2 * s;
   r3[1] -= m3 * s;
   s = r0[2];
   r1[2] -= m1 * s;
   r2[2] -= m2 * s;
   r3[2] -= m3 * s;
   s = r0[3];
   r1[3] -= m1 * s;
   r2[3] -= m2 * s;
   r3[3] -= m3 * s;
   s = r0[4];
   if (s != 0.0) {
      r1[4] -= m1 * s;
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r0[5];
   if (s != 0.0) {
      r1[5] -= m1 * s;
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r0[6];
   if (s != 0.0) {
      r1[6] -= m1 * s;
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r0[7];
   if (s != 0.0) {
      r1[7] -= m1 * s;
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[1]) > fabs(r2[1]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[1]) > fabs(r1[1]))
      SWAP_ROWS(r2, r1);
   if (0.0 == r1[1])
       return false;

   /* eliminate second variable */
   m2 = r2[1] / r1[1];
   m3 = r3[1] / r1[1];
   r2[2] -= m2 * r1[2];
   r3[2] -= m3 * r1[2];
   r2[3] -= m2 * r1[3];
   r3[3] -= m3 * r1[3];
   s = r1[4];
   if (0.0 != s) {
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r1[5];
   if (0.0 != s) {
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r1[6];
   if (0.0 != s) {
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r1[7];
   if (0.0 != s) {
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[2]) > fabs(r2[2]))
      SWAP_ROWS(r3, r2);
   if (0.0 == r2[2])
      return false;

   /* eliminate third variable */
   m3 = r3[2] / r2[2];
   r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
      r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6], r3[7] -= m3 * r2[7];

   /* last check */
   if (0.0 == r3[3])
      return false;

   s = 1.0 / r3[3];     /* now back substitute row 3 */
   r3[4] *= s;
   r3[5] *= s;
   r3[6] *= s;
   r3[7] *= s;

   m2 = r2[3];          /* now back substitute row 2 */
   s = 1.0 / r2[2];
   r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
      r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
   m1 = r1[3];
   r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
      r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
   m0 = r0[3];
   r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
      r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

   m1 = r1[2];          /* now back substitute row 1 */
   s = 1.0 / r1[1];
   r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
      r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
   m0 = r0[2];
   r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
      r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

   m0 = r0[1];          /* now back substitute row 0 */
   s = 1.0 / r0[0];
   r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
      r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

   MAT(out, 0, 0) = r0[4];
   MAT(out, 0, 1) = r0[5], MAT(out, 0, 2) = r0[6];
   MAT(out, 0, 3) = r0[7], MAT(out, 1, 0) = r1[4];
   MAT(out, 1, 1) = r1[5], MAT(out, 1, 2) = r1[6];
   MAT(out, 1, 3) = r1[7], MAT(out, 2, 0) = r2[4];
   MAT(out, 2, 1) = r2[5], MAT(out, 2, 2) = r2[6];
   MAT(out, 2, 3) = r2[7], MAT(out, 3, 0) = r3[4];
   MAT(out, 3, 1) = r3[5], MAT(out, 3, 2) = r3[6];
   MAT(out, 3, 3) = r3[7];

   return true;

#undef MAT
#undef SWAP_ROWS
}

/*
 * Compute the inverse of a 4x4 matrix.
 *
 * From an algorithm by V. Strassen, 1969, _Numerishe Mathematik_, vol. 13,
 * pp. 354-356.
 * 60 multiplies, 24 additions, 10 subtractions, 8 negations, 2 divisions,
 * 48 assignments, _0_ branches
 *
 * This implementation by Scott McCaskill
 */

typedef float Mat2[2][2];

enum {
    M00 = 0, M01 = 4, M02 = 8, M03 = 12,
    M10 = 1, M11 = 5, M12 = 9, M13 = 13,
    M20 = 2, M21 = 6, M22 = 10,M23 = 14,
    M30 = 3, M31 = 7, M32 = 11,M33 = 15
};

static void invert_matrix_general( const float *m, float *out )
{
   Mat2 r1, r2, r3, r4, r5, r6, r7;
   const float * A = m;
   float *       C = out;
   float one_over_det;

   /*
    * A is the 4x4 source matrix (to be inverted).
    * C is the 4x4 destination matrix
    * a11 is the 2x2 matrix in the upper left quadrant of A
    * a12 is the 2x2 matrix in the upper right quadrant of A
    * a21 is the 2x2 matrix in the lower left quadrant of A
    * a22 is the 2x2 matrix in the lower right quadrant of A
    * similarly, cXX are the 2x2 quadrants of the destination matrix
    */

   /* R1 = inverse( a11 ) */
   one_over_det = 1.0f / ( ( A[M00] * A[M11] ) - ( A[M10] * A[M01] ) );
   r1[0][0] = one_over_det * A[M11];
   r1[0][1] = one_over_det * -A[M01];
   r1[1][0] = one_over_det * -A[M10];
   r1[1][1] = one_over_det * A[M00];

   /* R2 = a21 x R1 */
   r2[0][0] = A[M20] * r1[0][0] + A[M21] * r1[1][0];
   r2[0][1] = A[M20] * r1[0][1] + A[M21] * r1[1][1];
   r2[1][0] = A[M30] * r1[0][0] + A[M31] * r1[1][0];
   r2[1][1] = A[M30] * r1[0][1] + A[M31] * r1[1][1];

   /* R3 = R1 x a12 */
   r3[0][0] = r1[0][0] * A[M02] + r1[0][1] * A[M12];
   r3[0][1] = r1[0][0] * A[M03] + r1[0][1] * A[M13];
   r3[1][0] = r1[1][0] * A[M02] + r1[1][1] * A[M12];
   r3[1][1] = r1[1][0] * A[M03] + r1[1][1] * A[M13];

   /* R4 = a21 x R3 */
   r4[0][0] = A[M20] * r3[0][0] + A[M21] * r3[1][0];
   r4[0][1] = A[M20] * r3[0][1] + A[M21] * r3[1][1];
   r4[1][0] = A[M30] * r3[0][0] + A[M31] * r3[1][0];
   r4[1][1] = A[M30] * r3[0][1] + A[M31] * r3[1][1];

   /* R5 = R4 - a22 */
   r5[0][0] = r4[0][0] - A[M22];
   r5[0][1] = r4[0][1] - A[M23];
   r5[1][0] = r4[1][0] - A[M32];
   r5[1][1] = r4[1][1] - A[M33];

   /* R6 = inverse( R5 ) */
   one_over_det = 1.0f / ( ( r5[0][0] * r5[1][1] ) - ( r5[1][0] * r5[0][1] ) );
   r6[0][0] = one_over_det * r5[1][1];
   r6[0][1] = one_over_det * -r5[0][1];
   r6[1][0] = one_over_det * -r5[1][0];
   r6[1][1] = one_over_det * r5[0][0];

   /* c12 = R3 x R6 */
   C[M02] = r3[0][0] * r6[0][0] + r3[0][1] * r6[1][0];
   C[M03] = r3[0][0] * r6[0][1] + r3[0][1] * r6[1][1];
   C[M12] = r3[1][0] * r6[0][0] + r3[1][1] * r6[1][0];
   C[M13] = r3[1][0] * r6[0][1] + r3[1][1] * r6[1][1];

   /* c21 = R6 x R2 */
   C[M20] = r6[0][0] * r2[0][0] + r6[0][1] * r2[1][0];
   C[M21] = r6[0][0] * r2[0][1] + r6[0][1] * r2[1][1];
   C[M30] = r6[1][0] * r2[0][0] + r6[1][1] * r2[1][0];
   C[M31] = r6[1][0] * r2[0][1] + r6[1][1] * r2[1][1];

   /* R7 = R3 x c21 */
   r7[0][0] = r3[0][0] * C[M20] + r3[0][1] * C[M30];
   r7[0][1] = r3[0][0] * C[M21] + r3[0][1] * C[M31];
   r7[1][0] = r3[1][0] * C[M20] + r3[1][1] * C[M30];
   r7[1][1] = r3[1][0] * C[M21] + r3[1][1] * C[M31];

   /* c11 = R1 - R7 */
   C[M00] = r1[0][0] - r7[0][0];
   C[M01] = r1[0][1] - r7[0][1];
   C[M10] = r1[1][0] - r7[1][0];
   C[M11] = r1[1][1] - r7[1][1];

   /* c22 = -R6 */
   C[M22] = -r6[0][0];
   C[M23] = -r6[0][1];
   C[M32] = -r6[1][0];
   C[M33] = -r6[1][1];
}




class MatrixfInverse
{
public:
    static const char* GetName() { return "Inverse"; }
    bool operator()(const Matrixf& a) const
    {
//      Matrixf m(Quatf(Vector3f(0.3f,1.5f,0.7f), -PI/2.0f), Vector3f(1.0, 2.0f, 3.0f));
//      m.SetScale(Vector3f(1,2,1));

        float fv[16];

        // convert the transpose of a into a float vec. (why? because invert_matrix expects an opengl style matrix)
        MatrixfToFloatVec(fv, a.Transpose());
        float rfv[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};  // initalize it to stop compiler warnings.
        bool invertable = invert_matrix(fv, rfv);

        float rfv2[16];
        invert_matrix_general(fv, rfv2);

        Matrixf r;
        bool matrix_invertable = FullInverse(a, r);

        if (invertable != matrix_invertable)
        {
            printf("Failure because invertable mismatch\n");
            return false;  // test failed
        }

        if (invertable)
        {
            bool rval;
            // compare r with the transpose of rfv.  (why? becaaue invert_matrix returns an opengl style matrix)
            const float epsilon = 0.05f;
            rval = FuzzyFloatTest(rfv[0], r.Elem(0,0), epsilon) &&
                FuzzyFloatTest(rfv[4], r.Elem(0,1), epsilon) &&
                FuzzyFloatTest(rfv[8], r.Elem(0,2), epsilon) &&
                FuzzyFloatTest(rfv[12], r.Elem(0,3), epsilon) &&
                FuzzyFloatTest(rfv[1], r.Elem(1,0), epsilon) &&
                FuzzyFloatTest(rfv[5], r.Elem(1,1), epsilon) &&
                FuzzyFloatTest(rfv[9], r.Elem(1,2), epsilon) &&
                FuzzyFloatTest(rfv[13], r.Elem(1,3), epsilon) &&
                FuzzyFloatTest(rfv[2], r.Elem(2,0), epsilon) &&
                FuzzyFloatTest(rfv[6], r.Elem(2,1), epsilon) &&
                FuzzyFloatTest(rfv[10], r.Elem(2,2), epsilon) &&
                FuzzyFloatTest(rfv[14], r.Elem(2,3), epsilon) &&
                FuzzyFloatTest(rfv[3], r.Elem(3,0), epsilon) &&
                FuzzyFloatTest(rfv[7], r.Elem(3,1), epsilon) &&
                FuzzyFloatTest(rfv[11], r.Elem(3,2), epsilon) &&
                FuzzyFloatTest(rfv[15], r.Elem(3,3), epsilon);

            if (!rval)
            {
                printf("Failure because matrices are not equal\n");

                printf("rfv = \n");
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[0], rfv[4], rfv[8], rfv[12]);
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[1], rfv[5], rfv[9], rfv[13]);
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[2], rfv[6], rfv[10], rfv[14]);
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[3], rfv[7], rfv[11], rfv[15]);

                printf("rfv2 = \n");
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[0], rfv2[4], rfv2[8], rfv2[12]);
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[1], rfv2[5], rfv2[9], rfv2[13]);
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[2], rfv2[6], rfv2[10], rfv2[14]);
                printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[3], rfv2[7], rfv2[11], rfv2[15]);


                printf("r = \n");
                PrintMatrix(r);
            }

            return rval;
        }
        else
        {
            return true;  // both matrices were determined to be non-invertable.
        }
    }
};

#ifdef CMD_LINE
int main(int argc, char* argv[])
#else
int AbaciUnitTest()
#endif
{
    InitTestData();

    TestSuite randomSuite("random");
    randomSuite.AddTest(new RandomIntTest());
    randomSuite.AddTest(new RandomFloatTest());
    randomSuite.RunTests();

    TestSuite floatSuite("float");
    floatSuite.AddTest(new FloatUnaryOpTest<FloatDegToRad>());
    floatSuite.AddTest(new FloatUnaryOpTest<FloatRadToDeg>());
    floatSuite.AddTest(new FloatTernaryOpTest<FloatClamp>());
    floatSuite.AddTest(new FloatUnaryOpTest<FloatLimitPi>());
    floatSuite.AddTest(new LimitPiTest());
    floatSuite.AddTest(new FloatUnaryOpTest<FloatMod2Pi>());
    floatSuite.AddTest(new FloatBinaryOpTest<FloatFuzzyEqual>());
    floatSuite.AddTest(new FloatBinaryOpTest<FloatLerp>());
    floatSuite.RunTests();

    TestSuite vector2Suite("Vector2f");
    vector2Suite.AddTest(new Vector2fUnaryOpTest<Vector2fNegation>());
    vector2Suite.AddTest(new Vector2fUnaryOpTest<Vector2fLength>());
    vector2Suite.AddTest(new Vector2fUnaryOpTest<Vector2fUnitVec>());
    vector2Suite.AddTest(new Vector2fUnaryOpTest<Vector2fScalarMultiplication>());
    vector2Suite.AddTest(new Vector2fUnaryOpTest<Vector2fScalarDivision>());
    vector2Suite.AddTest(new Vector2fBinaryOpTest<Vector2fAddition>());
    vector2Suite.AddTest(new Vector2fBinaryOpTest<Vector2fSubtraction>());
    vector2Suite.AddTest(new Vector2fBinaryOpTest<Vector2fDotProduct>());
    vector2Suite.AddTest(new Vector2fBinaryOpTest<Vector2fCompMul>());
    vector2Suite.AddTest(new Vector2fBinaryOpTest<Vector2fCompDiv>());
    vector2Suite.AddTest(new Vector2fSetZeroTest());
    vector2Suite.AddTest(new Vector2fUnaryOpTest<Vector2fSet>());
    vector2Suite.RunTests();
    // TODO: Construct from complex.
    // TODO: Construct from 2 scalars.
    // TODO: Set()
    // TODO: LenSq()
    // TODO: Lerp()
    // TODO: []
    // TODO: MinLen()
    // TODO: -=
    // TODO: +=
    // TODO: vector * scalar
    // TODO: scalar * vector
    // TODO: vector / scalar
    // TODO: scalar / vector

    TestSuite vector3Suite("Vector3f");
    vector3Suite.AddTest(new Vector3fUnaryOpTest<Vector3fNegation>());
    vector3Suite.AddTest(new Vector3fUnaryOpTest<Vector3fLength>());
    vector3Suite.AddTest(new Vector3fUnaryOpTest<Vector3fUnitVec>());
    vector3Suite.AddTest(new Vector3fUnaryOpTest<Vector3fScalarMultiplication>());
    vector3Suite.AddTest(new Vector3fUnaryOpTest<Vector3fScalarDivision>());
    vector3Suite.AddTest(new Vector3fBinaryOpTest<Vector3fAddition>());
    vector3Suite.AddTest(new Vector3fBinaryOpTest<Vector3fSubtraction>());
    vector3Suite.AddTest(new Vector3fBinaryOpTest<Vector3fDotProduct>());
    vector3Suite.AddTest(new Vector3fBinaryOpTest<Vector3fCompMul>());
    vector3Suite.AddTest(new Vector3fBinaryOpTest<Vector3fCompDiv>());
    vector3Suite.AddTest(new Vector3fBinaryOpTest<Vector3fCrossProduct>());
    vector3Suite.RunTests();
    // TODO: construct from 3 scalars.
    // TODO: SetZero()
    // TODO: Set()
    // TODO: LenSq()
    // TODO: Lerp()
    // TODO: []
    // TODO: MinLen()
    // TODO: -=
    // TODO: +=
    // TODO: vector * scalar
    // TODO: scalar * vector
    // TODO: vector / scalar
    // TODO: scalar / vector

    TestSuite vector4Suite("Vector4f");
    vector4Suite.AddTest(new Vector4fUnaryOpTest<Vector4fNegation>());
    vector4Suite.AddTest(new Vector4fUnaryOpTest<Vector4fLength>());
    vector4Suite.AddTest(new Vector4fUnaryOpTest<Vector4fUnitVec>());
    vector4Suite.AddTest(new Vector4fUnaryOpTest<Vector4fScalarMultiplication>());
    vector4Suite.AddTest(new Vector4fUnaryOpTest<Vector4fScalarDivision>());
    vector4Suite.AddTest(new Vector4fBinaryOpTest<Vector4fAddition>());
    vector4Suite.AddTest(new Vector4fBinaryOpTest<Vector4fSubtraction>());
    vector4Suite.AddTest(new Vector4fBinaryOpTest<Vector4fDotProduct>());
    vector4Suite.AddTest(new Vector4fBinaryOpTest<Vector4fCompMul>());
    vector4Suite.AddTest(new Vector4fBinaryOpTest<Vector4fCompDiv>());
    vector4Suite.RunTests();
    // TODO: construct from 4 scalars.
    // TODO: SetZero()
    // TODO: Set()
    // TODO: LenSq()
    // TODO: Lerp()
    // TODO: []
    // TODO: MinLen()
    // TODO: -=
    // TODO: +=
    // TODO: vector * scalar
    // TODO: scalar * vector
    // TODO: vector / scalar
    // TODO: scalar / vector

    TestSuite quatSuite("Quatf");
    quatSuite.AddTest(new QuatfUnaryOpTest<QuatfRotate>());
    quatSuite.AddTest(new QuatfUnaryOpTest<QuatfConjugate>());
    quatSuite.AddTest(new QuatfUnaryOpTest<QuatfNegation>());
    quatSuite.AddTest(new QuatfUnaryOpTest<QuatfLength>());
    quatSuite.AddTest(new QuatfUnaryOpTest<QuatfUnitVec>());
    quatSuite.AddTest(new QuatfUnaryOpTest<QuatfExponential>());
    quatSuite.AddTest(new QuatfUnaryOpTest<QuatfLogarithm>());
    quatSuite.AddTest(new QuatfBinaryOpTest<QuatfAddition>());
    quatSuite.AddTest(new QuatfBinaryOpTest<QuatfSubtraction>());
    quatSuite.RunTests();
    // TODO: construct from 4 scalars.
    // TODO: Set()
    // TODO: SetZero()
    // TODO: AxisAngle()
    // TODO: SetZero()
    // TODO: LenSq()
    // TODO: Dot()
    // TODO: quat

    TestSuite matrixSuite("Matrixf");
    matrixSuite.AddTest(new MatrixfUnaryOpTest<MatrixfTransform3x3>());
    matrixSuite.AddTest(new MatrixfUnaryOpTest<MatrixfTransform3x4>());
    matrixSuite.AddTest(new MatrixfUnaryOpTest<MatrixfTransform4x4>());
    matrixSuite.AddTest(new MatrixfBinaryOpTest<MatrixfMultiplication>());
    matrixSuite.AddTest(new MatrixfUnaryOpTest<MatrixfInverse>());

    // TODO: matrix from quat
    // TODO: matrix from quat & trans
    // TODO: matrix from quat, trans & scale
    // TODO: matrix from axis angle.
    // TODO: matrix from axes
    // TODO: matrix from rows
    // TODO: matrix set scale.
    // TODO: Frustum proj
    // TODO: Ortho proj
    // TODO: LookAt
    // TODO: Identity
    // TODO: matrix addition
    // TODO: matrix subtraction
    // TODO: OrthonormalInverse
    // TODO: Axis getters
    // TODO: Axis setters
    // TODO: SetScale uniform
    // TODO: SetScale non-uniform
    // TODO: Elem
    // TODO: GetCol
    // TODO: GetQuat
    // TODO: Transpose
    matrixSuite.RunTests();

    // TODO: Complex Tests
    // TODO: construct from 2 scalars
    // TODO: construct from vec2
    // TODO: ~
    // TODO: -
    // TODO: +
    // TODO: *
    // TODO: /
    // TODO: * scalar
    // TODO: ExpI
    // TODO: sqrt
    // TODO: exp
    // TODO: log

    return 0;
}
