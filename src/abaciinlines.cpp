#ifdef ABACI_H

#include <stdio.h>
#include <stdlib.h>

// Convert from degrees to radians
inline float DegToRad(float deg)
{
    return deg * (PI / 180.0);
}

// Convert from radians to degrees.
inline float RadToDeg(float rad)
{
    return rad * (180.0 / PI);
}

// Clamps value between min & max.
inline float Clamp(float value, float min, float max)
{
    float result = value;
    if (value > max)
        result = max;
    else if (value < min)
        result = min;
    return result;
}

// Limits angle between -PI & PI using modulus arithmetic
inline float LimitPi(float theta)
{
    return Mod2Pi(theta + PI) - PI;
}

// Limits angle between zero & 2 PI using modulus arithmetic
inline float Mod2Pi(float theta)
{
    return theta - 2*PI * floor(theta/(2*PI));
}

// Fuzzy comparison between two float values.
inline bool FuzzyEqual(float rhs, float lhs, float epsilon)
{
    return fabs(rhs - lhs) <= epsilon;
}

// Linear interpolate between two float values.
inline float Lerp(float a, float b, float t)
{
    return a + (b - a) * t;
}

inline int RandomInt(int min, int max)
{
    return min + (rand() % (max - min + 1));
}

float RandomFloat(float min, float max)
{
    float t = (float)((rand() % 16127) / 16126.0);
    return Lerp(min, max, t);
}

////////////////////////////////////////////////////////////////////////////////
// Vector2

// returns a random vector on the unit circle.
Vector2 Vector2::RandomUnitVector()
{
    float theta = RandomFloat(0.0f, (float)(2 * PI));
    return Vector2(cos(theta), sin(theta));
}

// Construct from two floats
inline Vector2::Vector2(float xIn, float yIn) : x(xIn), y(yIn) {}

// Construct from Complex
inline Vector2::Vector2(const Complex& complexIn) :
    x(complexIn.r), y(complexIn.i) {}

// Set from two floats
inline void Vector2::Set(float xIn, float yIn)
{
    x = xIn;
    y = yIn;
}

// Set from Complex
void Vector2::Set(const Complex& complexIn)
{
    x = complexIn.r;
    y = complexIn.i;
}

// Sets all elements to zero.
inline void Vector2::SetZero()
{
    x = y = 0;
}

// const array accessor
float Vector2::operator[](int i) const
{
    return *(&x + i);
}

// array accessor
float& Vector2::operator[](int i)
{
    return *(&x + i);
}

// Returns a vector with same direction but unit length.
inline Vector2 Vector2::Unit() const
{
    return *this / Len();
}

// Returns vector length.
inline float Vector2::Len() const
{
    return sqrt(Dot(*this, *this));
}

// Returns length squared.
inline float Vector2::LenSq() const
{
    return Dot(*this, *this);
}

// Returns a vector with the same direction, but with a Len() <= len.
inline Vector2 Vector2::MinLen(float len) const
{
    float l = Len();
    if (l > len)
        return (*this / l) * len;
    else
        return *this;
}

// Dot product of two vectors.
inline float Dot(const Vector2& a, const Vector2& b)
{
    return a.x * b.x + a.y * b.y;
}

// Linear interpolation between two vectors
inline Vector2 Lerp(const Vector2& a, const Vector2& b, float t)
{
    return a + (b - a) * t;
}

// Unary minus
inline Vector2 operator-(const Vector2& a)
{
    return Vector2(-a.x, -a.y);
}

// Vector subtraction.
inline Vector2 operator-(const Vector2& a, const Vector2& b)
{
    return Vector2(a.x - b.x, a.y - b.y);
}

// Vector addition.
inline Vector2 operator+(const Vector2& a, const Vector2& b)
{
    return Vector2(a.x + b.x, a.y + b.y);
}

// Vector decrement.
Vector2& operator-=(Vector2& a, const Vector2& b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

// Vector increment.
Vector2& operator+=(Vector2& a, const Vector2& b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

// Multplies all elements of a vector by a scalar.
inline Vector2 operator*(const Vector2& v, float scalar)
{
    return Vector2(scalar * v.x, scalar * v.y);
}

// Multplies all elements of a vector by a scalar.
inline Vector2 operator*(float scalar, const Vector2& v)
{
    return Vector2(scalar * v.x, scalar * v.y);
}

// Vector multiplication.
inline Vector2 operator*(const Vector2& a, const Vector2& b)
{
    return Vector2(a.x * b.x, a.y * b.y);
}

// Divides all elements of a vector by a scalar.
inline Vector2 operator/(const Vector2& v, float denominator)
{
    return Vector2(v.x / denominator, v.y / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
inline Vector2 operator/(float numerator, const Vector2& v)
{
    return Vector2(numerator / v.x, numerator / v.y);
}

// Vector division.
inline Vector2 operator/(const Vector2& a, const Vector2& b)
{
    return Vector2(a.x / b.x, a.y / b.y);
}

// Fuzzy comparison between two Vector2 values.
bool FuzzyEqual(const Vector2& rhs, const Vector2& lhs, float epsilon)
{
    return (FuzzyEqual(rhs.x, lhs.x, epsilon) &&
            FuzzyEqual(rhs.y, lhs.y, epsilon));
}

////////////////////////////////////////////////////////////////////////////////
// Vector3

// returns a random vector on the unit sphere.
Vector3 Vector3::RandomUnitVector()
{
    float theta = RandomFloat(0.0f, (float)(2 * PI));
    float phi = RandomFloat(0.0f, (float)(2 * PI));
    return Vector3(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi));
}

// Construct from three floats
inline Vector3::Vector3(float xIn, float yIn, float zIn) : x(xIn), y(yIn), z(zIn) {}

// Set from three floats
inline void Vector3::Set(float xIn, float yIn, float zIn)
{
    x = xIn;
    y = yIn;
    z = zIn;
}

// Set all elements to zero
inline void Vector3::SetZero()
{
    x = y = z = 0;
}

// const array accessor
float Vector3::operator[](int i) const
{
    return *(&x + i);
}

// array accessor
float& Vector3::operator[](int i)
{
    return *(&x + i);
}

// Returns a vector with same direction but unit length.
inline Vector3 Vector3::Unit() const
{
    return (*this) / Len();
}

// Returns vector length.
inline float Vector3::Len() const
{
    return sqrt(Dot(*this, *this));
}

// Returns vector length squared.
inline float Vector3::LenSq() const
{
    return Dot(*this, *this);
}

// Returns a vector with the same direction, but with a Len() <= len.
inline Vector3 Vector3::MinLen(float len) const
{
    float l = Len();
    if (l > len)
        return (*this / l) * len;
    else
        return *this;
}

// Generate vectors j and k which are orthognal to i and each other.
inline void Vector3::Basis(Vector3& iOut, Vector3& jOut, Vector3& kOut)
{
    iOut = Unit();
    float xabs = fabs(iOut.x);
    float yabs = fabs(iOut.y);
    float zabs = fabs(iOut.z);

    Vector3 j;
    if (xabs <= yabs && xabs <= zabs)
        j.Set(0, -iOut.z, iOut.y);
    else if (yabs <= xabs && yabs <= zabs)
        j.Set(-iOut.z, 0, iOut.x);
    else
        j.Set(-iOut.y, iOut.x, 0);

    jOut = j.Unit();
    kOut = Cross(iOut, jOut);
}

// Dot product of two vectors.
inline float Dot(const Vector3& a, const Vector3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Cross product of two vectors.
inline Vector3 Cross(const Vector3& a, const Vector3& b)
{
    return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

// Linear interpolation between two vectors
inline Vector3 Lerp(const Vector3& a, const Vector3& b, float t)
{
    return a + (b - a) * t;
}

// Unary minus
inline Vector3 operator-(const Vector3& a)
{
    return Vector3(-a.x, -a.y, -a.z);
}

// Vector subtraction.
inline Vector3 operator-(const Vector3& a, const Vector3& b)
{
    return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

// Vector addition.
inline Vector3 operator+(const Vector3& a, const Vector3& b)
{
    return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

// Vector decrement.
Vector3& operator-=(Vector3& a, const Vector3& b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

// Vector increment.
Vector3& operator+=(Vector3& a, const Vector3& b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

// Multplies all elements of a vector by a scalar.
inline Vector3 operator*(const Vector3& v, float factor)
{
    return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

// Multplies all elements of a vector by a scalar.
inline Vector3 operator*(float factor, const Vector3& v)
{
    return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

// Vector multiplication
inline Vector3 operator*(const Vector3& a, const Vector3& b)
{
    return Vector3(a.x * b.x, a.y * b.y, a.z * b.z);
}

// Divides all elements of a vector by a scalar.
inline Vector3 operator/(const Vector3& v, float denominator)
{
    return Vector3(v.x / denominator, v.y / denominator, v.z / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
inline Vector3 operator/(float numerator, const Vector3& v)
{
    return Vector3(numerator / v.x, numerator / v.y, numerator / v.z);
}

// Vector division.
inline Vector3 operator/(const Vector3& a, const Vector3& b)
{
    return Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
}

// Fuzzy comparison between two Vector3 values.
bool FuzzyEqual(const Vector3& rhs, const Vector3& lhs, float epsilon)
{
    return (FuzzyEqual(rhs.x, lhs.x, epsilon) &&
            FuzzyEqual(rhs.y, lhs.y, epsilon) &&
            FuzzyEqual(rhs.z, lhs.z, epsilon));
}

////////////////////////////////////////////////////////////////////////////////
// Vector4

// Construct from four floats.
inline Vector4::Vector4(float xIn, float yIn, float zIn, float wIn) : x(xIn), y(yIn), z(zIn), w(wIn) {}

// Set from four floats.
inline void Vector4::Set(float xIn, float yIn, float zIn, float wIn)
{
    x = xIn;
    y = yIn;
    z = zIn;
    w = wIn;
}

// Set all elements to zero.
inline void Vector4::SetZero()
{
    x = y = z = w = 0;
}

// Returns a vector with same direction but unit length.
inline Vector4 Vector4::Unit() const
{
    return (*this) / Len();
}

// Returns vector length.
inline float Vector4::Len() const
{
    return sqrt(Dot(*this, *this));
}

// Returns length squared.
inline float Vector4::LenSq() const
{
    return Dot(*this, *this);
}

// Returns a vector with the same direction, but with a Len() <= len.
inline Vector4 Vector4::MinLen(float len) const
{
    float l = Len();
    if (l > len)
        return (*this / l) * len;
    else
        return *this;
}

// const array accessor
inline float Vector4::operator[](int i) const
{
    return *(&x + i);
}

// array accessor
inline float& Vector4::operator[](int i)
{
    return *(&x + i);
}

// Dot product of two vectors.
inline float Dot(const Vector4& a, const Vector4& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

// Linear interpolation between two vectors
inline Vector4 Lerp(const Vector4& a, const Vector4& b, float t)
{
    return a + (b - a) * t;
}

// Unary minus.
inline Vector4 operator-(const Vector4& v)
{
    return Vector4(-v.x, -v.y, -v.z, -v.w);
}

// Vector subtraction.
inline Vector4 operator-(const Vector4& a, const Vector4& b)
{
    return Vector4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

// Vector addition.
inline Vector4 operator+(const Vector4& a, const Vector4& b)
{
    return Vector4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

// Vector decrement.
Vector4& operator-=(Vector4& a, const Vector4& b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    a.w -= b.w;
    return a;
}

// Vector increment.
Vector4& operator+=(Vector4& a, const Vector4& b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
    return a;
}

// Multplies all elements of a vector by a scalar.
inline Vector4 operator*(const Vector4& v, float factor)
{
    return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

// Multplies all elements of a vector by a scalar.
inline Vector4 operator*(float factor, const Vector4& v)
{
    return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

// Vector multiplication.
inline Vector4 operator*(const Vector4& a, const Vector4& b)
{
    return Vector4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

// Divides all elements of a vector by a scalar.
inline Vector4 operator/(const Vector4& v, float denominator)
{
    return Vector4(v.x / denominator, v.y / denominator, v.z / denominator, v.w / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
inline Vector4 operator/(float numerator, const Vector4& v)
{
    return Vector4(numerator / v.x, numerator / v.y, numerator / v.z, numerator / v.w);
}

// Vector division.
inline Vector4 operator/(const Vector4& a, const Vector4& b)
{
    return Vector4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

// Fuzzy comparison between two Vector4 values.
bool FuzzyEqual(const Vector4& rhs, const Vector4& lhs, float epsilon)
{
    return (FuzzyEqual(rhs.x, lhs.x, epsilon) &&
            FuzzyEqual(rhs.y, lhs.y, epsilon) &&
            FuzzyEqual(rhs.z, lhs.z, epsilon) &&
            FuzzyEqual(rhs.w, lhs.w, epsilon));
}

////////////////////////////////////////////////////////////////////////////////
// Quat

// Create from axis and angle
inline Quat Quat::AxisAngle(const Vector3& axis, float angle)
{
    Vector3 n = axis.Unit() * static_cast<float>(sin(angle/2.0f));
    float w = cos(angle/2.0f);
    return Quat(n.x, n.y, n.z, w);
}

// Create identity Quat
inline Quat Quat::Identity()
{
    return Quat(0, 0, 0, 1);
}

// Construct from four floats.
inline Quat::Quat(float iIn, float jIn, float kIn, float rIn) : i(iIn), j(jIn), k(kIn), r(rIn) {}

// Set from four floats.
inline void Quat::Set(float iIn, float jIn, float kIn, float rIn)
{
    i = iIn; j = jIn; k = kIn; r = rIn;
}

// const array accessor
float Quat::operator[](int index) const
{
    return *(&i + index);
}

// array accessor
float& Quat::operator[](int index)
{
    return *(&i + index);
}

// Set all elements to zero.
inline void Quat::SetZero()
{
    i = j = k = r = 0;
}

// Returns this quaternion normalized.
inline Quat Quat::Unit() const
{
    float len = Len();
    return Quat(i / len, j / len, k / len, r / len);
}

// Returns quat length.
inline float Quat::Len() const
{
    return sqrt(Dot(*this, *this));
}

// Returns quat length squared.
inline float Quat::LenSq() const
{
    return Dot(*this, *this);
}

// Rotate a vector.
inline Vector3 Quat::Rotate(const Vector3& v) const
{
    Quat r = (*this) * Quat(v.x, v.y, v.z, 0.0f) * ~(*this);
    return Vector3(r.i, r.j, r.k);
}

// Dot product
inline float Dot(const Quat& a, const Quat& b)
{
    return a.i * b.i + a.j * b.j + a.k * b.k + a.r * b.r;
}

// Linear interpolation between two quaternions
inline Quat Lerp(const Quat& a, const Quat& b, float t)
{
    Quatf diff(b - a);
    diff.Set(diff.i * t, diff.j * t, diff.k * t, diff.r * t);
    return a + diff;
}

// Quaternion conjugate
inline Quat operator~(const Quat& v)
{
    return Quat(-v.i, -v.j, -v.k, v.r);
}

// Unary minus.
inline Quat operator-(const Quat& v)
{
    return Quat(-v.i, -v.j, -v.k, -v.r);
}

// Quaternion subtraction.
inline Quat operator-(const Quat& a, const Quat& b)
{
    return Quat(a.i - b.i, a.j - b.j, a.k - b.k, a.r - b.r);
}

// Quaternion addition.
inline Quat operator+(const Quat& a, const Quat& b)
{
    return Quat(a.i + b.i, a.j + b.j, a.k + b.k, a.r + b.r);
}

// Quaternion multplication.
inline Quat operator*(const Quat& q1, const Quat& q2)
{
    return Quat(q1.i * q2.r + q1.j * q2.k - q1.k * q2.j + q1.r * q2.i,
               -q1.i * q2.k + q1.j * q2.r + q1.k * q2.i + q1.r * q2.j,
                q1.i * q2.j - q1.j * q2.i + q1.k * q2.r + q1.r * q2.k,
               -q1.i * q2.i - q1.j * q2.j - q1.k * q2.k + q1.r * q2.r);
}

// Quaternian exponential, e ^ x
Quat Exp(const Quat& x)
{
    float angle = Vector3(x.i, x.j, x.k).Len();
    Vector3 n;
    if (angle > 0.0001f)
        n = Vector3(x.i, x.j, x.k).Unit() * static_cast<float>(sin(angle / 2));
    else
        n.Set(0,0,0);

    return Quat(n.x, n.y, n.z, cos(angle / 2));
}

// Quaternian logarithm, ln(x)
Quat Log(const Quat& x)
{
    float cos_a = x.r;
    if (cos_a > 1) cos_a = 1;
    if (cos_a < -1) cos_a = -1;

    float sin_a = sqrt(1 - cos_a * cos_a);

    if (fabs(sin_a) < 0.0005)
        sin_a = 1;
    else
        sin_a = 1/sin_a;

    float angle = 2 * acos(cos_a);
    Quat log;
    log.i = x.i * sin_a * angle;
    log.j = x.j * sin_a * angle;
    log.k = x.k * sin_a * angle;
    log.r = 0;
    return log;
}

// Fuzzy comparison between two Quat values.
bool FuzzyEqual(const Quat& rhs, const Quat& lhs, float epsilon)
{
    return (FuzzyEqual(rhs.i, lhs.i, epsilon) &&
            FuzzyEqual(rhs.j, lhs.j, epsilon) &&
            FuzzyEqual(rhs.k, lhs.k, epsilon) &&
            FuzzyEqual(rhs.r, lhs.r, epsilon));
}

////////////////////////////////////////////////////////////////////////////////
// Matrix helpers

inline float Dot3(const Vector3& a, const Vector3& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline float Dot3(const Vector3& a, const Vector4& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline float Dot3(const Vector4& a, const Vector3& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline float Dot3(const Vector4& a, const Vector4& b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

////////////////////////////////////////////////////////////////////////////////
// Matrix

// Create Matrix from three principle axes and a translation.
inline Matrix Matrix::Axes(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis, const Vector3& trans)
{
    Matrix m;
    m.SetXAxis(xAxis);
    m.SetYAxis(yAxis);
    m.SetZAxis(zAxis);
    m.SetTrans(trans);
    return m;
}

// Create a Matrix from four row vectors.
inline Matrix Matrix::Rows(const Vector4& row0In, const Vector4& row1In, const Vector4& row2In, const Vector4& row3In)
{
    Matrix m;
    m.col0.x = row0In.x; m.col1.x = row0In.y; m.col2.x = row0In.z; m.col3.x = row0In.w;
    m.col0.y = row1In.x; m.col1.y = row1In.y; m.col2.y = row1In.z; m.col3.y = row1In.w;
    m.col0.z = row2In.x; m.col1.z = row2In.y; m.col2.z = row2In.z; m.col3.z = row2In.w;
    m.col0.w = row3In.x; m.col1.w = row3In.y; m.col2.w = row3In.z; m.col3.w = row3In.w;
    return m;
}


// Create a Matrix from a translation.
inline Matrix Matrix::Trans(const Vector3& trans)
{
    Matrix m;
    m.SetXAxis(Vector3(1, 0, 0));
    m.SetYAxis(Vector3(0, 1, 0));
    m.SetZAxis(Vector3(0, 0, 1));
    m.SetTrans(trans);
    return m;
}

// Create a Matrix from a quaternion.
inline Matrix Matrix::FromQuat(const Quat& q)
{
    Matrix m;
    m.SetXAxis(q.Rotate(Vector3(1, 0, 0)));
    m.SetYAxis(q.Rotate(Vector3(0, 1, 0)));
    m.SetZAxis(q.Rotate(Vector3(0, 0, 1)));
    m.SetTrans(Vector3(0, 0, 0));
    return m;
}

// Create a Matrix from a Quat and a translation.
Matrix Matrix::QuatTrans(const Quat& q, const Vector3& trans)
{
    Matrix m;
    m.SetXAxis(q.Rotate(Vector3(1, 0, 0)));
    m.SetYAxis(q.Rotate(Vector3(0, 1, 0)));
    m.SetZAxis(q.Rotate(Vector3(0, 0, 1)));
    m.SetTrans(trans);
    return m;
}

// Create a Matrix from a Scale vector, a Quat and a translation.
Matrix Matrix::ScaleQuatTrans(const Vector3& scale, const Quat& q, const Vector3& trans)
{
    Matrix m;
    m.SetXAxis(q.Rotate(Vector3(scale.x, 0, 0)));
    m.SetYAxis(q.Rotate(Vector3(0, scale.y, 0)));
    m.SetZAxis(q.Rotate(Vector3(0, 0, scale.z)));
    m.SetTrans(trans);
    return m;
}

// Create a Matrix from a rotation represented by an axis and an angle.
Matrix Matrix::AxisAngle(const Vector3& axis, float angle)
{
    Matrix m;
    Quat q = Quat::AxisAngle(axis, angle);
    m.SetXAxis(q.Rotate(Vector3(1, 0, 0)));
    m.SetYAxis(q.Rotate(Vector3(0, 1, 0)));
    m.SetZAxis(q.Rotate(Vector3(0, 0, 1)));
    m.SetTrans(Vector3(0, 0, 0));
    return m;
}

// Create an identity Matrix.
inline Matrix Matrix::Identity()
{
    Matrix m;
    m.col0.Set(1, 0, 0, 0);
    m.col1.Set(0, 1, 0, 0);
    m.col2.Set(0, 0, 1, 0);
    m.col3.Set(0, 0, 0, 1);
    return m;
}

// Create a persective transformation Matrix.
inline Matrix Matrix::Frustum(float fovy, float aspect, float nearVal, float farVal)
{
    Matrix m;
    float f = 1.0f / tan(fovy / 2);
    m.col0.Set(f / aspect, 0, 0, 0);
    m.col1.Set(0, f, 0, 0);
    m.col2.Set(0, 0, (farVal + nearVal) / (nearVal - farVal), -1);
    m.col3.Set(0, 0, (2 * farVal * nearVal) / (nearVal - farVal), 0);
    return m;
}

// Create an orthograpic projection Matrix.
inline Matrix Matrix::Ortho(float left, float right, float bottom, float top, float nearVal, float farVal)
{
    Matrix m;
    float tx = -(right + left) / (right - left);
    float ty = -(top + bottom) / (top - bottom);
    float tz = -(farVal + nearVal) / (farVal - nearVal);
    m.col0.Set(2.0f / (right - left), 0, 0, 0);
    m.col1.Set(0, 2.0f / (top - bottom), 0, 0);
    m.col2.Set(0, 0, -2.0f / (farVal - nearVal), 0);
    m.col3.Set(tx, ty, tz, 1);
    return m;
}

// Create a look at matrix. (x is forward)
inline Matrix Matrix::LookAt(const Vector3& eye, const Vector3& target, const Vector3& up)
{
    const float Epsilon = 0.001f;
    const float EpsilonSq = Epsilon * Epsilon;

    Vector3 x;
    if ((target - eye).LenSq() < EpsilonSq)
        x = Vector3(1, 0, 0);  // target == eye, so pick (1,0,0) for the x-axis
    else
        x = (target - eye).Unit();

    Vector3 u;
    if (up.LenSq() < EpsilonSq)
        u = Vector3(0, 1, 0);  // up is zero, so pick (0,1,0) for the y-axis.
    else
        u = up;

    if (fabs(Dot(u, x)) > (u.LenSq() - EpsilonSq)) // are u & x parallel?
        u = (u + Vector3(0.2,0.2,0.2)).Unit();  // nudge u so it's no longer parallel.

    Vector3 z = Cross(x, u).Unit();
    Vector3 y = Cross(z, x).Unit();

    Matrix m;
    m.SetXAxis(x);
    m.SetYAxis(y);
    m.SetZAxis(z);
    m.SetTrans(eye);
    return m;
}

// Axes accessors
inline Vector3 Matrix::GetXAxis() const
{
    return Vector3(col0.x, col0.y, col0.z);
}

inline Vector3 Matrix::GetYAxis() const
{
    return Vector3(col1.x, col1.y, col1.z);
}

inline Vector3 Matrix::GetZAxis() const
{
    return Vector3(col2.x, col2.y, col2.z);
}

inline Vector3 Matrix::GetTrans() const
{
    return Vector3(col3.x, col3.y, col3.z);
}

inline void Matrix::SetXAxis(const Vector3& xAxis)
{
    col0.Set(xAxis.x, xAxis.y, xAxis.z, 0);
}

inline void Matrix::SetYAxis(const Vector3& yAxis)
{
    col1.Set(yAxis.x, yAxis.y, yAxis.z, 0);
}

inline void Matrix::SetZAxis(const Vector3& zAxis)
{
    col2.Set(zAxis.x, zAxis.y, zAxis.z, 0);
}

inline void Matrix::SetTrans(const Vector3& trans)
{
    col3.Set(trans.x, trans.y, trans.z, 1);
}

// Row accessors
inline Vector4 Matrix::GetRow(int index) const
{
    float* p = (float*)&col0;
    return Vector4(p[0 + index], p[4 + index], p[8 + index], p[12 + index]);
}

inline void Matrix::SetRow(int index, const Vector4& row)
{
    float* p = (float*)&col0;
    p[0 + index] = row.x;
    p[4 + index] = row.y;
    p[8 + index] = row.z;
    p[12 + index] = row.w;
}

// Column accessors
inline Vector4 Matrix::GetCol(int index) const
{
    Vector4* v = (Vector4*)&col0;
    return v[index];
}

inline void Matrix::SetCol(int index, const Vector4& col)
{
    Vector4* v = (Vector4*)&col0;
    v[index] = col;
}

// Element accessors
inline const float& Matrix::Elem(int r, int c) const
{
    return ((float*)&col0)[c * 4 + r];
}

inline float& Matrix::Elem(int r, int c)
{
    return ((float*)&col0)[c * 4 + r];
}

// Multiplies by uniform scale.
inline void Matrix::SetScale(float scale)
{
    SetXAxis(GetXAxis() * scale);
    SetYAxis(GetYAxis() * scale);
    SetZAxis(GetZAxis() * scale);
}

// Multiplies by non-uniform scale.
inline void Matrix::SetScale(const Vector3& scale)
{
    SetXAxis(GetXAxis() * scale.x);
    SetYAxis(GetYAxis() * scale.y);
    SetZAxis(GetZAxis() * scale.z);
}

// Returns the rotation component of this Matrix.
Quat Matrix::GetQuat() const
{
    // TODO: add to unit test!

    int x = 0, y = 1, z = 2;
    // create a matrix with no scale
    Matrix m = Axes(GetXAxis().Unit(), GetYAxis().Unit(), GetZAxis().Unit(), Vector3(0,0,0));
    float trace = m.Elem(0,0) + m.Elem(1,1) + m.Elem(2,2);
    Vector4 q;
    if (trace < -1.0)
    {
        int i = x, j = y, k = z;
        if (m.Elem(y,y) > m.Elem(x,x)) { i = y; j = z; k = x; }
        if (m.Elem(z,z) > m.Elem(i,i)) { i = z; j = x; k = y; }
        float r = sqrt(m.Elem(i,i) - m.Elem(j,j) - m.Elem(k,k) + 1);
        float s = 0.5 / r;
        q[i] = 0.5 * r;
        q[j] = (m.Elem(i,j) + m.Elem(j,i)) * s;
        q[k] = (m.Elem(k,i) + m.Elem(i,k)) * s;
        q[3] = (m.Elem(k,j) + m.Elem(j,k)) * s;
    }
    else
    {
        float r = sqrt(trace + 1.0);
        float s = 0.5 / r;
        q[0] = (m.Elem(z,y) - m.Elem(y,z)) * s;
        q[1] = (m.Elem(x,z) - m.Elem(z,x)) * s;
        q[2] = (m.Elem(y,x) - m.Elem(x,y)) * s;
        q[3] = 0.5 * r;
    }
    return Quat(q.x, q.y, q.z, q.w);
}

// Multiply the 3x3 component of this Matrix with a column vector.
inline Vector3 Matrix::Mul3x3(const Vector3& v) const
{
    Vector4 r = col0 * v.x;
    r += col1 * v.y;
    r += col2 * v.z;
    return Vector3(r.x, r.y, r.z);
}

// Multiply the 3x4 component of this Matrix with a column vector. (w component of vector is 1.0)
inline Vector3 Matrix::Mul3x4(const Vector3& v) const
{
    Vector4 r = col0 * v.x;
    r += col1 * v.y;
    r += col2 * v.z;
    r += col3;
    return Vector3(r.x, r.y, r.z);
}

// Multiply this Matrix with a column vector.
inline Vector4 Matrix::Mul4x4(const Vector4& v) const
{
    Vector4 r = col0 * v.x;
    r += col1 * v.y;
    r += col2 * v.z;
    r += col3 * v.w;
    return r;
}

// Returns the transpose of this Matrix
inline Matrix Matrix::Transpose() const
{
    return Rows(col0, col1, col2, col3);
}

// If the 3x3 portion of this Matrix is Orthogonal (i.e. columns are orthogonal unit vectors)
// this will return the Inverse of that matrix.
Matrix Matrix::OrthoInverse() const
{
    Matrix r(*this);
    r.SetTrans(Vector3(0, 0, 0));
    r = r.Transpose();
    r.SetTrans(-r.Mul3x4(GetTrans()));
    return r;
}

// Full 4x4 Matrix Inverse, returns Identity if matrix has no inverse.
Matrix Matrix::FullInverse() const
{
    Matrix m;
    if (::FullInverse((*this), m))
        return m;
    else
        return Identity();
}

// Matrix addition
Matrix operator+(const Matrix& a, const Matrix& b)
{
    Matrix m;
    m.col0 = a.col0 + b.col0;
    m.col1 = a.col1 + b.col1;
    m.col2 = a.col2 + b.col2;
    m.col3 = a.col3 + b.col3;
    return m;
}

// Matrix subtraction
Matrix operator-(const Matrix& a, const Matrix& b)
{
    Matrix m;
    m.col0 = a.col0 - b.col0;
    m.col1 = a.col1 - b.col1;
    m.col2 = a.col2 - b.col2;
    m.col3 = a.col3 - b.col3;
    return m;
}

// Matrix multiplication
Matrix operator*(const Matrix& a, const Matrix& b)
{
    Vector4 a_row0 = a.GetRow(0);
    Vector4 a_row1 = a.GetRow(1);
    Vector4 a_row2 = a.GetRow(2);
    Vector4 a_row3 = a.GetRow(3);

    return Matrix::Rows(Vector4(Dot(a_row0, b.col0), Dot(a_row0, b.col1),
                                                Dot(a_row0, b.col2), Dot(a_row0, b.col3)),
                                Vector4(Dot(a_row1, b.col0), Dot(a_row1, b.col1),
                                                Dot(a_row1, b.col2), Dot(a_row1, b.col3)),
                                Vector4(Dot(a_row2, b.col0), Dot(a_row2, b.col1),
                                                Dot(a_row2, b.col2), Dot(a_row2, b.col3)),
                                Vector4(Dot(a_row3, b.col0), Dot(a_row3, b.col1),
                                                Dot(a_row3, b.col2), Dot(a_row3, b.col3)));
}

template <typename T>
void AbaciSwap(T& a, T& b)
{
    T temp = a;
    a = b;
    b = temp;
}

// Full 4x4 Matrix inverse, returns false if Matrix has no inverse.
bool FullInverse(const Matrix& m, Matrix& result)
{
    // Gaussian-Jordan Elimination, TODO: There are more numerically stable ways to do this...

    float temp[4][8];
    float* row[4];

    // initialize the r pointers
    for (int r = 0; r < 4; ++r)
        row[r] = temp[r];

    // initialize the augmented temp matrix.
    // the first four columns are from m
    temp[0][0] = m.col0.x; temp[0][1] = m.col1.x; temp[0][2] = m.col2.x; temp[0][3] = m.col3.x;
    temp[1][0] = m.col0.y; temp[1][1] = m.col1.y; temp[1][2] = m.col2.y; temp[1][3] = m.col3.y;
    temp[2][0] = m.col0.z; temp[2][1] = m.col1.z; temp[2][2] = m.col2.z; temp[2][3] = m.col3.z;
    temp[3][0] = m.col0.w; temp[3][1] = m.col1.w; temp[3][2] = m.col2.w; temp[3][3] = m.col3.w;

    // the second four are identity
    temp[0][4] = 1; temp[0][5] = 0; temp[0][6] = 0; temp[0][7] = 0;
    temp[1][4] = 0; temp[1][5] = 1; temp[1][6] = 0; temp[1][7] = 0;
    temp[2][4] = 0; temp[2][5] = 0; temp[2][6] = 1; temp[2][7] = 0;
    temp[3][4] = 0; temp[3][5] = 0; temp[3][6] = 0; temp[3][7] = 1;

    // bubble up row with largest leading number (partial pivot)
    if (fabs(row[0][0]) < fabs(row[1][0]))
        AbaciSwap(row[0], row[1]);
    if (fabs(row[0][0]) < fabs(row[2][0]))
        AbaciSwap(row[0], row[2]);
    if (fabs(row[0][0]) < fabs(row[3][0]))
        AbaciSwap(row[0], row[3]);

    if (fabs(row[0][0]) < 0.00001)  // column is all zeros, there is no inverse.
    {
        return false;
    }

    // mult row[0] by 1/row[0][0].  To introduce a leading 1.
    float s = 1 / row[0][0];
    for (int c = 0; c < 8; ++c) row[0][c] *= s;

    // add multiples of top row to lower rows so that all entries below leading 1 become zeros.
    for (int r = 1; r < 4; ++r)
    {
        float s = row[r][0];
        for (int c = 0; c < 8; ++c) row[r][c] -= s * row[0][c];
    }

    // move row with largest leading number
    if (fabs(row[1][1]) < fabs(row[2][1]))
        AbaciSwap(row[1], row[2]);
    if (fabs(row[1][1]) < fabs(row[3][1]))
        AbaciSwap(row[1], row[3]);

    if (fabs(row[1][1]) < 0.00001)  // column is all zeros, there is no inverse.
    {
        return false;
    }

    // mult row[1] by 1/row[1][1].  To introduce a leading 1.
    s = 1 / row[1][1];
    for (int c = 0; c < 8; ++c) row[1][c] *= s;

    // add multiples of top row to lower rows so that all entries below leading 1 become zeros.
    for (int r = 2; r < 4; ++r)
    {
        float s = row[r][1];
        for (int c = 0; c < 8; ++c) row[r][c] -= s * row[1][c];
    }

    // move row with largest leading number
    if (fabs(row[2][2]) < fabs(row[3][2]))
        AbaciSwap(row[2], row[3]);

    if (fabs(row[2][2]) < 0.00001)  // column is all zeros, there is no inverse.
    {
        return false;
    }

    // mult row[2] by 1/row[2][2].  To introduce a leading 1.
    s = 1.0f / row[2][2];
    for (int c = 0; c < 8; ++c) row[2][c] *= s;

    // add multiples of top row to lower rows so that all entries below leading 1 become zeros.
    for (int r = 3; r < 4; ++r)
    {
        float s = row[r][2];
        for (int c = 0; c < 8; ++c) row[r][c] -= s * row[2][c];
    }

    // mult row[3] by 1/row[3][3]. To introduce a leading 1.
    s = 1.0f / row[3][3];
    for (int c = 0; c < 8; ++c) row[3][c] *= s;

    // at this point row matrix should be in row-echelon form.

    // add multiples of row[3] to above rows, to zero out that column
    for (int r = 0; r < 3; ++r)
    {
        float s = row[r][3];
        for (int c = 0; c < 8; ++c) row[r][c] -= s * row[3][c];
    }

    // add multiples of row[2] to above rows, to zero out that column
    for (int r = 0; r < 2; ++r)
    {
        float s = row[r][2];
        for (int c = 0; c < 8; ++c) row[r][c] -= s * row[2][c];
    }

    // add multiples of row[1] to above row, to zero out that column
    for (int r = 0; r < 1; ++r)
    {
        float s = row[r][1];
        for (int c = 0; c < 8; ++c) row[r][c] -= s * row[1][c];
    }

    // init result
    result.col0.Set(row[0][4], row[1][4], row[2][4], row[3][4]);
    result.col1.Set(row[0][5], row[1][5], row[2][5], row[3][5]);
    result.col2.Set(row[0][6], row[1][6], row[2][6], row[3][6]);
    result.col3.Set(row[0][7], row[1][7], row[2][7], row[3][7]);

    return true;
}


// Print to stdout.
void PrintMatrix(const Matrix& m)
{
    // display in mathematical column-vector style, which is different then memory layout.
    printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.col0.x, m.col1.x, m.col2.x, m.col3.x);
    printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.col0.y, m.col1.y, m.col2.y, m.col3.y);
    printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.col0.z, m.col1.z, m.col2.z, m.col3.z);
    printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.col0.w, m.col1.w, m.col2.w, m.col3.w);
}

////////////////////////////////////////////////////////////////////////////////
// Complex

// Construct from two floats
inline Complex::Complex(float rIn, float iIn) : r(rIn), i(iIn) {}

// Construct from a Vector2
inline Complex::Complex(const Vector2& vector2In) :
    r(vector2In.x), i(vector2In.y) {}

// Set from two floats
inline void Complex::Set(float rIn, float iIn)
{
    r = rIn;
    i = iIn;
}

// Length
inline float Complex::Len() const
{
    return sqrt(Dot(*this, *this));
}

// Square of length
inline float Complex::LenSq() const
{
    return Dot(*this, *this);
}

// Unit
inline Complex Complex::Unit() const
{
    return *this / Complex(Len(), 0);
}

// Returns a vector with the same direction, but with a Len() <= len.
inline Complex Complex::MinLen(float len) const
{
    float l = Len();
    if (l > len)
        return (*this / Complex(l, 0)) * Complex(len, 0);
    else
        return *this;
}

// Dot product
inline float Dot(const Complex& a, const Complex& b)
{
    return a.r * b.r + a.i * b.i;
}

// Complex conjugate
inline Complex operator~(const Complex& a)
{
    return Complex(a.r, -a.i);
}

// Unary minus.
inline Complex operator-(const Complex& a)
{
    return Complex(-a.r, -a.i);
}

// Unary plus.
inline Complex operator+(const Complex& a)
{
    return Complex(a.r, a.i);
}

// Complex addition.
inline Complex operator+(const Complex& a, const Complex& b)
{
    return Complex(a.r + b.r, a.i + b.i);
}

// Complex subtraction.
inline Complex operator-(const Complex& a, const Complex& b)
{
    return Complex(a.r - b.r, a.i - b.i);
}

// Complex multiplication.
inline Complex operator*(const Complex& a, const Complex& b)
{
    float aa = a.r;
    float bb = a.i;
    float cc = b.r;
    float dd = b.i;

    // four muls & two adds
    // return Complex(aa * cc - (bb * dd), aa * dd + bb * cc);

    float a_c = aa * cc;
    float b_d = bb * dd;

    // three muls & five adds (faster?)
    return Complex(a_c - b_d, (aa + bb) * (cc + dd) - a_c - b_d);
}

// Multiplication by a real number.
inline Complex operator*(float scalar, const Complex& c)
{
    return Complex(scalar, 0) * c;
}

// Multiplication by a real number.
inline Complex operator*(const Complex& c, float scalar)
{
    return c * Complex(scalar, 0);
}

// Complex division.
inline Complex operator/(const Complex& a, const Complex& b)
{
    float aa = a.r;
    float bb = a.i;
    float cc = b.r;
    float dd = b.i;
    float denom = cc * cc + dd * dd;

    return Complex((aa * cc + bb * dd) / denom, (bb * cc - aa * dd) / denom);
}

// Complex division by a scalar denominator (x + 0i)
Complex operator/(const Complex& c, float denominator)
{
    return c / Complex(denominator, 0);
}

// Complex division with a scalar numerator (x + 0i)
Complex operator/(float numerator, const Complex& c)
{
    return Complex(numerator, 0) / c;
}

// Square root
Complex sqrt(const Complex& z)
{
    float x = z.r;
    float y = z.i;

    if (x == float())
    {
        float t = sqrt(fabs(y) / 2);
        return Complex(t, y < float() ? -t : t);
    }
    else
    {
        float t = sqrt(2 * (z.Len() + fabs(x)));
        float u = t / 2;
        return x > float() ? Complex(u, y / t) :
                              Complex(fabs(y) / t, y < float() ? -u : u);
    }
}

// Exponent e ^ z
Complex Exp(const Complex& z)
{
    float e = exp(z.r);
    return Complex(e * cos(z.i), e * sin(z.i));
}

// e ^ (0 + xi)
inline Complex ExpI(float x)
{
    return Complex(cos(x), sin(x));
}

// Natural Logarithm, ln(z)
Complex Log(const Complex& z)
{
    return Complex(log(z.Len()), atan2(z.i, z.r));
}

// Fuzzy comparison between two Complex values.
bool FuzzyEqual(const Complex& rhs, const Complex& lhs, float epsilon)
{
    return (FuzzyEqual(rhs.r, lhs.r, epsilon) &&
            FuzzyEqual(rhs.i, lhs.i, epsilon));
}

#endif
