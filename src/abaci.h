#ifndef ABACI_H
#define ABACI_H

#include <math.h>

#ifdef ABACI_NAMESPACE
namespace ABACI_NAMESPACE
{
#endif

#define PI M_PI

// Convert from degrees to radians
float DegToRad(float deg);

// Convert from radians to degrees.
float RadToDeg(float rad);

// Clamps value between min & max.
float Clamp(float value, float min, float max);

// Limits angle between -PI & PI using modulus arithmetic
float LimitPi(float theta);

// Limits angle between zero & 2 PI using modulus arithmetic
float Mod2Pi(float theta);

// Fuzzy comparison between two float values.
bool FuzzyEqual(float rhs, float lhs, float epsilon = 0.0001);

float Lerp(float a, float b, float t);

// forward declare
class Complex;

// returns a random integer between min and max.
// Note: closed interval. i.e. the values of min & max can be returned.
int RandomInt(int min, int max);

// returns a random scalar between min and max
// Note: closed interval. i.e. the values of min & max can be returned.
float RandomFloat(float min, float max);

//////////////////////////////////////////////////////

struct Vector2
{
    // Generates a random vector on the unit circle.
    static Vector2 RandomUnitVector();

    // Uninitialized by default.
    Vector2() {}

    // Construct from two floats
    Vector2(float xIn, float yIn);

    // Construct from complex
    Vector2(const Complex& complexIn);

    // Set from two floats
    void Set(float xIn, float yIn);

    // Set from complex
    void Set(const Complex& complexIn);

    // Sets all elements to zero.
    void SetZero();

    // const array accessor
    float operator[](int i) const;

    // array accessor
    float& operator[](int i);

    // Returns a vector with same direction but unit length.
    Vector2 Unit() const;

    // Returns vector length.
    float Len() const;

    // Returns length squared.
    float LenSq() const;

    // Returns a vector with the same direction, but with a Len() <= len.
    Vector2 MinLen(float len) const;

    float x;
    float y;
};

typedef Vector2 Vector2f;

// Dot product of two vectors.
float Dot(const Vector2& a, const Vector2& b);

// Linear interpolation between two vectors
Vector2 Lerp(const Vector2& a, const Vector2& b, float t);

// Unary minus
Vector2 operator-(const Vector2& a);

// Vector subtraction.
Vector2 operator-(const Vector2& a, const Vector2& b);

// Vector addition.
Vector2 operator+(const Vector2& a, const Vector2& b);

// Vector decrement.
Vector2& operator-=(Vector2& a, const Vector2& b);

// Vector increment.
Vector2& operator+=(Vector2& a, const Vector2& b);

// Multplies all elements of a vector by a scalar.
Vector2 operator*(const Vector2& v, float scalar);

// Multplies all elements of a vector by a scalar.
Vector2 operator*(float factor, const Vector2& v);

// Vector multiplication
Vector2 operator*(const Vector2& a, const Vector2& b);

// Divides all elements of a vector by a scalar.
Vector2 operator/(const Vector2& v, float denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
Vector2 operator/(float numerator, const Vector2& v);

// Vector division.
Vector2 operator/(const Vector2& a, const Vector2& b);

// Fuzzy comparison between two Vector2 values.
bool FuzzyEqual(const Vector2& rhs, const Vector2& lhs, float epsilon = 0.001f);

//////////////////////////////////////////////////////

struct Vector3
{
    // Generates a random vector on the unit sphere.
    static Vector3 RandomUnitVector();

    // Uninitialized by default.
    Vector3() {}

    // Construct from three floats
    Vector3(float xIn, float yIn, float zIn);

    // Set from three floats
    void Set(float xIn, float yIn, float zIn);

    // Set all elements to zero
    void SetZero();

    // const array accessor
    float operator[](int i) const;

    // array accessor
    float& operator[](int i);

    // Returns a vector with same direction but unit length.
    Vector3 Unit() const;

    // Returns vector length.
    float Len() const;

    // Returns length squared.
    float LenSq() const;

    // Returns a vector with the same direction, but with a Len() <= len.
    Vector3 MinLen(float len) const;

    // Returns basis vectors i, j & k.  Such that i is parallel to this.
    // j & k are orthonal to i and each other.
    void Basis(Vector3& iOut, Vector3& jOut, Vector3& kOut);

    float x;
    float y;
    float z;
};

typedef Vector3 Vector3f;

// Dot product of two vectors.
float Dot(const Vector3& a, const Vector3& b);

// Cross product of two vectors.
Vector3 Cross(const Vector3& a, const Vector3& b);

// Linear interpolation between two vectors
Vector3 Lerp(const Vector3& a, const Vector3& b, float t);

// Unary minus.
Vector3 operator-(const Vector3& a);

// Vector subtraction.
Vector3 operator-(const Vector3& a, const Vector3& b);

// Vector addition.
Vector3 operator+(const Vector3& a, const Vector3& b);

// Vector addition.
Vector3 operator+(const Vector3& a, const Vector3& b);

// Vector decrement.
Vector3& operator-=(Vector3& a, const Vector3& b);

// Multplies all elements of a vector by a scalar.
Vector3 operator*(const Vector3& v, float factor);

// Multplies all elements of a vector by a scalar.
Vector3 operator*(float factor, const Vector3& v);

// Vector multiplication
Vector3 operator*(const Vector3& a, const Vector3& b);

// Divides all elements of a vector by a scalar.
Vector3 operator/(const Vector3& v, float denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
Vector3 operator/(float numerator, const Vector3& v);

// Vector division.
Vector3 operator/(const Vector3& a, const Vector3& b);

// Fuzzy comparison between two Vector3 values.
bool FuzzyEqual(const Vector3& rhs, const Vector3& lhs, float epsilon = 0.001f);

//////////////////////////////////////////////////////

struct Vector4
{
    // Uninitialized by default.
    Vector4() {}

    // Construct from four floats.
    Vector4(float xIn, float yIn, float zIn, float wIn);

    // Set from four floats.
    void Set(float xIn, float yIn, float zIn, float wIn);

    // Set all elements to zero.
    void SetZero();

    // Returns a vector with same direction but unit length.
    Vector4 Unit() const;

    // Returns vector length.
    float Len() const;

    // Returns length squared.
    float LenSq() const;

    // Returns a vector with the same direction, but with a Len() <= len.
    Vector4 MinLen(float len) const;

    // const array accessor
    float operator[](int i) const;

    // array accessor
    float& operator[](int i);

    float x;
    float y;
    float z;
    float w;
};

typedef Vector4 Vector4f;

// Dot product of two vectors.
float Dot(const Vector4& a, const Vector4& b);

// Linear interpolation between two vectors
Vector4 Lerp(const Vector4& a, const Vector4& b, float t);

// Unary minus.
Vector4 operator-(const Vector4& v);

// Vector subtraction.
Vector4 operator-(const Vector4& a, const Vector4& b);

// Vector addition.
Vector4 operator+(const Vector4& a, const Vector4& b);

// Vector addition.
Vector4 operator+(const Vector4& a, const Vector4& b);

// Vector decrement.
Vector4& operator-=(Vector4& a, const Vector4& b);

// Multplies all elements of a vector by a scalar.
Vector4 operator*(const Vector4& v, float factor);

// Multplies all elements of a vector by a scalar.
Vector4 operator*(float factor, const Vector4& v);

// Vector multiplication.
Vector4 operator*(const Vector4& a, const Vector4& b);

// Divides all elements of a vector by a scalar.
Vector4 operator/(const Vector4& v, float denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
Vector4 operator/(float numerator, const Vector4& v);

// Vector division.
Vector4 operator/(const Vector4& a, const Vector4& b);

// Fuzzy comparison between two Vector4 values.
bool FuzzyEqual(const Vector4& rhs, const Vector4& lhs, float epsilon = 0.001f);

//////////////////////////////////////////////////////

struct Quat
{
    // Create from axis and angle
    static Quat AxisAngle(const Vector3& axis, float angle);

    // Create identity Quat
    static Quat Identity();

    // Uninitialized by default.
    Quat() {}

    // Construct from four floats.
    Quat(float iIn, float jIn, float kIn, float rIn);

    // const array accessor
    float operator[](int i) const;

    // array accessor
    float& operator[](int i);

    // Set from four floats.
    void Set(float iIn, float jIn, float kIn, float rIn);

    // Set all elements to zero.
    void SetZero();

    // Returns this quaternion normalized.
    Quat Unit() const;

    // Returns quat length.
    float Len() const;

    // Returns quat length squared.
    float LenSq() const;

    // Rotate a vector.
    Vector3 Rotate(const Vector3& v) const;

    float i;
    float j;
    float k;
    float r;  // real part
};

typedef Quat Quatf;

// Dot product
float Dot(const Quat& a, const Quat& b);

// Linear interpolation between two vectors
Quat Lerp(const Quat& a, const Quat& b, float t);

// Quaternion conjugate
Quat operator~(const Quat& v);

// Unary minus.
Quat operator-(const Quat& v);

// Quaternion subtraction.
Quat operator-(const Quat& a, const Quat& b);

// Quaternion addition.
Quat operator+(const Quat& a, const Quat& b);

// Quaternion multplication.
Quat operator*(const Quat& a, const Quat& b);

// Quaternian exponential, e ^ x
Quat Exp(const Quat& x);

// Quaternian logarithm, ln(x)
Quat Log(const Quat& x);

// Fuzzy comparison between two Quat values.
bool FuzzyEqual(const Quat& rhs, const Quat& lhs, float epsilon = 0.001f);

//////////////////////////////////////////////////////

// Memory layout is column-major for OpenGL.
// However, the API uses row-major notation. Vector4f vv = (A * B * C).Mul4x4(v), means
// apply transformation C to v then B then A.
// This follows the mathematical convention of using column vectors.
struct Matrix
{
    // Create a Matrix from three principle axes and a translation.
    static Matrix Axes(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis, const Vector3& trans = Vector3(0,0,0));

    // Create a Matrix from four row vectors.
    static Matrix Rows(const Vector4& row0In, const Vector4& row1In, const Vector4& row2In, const Vector4& row3In);

    // Create a Matrix from a translation.
    static Matrix Trans(const Vector3& trans);

    // Create a Matrix from a quaternion.
    static Matrix FromQuat(const Quat& q);

    // Create a Matrix from a Quat and a translation.
    static Matrix QuatTrans(const Quat& q, const Vector3& trans);

    // Create a Matrix from a Scale vector, a Quat and a translation.
    static Matrix ScaleQuatTrans(const Vector3& scale, const Quat& rot, const Vector3& trans);

    // Create a Matrix from a rotation represented by an axis and an angle.
    static Matrix AxisAngle(const Vector3& axis, float angle);

    // Create an identity Matrix.
    static Matrix Identity();

    // Create a persective projection Matrix.
    static Matrix Frustum(float fovy, float aspect, float nearVal, float farVal);

    // Create an orthograpic projection Matrix.
    static Matrix Ortho(float left, float right, float bottom, float top, float nearVal, float farVal);

    // Create a look at matrix. (x is forward)
    static Matrix LookAt(const Vector3& eye, const Vector3& target, const Vector3& up);

    // Uninitialized by default.
    Matrix() {}

    // Axes accessors
    Vector3 GetXAxis() const;
    Vector3 GetYAxis() const;
    Vector3 GetZAxis() const;
    Vector3 GetTrans() const;

    void SetXAxis(const Vector3& xAxis);
    void SetYAxis(const Vector3& yAxis);
    void SetZAxis(const Vector3& zAxis);
    void SetTrans(const Vector3& trans);

    // Row accessors
    Vector4 GetRow(int index) const;
    void SetRow(int index, const Vector4& row);

    // Column accessors
    Vector4 GetCol(int index) const;
    void SetCol(int index, const Vector4& col);

    // Element accessors
    const float& Elem(int r, int c) const;
    float& Elem(int r, int c);

    // Multiplies by uniform scale.
    void SetScale(float scale);

    // Multiplies by non-uniform scale.
    void SetScale(const Vector3& scale);

    // Returns the rotation component of this Matrix.
    Quat GetQuat() const;

    // Multiply the 3x3 component of this Matrix with a column vector.
    Vector3 Mul3x3(const Vector3& v) const;

    // Multiply the 3x4 component of this Matrix with a column vector. (w component of vector is 1.0)
    Vector3 Mul3x4(const Vector3& v) const;

    // Multiply this Matrix with a column vector.
    Vector4 Mul4x4(const Vector4& v) const;

    // Returns the transpose of this Matrix
    Matrix Transpose() const;

    // If the 3x3 portion of this Matrix is Orthogonal (i.e. columns are orthogonal unit vectors)
    // this will return the Inverse of that matrix.
    Matrix OrthoInverse() const;

    // Full 4x4 Matrix Inverse, returns Identity if matrix has no inverse.
    Matrix FullInverse() const;

    Vector4 col0;
    Vector4 col1;
    Vector4 col2;
    Vector4 col3;
};

typedef Matrix Matrixf;

// Matrix addition
Matrix operator+(const Matrix& a, const Matrix& b);

// Matrix subtraction
Matrix operator-(const Matrix& a, const Matrix& b);

// Matrix multiplication
Matrix operator*(const Matrix& a, const Matrix& b);

// Full 4x4 Matrix inverse, returns false if Matrix has no inverse.
bool FullInverse(const Matrix& m, Matrix& result);

// Print to stdout.
void PrintMatrix(const Matrix& m);

//////////////////////////////////////////////////////

struct Complex
{
    // Uninitialized by default.
    Complex() {}

    // Construct from two floats
    Complex(float rIn, float iIn);

    // Construct from a Vector2
    Complex(const Vector2& vector2In);

    // Set from two floats
    void Set(float rIn, float iIn);

    // Length
    float Len() const;

    // Square of length
    float LenSq() const;

    // Unit
    Complex Unit() const;

    // Returns a vector with the same direction, but with a Len() <= len.
    Complex MinLen(float len) const;

    float r;
    float i;
};

typedef Complex Complexf;

// Dot product
float Dot(const Complex& a, const Complex& b);

// Complex conjugate
Complex operator~(const Complex& a);

// Unary minus.
Complex operator-(const Complex& a);

// Unary plus.
Complex operator+(const Complex& a);

// Complex addition.
Complex operator+(const Complex& a, const Complex& b);

// Complex subtraction.
Complex operator-(const Complex& a, const Complex& b);

// Complex multiplication.
Complex operator*(const Complex& a, const Complex& b);

// Multiplication by a real number.
Complex operator*(float scalar, const Complex& c);

// Multiplication by a real number.
Complex operator*(const Complex& c, float scalar);

// Complex division.
Complex operator/(const Complex& a, const Complex& b);

// Complex division by a scalar denominator (x + 0i)
Complex operator/(const Complex& c, float denominator);

// Complex division with a scalar numerator (x + 0i)
Complex operator/(float numerator, const Complex& c);

// Square root
Complex Sqrt(const Complex& z);

// Exponent e ^ z
Complex Exp(const Complex& z);

// e ^ (0 + xi)
Complex ExpI(float x);

// Natural Logarithm, ln(z)
Complex Log(const Complex& z);

// Fuzzy comparison between two Complex values.
bool FuzzyEqual(const Complex& rhs, const Complex& lhs, float epsilon = 0.001f);

// inlines
#include "abaciinlines.cpp"

#ifdef ABACI_NAMESPACE
} // namespace
#endif

#endif
