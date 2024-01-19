#pragma once

#include <math.h>
#include <ostream>

class Math
{
public:
  static constexpr float Clamp(double, double = 0.0, double = 1.0);
  static constexpr float Mix(float, float, float);
  static constexpr float Sign(float);
  static constexpr float Abs(float);
  static constexpr float Unit(float, float);
  static constexpr float Binom(float, float);

  // Minimum and maximum
  static constexpr float Min(float, float);
  static constexpr float Max(float, float);
  static constexpr float Min(float, float, float);
  static constexpr float Max(float, float, float);

  static constexpr double DegreeToRadian(double);
  static constexpr double RadianToDegree(double);
};

/*!
\brief Clamp a double value between two bounds.
\param x Input value.
\param a, b Lower and upper bounds.
*/
inline constexpr float Math::Clamp(double x, double a, double b)
{
  return (x < a ? a : (x > b ? b : x));
}

inline constexpr float Math::Mix(float x, float a, float b)
{
  return x + b * (a - x);
}

inline constexpr float Math::Sign(float a)
{
  if (a < 0) {
    return -1.0;
  }
  else {
    return 1.0;
  }
}

inline constexpr float Math::Abs(float a)
{
  return a > 0 ? a : -a;
}

inline constexpr float Math::Unit(float a, float b)
{
  return a/b ;
}

inline constexpr float Math::Binom(float n, float k)
{
  if (k == 0 || k == n) {
        return 1;
    }
   else {
        return Binom(n - 1, k - 1) + Binom(n - 1, k);
    }
}

/*!
\brief Minimum of two reals.
\param a, b Real values.
*/
inline constexpr float Math::Min(float a, float b)
{
  return (a < b ? a : b);
}

/*!
\brief Maximum of two reals.
\param a, b Real values.
*/
inline constexpr float Math::Max(float a, float b)
{
  return (a > b ? a : b);
}

/*!
\brief Maximum of three reals.
\param a, b, c Real values.
*/
inline constexpr float Math::Max(float a, float b, float c)
{
  return Math::Max(Math::Max(a, b), c);
}

/*!
\brief Minimum of three reals.
\param a, b, c Real values.
*/
inline constexpr float Math::Min(float a, float b, float c)
{
  return Math::Min(Math::Min(a, b), c);
}

/*!
\brief Convert degrees to randians.
\param a Angle in degrees.
*/
inline constexpr double Math::DegreeToRadian(double a)
{
  return a * 3.14159265358979323846 / 180.0;
}

/*!
\brief Convert radian to degrees.
\param a Angle in radian.
*/
inline constexpr double Math::RadianToDegree(double a)
{
  return a * 180.0 / 3.14159265358979323846;
}

// Class
class Vector
{
protected:
  double c[3]; //!< Components.
public:
  //! Empty 
  Vector() {}

  explicit Vector(double);
  explicit Vector(double, double, double);

  // Access members
  double& operator[] (int);
  double operator[] (int) const;

  // Unary operators
  Vector operator+ () const;
  Vector operator- () const;

  // Assignment operators
  Vector& operator+= (const Vector&);
  Vector& operator-= (const Vector&);
  Vector& operator*= (const Vector&);
  Vector& operator/= (const Vector&);
  Vector& operator*= (double);
  Vector& operator/= (double);

  // Binary operators
  friend int operator> (const Vector&, const Vector&);
  friend int operator< (const Vector&, const Vector&);

  friend int operator>= (const Vector&, const Vector&);
  friend int operator<= (const Vector&, const Vector&);

  // Binary operators
  friend Vector operator+ (const Vector&, const Vector&);
  friend Vector operator- (const Vector&, const Vector&);

  friend constexpr double operator* (const Vector&, const Vector&);

  friend Vector operator* (const Vector&, double);
  friend Vector operator* (double, const Vector&);
  friend Vector operator/ (const Vector&, double);

  friend Vector operator/ (const Vector&, const Vector&);

  // Boolean functions
  friend int operator==(const Vector&, const Vector&);
  friend int operator!=(const Vector&, const Vector&);

  // Norm
  friend double Norm(const Vector&);
  friend double Norm2(const double&, const double&);
  friend double SquaredNorm(const Vector&);
  friend double SquaredNorm2(const Vector&, const Vector&);

  friend void Normalize(Vector&);
  friend Vector Normalized(const Vector&);

  // Compare functions
  static Vector Min(const Vector&, const Vector&);
  static Vector Max(const Vector&, const Vector&);

  // Abs
  friend Vector Abs(const Vector&);

  // Orthogonal and orthonormal vectors
  Vector Orthogonal() const;
  void Orthonormal(Vector&, Vector&) const;

  friend Vector Lerp(const Vector&, const Vector&, double);
  friend float Clamp(float, float, float);
  static Vector Bilinear(const Vector&, const Vector&, const Vector&, const Vector&, double, double);

  // Scale
  Vector Scaled(const Vector&) const;
  Vector Inverse() const;

  friend std::ostream& operator<<(std::ostream&, const Vector&);

public:
  static const Vector Null; //!< Null vector.
  static const Vector X; //!< Vector(1,0,0).
  static const Vector Y; //!< Vector(0,1,0).
  static const Vector Z; //!< Vector(0,0,1).
};

/*!
\brief Create a vector with the same coordinates.
\param a Real.
*/
inline Vector::Vector(double a)
{
  c[0] = c[1] = c[2] = a;
}

/*!
\brief Create a vector with argument coordinates.
\param a,b,c Coordinates.
*/
inline Vector::Vector(double a, double b, double c)
{
  Vector::c[0] = a;
  Vector::c[1] = b;
  Vector::c[2] = c;
}

//! Gets the i-th coordinate of vector.
inline double& Vector::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of vector.
inline double Vector::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vector Vector::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vector Vector::operator- () const
{
  return Vector(-c[0], -c[1], -c[2]);
}

// Assignment unary operators

//! Destructive addition.
inline Vector& Vector::operator+= (const Vector& u)
{
  c[0] += u.c[0]; c[1] += u.c[1]; c[2] += u.c[2];
  return *this;
}

//! Destructive subtraction.
inline Vector& Vector::operator-= (const Vector& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1]; c[2] -= u.c[2];
  return *this;
}

//! Destructive scalar multiply.
inline Vector& Vector::operator*= (double a)
{
  c[0] *= a; c[1] *= a; c[2] *= a;
  return *this;
}

/*!
\brief Scale a vector.
\param a Scaling vector.
*/
inline Vector Vector::Scaled(const Vector& a) const
{
  return Vector(c[0] * a[0], c[1] * a[1], c[2] * a[2]);
}

/*!
\brief Inverse of a vector.

This function inverses the components of the vector. This is the same as:
\code
Vector v=Vector(1.0/u[0],1.0/u[1],1.0/u[2]);
\endcode
*/
inline Vector Vector::Inverse() const
{
  return Vector(1.0 / c[0], 1.0 / c[1], 1.0 / c[2]);
}

//! Destructive division by a scalar.
inline Vector& Vector::operator/= (double a)
{
  c[0] /= a; c[1] /= a; c[2] /= a;
  return *this;
}

/*!
\brief Destructively scale a vector by another vector.

This is the same as Scale:
\code
Vector u(2.0,-1.0,1.0);
u=u.Scaled(Vector(3.0,1.0,2.0)); // u*=Vector(3.0,1.0,2.0);
\endcode
*/
inline Vector& Vector::operator*= (const Vector& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1]; c[2] *= u.c[2];
  return *this;
}

//! Destructively divide the components of a vector by another vector.
inline Vector& Vector::operator/= (const Vector& u)
{
  c[0] /= u.c[0]; c[1] /= u.c[1]; c[2] /= u.c[2];
  return *this;
}

//! Compare two vectors.
inline int operator> (const Vector& u, const Vector& v)
{
  return ((u.c[0] > v.c[0]) && (u.c[1] > v.c[1]) && (u.c[2] > v.c[2]));
}

//! Compare two vectors.
inline int operator< (const Vector& u, const Vector& v)
{
  return ((u.c[0] < v.c[0]) && (u.c[1] < v.c[1]) && (u.c[2] < v.c[2]));
}

//! Overloaded
inline int operator>= (const Vector& u, const Vector& v)
{
  return ((u.c[0] >= v.c[0]) && (u.c[1] >= v.c[1]) && (u.c[2] >= v.c[2]));
}

//! Overloaded
inline int operator<= (const Vector& u, const Vector& v)
{
  return ((u.c[0] <= v.c[0]) && (u.c[1] <= v.c[1]) && (u.c[2] <= v.c[2]));
}

//! Adds up two vectors.
inline Vector operator+ (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] + v.c[0], u.c[1] + v.c[1], u.c[2] + v.c[2]);
}

//! Difference between two vectors.
inline Vector operator- (const Vector& u, const Vector& v)
{
  return Vector(u.c[0] - v.c[0], u.c[1] - v.c[1], u.c[2] - v.c[2]);
}

//! Scalar product.
inline constexpr double operator* (const Vector& u, const Vector& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1] + u.c[2] * v.c[2]);
}

//! Right multiply by a scalar.
inline Vector operator* (const Vector& u, double a)
{
  return Vector(u.c[0] * a, u.c[1] * a, u.c[2] * a);
}

//! Left multiply by a scalar.
inline Vector operator* (double a, const Vector& v)
{
  return v * a;
}

//! Cross product.
inline Vector operator/ (const Vector& u, const Vector& v)
{
  return Vector(u.c[1] * v.c[2] - u.c[2] * v.c[1], u.c[2] * v.c[0] - u.c[0] * v.c[2], u.c[0] * v.c[1] - u.c[1] * v.c[0]);
}

//! Left multiply by a scalar
inline Vector operator/ (const Vector& u, double a)
{
  return Vector(u.c[0] / a, u.c[1] / a, u.c[2] / a);
}

// Boolean functions

//! Strong equality test.
inline int operator== (const Vector& u, const Vector& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]) && (u.c[2] == v.c[2]));
}

//! Strong difference test.
inline int operator!= (const Vector& u, const Vector& v)
{
  return (!(u == v));
}

/*!
\brief Compute the Euclidean norm of a vector.

This function involves a square root computation, it is in general more efficient to rely on
the squared norm of a vector instead.
\param u %Vector.
\sa SquaredNorm
*/
inline double Norm(const Vector& u)
{
  return sqrt(u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

inline double Norm2(const double& c1, const double& c2)
{
  return sqrt(c1 * c1 + c2 * c2);
}

/*!
\brief Compute the squared Euclidean norm of a vector.
\param u %Vector.
\sa Norm
*/
inline double SquaredNorm(const Vector& u)
{
  return (u.c[0] * u.c[0] + u.c[1] * u.c[1] + u.c[2] * u.c[2]);
}

inline double SquaredNorm2(const Vector& u, const Vector& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1] + u.c[2] * v.c[2]);
}
/*!
\brief Return a normalized vector.

Compute the inverse of its norm and scale the components.

This function does not check if the vector is null.
\param u %Vector.
*/
  inline Vector Normalized(const Vector& u)
  {
    return u * (1.0 / Norm(u));
  }

/*!
\brief Computes the absolute value of a vector.
\param u %Vector.
*/
inline Vector Abs(const Vector& u)
{
  return Vector(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1], u[2] > 0.0 ? u[2] : -u[2]);
}

/*!
\brief Return a vector with coordinates set to the minimum coordinates
of the two argument vectors.
*/
inline Vector Vector::Min(const Vector& a, const Vector& b)
{
  return Vector(a[0] < b[0] ? a[0] : b[0], a[1] < b[1] ? a[1] : b[1], a[2] < b[2] ? a[2] : b[2]);
}

/*!
\brief Return a vector with coordinates set to the maximum coordinates
of the two argument vectors.
*/
inline Vector Vector::Max(const Vector& a, const Vector& b)
{
  return Vector(a[0] > b[0] ? a[0] : b[0], a[1] > b[1] ? a[1] : b[1], a[2] > b[2] ? a[2] : b[2]);
}

/*!
\brief Linear interpolation between two vectors.
\param a,b Interpolated points.
\param t Interpolant.
*/
inline Vector Lerp(const Vector& a, const Vector& b, double t)
{
  return a + t * (b - a);
}

inline float Clamp(float v, float a, float b)
{
  return (v<a)?a:(v>b)?b:v;
}



/*!
\brief Bi-linear interpolation between four vectors.

The values are given in trigonometric order.

\param a00,a10,a11,a01 Interpolated vectors.
\param u,v Interpolation coefficients.

\sa Math::Bilinear
*/
inline Vector Vector::Bilinear(const Vector& a00, const Vector& a10, const Vector& a11, const Vector& a01, double u, double v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}
class Vec2
{
protected:
  double c[2]; //!< Components.
public:
  //! Empty 
  Vec2() {}

  explicit Vec2(double);
  explicit Vec2(double, double);

  // Access members
  double& operator[] (int);
  double operator[] (int) const;

  // Unary operators
  Vec2 operator+ () const;
  Vec2 operator- () const;

  // Assignment operators
  Vec2& operator+= (const Vec2&);
  Vec2& operator-= (const Vec2&);
  Vec2& operator*= (const Vec2&);
  Vec2& operator/= (const Vec2&);
  Vec2& operator*= (double);
  Vec2& operator/= (double);

  // Binary operators
  friend int operator> (const Vec2&, const Vec2&);
  friend int operator< (const Vec2&, const Vec2&);

  friend int operator>= (const Vec2&, const Vec2&);
  friend int operator<= (const Vec2&, const Vec2&);

  // Binary operators
  friend Vec2 operator+ (const Vec2&, const Vec2&);
  friend Vec2 operator- (const Vec2&, const Vec2&);

  friend constexpr double operator* (const Vec2&, const Vec2&);

  friend Vec2 operator* (const Vec2&, double);
  friend Vec2 operator* (double, const Vec2&);
  friend Vec2 operator/ (const Vec2&, double);

  // Boolean functions
  friend int operator==(const Vec2&, const Vec2&);
  friend int operator!=(const Vec2&, const Vec2&);

  // Norm
  friend double Norm(const Vec2&);
  friend double SquaredNorm(const Vec2&);
  friend double SquaredNorm2(const Vec2&, const Vec2&);

  friend void Normalize(Vec2&);
  friend Vec2 Normalized(const Vec2&);

  // Compare functions
  static Vec2 Min(const Vec2&, const Vec2&);
  static Vec2 Max(const Vec2&, const Vec2&);

  // Abs
  friend Vec2 Abs(const Vec2&);

  // Orthogonal and orthonormal Vec2s
  Vec2 Orthogonal() const;
  void Orthonormal(Vec2&, Vec2&) const;

  friend Vec2 Lerp(const Vec2&, const Vec2&, double);
  friend Vec2 Clamp(const Vec2&,double minVal, double maxVal);
  static Vec2 Bilinear(const Vec2&, const Vec2&, const Vec2&, const Vec2&, double, double);

  // Scale
  Vec2 Scaled(const Vec2&) const;
  Vec2 Inverse() const;

  friend std::ostream& operator<<(std::ostream&, const Vec2&);

public:
  static const Vec2 Null; //!< Null Vec2.
  static const Vec2 X; //!< Vec2(1,0,0).
  static const Vec2 Y; //!< Vec2(0,1,0).

};

/*!
\brief Create a vector with the same coordinates.
\param a Real.
*/
inline Vec2::Vec2(double a)
{
  c[0] = c[1] = a;
}

/*!
\brief Create a Vec2 with argument coordinates.
\param a,b,c Coordinates.
*/
inline Vec2::Vec2(double a, double b)
{
  Vec2::c[0] = a;
  Vec2::c[1] = b;
}

//! Gets the i-th coordinate of Vec2.
inline double& Vec2::operator[] (int i)
{
  return c[i];
}

//! Returns the i-th coordinate of Vec2.
inline double Vec2::operator[] (int i) const
{
  return c[i];
}

// Unary operators

//! Overloaded.
inline Vec2 Vec2::operator+ () const
{
  return *this;
}

//! Overloaded.
inline Vec2 Vec2::operator- () const
{
  return Vec2(-c[0], -c[1]);
}

// Assignment unary operators

//! Destructive addition.
inline Vec2& Vec2::operator+= (const Vec2& u)
{
  c[0] += u.c[0]; c[1] += u.c[1];
  return *this;
}

//! Destructive subtraction.
inline Vec2& Vec2::operator-= (const Vec2& u)
{
  c[0] -= u.c[0]; c[1] -= u.c[1];
  return *this;
}

//! Destructive scalar multiply.
inline Vec2& Vec2::operator*= (double a)
{
  c[0] *= a; c[1] *= a;
  return *this;
}

/*!
\brief Scale a Vec2.
\param a Scaling Vec2.
*/
inline Vec2 Vec2::Scaled(const Vec2& a) const
{
  return Vec2(c[0] * a[0], c[1] * a[1]);
}

/*!
\brief Inverse of a Vec2.

This function inverses the components of the Vec2. This is the same as:
\code
Vec2 v=Vec2(1.0/u[0],1.0/u[1],1.0/u[2]);
\endcode
*/
inline Vec2 Vec2::Inverse() const
{
  return Vec2(1.0 / c[0], 1.0 / c[1]);
}

//! Destructive division by a scalar.
inline Vec2& Vec2::operator/= (double a)
{
  c[0] /= a; c[1] /= a;
  return *this;
}

/*!
\brief Destructively scale a Vec2 by another Vec2.

This is the same as Scale:
\code
Vec2 u(2.0,-1.0,1.0);
u=u.Scaled(Vec2(3.0,1.0,2.0)); // u*=Vec2(3.0,1.0,2.0);
\endcode
*/
inline Vec2& Vec2::operator*= (const Vec2& u)
{
  c[0] *= u.c[0]; c[1] *= u.c[1];
  return *this;
}

//! Destructively divide the components of a Vec2 by another Vec2.
inline Vec2& Vec2::operator/= (const Vec2& u)
{
  c[0] /= u.c[0]; c[1] /= u.c[1];
  return *this;
}

//! Compare two Vec2s.
inline int operator> (const Vec2& u, const Vec2& v)
{
  return ((u.c[0] > v.c[0]) && (u.c[1] > v.c[1]));
}

//! Compare two Vec2s.
inline int operator< (const Vec2& u, const Vec2& v)
{
  return ((u.c[0] < v.c[0]) && (u.c[1] < v.c[1]));
}

//! Overloaded
inline int operator>= (const Vec2& u, const Vec2& v)
{
  return ((u.c[0] >= v.c[0]) && (u.c[1] >= v.c[1]));
}

//! Overloaded
inline int operator<= (const Vec2& u, const Vec2& v)
{
  return ((u.c[0] <= v.c[0]) && (u.c[1] <= v.c[1]));
}

//! Adds up two Vec2s.
inline Vec2 operator+ (const Vec2& u, const Vec2& v)
{
  return Vec2(u.c[0] + v.c[0], u.c[1] + v.c[1]);
}

//! Difference between two Vec2s.
inline Vec2 operator- (const Vec2& u, const Vec2& v)
{
  return Vec2(u.c[0] - v.c[0], u.c[1] - v.c[1]);
}

//! Scalar product.
inline constexpr double operator* (const Vec2& u, const Vec2& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1]);
}

//! Right multiply by a scalar.
inline Vec2 operator* (const Vec2& u, double a)
{
  return Vec2(u.c[0] * a, u.c[1] * a);
}

//! Left multiply by a scalar.
inline Vec2 operator* (double a, const Vec2& v)
{
  return v * a;
}

//! Left multiply by a scalar
inline Vec2 operator/ (const Vec2& u, double a)
{
  return Vec2(u.c[0] / a, u.c[1] / a);
}

// Boolean functions

//! Strong equality test.
inline int operator== (const Vec2& u, const Vec2& v)
{
  return ((u.c[0] == v.c[0]) && (u.c[1] == v.c[1]));
}

//! Strong difference test.
inline int operator!= (const Vec2& u, const Vec2& v)
{
  return (!(u == v));
}

/*!
\brief Compute the Euclidean norm of a Vec2.

This function involves a square root computation, it is in general more efficient to rely on
the squared norm of a Vec2 instead.
\param u %Vec2.
\sa SquaredNorm
*/
inline double Norm(const Vec2& u)
{
  return sqrt(u.c[0] * u.c[0] + u.c[1] * u.c[1]);
}

/*!
\brief Compute the squared Euclidean norm of a vector.
\param u %Vec2.
\sa Norm
*/
inline double SquaredNorm(const Vec2& u)
{
  return (u.c[0] * u.c[0] + u.c[1] * u.c[1]);
}

inline double SquaredNorm2(const Vec2& u, const Vec2& v)
{
  return (u.c[0] * v.c[0] + u.c[1] * v.c[1]);
}
/*!
\brief Return a normalized Vec2.

Compute the inverse of its norm and scale the components.

This function does not check if the Vec2 is null.
\param u %Vec2.
*/
inline Vec2 Normalized(const Vec2& u)
{
  return u * (1.0 / Norm(u));
}

/*!
\brief Computes the absolute value of a Vec2.
\param u %Vec2.
*/
inline Vec2 Abs(const Vec2& u)
{
  return Vec2(u[0] > 0.0 ? u[0] : -u[0], u[1] > 0.0 ? u[1] : -u[1]);
}

/*!
\brief Return a Vec2 with coordinates set to the minimum coordinates
of the two argument Vec2s.
*/
inline Vec2 Vec2::Min(const Vec2& a, const Vec2& b)
{
  return Vec2(a[0] < b[0] ? a[0] : b[0], a[1] < b[1] ? a[1] : b[1]);
}

/*!
\brief Return a Vec2 with coordinates set to the maximum coordinates
of the two argument Vec2s.
*/
inline Vec2 Vec2::Max(const Vec2& a, const Vec2& b)
{
  return Vec2(a[0] > b[0] ? a[0] : b[0], a[1] > b[1] ? a[1] : b[1]);
}

/*!
\brief Linear interpolation between two Vec2s.
\param a,b Interpolated points.
\param t Interpolant.
*/
inline Vec2 Lerp(const Vec2& a, const Vec2& b, double t)
{
  return a + t * (b - a);
}

inline Vec2 Clamp(const Vec2& a,double minVal, double maxVal)
{
  return Vec2(Math::Clamp(a[0], minVal, maxVal), Math::Clamp(a[1], minVal, maxVal));
}



/*!
\brief Bi-linear interpolation between four Vec2s.

The values are given in trigonometric order.

\param a00,a10,a11,a01 Interpolated Vec2s.
\param u,v Interpolation coefficients.

\sa Math::Bilinear
*/
inline Vec2 Vec2::Bilinear(const Vec2& a00, const Vec2& a10, const Vec2& a11, const Vec2& a01, double u, double v)
{
  return (1 - u) * (1 - v) * a00 + (1 - u) * (v)*a01 + (u) * (1 - v) * a10 + (u) * (v)*a11;
}
