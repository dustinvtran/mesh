#ifndef CS207_METRICS_HPP
#define CS207_METRICS_HPP

/**
 * @file Metrics.hpp
 * @brief Template class to keep track of metric units
 */

#include "Point.hpp"

namespace Metrics {
/**
 * @class Unit
 * @brief Template class representing a metric unit
 *
 * The template parameters i, j, and k represent the exponents of the mass,
 * distance, and time of the unit, respectively.
 */
template<int i, int j, int k, typename T>
struct Unit {
  typedef T value_type;
  T value;

  // CONSTRUCTORS

  /* Unit with given value */
  constexpr Unit(const T& v = T())
      : value(v) {
  }

  /* Output stream operator */
  friend std::ostream& operator<<(std::ostream& s, const Unit& u) {
    return s << u.value << " [g^" << i << "m^" << j << "s^" << k << "]";
  }
};

// TYPEDEFS - Common units
typedef Unit<0, 0,  0, double>  Scalar;
typedef Unit<1, 0,  0, double>  Mass;
typedef Unit<0, 1,  0, double>  Distance;
typedef Unit<0, 1,  0, Point>   Distance3;
typedef Unit<0, 0,  1, double>  Time;
typedef Unit<0, 1, -1, double>  Velocity;
typedef Unit<0, 1, -1, Point>   Velocity3;
typedef Unit<0, 1, -2, double>  Acceleration;
typedef Unit<0, 1, -2, Point>   Acceleration3;
typedef Unit<1, 1, -2, double>  Force;
typedef Unit<1, 1, -2, Point>   Force3;
typedef Unit<1, 0, -2, double>  SpringConstant;
typedef Unit<1, 0, -1, double>  DampingConstant;

// USER DEFINED LITERALS -- Sugar for constructing units
constexpr Mass operator"" _g(long double v) {
  return Mass(v);
}
constexpr Time operator"" _s(long double v) {
  return Time(v);
}
constexpr Time operator"" _s(unsigned long long v) {
  return Time(v);
}
constexpr Distance operator"" _m(long double v) {
  return Distance(v);
}

// COMPARATORS
template <int i, int j, int k, typename T>
inline bool operator==(const Unit<i, j, k, T>& a,
                       const Unit<i, j, k, T>& b) {
  return a.value == b.value;
}


// ADDITION AND SUBTRACTION
/* Addition and subtraction of units should not change their type */
template <int i, int j, int k, typename T>
inline Unit<i, j, k, T>& operator+=(Unit<i, j, k, T>& a,
                                    const Unit<i, j, k, T>& b) {
  a.value += b.value;
  return a;
}

template <int i, int j, int k, typename T>
inline Unit<i, j, k, T> operator+(const Unit<i, j, k, T>& a,
                                  const Unit<i, j, k, T>& b) {
  return Unit<i, j, k, T>(a.value + b.value);
}

template <int i, int j, int k, typename T>
inline Unit<i, j, k, T>& operator-=(Unit<i, j, k, T>& a,
                                    const Unit<i, j, k, T>& b) {
  a.value -= b.value;
  return a;
}

template <int i, int j, int k, typename T>
inline Unit<i, j, k, T> operator-(const Unit<i, j, k, T>& a,
                                  const Unit<i, j, k, T>& b) {
  return Unit<i, j, k, T>(a.value - b.value);
}

// UNARY FUNCTIONS

/* Return @a -a */
template <int i, int j, int k, typename T>
inline Unit<i, j, k, T> operator-(const Unit<i, j, k, T>& a) {
  return Unit<i, j, k, T>(-a.value);
}

/* Return @a a */
template <int i, int j, int k, typename T>
inline Unit<i, j, k, T> operator+(const Unit<i, j, k, T>& a) {
  return a;
}

// MULTIPLICATION AND DIVISION

/* Only scalars can scale a Unit */
template <int i, int j, int k, typename T>
inline Unit<i, j, k, T>& operator*=(Unit<i, j, k, T>& a,
                                    const Scalar& b) {
  a.value *= b.value;
  return a;
}

/* Only scalars can scale a Unit */
template <int i, int j, int k, typename T>
inline Unit<i, j, k, T>& operator/=(Unit<i, j, k, T>& a,
                                    const Scalar& b) {
  a.value /= b.value;
  return a;
}

/* Multiplication and division of units changes the type of the result since
 * the exponents are added or subtracted, respectively.  The type of the value
 * is determined as the type of the products and divisions of those values.
 */

/** Define a helper type that is the product of two types */
template <typename T, typename U>
using prod_type = decltype(std::declval<T>() * std::declval<U>());

template <int i, int j, int k, int l, int m, int n, typename T, typename C>
inline Unit<i+l, j+m, k+n, prod_type<T,C>> operator*(const Unit<i, j, k, T>& a,
                                                     const Unit<l, m, n, C>& b) {
  return Unit<i+l, j+m, k+n, prod_type<T,C>>(a.value * b.value);
}

/** Define a helper type that is the division of two types */
template <typename T, typename U>
using div_type = decltype(std::declval<T>() / std::declval<U>());

template <int i, int j, int k, int l, int m, int n, typename T, typename C>
inline Unit<i-l, j-m, k-n, div_type<T,C>> operator/(const Unit<i, j, k, T>& a,
                                                    const Unit<l, m, n, C>& b) {
  return Unit<i-l, j-m, k-n, div_type<T,C>>(a.value / b.value);
}

};

#endif // CS207_METRICS_HPP
