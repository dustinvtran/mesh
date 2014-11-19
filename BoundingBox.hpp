#pragma once
/** @file BoundingBox.hpp
 * @brief Define the BoundingBox class for 3D bounding boxes. */

#include <iostream>
#include <algorithm>
#include <cmath>

#include "Point.hpp"

/** @class BoundingBox
 * @brief Class representing 3D bounding boxes.
 *
 * A BoundingBox is a 3D volume. Its fundamental operations are contains(),
 * which tests whether a point is in the volume, and operator+=(), which
 * extends the volume as necessary to ensure that the volume contains a point.
 *
 * BoundingBoxes are implemented as boxes -- 3D rectangular cuboids -- whose
 * sides are aligned with the principal axes.
 */
class BoundingBox {
 public:
  static constexpr unsigned DIM = 3;
  typedef Point point_type;

  /** Construct an empty bounding box. */
  BoundingBox()
      : empty_(true), min_(), max_() {
  }
  /** Construct the minimal bounding box containing @a p.
   * @post contains(@a p) && min() == @a p && max() == @a p */
  explicit BoundingBox(const point_type& p)
      : empty_(false), min_(p), max_(p) {
  }
  /** Construct the minimal bounding box containing @a p1 and @a p2.
   * @post contains(@a p1) && contains(@a p2) */
  BoundingBox(const point_type& p1, const point_type& p2)
      : empty_(false), min_(p1), max_(p1) {
    *this |= p2;
  }
  /** Construct a bounding box containing the points in [first, last). */
  template <typename PointIter>
  BoundingBox(PointIter first, PointIter last)
      : empty_(true), min_(), max_() {
    insert(first, last);
  }

  /** Test if the bounding box is empty (contains no points). */
  bool empty() const {
    return empty_;
  }

  /** Test if the bounding box is nonempty.
   *
   * This function lets you write code such as "if (b) { ... }" or
   * "if (box1 & box2) std::cout << "box1 and box2 intersect\n". */
  operator bool() const {
    return empty();
  }

  /** Return the minimum corner of the bounding box.
   * @post empty() || contains(min())
   * @note An empty box has min() == point_type(). */
  const point_type& min() const {
    return min_;
  }

  /** Return the maximum corner of the bounding box.
   * @post empty() || contains(max())
   * @note An empty box has max() == point_type(). */
  const point_type& max() const {
    return max_;
  }

  /** Test if point @a p is in the bounding box. */
  bool contains(const point_type& p) const {
    if (empty())
      return false;
    for (unsigned i = 0; i != DIM; ++i)
      if (p[i] < min_[i] || p[i] > max_[i])
        return false;
    return true;
  }

  /** Test if @a b is entirely within this bounding box.
   * @returns true if all @a p with @a b.contains(@a p) implies contains(@a p) */
  bool contains(const BoundingBox& b) const {
    if (empty() || b.empty())
      return false;
    for (unsigned i = 0; i != DIM; ++i)
      if (b.min_[i] < min_[i] || b.min_[i] > max_[i] ||
          b.max_[i] < min_[i] || b.max_[i] > max_[i])
        return false;
    return true;
  }

  /** Test if @a b intersects this bounding box.
   * @returns true if there exists @a p such that
   *            contains(@a p) && b.contains(@a p) */
  bool intersects(const BoundingBox& b) const {
    if (empty() || b.empty())
      return false;
    for (unsigned i = 0; i != DIM; ++i)
      if (b.min_[i] > max_[i] || b.max_[i] < min_[i])
        return false;
    return true;
  }

  /** Extend the bounding box to contain @a p.
   * @post contains(@a p) is true
   * @post For all @a x with old contains(@a x),
             then new contains(@a x) is true. */
  BoundingBox& operator|=(const point_type& p) {
    if (empty()) {
      empty_ = false;
      min_ = max_ = p;
    } else {
      for (unsigned i = 0; i != DIM; ++i) {
        if (p[i] < min_[i])  min_[i] = p[i];
        if (p[i] > max_[i])  max_[i] = p[i];
      }
    }
    return *this;
  }

  /** Extend the bounding box to contain @a b.
   * @post contains(@a b) || @a b.empty()
   * @post For all @a x with old contains(@a x) or @a b.contains(@a x),
   *         then new contains(@a x) is true. */
  BoundingBox& operator|=(const BoundingBox& b) {
    if (!b.empty())
      (*this |= b.min()) |= b.max();
    return *this;
  }

  /** Extend the bounding box to contain the points in [first, last).
   * @post For all @a p in [@a first, @a last), contains(@a p) is true.
   * @post For all @a x with old contains(@a x),
   *         then new contains(@a x) is true. */
  template <typename PointIter>
  BoundingBox& insert(PointIter first, PointIter last) {
    for ( ; first != last; ++first)
      *this |= *first;
    return *this;
  }

  /** Intersect this bounding box with another bounding box @a b.
   * @post For all @a x with old contains(@a x) and @a b.contains(@a x),
   *         then new contains(@a x) is true. */
  BoundingBox& operator&=(const BoundingBox& b) {
    if (!intersects(b))
      return clear();
    for (unsigned i = 0; i != DIM; ++i) {
      if (min_[i] < b.min_[i])  min_[i] = b.min_[i];
      if (max_[i] > b.max_[i])  max_[i] = b.max_[i];
    }
    return *this;
  }

  /** Clear the bounding box.
   * @post empty() */
  BoundingBox& clear() {
    empty_ = true;
    min_ = max_ = point_type();
    return *this;
  }

 private:
  bool empty_;
  point_type min_;
  point_type max_;
};

/** Write a BoundingBox to an output stream.
 *
 * An empty BoundingBox is written as "[]". A nonempty BoundingBox is
 * written as "[minx miny minz : maxx maxy maxz]", where all coordinates
 * are double-precision numbers.
 */
std::ostream& operator<<(std::ostream& s, const BoundingBox& box) {
  if (box.empty())
    return s << '[' << ']';
  return s << '[' << box.min() << " : " << box.max() << ']';
}

/** Return a bounding box that contains @a b and @a p. */
BoundingBox operator|(BoundingBox b, const Point& p) {
  return b |= p;
}
/** Return the union of @a b1 and @a b2. */
BoundingBox operator|(BoundingBox b1, const BoundingBox& b2) {
  return b1 |= b2;
}
/** Return a bounding box that contains @a p1 and @a p2. */
BoundingBox operator|(const Point& p1, const Point& p2) {
  return BoundingBox(p1, p2);
}

/** Return the intersection of @a b1 and @a b2. */
BoundingBox operator&(BoundingBox b1, const BoundingBox& b2) {
  return b1 &= b2;
}
