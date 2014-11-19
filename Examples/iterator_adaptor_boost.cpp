/** @file iterator_adaptor_boost.cpp
 * @brief A simple example of using Boost's iterator adaptors.
 * Compare these implementations with those of hand-written iterators and/or
 * hand-written iterator-adpators!
 *
 For more information and examples, see
 * http://www.boost.org/doc/libs/1_55_0/libs/iterator/doc/index.html
 * http://www.boost.org/doc/libs/1_55_0/libs/iterator/doc/iterator_adaptor.html
 */

#include <boost/iterator/iterator_adaptor.hpp>
using boost::iterator_adaptor;

#include <vector>
#include <numeric>


/** Derived class implementing a strided iterator (random access).
 */
struct my_stride_iterator
    : public iterator_adaptor<my_stride_iterator,                // Derived type
                              std::vector<int>::iterator> {      // Base type
  typedef my_stride_iterator::iterator_adaptor_ super_type;
  my_stride_iterator(const std::vector<int>::iterator& _m)
      : super_type(_m) {
  }
 private:
  friend class boost::iterator_core_access;
  void increment() {
    base_reference() += 2;
  }
  void advance(typename super_type::difference_type n) {
    base_reference() += 2*n;
  }
};

/** Derived class implementing a transformed return value (random access).
 * NOTE: Because the elements are implicit, the dereference returns
 * a value rather a reference. Thus, this is a read-only iterator!
 */
struct my_transform_iterator
    : public iterator_adaptor<my_transform_iterator,             // Derived type
                              std::vector<int>::iterator,        // Base type
                              int,                               // Value type
                              std::random_access_iterator_tag,   // Category
                              int> {                             // Reference
  typedef my_transform_iterator::iterator_adaptor_ super_type;
  my_transform_iterator(const std::vector<int>::iterator& _m)
      : super_type(_m) {
  }
 private:
  friend class boost::iterator_core_access;
  int dereference() const {
    return (*base()) * (*base());
  }
};

/** Derived class implementing an iterator with an index into container.
 * NOTE: Because the base type is an int (index), the difference_type must be
 * specified as well since the base type has no default in iterator_traits!
 */
struct my_position_iterator
    : public iterator_adaptor<my_position_iterator,              // Derived type
                              int,                               // Base type
                              int,                               // Value type
                              std::random_access_iterator_tag,   // Category
                              int&,                              // Reference
                              std::ptrdiff_t> {                  // Difference
  typedef my_position_iterator::iterator_adaptor_ super_type;
  my_position_iterator(std::vector<int>& _c, int _pos)
      : super_type(_pos), c(_c) {
  }
 private:
  friend class boost::iterator_core_access;
  int& dereference() const {
    return c[base()];
  }

  std::vector<int>& c;
};


#include <algorithm>
#include <iostream>

int main() {
  std::vector<int> a = {0,2,4,6,1,2,3,4};
  std::copy(a.begin(), a.end(), std::ostream_iterator<int>(std::cout, ", "));
  std::cout << "\b\b" << std::endl;

  auto begin = a.begin();
  auto end   = a.end();

  std::cout << "Normal sum: " << std::accumulate(begin, end, 0) << std::endl;

  auto pbegin = my_position_iterator(a, 0);
  auto pend   = my_position_iterator(a, 8);

  std::cout << "By position: " << std::accumulate(pbegin, pend, 0) << std::endl;

  auto tbegin = my_transform_iterator(a.begin());
  auto tend   = my_transform_iterator(a.end());

  std::cout << "Squared sum: " << std::accumulate(tbegin, tend, 0) << std::endl;

  auto sbegin = my_stride_iterator(a.begin());
  auto send   = my_stride_iterator(a.end());

  std::cout << "Strided sum: " << std::accumulate(sbegin, send, 0) << std::endl;

  return 0;
}
