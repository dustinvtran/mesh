#include <iterator>
#include <numeric>


/** @class iterator_adaptor
 * @brief A simple adaptor class for forward iterators.
 * Classes may inherit and override appropriate operators.
 * This is primarily an example and should not be used in practice!
 */
template <typename Derived, typename Base, typename Value>
class iterator_adaptor {
  using this_type    = iterator_adaptor<Derived,Base,Value>;
  using derived_type = Derived;
  using base_type    = Base;

 public:
  // Typedefs for std::iterator_traits
  using value_type        = Value;
  using pointer           = value_type*;
  using reference         = value_type&;
  using iterator_category = std::forward_iterator_tag;
  using difference_type   = std::ptrdiff_t;

  // For convenience in derived classes
  typedef this_type                             super_type;

  // Constructors
  explicit iterator_adaptor(const base_type& b)
      : b_(b) {
  }
  iterator_adaptor(const this_type& m)
      : b_(m.base()) {
  }

  /** Overloadable forward iterator operators */

  derived_type& operator++() {                      // Prefix increment (++p)
    ++base();
    return derived();
  }
  derived_type operator++(int) {                    // Postfix increment (p++)
    derived_type tmp = derived();
    ++base();
    return tmp;
  }
  reference operator*() {                           // Dereference (*p)
    return *base();
  }
  const reference operator*() const {               // Const Dereference (*p)
    return *base();
  }
  bool operator==(const derived_type& m) const {    // Equality (p==q)
    return base() == m.base();
  }
  bool operator!=(const derived_type& m) const {    // Not Equality (p!=q)
    return !(derived() == m);
  }
 private:
  derived_type& derived() {
    return static_cast<derived_type&>(*this);
  }
  const derived_type& derived() const {
    return static_cast<const derived_type&>(*this);
  }

 protected:
  /** Access methods for the iterable base object */
  base_type& base() {
    return b_;
  }
  const base_type& base() const {
    return b_;
  }
 private:
  base_type b_;
};


#include <vector>

/** Derived class implementing a strided iterator.
 */
struct my_stride_iterator
    : public iterator_adaptor<my_stride_iterator,           // Derived type
                              std::vector<int>::iterator,   // Base type
                              int> {                        // Value type
  my_stride_iterator(const std::vector<int>::iterator& _m)
      : super_type(_m) {
  }
  my_stride_iterator& operator++() {
    base() += 2;
    return *this;
  }
};

/** Derived class implementing a transformed return value.
 */
struct my_transform_iterator
    : public iterator_adaptor<my_transform_iterator,        // Derived type
                              std::vector<int>::iterator,   // Base type
                              int> {                        // Value type
  my_transform_iterator(const std::vector<int>::iterator& _m)
      : super_type(_m) {
  }
  int operator*() const {
    return (*base()) * (*base());
  }
};

/** Derived class implementing an iterator with an index into container.
 */
struct my_position_iterator
    : public iterator_adaptor<my_position_iterator,         // Derived type
                              int,                          // Base type
                              int> {                        // Value_type
  my_position_iterator(std::vector<int>& _c, int _pos)
      : super_type(_pos), c(_c) {
  }
  const int& operator*() const {
    return c[base()];
  }
  int& operator*() {
    return c[base()];
  }
 private:
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
