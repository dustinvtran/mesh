#include <iostream>
#include <utility>
#include <type_traits>
#include <cstring>   // C-style memmove

#include "CS207/Util.hpp"

namespace cris {

// Inner namespace that users shouldn't use directly
namespace detail {

/** Generic copy for any type */
template <typename InputIter, typename OutputIter, typename Trait>
OutputIter copy(InputIter first, InputIter last, OutputIter out, Trait) {
  std::cout << "Generic Copy! " << std::endl;
  for ( ; first != last; ++first, ++out)
    *out = *first;
  return out;
}

/** Specialized copy optimized for types with trivial copy
 * Selects out the case where InputIter and OutputIter are T*
 *   and T is trivially copyable.
 */
template <typename T>
T* copy(T* first, T* last, T* out, std::true_type) {
  std::cout << "Memcpy Copy!" << std::endl;
  memcpy(out, first, (last-first)*sizeof(T));
  return out + (last-first);
}

} // end namespace detail

/** Publically available copy interface.
 * Dispatches to implementations depending on Iterator type and Value type
 * @pre ...
 * @post ...
 */
template <typename InputIter, typename OutputIter>
OutputIter copy(InputIter first, InputIter last, OutputIter out) {
  using T     = typename std::iterator_traits<InputIter>::value_type;
  using trait = typename std::is_trivial<T>::type;  // true_type or false_type
  return detail::copy(first, last, out, trait());
}

}  // end namespace cris


/** A quickie vector class to show off cris::copy */
template <typename T>
class Vector {
 public:
  using value_type = T;
  using size_type = unsigned;

  Vector()
    : v_(), n_(0), cap_(0) {
  }
  explicit Vector(size_type n)
    : v_(new T[n]), n_(n), cap_(n) {
  }
  ~Vector() {
    delete[] v_;
  }
  size_type size() const {
    return n_;
  }
  void push_back(const T& x) {
    if (n_ >= cap_) {
      cap_ = (cap_ == 0 ? 16 : 2*cap_);
      T* new_v = new T[cap_];
      cris::copy(v_, v_ + n_, new_v);
      delete[] v_;
      v_ = new_v;
    }
    v_[n_] = x;
    ++n_;
  }
  const T& operator[](size_type i) const {
    return v_[i];
  }

 private:
  // RI: size(v_) == cap_
  // RI: n_ <= cap_
  T* v_;
  size_t n_;
  size_t cap_;
};



int main()
{
  unsigned N = 30000000;
  CS207::Clock timer;


  //** TEST std::vector<int> **//
  {
  // A trivial vector of ints
  std::vector<int> v(N);

  std::cout << "std::is_trivial<int> == " << std::boolalpha
            << std::is_trivial<int>() << std::endl;

  timer.start();
  v.push_back(int());
  double time = timer.seconds();
  std::cout << "std::vector<int>:\n " << time << " secs\n"<< std::endl;
  }


  //** TEST Vector<int> **//
  {
  // A trivial vector of ints
  Vector<int> v(N);

  std::cout << "std::is_trivial<int> == " << std::boolalpha
            << std::is_trivial<int>() << std::endl;

  timer.start();
  v.push_back(int());
  double time = timer.seconds();
  std::cout << "Vector<int>:\n " << time << " secs\n" << std::endl;
  }


  //** TEST Vector<my_trivial_type> **//
  {
  // A custom trivial type
  struct my_trivial_type {
    int x;
    my_trivial_type() = default;
  };
  Vector<my_trivial_type> v(N);

  std::cout << "std::is_trivial<my_trivial_type> == " << std::boolalpha
            << std::is_trivial<my_trivial_type>() << std::endl;

  timer.start();
  v.push_back(my_trivial_type());
  double time = timer.seconds();
  std::cout << "Vector<my_trivial_type>:\n " << time << " secs\n" << std::endl;
  }


  //** TEST Vector<my_nontrivial_type> **//
  {
  // A non-trivial type
  struct my_nontrivial_type {
    int x;
    my_nontrivial_type() = default;
    // Because it has custom copy/assign!
    my_nontrivial_type(const my_nontrivial_type& other) { *this = other; }
    my_nontrivial_type& operator=(const my_nontrivial_type& other) {
      std::copy(&(other.x), &(other.x) + 1, &x);
      return *this;
    }
  };
  Vector<my_nontrivial_type> v(N);

  std::cout << "std::is_trivial<my_nontrivial_type> == " << std::boolalpha
            << std::is_trivial<my_nontrivial_type>() << std::endl;

  // gcc does very poorly here, but clang does very well!!
  timer.start();
  v.push_back(my_nontrivial_type());
  double time = timer.seconds();
  std::cout << "Vector<my_nontrivial_type>:\n " << time << " secs\n" << std::endl;
  }
}
