#include <iostream>
#include <utility>
#include <type_traits>
#include <iterator>

#include <vector>
#include <list>

#include "CS207/Util.hpp"

namespace cris {

// Inner namespace that users shouldn't use directly
namespace detail {

/** Generic distance for any iterator that is not random access.
 * Iterator categories inherit from weaker categories. Any iterator category
 * weaker than random access will be cast down to input_iterator_tag.
 */
template <typename Iter>
typename std::iterator_traits<Iter>::difference_type
distance(Iter first, Iter last, std::input_iterator_tag) {
  std::cout << "Generic Distance! " << std::endl;
  using size_type = typename std::iterator_traits<Iter>::difference_type;
  size_type n = 0;
  for ( ; first != last; ++first, ++n);
  return n;
}

/** Specialized distance optimized for random access iterators. */
template <typename Iter>
typename std::iterator_traits<Iter>::difference_type
distance(Iter first, Iter last, std::random_access_iterator_tag) {
  std::cout << "RAI Distance!" << std::endl;
  return {last - first};
}

} // end namespace detail

/** Publically available distance interface.
 * Dispatches to implementations depending on Iterator category.
 */
template <typename Iter>
typename std::iterator_traits<Iter>::difference_type
distance(Iter first, Iter last) {
  using trait = typename std::iterator_traits<Iter>::iterator_category;
  return detail::distance(first, last, trait());
}

}  // end namespace cris


/** A dumb iterator that wraps any iter, but only presents a forward_iterator */
template <typename Iter>
struct dumb_iterator {
  using value_type        = typename std::iterator_traits<Iter>::value_type;
  using pointer           = typename std::iterator_traits<Iter>::pointer;
  using reference         = typename std::iterator_traits<Iter>::reference;
  using difference_type   = typename std::iterator_traits<Iter>::difference_type;
  using iterator_category = std::forward_iterator_tag;

  dumb_iterator() {}
  dumb_iterator(const Iter& it) : it_(it) {}

  dumb_iterator& operator++() { it_++; return *this; }
  value_type& operator*() { return *it_; }

  bool operator==(const dumb_iterator& other) { return it_ == other.it_; }
  bool operator!=(const dumb_iterator& other) { return it_ != other.it_; }
 private:
  Iter it_;
};

// Make a dumb_iterator<Iter>
template <typename Iter>
dumb_iterator<Iter> make_dumb(const Iter& it) {
  return {it};
}



int main()
{
  unsigned N = 30000000;
  CS207::Clock timer;

  //** TEST std::list<int> **//
  {
  // A trivial list of ints
  std::list<int> v(N);

  timer.start();
  auto n = cris::distance(v.begin(), v.end());
  double time = timer.seconds();
  std::cout << "std::list<int>: " << n
            << "\n " << time << " secs\n" << std::endl;
  }

  //** TEST std::vector<int> **//
  {
  // A trivial vector of ints
  std::vector<int> v(N);

  timer.start();
  auto n = cris::distance(v.begin(), v.end());
  double time = timer.seconds();
  std::cout << "std::vector<int>: " << n
            << "\n " << time << " secs\n"<< std::endl;
  }

  //** TEST wrapped RAI iterator **//
  {
  // A trivial vector of ints
  std::vector<int> v(N);

  auto first = make_dumb(v.begin());
  auto last  = make_dumb(v.end());

  // The compiler is doing a very good job of recognizing these are pointers!!
  timer.start();
  auto n = cris::distance(first, last);
  double time = timer.seconds();
  std::cout << "std::vector<int>: " << n
            << "\n " << time << " secs\n"<< std::endl;
  }
}
