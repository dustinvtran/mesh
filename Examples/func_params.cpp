#include <functional>   // std::greater<T>
#include <vector>       // std::vector<T>
#include <iostream>     // std::cout, std::endl


/** Find the first element in a sorted array that is not less than a value.
 * @param[in] a Array to search.
 * @param[in] low,high The open range into @a a to search, [@a low, @a high).
 * @param[in] v Value to search for.
 * @return An index into array @a a or @a high.
 *
 * @tparam T    Type of the array.
 * @tparam Comp Type of the comparator object with signature:
 *                bool operator()(const T&, const T&)
 *              and returns true when the first argument is considered *less*.
 *
 * @pre 0 <= @a low <= @a high <= Size of the array @a a.
 * @pre For all i,j with @a low <= i < j < @a high, !(@a comp(@a a[j], @a a[i])).
 *
 * @post For all i with @a low <= i < result,    @a comp(@a a[i], @a v)
 *   and for all i with result <= i < @a high, !(@a comp(@a a[i], @a v)).
 *
 * Performs at most O(log(@a high - @a low)) operations.
 */
template <typename T, typename Comp = std::less<T> >
int lower_bound(const T* a, int low, int high, const T& v, Comp comp = Comp()) {
  while (low < high) {
    int mid = low + (high - low) / 2;
    if (comp(a[mid],v))
      low = mid + 1;
    else
      high = mid;
  }
  return low;
}

/** Convenience form of lower_bound to search that range [0,n)
 */
template <typename T, typename Comp = std::less<T> >
int lower_bound(const T* a, int n, const T& v, Comp comp = Comp()) {
  return lower_bound(a, 0, n, v, comp);
}


/** A simple function that returns @a a > @a b
 * Usage:
 * bool result = f_greater(4, 7);
 */
bool f_greater(int a, int b) {
  return a > b;
}

/** A function object (or functor) that can be called to return @a a > @a b
 * Usage:
 * s_greater gtr;   // Declare and construct an object
 * bool result = gtr(4, 7);
 *
 * Note that gtr(4, 7) is equivalent to gtr.operator()(4, 7)
 * This is just like (a < b) is equivalent to operator<(a, b)
 **/
struct s_greater {
  bool operator()(int a, int b) const {
    return a > b;
  }
};

int main()
{
  // Initialize a vector sorted in increasing order
  std::vector<int> v1 = {8, 18, 28, 36, 46, 65, 73, 83, 91};
  // Get the raw array
  int* v1_ptr = v1.data();

  // Default implementation uses operator< and only works with ascending arrays
  int idx = lower_bound(v1_ptr, v1.size(), 18);

  // Abstracted implementation uses a comparison function rather than op<
  // lower_bound will work with descending arrays if we use > as our comparison

  // Initialize a vector sorted in decreasing order
  std::vector<int> v2 = {91, 83, 73, 65, 56, 36, 28, 18, 8};
  // Get the raw array
  int* v2_ptr = v2.data();

  // Below are four (nearly) equivalent forms of passing a callable function
  // that returns a > b for two ints, a and b.

  // 1. Plain function parameter
  int f_idx = lower_bound(v2_ptr, v2.size(), 18,
                          f_greater);
  // 2. Function object (functor) parameter -- must construct an instance
  int s_idx = lower_bound(v2_ptr, v2.size(), 18,
                          s_greater());
  // 3. Standard library function object parameter
  int g_idx = lower_bound(v2_ptr, v2.size(), 18,
                          std::greater<int>());
  // 4. Lambda (unnamed function object) function parameter
  int l_idx = lower_bound(v2_ptr, v2.size(), 18,
                          [](int a, int b) { return a > b; });

  // Print results
  std::cout << "v1[" <<   idx << "] = " << v1[  idx] << std::endl;
  std::cout << "v2[" << f_idx << "] = " << v2[f_idx] << std::endl;
  std::cout << "v2[" << s_idx << "] = " << v2[s_idx] << std::endl;
  std::cout << "v2[" << g_idx << "] = " << v2[g_idx] << std::endl;
  std::cout << "v2[" << l_idx << "] = " << v2[l_idx] << std::endl;

  // Fancy stuff!

  // v2 is also "sorted" by the last digit of the elements (increasing)
  // Let's find one that ends in 5
  int k1 = lower_bound(v2_ptr, v2.size(), 05,
                       [](int a, int b) { return (a % 10) < (b % 10); });

  // v2 is also sorted by the first digit of the elements (decreasing)
  // Let's find one that begins in 5
  int k2 = lower_bound(v2_ptr, 0, v2.size(), 50,
                       [](int a, int b) { return (a / 10) > (b / 10); });

  // Print results
  std::cout << "v2[" << k1 << "] = " << v2[k1] << std::endl;
  std::cout << "v2[" << k2 << "] = " << v2[k2] << std::endl;

  // Note that this modularity comes at no cost.
  // We never modified or copied the data!
  // All operations are O(log(N)) and are as fast as if we had written
  //   a custom search function to accomplish each task above.

  // It even works with strings!
  std::vector<std::string> v3 = {"AlphaNumeric", "BetaGamma", "Omega", "Zeus"};
  // Get the raw array
  std::string* v3_ptr = v3.data();

  // v3 is sorted alphanumerically, we can use the default operator<
  int str_idx = lower_bound(v3_ptr, v3.size(), std::string("Omega"));

  // v3 is also sorted by string length (decreasing)
  // Construct a custom comparator
  auto comp = [](const std::string& a, const std::string& b) {
    return a.size() > b.size();
  };
  int str2_idx = lower_bound(v3_ptr, v3.size(), std::string("XXXXX"), comp);

  // Print results
  std::cout << "v3[" << str_idx  << "] = " << v3[str_idx]  << std::endl;
  std::cout << "v3[" << str2_idx << "] = " << v3[str2_idx] << std::endl;

  // Wow, lower_bound is already pretty powerful... Can it be improved?
  // Hint: Yes.
}
