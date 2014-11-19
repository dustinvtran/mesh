#include "CS207/Util.hpp"
#include <iostream>
using namespace std;

#include <vector>

/** Return true iff @a n is prime.
 * @pre @a n >= 0
 */

bool is_prime(int n)
{
  /**
  * One important implementation detail is DO NOT ADD @a n to @a primes
  * if @a n is a prime! Instead, just add the last prime that is no
  * greater than sqrt of @a n
  *
  * Invariant: @a primes is always sorted from 2,3,5,7,11,...
  *            and primes.last is always no greater than the sqrt of @a n!!
  */
  assert(n >= 0);
  assert(n != 1);

  static std::vector<int> primes;
  int root = (int)(sqrt(n));
  //std::cout << "n is: " << n <<std::endl;
  //std::cout << "root of n is:" << root << std::endl;

  /* First time enter*/
  if (primes.size() == 0)
  {
    /**
    * Find out all the primes between 2 and sqrt(n).
    * Complexity: O(N^(1/2+1/4)) = O(N^(3/4))
    */

    for (int i = 2; i <= root; ++i)
    {
      //std::cout << "first time --- i is:" << i << std::endl;
      int sqrt_i = (int)(sqrt(i));
      //std::cout << "first time --- root of i is:" << sqrt_i << std::endl;
      bool i_is_prime = true;
      for (int j = 2; j <= sqrt_i; ++j)
      {
        if (i%j == 0)
        {
          i_is_prime = false;
          break;
        }
      }
      if(i_is_prime)
      {
        //std::cout << "first time --- " << i << " is prime" << std::endl;
        primes.push_back(i);
      }
    }

    for (auto it = primes.cbegin(); (it!=primes.cend()) && ((*it) <= root); ++it)
    {
      if (n % (*it) == 0)
        return false;
    }
    //if (((primes.size() > 0) && (n > (*primes.cend()))) || (primes.size() == 0))
      //primes.push_back(n);
    return true;
  }
  /**
  * Most easy case, just scan through @a primes to see if @a n can be divided by
  * any of them. No need to add any primes at all~~~~
  * Complexity: O(size(primes))
  */
  else if (root < primes[primes.size()-1])
  {
    //std::cout << "root of n < last prime" << std::endl;
    //std::cout << "root of n < last prime --- last prime is:"  << primes[primes.size()-1]  << std::endl;
    for (auto it = primes.cbegin(); (it!=primes.cend()) && ((*it) <= root); ++it)
    {
      if (n % (*it) == 0)
        return false;
    }
    //if (n > (*primes.cend()))
      //primes.push_back(n);
    return true;
  }
  /*
  * @a primes holds its invariant currently. However, the primes between primes.last
  * and sqrt(n) are missing. So we should first add these primes back.
  */
  else
  {
    /**
    * Find out all the primes between primes.last and sqrt(n).
    * Complexity: O((N^(1/2)-size(primes))^(3/4))
    */
    for (int i = primes[primes.size()-1]+1; i <= root; ++i)
    {
      //std::cout << "root of n > last prime --- i is:" << i << std::endl;
      int sqrt_i = (int)(sqrt(i));
      //std::cout << "root of n > last prime --- root of i is:" << sqrt_i << std::endl;
      bool i_is_prime = true;
      for (int j = 2; j <= sqrt_i; ++j)
      {
        if (i%j == 0)
        {
          i_is_prime = false;
          break;
        }
      }
      if(i_is_prime)
      {
        //std::cout << "root of n > last prime --- " << i << " is prime" << std::endl;
        primes.push_back(i);
      }
    }
    for (auto it = primes.cbegin(); (it!=primes.cend()) && ((*it) <= root); ++it)
    {
      if (n % (*it) == 0)
        return false;
    }

    return true;
  }
}

/*bool is_prime(int n)
{
  assert(n >= 0);
  // special case checking
  if (n == 2) return true;
  // calculate the root of the input number
  int root = (int)(sqrt(n)) + 1;
  // should not test case: i = 1
  for (int i = root; i > 1; i--)
  {
    if (n % i == 0)
      return false;
  }
  return true;
  for (int i = n-1; i > 0; --i)
    if (n % i == 0)
      return true;
  return false;
}*/

int main()
{
  while (!std::cin.eof()) {
    // How many primes to test? And should we print them?
    std::cerr << "Input Number: ";
    int n = 0;
    CS207::getline_parsed(std::cin, n);
    if (n <= 0)
      break;

    std::cerr << "Print Primes (y/n): ";
    char confirm = 'n';
    CS207::getline_parsed(std::cin, confirm);
    bool print_primes = (confirm == 'y' || confirm == 'Y');

    CS207::Clock timer;

    // Loop and count primes from 2 up to n
    int num_primes = 0;
    for (int i = 2; i <= n; ++i) {
      if (is_prime(i)) {
        ++num_primes;
        if (print_primes)
          std::cout << i << std::endl;
      }
    }

    double elapsed_time = timer.seconds();

    std::cout << "There are " << num_primes
              << " primes less than or equal to " << n << ".\n"
              << "Found in " << (1000 * elapsed_time) << " milliseconds.\n\n";
  }

  return 0;
}
