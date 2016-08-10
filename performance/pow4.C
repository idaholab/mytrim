#include <iostream>
#include <chrono>
#include <cmath>

//////////////////

template <int N, typename T>
struct do_pow {
  static inline T apply (const T & x)
  {
    if (N%2) // odd exponent
      return x * do_pow<N-1,T>::apply(x);

    const T xNover2 = do_pow<N/2,T>::apply(x);

    return xNover2*xNover2;
  }
};

// An efficient compiler would distill N=6 down to 3
// multiplications, but an inefficient one (or a complicated
// T::operator*) might do worse, so we'll specialize here.
template <typename T>
struct do_pow<6,T> {
  static inline T apply (const T & x)
  {
    const T x2 = x*x,
      x4 = x2*x2;

    return x4*x2;
  }
};

template <typename T>
struct do_pow<1,T> {
  static inline T apply (const T & x) { return x; }
};

template <typename T>
struct do_pow<0,T> {
  static inline T apply (const T &) { return 1; }
};


template <int N, typename T>
inline
T pow(const T & x)
{
  return do_pow<N,T>::apply(x);
}

//////////////////

int main()
{
  double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;

  {
    auto start = std::chrono::high_resolution_clock::now();
    for (double x = 0.0; x < 10000000.0; x += 1.0)
      sum1 += std::pow(x, 4.0);

    std::cout << "std::pow(x, 4.0) " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "ms\n";
  }

  {
    auto start = std::chrono::high_resolution_clock::now();
    for (double x = 0.0; x < 10000000.0; x += 1.0)
      sum2 += x * x * x * x;

    std::cout << "x * x * x * x " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "ms\n";
  }

  {
    auto start = std::chrono::high_resolution_clock::now();
    for (double x = 0.0; x < 10000000.0; x += 1.0)
      sum3 += pow<4>(x);

    std::cout << "pow<4>(x) " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "ms\n";
  }

  std::cout << "Result ratio sum1/sum2 = " << sum1/sum2 << '\n';
  std::cout << "Result ratio sum2/sum3 = " << sum2/sum3 << '\n';

  return 0;
}
