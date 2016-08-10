#include <iostream>
#include <chrono>
#include <cmath>

int main()
{
  double sum1 = 0.0, sum2 = 0.0;

  {
    auto start = std::chrono::high_resolution_clock::now();
    for (double x = 0.0; x < 1000000.0; x += 1.0)
      sum1 += std::pow(x, 0.3333333);

    std::cout << "std::pow " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "ms\n";
  }

  {
    double sum = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
    for (double x = 0.0; x < 1000000.0; x += 1.0)
      sum2 += std::cbrt(x);

    std::cout << "std::cbrt " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "ms\n";
  }

  std::cout << "Result ratio = " << sum1/sum2 << '\n';

  return 0;
}
