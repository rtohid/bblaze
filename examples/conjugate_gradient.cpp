//   Copyright (c) 2019 R. Tohid
//
//   Distributed under the Boost Software License, Version 1.0. (See
//   accompanying file LICENSE_1_0.txt or copy at
//   http://www.boost.org/LICENSE_1_0.txt)

// Source: https://bitbucket.org/blaze-lib/blaze/wiki/Getting%20Started

#include <blaze/Math.h>
#include <iostream>

int main() {
  const size_t N(100UL);
  const size_t iterations(10UL);

  const size_t NN(N * N);

  blaze::CompressedMatrix<double, blaze::rowMajor> A(NN, NN);
  blaze::DynamicVector<double, blaze::columnVector> x(NN, 1.0), b(NN, 0.0),
      r(NN), p(NN), Ap(NN);
  double alpha, beta, delta;

  // ... Initializing the sparse matrix A

  // Performing the CG algorithm
  r = b - A * x;
  p = r;
  delta = (r, r);

  //   std::cout << "r = " << r << "\n";
  //   std::cout << "p = " << p << "\n";
  //   std::cout << "delta = " << delta << "\n";
  //
  //   std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<" << "\n";
  for (size_t iteration = 0UL; iteration < iterations; ++iteration) {
    Ap = A * p;
    alpha = delta / (p, Ap);
    x += alpha * p;
    r -= alpha * Ap;
    beta = (r, r);
    if (std::sqrt(beta) < 1E-8)
      break;
    p = r + (beta / delta) * p;
    delta = beta;
  }
  //   std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>" << "\n";
  //
  //   std::cout << "r = " << r << "\n";
  //   std::cout << "p = " << p << "\n";
  //   std::cout << "delta = " << delta << "\n";
}
