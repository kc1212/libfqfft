#include <cstdio>
#include <vector>

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libfqfft/polynomial_arithmetic/basic_operations.hpp>

using namespace libfqfft;

template<typename FieldT>
void _naive_multiplication(std::vector<FieldT> &c, const std::vector<FieldT> &a, const std::vector<FieldT> &b) {
  c.clear();
  c.resize(a.size() + b.size() - 1, FieldT(0, 0));

  for (unsigned i = 0; i < a.size(); i++) {
    for (unsigned j = 0; j < b.size(); j++) {
      c[i+j] += a[i]*b[j];
    }
  }
}

template<typename FieldT>
void print_poly(const std::vector<FieldT> &c) {
  for (size_t i = 0; i < c.size(); i++)
  {
    auto coefficient = c[i];

    if (i == 0) {
      std::cout << "(" <<  coefficient.c0.as_ulong() << " + " << coefficient.c1.as_ulong() << "i)" << " + ";
    }  else if (i < c.size()-1) {
      std::cout << "(" <<  coefficient.c0.as_ulong() << " + " << coefficient.c1.as_ulong() << "i)" << "x^" << i << " + ";
    } else {
      std::cout << "(" <<  coefficient.c0.as_ulong() << " + " << coefficient.c1.as_ulong() << "i)" << "x^" << i << std::endl;
    }
  }
}

/* Polynomial Multiplication via FFT */
template <typename FieldT>
void polynomial_multiplication_on_FFT_example()
{

  /* Polynomial a = 1 + 2x + 3x^2 + x^3 */
  std::vector<FieldT> a = { FieldT(1,0), FieldT(2,1), FieldT(3,0), FieldT(1,0) };

  /* Polynomial b = 1 + 2x + x^2 + x^3 */
  std::vector<FieldT> b = { FieldT(1,0), FieldT(2,0), FieldT(1,1), FieldT(1,0) };

  /*
   * c = a * b
   *   = (1 + 2x + 3x^2 + x^3) * (1 + 2x + x^2 + x^3)
   *   = 1 + 4x + 8x^2 + 10x^3 + 7x^4 + 4x^5 + x^6
   */
  std::vector<FieldT> c(1, FieldT::zero());
  _polynomial_multiplication(c, a, b);

  std::vector<FieldT> expected;
  _naive_multiplication(expected, a, b);

  /* Print out the polynomial in human-readable form */
  std::cout << "actual:" << std::endl;
  print_poly(c);

  std::cout << "expected:" << std::endl;
  print_poly(expected);
}

template <typename FieldT>
void polynomial_multiplication_on_FFT_example_random() {
  std::vector<FieldT> a = { FieldT::random_element(), FieldT::random_element(), FieldT::random_element(), FieldT::random_element() };
  std::vector<FieldT> b = { FieldT::random_element(), FieldT::random_element(), FieldT::random_element(), FieldT::random_element() };

  std::vector<FieldT> c(1, FieldT::zero());
  _polynomial_multiplication(c, a, b);

  std::vector<FieldT> expected;
  _naive_multiplication(expected, a, b);

  std::cout << "actual:" << std::endl;
  print_poly(c);

  std::cout << "expected:" << std::endl;
  print_poly(expected);
}

int main()
{
  // polynomial_multiplication_on_FFT_example<libff::Double> ();
  libff::alt_bn128_pp::init_public_params();
  polynomial_multiplication_on_FFT_example<libff::alt_bn128_pp::Fqe_type>();
polynomial_multiplication_on_FFT_example_random<libff::alt_bn128_pp::Fqe_type>();
}
