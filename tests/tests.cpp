#include "chebyshev.h"

#include "Random.h"

#define BOOST_TEST_MODULE Legendre
#include <boost/test/included/unit_test.hpp>

using namespace Storage_B;

namespace
{
  Random<double> dRan(-1.0, 1.0);
  Random<float> fRan(-1.0f, 1.0f);
  Random<unsigned int> iRan(3, 10);

  auto tol = boost::test_tools::tolerance(0.00049);
}

// T0
BOOST_AUTO_TEST_CASE(T0)
{
  auto val1 = Chebyshev::T0<double>(dRan());
  auto val2 = Chebyshev::Tn<double>(0, dRan());

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);
}

// T1
BOOST_AUTO_TEST_CASE(T1)
{
  auto x = fRan();

  auto val1 = Chebyshev::T1<float>(x);
  auto val2 = Chebyshev::Tn<float>(1, x);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == x);
}

// Tn when x = 1
BOOST_AUTO_TEST_CASE(TnXequals1)
{
  auto val1 = Chebyshev::T0<double>(1.0);
  auto val2 = Chebyshev::Tn<double>(0, 1.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);

  val1 = Chebyshev::T1<double>(1.0);
  val2 = Chebyshev::Tn<double>(1, 1.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);

  val1 = Chebyshev::T2<double>(1.0);
  val2 = Chebyshev::Tn<double>(2, 1.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);

  auto n = iRan();

  val1 = Chebyshev::Tn<double>(n, 1.0);

  BOOST_TEST(val1 == 1.0);
}

// Tn when x = -1
BOOST_AUTO_TEST_CASE(TnXequalsMinus1)
{
  auto val1 = Chebyshev::T0<double>(-1.0);
  auto val2 = Chebyshev::Tn<double>(0, -1.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);

  val1 = Chebyshev::T1<double>(-1.0);
  val2 = Chebyshev::Tn<double>(1, -1.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == -1.0);

  val1 = Chebyshev::T2<double>(-1.0);
  val2 = Chebyshev::Tn<double>(2, -1.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);

  auto n1 = iRan();
  auto n2 = n1 + 1;

  val1 = Chebyshev::Tn<double>(n1, -1.0);
  val2 = Chebyshev::Tn<double>(n2, -1.0);

  BOOST_TEST(val1 == (n1 % 2 == 0 ? 1.0 : -1.0));
  BOOST_TEST(val2 == (n2 % 2 == 0 ? 1.0 : -1.0));
}

// Tn when x = 0
BOOST_AUTO_TEST_CASE(TnXequalsZero)
{
  auto val1 = Chebyshev::T0<double>(0.0);
  auto val2 = Chebyshev::Tn<double>(0, 0.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 1.0);

  val1 = Chebyshev::T1<double>(0.0);
  val2 = Chebyshev::Tn<double>(1, 0.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == 0.0);

  val1 = Chebyshev::T2<double>(0.0);
  val2 = Chebyshev::Tn<double>(2, 0.0);

  BOOST_TEST(val1 == val2);
  BOOST_TEST(val1 == -1.0);

  auto n1 = iRan();
  auto n2 = (2 * n1) + 1;

  val1 = Chebyshev::Tn<double>(n1 * 2, 0.0);
  val2 = Chebyshev::Tn<double>(n2, 0.0);

  BOOST_TEST(val1 == (n1 % 2 == 0 ? 1.0 : -1.0));
  BOOST_TEST(val2 == 0.0);
}

BOOST_AUTO_TEST_CASE(TnMinusPlus)
{
  auto n1 = iRan();
  auto n2 = n1 + 1;

  auto x = dRan();

  auto val1 = Chebyshev::Tn<double>(n1, x);
  auto val2 = Chebyshev::Tn<double>(n1, -x);

  BOOST_TEST(val1 == (n1 % 2 == 0 ? val2 : -val2));

  val1 = Chebyshev::Tn<double>(n2, x);
  val2 = Chebyshev::Tn<double>(n2, -x);

  BOOST_TEST(val1 == (n2 % 2 == 0 ? val2 : -val2));
}

BOOST_AUTO_TEST_CASE(Tn05)
{
  BOOST_TEST(Chebyshev::Tn<float>(2, -0.05) == -0.995, tol);
  BOOST_TEST(Chebyshev::Tn<float>(3, -0.05) == 0.1495, tol);
  BOOST_TEST(Chebyshev::Tn<float>(4, -0.05) == 0.98005, tol);
  BOOST_TEST(Chebyshev::Tn<float>(5, -0.05) == -0.247595, tol);

  BOOST_TEST(Chebyshev::Tn<float>(2, 0.05) == -0.995, tol);
  BOOST_TEST(Chebyshev::Tn<float>(3, 0.05) == -0.1495, tol);
  BOOST_TEST(Chebyshev::Tn<float>(4, 0.05) == 0.98005, tol);
  BOOST_TEST(Chebyshev::Tn<float>(5, 0.05) == 0.247595, tol);
}

BOOST_AUTO_TEST_CASE(Tn50)
{
  BOOST_TEST(Chebyshev::Tn<float>(2, -0.5) == -0.5);
  BOOST_TEST(Chebyshev::Tn<float>(3, -0.5) == 1.0);
  BOOST_TEST(Chebyshev::Tn<float>(4, -0.5) == -0.5);
  BOOST_TEST(Chebyshev::Tn<float>(5, -0.5) == -0.5);

  BOOST_TEST(Chebyshev::Tn<float>(2, 0.5) == -0.5);
  BOOST_TEST(Chebyshev::Tn<float>(3, 0.5) == -1.0);
  BOOST_TEST(Chebyshev::Tn<float>(4, 0.5) == -0.5);
  BOOST_TEST(Chebyshev::Tn<float>(5, 0.5) == 0.5);
}

BOOST_AUTO_TEST_CASE(Pn95)
{
  BOOST_TEST(Chebyshev::Tn<float>(2, -0.95) == 0.805, tol);
  BOOST_TEST(Chebyshev::Tn<float>(3, -0.95) == -0.5795, tol);
  BOOST_TEST(Chebyshev::Tn<float>(4, -0.95) == 0.2961, tol);
  BOOST_TEST(Chebyshev::Tn<float>(5, -0.95) == 0.017, tol);

  BOOST_TEST(Chebyshev::Tn<float>(2, 0.95) == 0.805, tol);
  BOOST_TEST(Chebyshev::Tn<float>(3, 0.95) == 0.5795, tol);
  BOOST_TEST(Chebyshev::Tn<float>(4, 0.95) == 0.2961, tol);
  BOOST_TEST(Chebyshev::Tn<float>(5, 0.95) == -0.017, tol);
}
