#include "chebyshev.h"

#include <iostream>
#include <cmath>

using namespace std;
using namespace Storage_B::Chebyshev;

int main()
{
  double val ;

  cout.precision(5) ;
  for (auto n = 0u ; n <= 5 ; ++n)
  {
    for (auto x = -1.0 ; x <= 1.0 ; x += 0.1)
    { 
      val = Tn<double>(n, x) ;
      cout << "T" << n << "(" << x << ") = " << val << endl ;
    }
    cout << endl ;
  }
  cout << endl ;

  for (auto n = 0u ; n <= 5 ; ++n)
  {
    for (auto x = -1.0 ; x <= 1.0 ; x += 0.1)
    { 
      val = Un<double>(n, x) ;
      cout << "U" << n << "(" << x << ") = " << val << endl ;
    }
    cout << endl ;
  }

  return 0 ;
}
