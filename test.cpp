#include <iostream>
#include "chebyshev.h"

using namespace std ;
using namespace Chebyshev ;

int main()
{
  double tn ;

  cout.precision(5) ;
  for (unsigned int n = 0 ; n <= 5 ; n++)
  {
    for (double x = -1.0 ; x <= 1.0 ; x = x + 0.1)
    { 
      tn = Tn(n, x) ;
      cout << "T" << n << "(" << x << ") = " << tn << endl ;
    }
    cout << endl ;
  }

  return 0 ;
}
