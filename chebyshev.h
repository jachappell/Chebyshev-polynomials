//==================================================================
/**
 *  chebyshev.h -- C++ functions to evaluate Chebyshev polynomials
 *
 *  Copyright (C) 2019 by James A. Chappell (rlrrlrll@gmail.com)
 *
 *  Permission is hereby granted, free of charge, to any person
 *  obtaining a copy of this software and associated documentation
 *  files (the "Software"), to deal in the Software without
 *  restriction, including without limitation the rights to use,
 *  copy, modify, merge, publish, distribute, sublicense, and/or
 *  sell copies of the Software, and to permit persons to whom the
 *  Software is furnished to do so, subject to the following
 *  condition:
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *  HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *  WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *  OTHER DEALINGS IN THE SOFTWARE.
 */
//=================================================================
/*
 * chebyshev.h:  Version 0.02
 * Created by James A. Chappell <rlrrlrll@gmail.com>
 * http://www.storage-b.com/math-numerical-analysis/27
 * Created 28 July 2007
 *
 * History:
 * 28-jul-2007  created
 * 14-nov-2014  templates
 * 15-nov-2019  deduced types
 */
//==============

#ifndef __CHEBYSHEV_H__
#define __CHEBYSHEV_H__
/*
 *	Function calculates Chebyshev Polynomials Tn(x)
 */

namespace Chebyshev
{
  // n = 0
  template <class T> inline auto T0(const T& x)
  {
    return static_cast<T>(1);
  }

  // n = 1
  template <class T> inline auto T1(const T& x)
  {
    return x;
  }

  // n = 2
  template <class T> inline auto T2(const T& x)
  {
    return (static_cast<T>(2) * x*x) - static_cast<T>(1);
  }

/*
 *	Tn(x)
 */
  template <class T> inline auto Tn(unsigned int n, const T& x)
  {
    switch(n)
    {
      case 0:
        return T0<T>(x);

      case 1: 
        return T1<T>(x);

      case 2: 
        return T2<T>(x);
   
      default:
        break;
    }

/* We could simply do this:
    return (static_cast<T>(2) * x * Tn(n - 1, x)) - Tn(n - 2, x);
   but it could be slow for large n */
 
    auto tnm1(T2<T>(x));
    auto tnm2(T1<T>(x));
    auto tn(tnm1);

    for (auto l = 3u ; l <= n ; ++l)
    { 
      tn = (static_cast<T>(2) * x * tnm1) - tnm2;
      tnm2 = tnm1;
      tnm1 = tn;
    }

    return tn;
  }
}
#endif
