//==================================================================
/**
 *  chebyshev.h -- C++ functions to evaluate Chebyshev polynomials
 *
 *  Copyright (C) 2014 by James A. Chappell (rlrrlrll@gmail.com)
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
  template <class T> inline T T0(const T& x)
  {
    return static_cast<T>(1.0) ;
  }

  // n = 1
  template <class T> inline T T1(const T& x)
  {
    return x ;
  }

  // n = 2
  template <class T> inline T T2(const T& x)
  {
    return (static_cast<T>(2.0) * x*x) - static_cast<T>(1.0) ;
  }

/*
 *	Tn(x)
 */
  template <class T> inline T Tn(unsigned int n, const T& x)
  {
    if (n == 0)
    {
      return T0<T>(x) ;
    }
    else if (n == 1)
    {
      return T1<T>(x) ;
    }
    else if (n == 2)
    {
      return T2<T>(x) ;
    }

/* We could simply do this:
    return (2.0 * x * Tn(n - 1, x)) - Tn(n - 2, x) ;
   but it could be slow for large n */
 
    T tnm1(T2<T>(x)) ;
    T tnm2(T1<T>(x)) ;
    T tn(tnm1) ;

    for (unsigned int l = 3 ; l <= n ; l++)
    { 
      tn = (static_cast<T>(2.0) * x * tnm1) - tnm2 ;
      tnm2 = tnm1;
      tnm1 = tn;
    }

    return tn ;
  }
}
#endif
