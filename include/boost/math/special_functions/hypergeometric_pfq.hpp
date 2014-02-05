//  Copyright (c) 2014 Anton Bikineev
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_HYPERGEOMETRIC_PFQ_HPP
#define BOOST_MATH_HYPERGEOMETRIC_PFQ_HPP

#include <boost/math/tools/promotion.hpp>
#include <boost/math/tools/series.hpp>
#include <boost/math/tools/fraction.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/detail/hypergeometric_pfq_series.hpp>

#ifdef _MSC_VER
#pragma once
#endif

namespace boost{ namespace math{ namespace detail{

//template <class T>
//struct hypergeometric_1F1_continued_fraction_term
//{
   //typedef std::pair<T,T> result_type;
   //hypergeometric_1F1_continued_fraction_term(T a_, T b_, T z_):
      //N(1), a(a_), b(b_), z(z_), numer((a_*b_)/z_),
      //term(std::make_pair(T(0), T(1)))
   //{
   //}
   //result_type operator()()
   //{
      //result_type result = term;
      //++a;
      //++b;
      //++N;
      //numer = -((a * z) / (b * N));
      //term = std::make_pair<T,T>(numer, 1-numer);
      //return result;
   //}
//private:
   //unsigned N;
   //T a, b, z;
   //T numer;
   //result_type term;
//};

//template <class T, class Policy>
//inline T hypergeometric_1F1_continued_fraction(T a, T b, T z, const Policy& pol)
//{
   //BOOST_MATH_STD_USING
   //boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
   //hypergeometric_1F1_continued_fraction_term<T> f(a, b, z);
   //if (a < 0)
      //max_iter = itrunc(-a);
   //T result = tools::continued_fraction_b(
      //f,
      //boost::math::policies::get_epsilon<T, Policy>(),
      //max_iter);
   //boost::math::policies::check_series_iterations<T>(
      //"boost::math::hypergeometric_1F1_continued_fraction<%1%>(%1%,%1%,%1%)",
      //max_iter,
      //pol);
   //result = (((a * z) / b) / result) + 1;
   //return result;
//}

//template <class T, class Policy>
//T hypergeometric_1F1_imp(T a, T b, T z, const Policy& pol)
//{
   //BOOST_MATH_STD_USING
   ////
   //// Special case: complex infinity:
   ////
   //if (floor(a) == a && floor(b) == b && a < 0 && b > a)
      //return policies::raise_domain_error<T>("boost::math::hypergeometric_1F1<%1%>(%1%,%1%,%1%)",
         //"Got z = %1%, but function is not finite for these a and b", z, pol);
   ////
   //// Special case: z -> 0:
   ////
   //if (z == 0 || abs(z) < tools::epsilon<T>())
      //return T(1);
   ////
   //// Special cases for specific values:
   ////
   //if (a == 0)
      //return T(1);
   //if (b == 0)
      //return policies::raise_domain_error<T>("boost::math::hypergeometric_1F1<%1%>(%1%,%1%,%1%)",
         //"Got z = %1%, but function is not finite for and b == 0", z, pol);
   //if (a == b)
      //return a >= 0 ? exp(z) : detail::hypergeometric_1F1_truncated_exp_series(z, itrunc(-a) + 1, pol);
   //if (a == -1)
      //return 1 - (z / b);
   //if (a == -2)
      //return 1 - ((2 * z) / b) + ((z * z) / (b * (b + 1)));
   ////
   //// Default series:
   ////
   //return detail::hypergeometric_1F1_continued_fraction(a, b, z, pol);
//}

template <class T, class Policy>
inline T hypergeometric_0f1_imp(T b, T z, const Policy& pol)
{
   // some special cases
   // ...

   return detail::hypergeometric_0f1_generic_series(b, z, pol);
}

template <class T, class Policy>
inline T hypergeometric_1f0_imp(T a, T z, const Policy& pol)
{
   // some special cases
   // ...

   return detail::hypergeometric_1f0_generic_series(a, z, pol);
}

template <class T, class Policy>
inline T hypergeometric_1f1_imp(T a, T b, T z, const Policy& pol)
{
   // some special cases
   // ...

   return detail::hypergeometric_1f1_generic_series(a, b, z, pol);
}

template <class T, class Policy>
inline T hypergeometric_1f2_imp(T a, T b1, T b2, T z, const Policy& pol)
{
   // some special cases
   // ...

   return detail::hypergeometric_1f2_generic_series(a, b1, b2, z, pol);
}

template <class T, class Policy>
inline T hypergeometric_2f1_imp(T a1, T a2, T b, T z, const Policy& pol)
{
   // some special cases
   // ...

   return detail::hypergeometric_2f1_generic_series(a1, a2, b, z, pol);
}

} // namespace detail

template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type hypergeometric_0f1(T1 b, T2 z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, Policy>(
         detail::hypergeometric_0f1_imp<value_type>(
               static_cast<value_type>(b),
               static_cast<value_type>(z),
               forwarding_policy()),
         "boost::math::hypergeometric_0f1<%1%>(%1%,%1%)");
}

template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type hypergeometric_0f1(T1 b, T2 z)
{
   return hypergeometric_0f1(b, z, policies::policy<>());
}

template <class T1, class T2, class Policy>
inline typename tools::promote_args<T1, T2>::type hypergeometric_1f0(T1 a, T2 z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, Policy>(
         detail::hypergeometric_1f0_imp<value_type>(
               static_cast<value_type>(a),
               static_cast<value_type>(z),
               forwarding_policy()),
         "boost::math::hypergeometric_1f0<%1%>(%1%,%1%)");
}

template <class T1, class T2>
inline typename tools::promote_args<T1, T2>::type hypergeometric_1f0(T1 a, T2 z)
{
   return hypergeometric_1f0(a, z, policies::policy<>());
}

template <class T1, class T2, class T3, class Policy>
inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_1f1(T1 a, T2 b, T3 z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2, T3>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, Policy>(
         detail::hypergeometric_1f1_imp<value_type>(
               static_cast<value_type>(a),
               static_cast<value_type>(b),
               static_cast<value_type>(z),
               forwarding_policy()),
         "boost::math::hypergeometric_1f1<%1%>(%1%,%1%,%1%)");
}

template <class T1, class T2, class T3>
inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_1f1(T1 a, T2 b, T3 z)
{
   return hypergeometric_1f1(a, b, z, policies::policy<>());
}

template <class T1, class T2, class T3, class T4, class Policy>
inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_1f2(T1 a, T2 b1, T3 b2, T4 z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2, T3, T4>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, Policy>(
         detail::hypergeometric_1f2_imp<value_type>(
               static_cast<value_type>(a),
               static_cast<value_type>(b1),
               static_cast<value_type>(b2),
               static_cast<value_type>(z),
               forwarding_policy()),
         "boost::math::hypergeometric_1f2<%1%>(%1%,%1%,%1%,%1%)");
}

template <class T1, class T2, class T3, class T4>
inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_1f2(T1 a, T2 b1, T3 b2, T4 z)
{
   return hypergeometric_1f2(a, b1, b2, z, policies::policy<>());
}

template <class T1, class T2, class T3, class T4, class Policy>
inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_2f1(T1 a1, T2 a2, T3 b, T4 z, const Policy& /* pol */)
{
   BOOST_FPU_EXCEPTION_GUARD
   typedef typename tools::promote_args<T1, T2, T3, T4>::type result_type;
   typedef typename policies::evaluation<result_type, Policy>::type value_type;
   typedef typename policies::normalise<
      Policy,
      policies::promote_float<false>,
      policies::promote_double<false>,
      policies::discrete_quantile<>,
      policies::assert_undefined<> >::type forwarding_policy;
   return policies::checked_narrowing_cast<result_type, Policy>(
         detail::hypergeometric_2f1_imp<value_type>(
               static_cast<value_type>(a1),
               static_cast<value_type>(a2),
               static_cast<value_type>(b),
               static_cast<value_type>(z),
               forwarding_policy()),
         "boost::math::hypergeometric_2f1<%1%>(%1%,%1%,%1%,%1%)");
}

template <class T1, class T2, class T3, class T4>
inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_2f1(T1 a1, T2 a2, T3 b, T4 z)
{
   return hypergeometric_2f1(a1, a2, b, z, policies::policy<>());
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_HYPERGEOMETRIC_PFQ_HPP
