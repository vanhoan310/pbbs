#ifndef _ITEMGEN_INCLUDED
#define _ITEMGEN_INCLUDED

#include <iostream>
#include <algorithm>
#include "utils.h"

namespace dataGen {

  using namespace std;

#define HASH_MAX_INT ((unsigned) 1 << 31)

  //#define HASH_MAX_LONG ((unsigned long) 1 << 63)

  template <class T> T myhash(intT i);
  
  template <>
  intT myhash<intT>(intT i) {
    return utils::myhash(i) & (HASH_MAX_INT-1);}

  template <>
  uintT myhash<uintT>(intT i) {
    return utils::myhash(i);}

  template <>
  double myhash<double>(intT i) {
    return ((double) myhash<intT>(i)/((double) HASH_MAX_INT));}

  /* template <class T> T myhash(long i); */

  /* template <> */
  /* long myhash<long>(long i) { */
  /*   return utils::myhash(i) & (HASH_MAX_INT-1);} */
  
  /* template <> */
  /* int myhash<int>(long i) { */
  /*   return utils::myhash(i) & (HASH_MAX_INT-1);} */

  /* template <> */
  /* unsigned int myhash<unsigned int>(long i) { */
  /*   return utils::myhash(i);} */

  /* template <> */
  /* double myhash<double>(long i) { */
  /*   return ((double) myhash<long>(i)/((double) HASH_MAX_INT));} */

};

#endif // _ITEMGEN_INCLUDED
