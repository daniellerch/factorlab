#ifndef NTLX_svec_ZZ_p__H
#define NTLX_svec_ZZ_p__H

#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_long.h>
#include "svector.h"

NTL_OPEN_NNS;

NTL_svector_decl(ZZ_p,vec_ZZ_p,svec_ZZ_p,long,vec_long);
NTL_math_svector_decl(ZZ_p,vec_ZZ_p,svec_ZZ_p);
NTL_eq_svector_decl(ZZ_p,svec_ZZ_p);
NTL_io_svector_decl(ZZ_p,svec_ZZ_p);

NTL_CLOSE_NNS;

#endif
