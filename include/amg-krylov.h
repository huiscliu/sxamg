
#ifndef SX_AMG_KRYLOV
#define SX_AMG_KRYLOV_

#include "amg-type-def.h"

#ifdef __cplusplus
extern "C" {
#endif


SX_RTN sx_solver_gmres_itnl(SX_KRYLOV *ks, SX_AMG *mg);
SX_RTN sx_solver_cg_itnl(SX_KRYLOV *ks, SX_AMG *mg);
SX_RTN sx_solver_bicgstab_itnl(SX_KRYLOV *ks, SX_AMG *mg);

#ifdef __cplusplus
}
#endif

#endif
