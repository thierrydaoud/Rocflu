#ifndef PPICLF_LRP
#define PPICLF_LRP 1
#endif

#ifndef PPICLF_LIP
#define PPICLF_LIP 11
#endif

#ifndef PPICLF_INTERP
#define PPICLF_INTERP 1
#endif

#ifndef PPICLF_LRP_INT
#undef PPICLF_INTERP
#define PPICLF_INTERP 0
#define PPICLF_LRP_INT 1
#endif

#ifndef PPICLF_PROJECT
#define PPICLF_PROJECT 1
#endif

#ifndef PPICLF_LRP_PRO
#undef PPICLF_PROJECT
#define PPICLF_PROJECT 0
#define PPICLF_LRP_PRO 1
#endif

! Ghosts carry y and rprop, which is all the neighbour search reads.
! Projection does not use ghosts (it goes through the overlap grid).
#ifndef PPICLF_LRP_GP
#define PPICLF_LRP_GP (PPICLF_LRS+PPICLF_LRP)
#endif

#ifndef PPICLF_LIP_GP
#define PPICLF_LIP_GP 8
#endif

#ifndef PPICLF_BMAX
#define PPICLF_BMAX 1
#endif

#ifndef PPICLF_BX1
#define PPICLF_BX1 1
#endif

#ifndef PPICLF_BY1
#define PPICLF_BY1 1
#endif

#ifndef PPICLF_BZ1
#define PPICLF_BZ1 1
#endif

#ifndef PPICLF_LRMAX
#define PPICLF_LRMAX 6
#endif

