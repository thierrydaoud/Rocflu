#define PPICLF_LRS 10
#define PPICLF_LRP 53
#define PPICLF_LRP_INT 29
#define PPICLF_LRP_PRO 31

! The viscous-unsteady history (PPICLF_VU, PPICLF_LRP3) is sized at run
! time from the input (VISCOUSUNSTEADY, VUSTEPS); see ppiclf_unsteady.f

! Particle, ghost, cell and wall-triangle counts (formerly PPICLF_LPART,
! PPICLF_LPART_GP, PPICLF_LEE, PPICLF_LWALL) are not set here: ppiclF
! sizes those arrays at run time and grows them (ppiclf_data.f)

! y, y1, ydot, ydotc: PPICLF_LRS
#define PPICLF_JX 1
#define PPICLF_JY 2
#define PPICLF_JZ 3
#define PPICLF_JVX 4
#define PPICLF_JVY 5
#define PPICLF_JVZ 6
#define PPICLF_JT 7
#define PPICLF_JOX 8
#define PPICLF_JOY 9
#define PPICLF_JOZ 10

! rprop: PPICLF_LRP, grouped by role
!--- Particle properties: set at start-up, used every step
#define PPICLF_R_JRHOP 1
#define PPICLF_R_JDP 2
#define PPICLF_R_JVOLP 3
!    superparticle loading (always 1 for now)
#define PPICLF_R_JSPL 4
!--- Model state carried from step to step
!    fluctuating quasi-steady force (QSFLUCT)
#define PPICLF_R_FLUCTFX 5
#define PPICLF_R_FLUCTFY 6
#define PPICLF_R_FLUCTFZ 7
!    pseudo-turbulence stochastic variables (PSEUDOTURB)
#define PPICLF_R_XIPAR 8
#define PPICLF_R_XIPERP 9
#define PPICLF_R_XIT 10
!    density-weighted relative acceleration, read from neighbours
!    (AMFORCE 2); must stay contiguous (x, y, z)
#define PPICLF_R_WDOTX 11
#define PPICLF_R_WDOTY 12
#define PPICLF_R_WDOTZ 13
!--- Gas at the particle: interpolated from Rocflu every stage
!    (PPICLF_LRP_INT of them, PICL_TEMP_Runge InterpFieldUser calls)
#define PPICLF_R_JRHOF 14
!    particle volume fraction: the projected particle volume,
!    interpolated back to the particle
#define PPICLF_R_JPHIP 15
#define PPICLF_R_JUX 16
#define PPICLF_R_JUY 17
#define PPICLF_R_JUZ 18
#define PPICLF_R_JCS 19
#define PPICLF_R_JT 20
#define PPICLF_R_JDPDX 21
#define PPICLF_R_JDPDY 22
#define PPICLF_R_JDPDZ 23
!    substantial derivative of the gas velocity, Du/Dt
#define PPICLF_R_JSDRX 24
#define PPICLF_R_JSDRY 25
#define PPICLF_R_JSDRZ 26
!    d(rho phi_g)/dt and grad(rho phi_g)
#define PPICLF_R_JRHSR 27
#define PPICLF_R_JPGCX 28
#define PPICLF_R_JPGCY 29
#define PPICLF_R_JPGCZ 30
!    vorticity
#define PPICLF_R_JXVOR 31
#define PPICLF_R_JYVOR 32
#define PPICLF_R_JZVOR 33
!    gas density gradient
#define PPICLF_R_JRHOGX 34
#define PPICLF_R_JRHOGY 35
#define PPICLF_R_JRHOGZ 36
!    divergence of the viscous stress
#define PPICLF_R_JDPVDX 37
#define PPICLF_R_JDPVDY 38
#define PPICLF_R_JDPVDZ 39
!    substantial derivative of the vorticity
#define PPICLF_R_JSDOX 40
#define PPICLF_R_JSDOY 41
#define PPICLF_R_JSDOZ 42
!--- Pseudo-turbulence diagnostics (PSEUDOTURB), output only
#define PPICLF_R_RTildePar 43
#define PPICLF_R_RTildePerp 44
#define PPICLF_R_R11 45
#define PPICLF_R_R22 46
#define PPICLF_R_R33 47
#define PPICLF_R_KMehra 48
#define PPICLF_R_KOsnes 49
#define PPICLF_R_RparMehra 50
#define PPICLF_R_RperpMehra 51
#define PPICLF_R_RTerm2Par 52
#define PPICLF_R_RTerm2Perp 53

! rprop slots sent with ghost particles: the ones the neighbour routine
! (ppiclf_user_EvalNearestNeighbor) reads from its neighbour (rpropj).
! A ghost carries y (PPICLF_LRS) and these slots only; the other slots
! of rpropj are zero for a ghost. Add a slot here if that routine
! starts reading it from rpropj.
#define PPICLF_LRP_GHOST 6
! Module ppiclf_data turns the list below into ppiclf_ghost_slots.
! The macro keeps a different name on purpose: cpp would replace an
! upper-case PPICLF_GHOST_SLOTS written for the array (Fortran ignores
! case) with the list.
#define PPICLF_GHOST_SLOTS_INDEX PPICLF_R_JDP, PPICLF_R_JRHOP, PPICLF_R_JVOLP, PPICLF_R_WDOTX, PPICLF_R_WDOTY, PPICLF_R_WDOTZ

! map: PPICLF_LRP_PRO
!--- Particle Volume Fraction Feedback
#define PPICLF_P_JPHIP 1
!--- x,y,z Forces Feedback
#define PPICLF_P_JFX 2
#define PPICLF_P_JFY 3
#define PPICLF_P_JFZ 4
!---Energy Feedback
#define PPICLF_P_JE 5
!--- More VF quanities. ***NEED TO CONFIRM THEY ARE USED ***
#define PPICLF_P_JPHIPD 6
#define PPICLF_P_JPHIPU 7
#define PPICLF_P_JPHIPV 8
#define PPICLF_P_JPHIPW 9
#define PPICLF_P_JPHIPT 10
!--- Reynolds Subgrid Stress Tensor
#define PPICLF_P_JRSG11 11
#define PPICLF_P_JRSG12 12
#define PPICLF_P_JRSG13 13
#define PPICLF_P_JRSG21 14
#define PPICLF_P_JRSG22 15
#define PPICLF_P_JRSG23 16
#define PPICLF_P_JRSG31 17
#define PPICLF_P_JRSG32 18
#define PPICLF_P_JRSG33 19
!--- Pseudo Turbulent Kinetic Energy
#define PPICLF_P_JTSG1 20
#define PPICLF_P_JTSG2 21
#define PPICLF_P_JTSG3 22
!--- Alpha Term Pseudo Turbulent Heat Flux
#define PPICLF_P_JAlphaPT11 23
#define PPICLF_P_JAlphaPT12 24
#define PPICLF_P_JAlphaPT13 25
#define PPICLF_P_JAlphaPT21 26
#define PPICLF_P_JAlphaPT22 27
#define PPICLF_P_JAlphaPT23 28
#define PPICLF_P_JAlphaPT31 29
#define PPICLF_P_JAlphaPT32 30
#define PPICLF_P_JAlphaPT33 31
