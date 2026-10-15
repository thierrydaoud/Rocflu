!-----------------------------------------------------------------------
!
! Particle, ghost, overlap-grid and wall arrays of ppiclF. They used to
! be fixed-size COMMON arrays (PPICLF_LPART, PPICLF_LPART_GP, PPICLF_LEE,
! PPICLF_LWALL); they are allocated by ppiclf_data_Allocate.
!
! Capacities (number of columns allocated). They start at 1 and only
! grow (ppiclf_data_Grow*): before a transfer, a rank learns how many
! columns it will receive (ppiclf_comm_NumIncoming) and grows first,
! because gslib's crystal router drops what does not fit.
!
! Capacities (number of columns allocated):
!   ppiclf_maxpart  : real particles on this rank
!   ppiclf_maxghost : ghost particles on this rank
!   ppiclf_maxcell  : overlap-grid cells (cells sent to, or received
!                     from, the bins; cell_map, picl_grid, int_fld, ...)
!   ppiclf_maxfv    : fluid-solver cells of this rank (fluid_grid,
!                     int_fld_input, pro_fld)
!   ppiclf_maxwall  : particle-wall triangles
!
! Comes right after ppiclf_unsteady.f in ppiclf.f: the routines below
! USE it.
!
!-----------------------------------------------------------------------
      module ppiclf_data
      implicit none
      integer*4 ppiclf_maxpart, ppiclf_maxghost, ppiclf_maxcell,
     >          ppiclf_maxfv, ppiclf_maxwall
      data ppiclf_maxpart, ppiclf_maxghost, ppiclf_maxcell,
     >     ppiclf_maxfv, ppiclf_maxwall /5*0/

! Real particles
      real*8, allocatable :: ppiclf_y(:,:), ppiclf_y1(:,:),
     >   ppiclf_ydot(:,:), ppiclf_ydotc(:,:), ppiclf_rprop(:,:),
     >   ppiclf_feedbk(:,:)
      integer*4, allocatable :: ppiclf_iprop(:,:)

! Ghost particles, and the real particles' copy for the neighbour search.
! cp_map(:,i): y then all of rprop (PPICLF_LRS+PPICLF_LRP).
! rprop_gp(:,j): y then the rprop slots ppiclf_ghost_slots
! (PPICLF_LRP_GP = PPICLF_LRS+PPICLF_LRP_GHOST).
      real*8, allocatable :: ppiclf_rprop_gp(:,:), ppiclf_cp_map(:,:)
      integer*4, allocatable :: ppiclf_iprop_gp(:,:)
      integer*4, parameter :: ppiclf_ghost_slots(PPICLF_LRP_GHOST) =
     >   (/ PPICLF_GHOST_SLOTS_INDEX /)

! Overlap grid. PPICLF_PART2CELL_DIST(ip,i): squared particle-to-cell-
! center offset, each direction scaled by that cell's length,
! sum_l (dx_l/L_l)**2
      real*8, allocatable :: ppiclf_fluid_grid(:,:),
     >   ppiclf_int_fld_input(:,:), ppiclf_pro_fld(:,:)
      real*8, allocatable :: ppiclf_picl_grid(:,:),
     >   ppiclf_int_fld(:,:), ppiclf_pro_fld_picl(:,:)
      integer*4, allocatable :: ppiclf_cell_map(:,:),
     >   ppiclf_cell_map_orig(:,:), ppiclf_cell_map_interp(:,:),
     >   ppiclf_cell_map_proj(:,:)
      real*8, allocatable :: ppiclf_part2cell_dist(:,:)
      integer*4, allocatable :: ppiclf_part2cell_map(:,:),
     >   ppiclf_npart2cell(:)

! Particle walls
      real*8, allocatable :: ppiclf_wall_c(:,:), ppiclf_wall_n(:,:)

      contains
!-----------------------------------------------------------------------
! Resize an allocatable array along its particle/cell dimension to n,
! keeping the existing entries and zeroing the new ones.
!-----------------------------------------------------------------------
      subroutine ppiclf_resize_r2(a, n, idim)
      real*8, allocatable, intent(inout) :: a(:,:)
      integer*4 n, idim, m, ierr
      real*8, allocatable :: t(:,:)
      if (idim .eq. 2) then
        allocate(t(size(a,1),n), stat=ierr)
      else
        allocate(t(n,size(a,2)), stat=ierr)
      end if
      if (ierr .ne. 0) call ppiclf_exittr(
     >   'Cannot grow a ppiclF array$',0.0d0,n)
      t = 0.0d0
      m = min(size(a,idim), n)
      if (m .gt. 0) then
        if (idim .eq. 2) then
          t(:,1:m) = a(:,1:m)
        else
          t(1:m,:) = a(1:m,:)
        end if
      end if
      call move_alloc(t, a)
      end subroutine ppiclf_resize_r2
!-----------------------------------------------------------------------
      subroutine ppiclf_resize_i2(a, n, idim)
      integer*4, allocatable, intent(inout) :: a(:,:)
      integer*4 n, idim, m, ierr
      integer*4, allocatable :: t(:,:)
      if (idim .eq. 2) then
        allocate(t(size(a,1),n), stat=ierr)
      else
        allocate(t(n,size(a,2)), stat=ierr)
      end if
      if (ierr .ne. 0) call ppiclf_exittr(
     >   'Cannot grow a ppiclF array$',0.0d0,n)
      t = 0
      m = min(size(a,idim), n)
      if (m .gt. 0) then
        if (idim .eq. 2) then
          t(:,1:m) = a(:,1:m)
        else
          t(1:m,:) = a(1:m,:)
        end if
      end if
      call move_alloc(t, a)
      end subroutine ppiclf_resize_i2
!-----------------------------------------------------------------------
      subroutine ppiclf_resize_i1(a, n)
      integer*4, allocatable, intent(inout) :: a(:)
      integer*4 n, m, ierr
      integer*4, allocatable :: t(:)
      allocate(t(n), stat=ierr)
      if (ierr .ne. 0) call ppiclf_exittr(
     >   'Cannot grow a ppiclF array$',0.0d0,n)
      t = 0
      m = min(size(a), n)
      if (m .gt. 0) t(1:m) = a(1:m)
      call move_alloc(t, a)
      end subroutine ppiclf_resize_i1
!-----------------------------------------------------------------------
      subroutine ppiclf_resize_r3(a, n)
      ! along the third (particle) dimension
      real*8, allocatable, intent(inout) :: a(:,:,:)
      integer*4 n, m, ierr
      real*8, allocatable :: t(:,:,:)
      allocate(t(size(a,1),size(a,2),n), stat=ierr)
      if (ierr .ne. 0) call ppiclf_exittr(
     >   'Cannot grow a ppiclF array$',0.0d0,n)
      t = 0.0d0
      m = min(size(a,3), n)
      if (m .gt. 0) t(:,:,1:m) = a(:,:,1:m)
      call move_alloc(t, a)
      end subroutine ppiclf_resize_r3
!-----------------------------------------------------------------------
      ! New capacity: at least nneed, and at least 1.5 times the old one
      ! so that a slowly growing count reallocates rarely
      integer*4 function ppiclf_newcap(nold, nneed)
      integer*4 nold, nneed
      ppiclf_newcap = max(nneed, nold + nold/2, 1)
      end function ppiclf_newcap

      end module ppiclf_data
!-----------------------------------------------------------------------
      subroutine ppiclf_data_Allocate
!
! Allocates every array of modules ppiclf_data and ppiclf_unsteady that
! is not allocated yet, with the current capacities, and zeroes it.
!
      use ppiclf_data
      use ppiclf_unsteady
      implicit none
!
      include "PPICLF"
!
      integer*4 n, ierr
!
      ierr = 0
      if (.not. allocated(ppiclf_y)) then
        n = ppiclf_maxpart
        allocate(ppiclf_y(PPICLF_LRS,n), ppiclf_y1(PPICLF_LRS,n),
     >     ppiclf_ydot(PPICLF_LRS,n), ppiclf_ydotc(PPICLF_LRS,n),
     >     ppiclf_rprop(PPICLF_LRP,n),
     >     ppiclf_feedbk(PPICLF_LRP_PRO,n), ppiclf_iprop(PPICLF_LIP,n),
     >     ppiclf_cp_map(PPICLF_LRS+PPICLF_LRP,n),
     >     ppiclf_part2cell_dist(n,27), ppiclf_part2cell_map(n,27),
     >     ppiclf_npart2cell(n), stat=ierr)
        if (ierr .ne. 0) call ppiclf_exittr(
     >     'Cannot allocate the particle arrays$',0.0d0,n)
        ppiclf_y = 0.0d0
        ppiclf_y1 = 0.0d0
        ppiclf_ydot = 0.0d0
        ppiclf_ydotc = 0.0d0
        ppiclf_rprop = 0.0d0
        ppiclf_feedbk = 0.0d0
        ppiclf_iprop = 0
        ppiclf_cp_map = 0.0d0
        ppiclf_part2cell_dist = 0.0d0
        ppiclf_part2cell_map = 0
        ppiclf_npart2cell = 0
      end if

      if (.not. allocated(ppiclf_rprop_gp)) then
        n = ppiclf_maxghost
        allocate(ppiclf_rprop_gp(PPICLF_LRP_GP,n),
     >     ppiclf_iprop_gp(PPICLF_LIP_GP,n), stat=ierr)
        if (ierr .ne. 0) call ppiclf_exittr(
     >     'Cannot allocate the ghost arrays$',0.0d0,n)
        ppiclf_rprop_gp = 0.0d0
        ppiclf_iprop_gp = 0
      end if

      if (.not. allocated(ppiclf_fluid_grid)) then
        n = ppiclf_maxfv
        allocate(ppiclf_fluid_grid(7,n),
     >     ppiclf_int_fld_input(n,PPICLF_LRP_INT),
     >     ppiclf_pro_fld(n,PPICLF_LRP_PRO), stat=ierr)
        if (ierr .ne. 0) call ppiclf_exittr(
     >     'Cannot allocate the fluid-cell arrays$',0.0d0,n)
        ppiclf_fluid_grid = 0.0d0
        ppiclf_int_fld_input = 0.0d0
        ppiclf_pro_fld = 0.0d0
      end if

      if (.not. allocated(ppiclf_picl_grid)) then
        n = ppiclf_maxcell
        allocate(ppiclf_picl_grid(7,n),
     >     ppiclf_int_fld(PPICLF_LRP_INT,n),
     >     ppiclf_pro_fld_picl(PPICLF_LRP_PRO,n),
     >     ppiclf_cell_map(PPICLF_LRMAX,n),
     >     ppiclf_cell_map_orig(PPICLF_LRMAX,n),
     >     ppiclf_cell_map_interp(PPICLF_LRMAX,n),
     >     ppiclf_cell_map_proj(PPICLF_LRMAX,n), stat=ierr)
        if (ierr .ne. 0) call ppiclf_exittr(
     >     'Cannot allocate the overlap-grid arrays$',0.0d0,n)
        ppiclf_picl_grid = 0.0d0
        ppiclf_int_fld = 0.0d0
        ppiclf_pro_fld_picl = 0.0d0
        ppiclf_cell_map = 0
        ppiclf_cell_map_orig = 0
        ppiclf_cell_map_interp = 0
        ppiclf_cell_map_proj = 0
      end if

      if (.not. allocated(ppiclf_wall_c)) then
        n = ppiclf_maxwall
        allocate(ppiclf_wall_c(9,n), ppiclf_wall_n(4,n), stat=ierr)
        if (ierr .ne. 0) call ppiclf_exittr(
     >     'Cannot allocate the wall arrays$',0.0d0,n)
        ppiclf_wall_c = 0.0d0
        ppiclf_wall_n = 0.0d0
      end if

      ! Viscous-unsteady history: ppiclf_nUnsteadyData steps, set by
      ! Rocflu from the input (0 unless VISCOUSUNSTEADY >= 1)
      if (.not. allocated(ppiclf_rprop3)) then
        n = max(ppiclf_nUnsteadyData, 0)
        ! history packed per particle: drudt Mixt/Plag (3 each),
        ! dTdt Mixt/Plag, T Mixt/Plag
        ppiclf_lrp3 = 10*n
        allocate(ppiclf_rprop3(ppiclf_lrp3,ppiclf_maxpart),
     >     ppiclf_TimeBH(n),
     >     ppiclf_drudtPlag(3,n,ppiclf_maxpart),
     >     ppiclf_drudtMixt(3,n,ppiclf_maxpart),
     >     ppiclf_dTdtPlag(n,ppiclf_maxpart),
     >     ppiclf_dTdtMixt(n,ppiclf_maxpart),
     >     ppiclf_TMixt(n,ppiclf_maxpart),
     >     ppiclf_TPlag(n,ppiclf_maxpart), stat=ierr)
        if (ierr .ne. 0) call ppiclf_exittr(
     >     'Cannot allocate the unsteady history$',0.0d0,n)
      end if

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_data_GrowPart(nneed)
!
! Makes room for nneed real particles on this rank: the particle arrays
! of ppiclf_data and the viscous-unsteady history of ppiclf_unsteady.
!
      use ppiclf_data
      use ppiclf_unsteady
      implicit none
      integer*4 nneed, n
!
      if (nneed .le. ppiclf_maxpart) return
      n = ppiclf_newcap(ppiclf_maxpart, nneed)
      call ppiclf_resize_r2(ppiclf_y,      n, 2)
      call ppiclf_resize_r2(ppiclf_y1,     n, 2)
      call ppiclf_resize_r2(ppiclf_ydot,   n, 2)
      call ppiclf_resize_r2(ppiclf_ydotc,  n, 2)
      call ppiclf_resize_r2(ppiclf_rprop,  n, 2)
      call ppiclf_resize_r2(ppiclf_feedbk, n, 2)
      call ppiclf_resize_i2(ppiclf_iprop,  n, 2)
      call ppiclf_resize_r2(ppiclf_cp_map, n, 2)
      call ppiclf_resize_r2(ppiclf_part2cell_dist, n, 1)
      call ppiclf_resize_i2(ppiclf_part2cell_map,  n, 1)
      call ppiclf_resize_i1(ppiclf_npart2cell, n)
      if (allocated(ppiclf_rprop3)) then
        call ppiclf_resize_r2(ppiclf_rprop3,   n, 2)
        call ppiclf_resize_r3(ppiclf_drudtPlag, n)
        call ppiclf_resize_r3(ppiclf_drudtMixt, n)
        call ppiclf_resize_r2(ppiclf_dTdtPlag, n, 2)
        call ppiclf_resize_r2(ppiclf_dTdtMixt, n, 2)
        call ppiclf_resize_r2(ppiclf_TMixt,    n, 2)
        call ppiclf_resize_r2(ppiclf_TPlag,    n, 2)
      end if
      ppiclf_maxpart = n

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_data_GrowGhost(nneed)
!
! Makes room for nneed ghost particles on this rank
!
      use ppiclf_data
      implicit none
      integer*4 nneed, n
!
      if (nneed .le. ppiclf_maxghost) return
      n = ppiclf_newcap(ppiclf_maxghost, nneed)
      call ppiclf_resize_r2(ppiclf_rprop_gp, n, 2)
      call ppiclf_resize_i2(ppiclf_iprop_gp, n, 2)
      ppiclf_maxghost = n

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_data_GrowCell(nneed)
!
! Makes room for nneed overlap-grid cells on this rank
!
      use ppiclf_data
      implicit none
      integer*4 nneed, n
!
      if (nneed .le. ppiclf_maxcell) return
      n = ppiclf_newcap(ppiclf_maxcell, nneed)
      call ppiclf_resize_r2(ppiclf_picl_grid,       n, 2)
      call ppiclf_resize_r2(ppiclf_int_fld,         n, 2)
      call ppiclf_resize_r2(ppiclf_pro_fld_picl,    n, 2)
      call ppiclf_resize_i2(ppiclf_cell_map,        n, 2)
      call ppiclf_resize_i2(ppiclf_cell_map_orig,   n, 2)
      call ppiclf_resize_i2(ppiclf_cell_map_interp, n, 2)
      call ppiclf_resize_i2(ppiclf_cell_map_proj,   n, 2)
      ppiclf_maxcell = n

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_data_GrowFV(nneed)
!
! Makes room for nneed fluid-solver cells (set once, from Rocflu)
!
      use ppiclf_data
      implicit none
      integer*4 nneed
!
      if (nneed .le. ppiclf_maxfv) return
      call ppiclf_resize_r2(ppiclf_fluid_grid,    nneed, 2)
      call ppiclf_resize_r2(ppiclf_int_fld_input, nneed, 1)
      call ppiclf_resize_r2(ppiclf_pro_fld,       nneed, 1)
      ppiclf_maxfv = nneed

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_data_GrowWall(nneed)
!
! Makes room for nneed particle-wall triangles
!
      use ppiclf_data
      implicit none
      integer*4 nneed, n
!
      if (nneed .le. ppiclf_maxwall) return
      n = ppiclf_newcap(ppiclf_maxwall, nneed)
      call ppiclf_resize_r2(ppiclf_wall_c, n, 2)
      call ppiclf_resize_r2(ppiclf_wall_n, n, 2)
      ppiclf_maxwall = n

      return
      end
