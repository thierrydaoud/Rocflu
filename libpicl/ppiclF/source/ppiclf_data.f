!-----------------------------------------------------------------------
!
! Particle, ghost, overlap-grid and wall arrays of ppiclF. They used to
! be fixed-size COMMON arrays (PPICLF_LPART, PPICLF_LPART_GP, PPICLF_LEE,
! PPICLF_LWALL); they are allocated by ppiclf_data_Allocate.
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
! Must come first in ppiclf.f: the routines below USE it.
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
     >   ppiclf_rprop2(:,:), ppiclf_rprop4(:,:), ppiclf_rprop5(:,:),
     >   ppiclf_feedbk(:,:)
      integer*4, allocatable :: ppiclf_iprop(:,:)

! Ghost particles, and the real particles' copy for the neighbour search
      real*8, allocatable :: ppiclf_rprop_gp(:,:), ppiclf_cp_map(:,:)
      integer*4, allocatable :: ppiclf_iprop_gp(:,:)

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

      end module ppiclf_data
!-----------------------------------------------------------------------
      subroutine ppiclf_data_Allocate
!
! Allocates every array of module ppiclf_data that is not allocated
! yet, with the current capacities, and zeroes it.
!
      use ppiclf_data
      implicit none
!
      integer*4 n, ierr
!
      ierr = 0
      if (.not. allocated(ppiclf_y)) then
        n = ppiclf_maxpart
        allocate(ppiclf_y(PPICLF_LRS,n), ppiclf_y1(PPICLF_LRS,n),
     >     ppiclf_ydot(PPICLF_LRS,n), ppiclf_ydotc(PPICLF_LRS,n),
     >     ppiclf_rprop(PPICLF_LRP,n), ppiclf_rprop2(PPICLF_LRP2,n),
     >     ppiclf_rprop4(PPICLF_LRP4,n), ppiclf_rprop5(PPICLF_LRP5,n),
     >     ppiclf_feedbk(PPICLF_LRP_PRO,n), ppiclf_iprop(PPICLF_LIP,n),
     >     ppiclf_cp_map(PPICLF_LRP_GP,n),
     >     ppiclf_part2cell_dist(n,27), ppiclf_part2cell_map(n,27),
     >     ppiclf_npart2cell(n), stat=ierr)
        if (ierr .ne. 0) call ppiclf_exittr(
     >     'Cannot allocate the particle arrays$',0.0d0,n)
        ppiclf_y = 0.0d0
        ppiclf_y1 = 0.0d0
        ppiclf_ydot = 0.0d0
        ppiclf_ydotc = 0.0d0
        ppiclf_rprop = 0.0d0
        ppiclf_rprop2 = 0.0d0
        ppiclf_rprop4 = 0.0d0
        ppiclf_rprop5 = 0.0d0
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

      return
      end
