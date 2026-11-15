!-----------------------------------------------------------------------
#ifdef PPICLC
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin,
     > npart,dp_max)
     > bind(C, name="ppiclc_io_ReadParticleVTU")
#else
      subroutine ppiclf_io_ReadParticleVTU(filein1,istartoutin,
     > npart,dp_max)
#endif
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      character*1 filein1(132)
!
! Internal:
!
      ! sized by this rank's particle count
      real*4, allocatable :: rout_pos(:), rout_sln(:), rout_lrp(:),
     >                       rout_lip(:)
      real*8 dp_max
      character*1 dum_read
      character*132 filein2
      integer*8 idisp_pos,idisp_sln,idisp_lrp,idisp_lip,stride_len
      integer*4 vtu, isize, jx, jy, jz, ivtu_size, ifound, i, npt_total,
     >          npart_min, npmax, npart, ndiff, iorank, icount_pos,
     >          icount_sln, icount_lrp, icount_lip, j, ic_pos, ic_sln,
     >          ic_lrp, ic_lip, pth, ierr, indx1
      integer*4 ppiclf_indx1
      external ppiclf_indx1
      character*132 PPICLF_CHSTR
      EXTERNAL PPICLF_CHSTR
      ! Sam - this modifies the interface for ppiclC. Will now need to
      ! include NULL for optional arguments
      integer*4, optional :: istartoutin
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!
      call ppiclf_prints(' *Begin ReadParticleVTU$')
      
      call ppiclf_solve_InitZero
      PPICLF_RESTART = .true.

      if (present(istartoutin)) then
        istartout = istartoutin
      else
        istartout = 0
      end if


      indx1 = ppiclf_indx1(filein1,'.')
      indx1 = indx1 + 3 ! v (1) t (2) u (3)
      filein2 = ppiclf_chstr(filein1(1:indx1),indx1)
      if (ppiclf_nid == 0 .and. ppiclf_verbose > 0) then
         print*,'   * ParticleVTU filein2 ',trim(filein2)
      endif

      isize = 4
      jx    = 1
      jy    = 2
      jz    = 3

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=trim(filein2)
     >    ,access='stream',form="unformatted")

      ivtu_size = -1
      ifound = 0
      do i=1,1000000
      read(vtu) dum_read
      if (dum_read == '_') ifound = ifound + 1
      if (ifound .eq. 2) then
         ivtu_size = i
         exit
      endif
      enddo
      read(vtu) npt_total
      close(vtu)
      npt_total = npt_total/isize/3
      endif

      call ppiclf_bcast(npt_total,isize)
      call ppiclf_bcast(ivtu_size,isize)


      ! Each rank reads a block of ceiling(npt_total/np) particles, the
      ! last ranks fewer (as Rocflu splits points.dat); ppiclF then
      ! moves them to the ranks that own their bins
      npart_min = (npt_total + ppiclf_np - 1)/ppiclf_np
      stride_len = int(ppiclf_nid,8)*npart_min
      npart = max(0, min(npart_min, npt_total - ppiclf_nid*npart_min))
      call ppiclf_data_GrowPart(npart)
      allocate(rout_pos(3*max(npart,1)),
     >         rout_sln(ppiclf_lrs*max(npart,1)),
     >         rout_lrp(ppiclf_lrp*max(npart,1)),
     >         rout_lip(3*max(npart,1)))

      iorank = -1

      call ppiclf_byte_open_mpi(trim(filein2),pth,.true.,ierr)

      idisp_pos = ivtu_size + isize*(3*stride_len + 1)
      icount_pos = npart*3   
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_read_mpi(rout_pos,icount_pos,pth,ierr)

      do i=1,ppiclf_lrs
         idisp_sln = ivtu_size + isize*(3*npt_total 
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + i)
         j = 1 + (i-1)*npart
         icount_sln = npart

         call ppiclf_byte_set_view(idisp_sln,pth)
         call ppiclf_byte_read_mpi(rout_sln(j),icount_sln
     >                            ,pth,ierr)
      enddo

      do i=1,ppiclf_lrp
         idisp_lrp = ivtu_size + isize*(3*npt_total  
     >                         + ppiclf_lrs*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + ppiclf_lrs + i)

         j = 1 + (i-1)*npart
         icount_lrp = npart

         call ppiclf_byte_set_view(idisp_lrp,pth)
         call ppiclf_byte_read_mpi(rout_lrp(j),icount_lrp
     >                            ,pth,ierr)
      enddo

      do i=1,3
         idisp_lip = ivtu_size + isize*(3*npt_total  
     >                         + ppiclf_lrs*npt_total
     >                         + ppiclf_lrp*npt_total
     >                         + (i-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + ppiclf_lrs + ppiclf_lrp + i)

         j = 1 + (i-1)*npart
         icount_lip = npart

         call ppiclf_byte_set_view(idisp_lip,pth)
         call ppiclf_byte_read_mpi(rout_lip(j),icount_lip
     >                            ,pth,ierr)
      enddo

      call ppiclf_byte_close_mpi(pth,ierr)

      ic_pos = 0
      ic_sln = 0
      ic_lrp = 0
      ic_lip = 0
      do i=1,npart
         ic_pos = ic_pos + 1
         ppiclf_y(jx,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         ppiclf_y(jy,i) = rout_pos(ic_pos)
         ic_pos = ic_pos + 1
         if (ppiclf_ndim .eq. 3) then
         ppiclf_y(jz,i) = rout_pos(ic_pos)
         endif
      enddo
      do j=1,ppiclf_lrs
      do i=1,npart
         ic_sln = ic_sln + 1
         ppiclf_y(j,i) = rout_sln(ic_sln)
      enddo
      enddo
      do j=1,ppiclf_lrp
      do i=1,npart
         ic_lrp = ic_lrp + 1
         ppiclf_rprop(j,i) = rout_lrp(ic_lrp)
      enddo
      enddo
      ! This reads the particle tag infomation.
      do j=1,3
      do i=1,npart
         ic_lip = ic_lip + 1
         ppiclf_iprop(j,i) = int(rout_lip(ic_lip))
      enddo
      enddo

      ppiclf_npart = npart

      dp_max = MAXVAL(ppiclf_rprop(PPICLF_R_JDP,:))
      call ppiclf_printsi('  End ReadParticleVTU$',npt_total)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteBinVTU(filein1)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      character*3 filein
      character*12 vtufile
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth, nvtx_total, ncll_total
      integer*8 idisp_pos,idisp_cll,stride_lenv(8),stride_lenc
      integer*4 iint, nnp, nxx, ndxgpp1, ndygpp1, ndxygpp1, if_sz, ibin,
     >          jbin, kbin, il, ir, jl, jr, kl, kr, nbinpa, nbinpb,
     >          nbinpc, nbinpd, nbinpe, nbinpf, nbinpg, nbinph, i, j, k,
     >          ii, npa, npb, npc, npd, npe, npf, npg, nph, 
     >          ioff_dum, itype, iorank, if_cll, if_pos, icount_pos, 
     >          icount_cll, ierr, isize, ivtu_size
      real*4 rpoint(3)
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!


      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

      !ppiclf_printbinvtu set true in ppiclf_comm_CreateBin
      IF(ppiclf_printbinvtu .OR. ppiclf_time .EQ. 0.0) THEN
        call ppiclf_printsi(' *Begin WriteBinVTU$',ppiclf_cycle)
        ! Set false for next iteration
        ppiclf_printbinvtu = .FALSE. 
      ELSE
        RETURN
      END IF

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART

      nvtx_total = (ppiclf_n_bins(1)+1)*(ppiclf_n_bins(2)+1)
      if (ppiclf_ndim .gt. 2) 
     >    nvtx_total = nvtx_total*(ppiclf_n_bins(3)+1)
      ncll_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)
      if (ppiclf_ndim .gt. 2) ncll_total = ncll_total*ppiclf_n_bins(3)

      ndxgpp1 = ppiclf_n_bins(1) + 1
      ndygpp1 = ppiclf_n_bins(2) + 1
      ndxygpp1 = ndxgpp1*ndygpp1

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'bin'
      else 
         write(filein,'(A3)') filein1
      endif

      isize  = 4

      ! get which bin this processor holds
      ibin = modulo(ppiclf_nid,ppiclf_n_bins(1))
      jbin = modulo(ppiclf_nid/ppiclf_n_bins(1),ppiclf_n_bins(2))
      kbin = 0
      if (ppiclf_ndim .eq. 3)
     >kbin = ppiclf_nid/(ppiclf_n_bins(1)*ppiclf_n_bins(2))

      il = ibin
      ir = ibin+1
      jl = jbin
      jr = jbin+1
      kl = kbin
      kr = kbin
      if (ppiclf_ndim .eq. 3) then
         kl = kbin
         kr = kbin+1
      endif

      nbinpa = il + ndxgpp1*jl + ndxygpp1*kl
      nbinpb = ir + ndxgpp1*jl + ndxygpp1*kl
      nbinpc = il + ndxgpp1*jr + ndxygpp1*kl
      nbinpd = ir + ndxgpp1*jr + ndxygpp1*kl
      if (ppiclf_ndim .eq. 3) then
         nbinpe = il + ndxgpp1*jl + ndxygpp1*kr
         nbinpf = ir + ndxgpp1*jl + ndxygpp1*kr
         nbinpg = il + ndxgpp1*jr + ndxygpp1*kr
         nbinph = ir + ndxgpp1*jr + ndxygpp1*kr
      endif

      stride_lenv(1) = 0
      stride_lenv(2) = 0
      stride_lenv(3) = 0
      stride_lenv(4) = 0
      stride_lenv(5) = 0
      stride_lenv(6) = 0
      stride_lenv(7) = 0
      stride_lenv(8) = 0
 
      stride_lenc = 0
      if (ppiclf_nid .le. 
     >      ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         stride_lenv(1) = nbinpa
         stride_lenv(2) = nbinpb
         stride_lenv(3) = nbinpc
         stride_lenv(4) = nbinpd
         if (ppiclf_ndim .eq. 3) then
            stride_lenv(5) = nbinpe
            stride_lenv(6) = nbinpf
            stride_lenv(7) = nbinpg
            stride_lenv(8) = nbinph
         endif
         
         stride_lenc    = ppiclf_nid
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') nvtx_total
      write(vtu,'(A)',advance='no') '" NumberOfCells="'
      write(vtu,'(I0)',advance='no') ncll_total
      write(vtu,'(A)',advance='yes') '"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*nvtx_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'
      write(vtu,'(A)',advance='yes') '   </PointData> '
      write(vtu,'(A)',advance='yes') '   <CellData>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"PPR",1,iint)
      iint = iint + 1   *isize*ncll_total + isize
      write(vtu,'(A)',advance='yes') '   </CellData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write connectivity here
      do ii=0,ncll_total-1
         i = modulo(ii,ppiclf_n_bins(1))
         j = modulo(ii/ppiclf_n_bins(1),ppiclf_n_bins(2))
         k = 0
         if (ppiclf_ndim .eq. 3)
     >   k = ii/(ppiclf_n_bins(1)*ppiclf_n_bins(2))
          
c     do K=0,ppiclf_n_bins(3)-1
         kl = K
         kr = K
         if (ppiclf_ndim .eq. 3) then
            kl = K
            kr = K+1
         endif
c     do J=0,ppiclf_n_bins(2)-1
         jl = J
         jr = J+1
c     do I=0,ppiclf_n_bins(1)-1
         il = I
         ir = I+1

         if (ppiclf_ndim .eq. 3) then
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            npe = il + ndxgpp1*jl + ndxygpp1*kr
            npf = ir + ndxgpp1*jl + ndxygpp1*kr
            npg = il + ndxgpp1*jr + ndxygpp1*kr
            nph = ir + ndxgpp1*jr + ndxygpp1*kr
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0,A)',advance='no')  npd, ' '
            write(vtu,'(I0,A)',advance='no')  npe, ' '
            write(vtu,'(I0,A)',advance='no')  npf, ' '
            write(vtu,'(I0,A)',advance='no')  npg, ' '
            write(vtu,'(I0)'  ,advance='yes') nph
         else
            npa = il + ndxgpp1*jl + ndxygpp1*kl
            npb = ir + ndxgpp1*jl + ndxygpp1*kl
            npc = il + ndxgpp1*jr + ndxygpp1*kl
            npd = ir + ndxgpp1*jr + ndxygpp1*kl
            write(vtu,'(I0,A)',advance='no')  npa, ' '
            write(vtu,'(I0,A)',advance='no')  npb, ' '
            write(vtu,'(I0,A)',advance='no')  npc, ' '
            write(vtu,'(I0)'  ,advance='yes') npd
         endif
c     enddo
c     enddo
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      ! write offsetts here
      ioff_dum = 4
      if (ppiclf_ndim .eq. 3) ioff_dum = 8
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') ioff_dum*i
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="UInt8" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"> '
      itype = 8
      if (ppiclf_ndim .eq. 3) itype = 11
      ! write types here
      do i=1,ncll_total
         write(vtu,'(I0)',advance='yes') itype
      enddo
      write(vtu,'(A)',advance='yes')  '    </DataArray>'

      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

c1511 continue

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call ppiclf_bcast(ivtu_size, isize)

      iorank = -1

      if_pos = 3*isize*nvtx_total


      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif


      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      ! point A
      icount_pos = 0
      if (ppiclf_nid .le. 
     >    ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         icount_pos = 3
      endif
      idisp_pos  = ivtu_size + isize*(3*stride_lenv(1) + 1)
      rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
      rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
      rpoint(3)  = 0.0
      if (ppiclf_ndim .eq. 3)
     >rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      ! 3d
      if (ppiclf_ndim .gt. 2) then

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point E
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(5) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(1,3))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point F
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(6) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point G
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(7) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point H
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3   *stride_lenv(8) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
         if (kbin .eq. ppiclf_n_bins(3)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = sngl(ppiclf_bin_pos(2,3))
         endif
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)


      ! 2d
      else

         ! point B
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(2) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(1,2))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point C
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(3) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(1,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = 0.0
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)
         
         ! point D
         icount_pos = 0
         idisp_pos  = ivtu_size + isize*(3*stride_lenv(4) + 1)
         if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         if (ibin .eq. ppiclf_n_bins(1)-1) then
         if (jbin .eq. ppiclf_n_bins(2)-1) then
            icount_pos = 3
            rpoint(1)  = sngl(ppiclf_bin_pos(2,1))
            rpoint(2)  = sngl(ppiclf_bin_pos(2,2))
            rpoint(3)  = 0.0
         endif
         endif
         endif
         call ppiclf_byte_set_view(idisp_pos,pth)
         call ppiclf_byte_write_mpi(rpoint,icount_pos,iorank,pth,ierr)

      endif

      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      if_cll = 1*isize*ncll_total

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_cll
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write points first
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)

      !
      ! cell values
      idisp_cll = ivtu_size + isize*(3*(nvtx_total) 
     >     + 1*stride_lenc + 2)
      icount_cll = 0
      if (ppiclf_nid .le. 
     >       ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)-1) then
         icount_cll = 1
      endif
      rpoint(1)  = real(nxx)
      call ppiclf_byte_set_view(idisp_cll,pth)
      call ppiclf_byte_write_mpi(rpoint,icount_cll,iorank,pth,ierr)

      call ppiclf_byte_close_mpi(pth,ierr)

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi('  End WriteBinVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteParticleVTU(filein1)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) filein1
!
! Internal:
!
      ! sized by this rank's particle count
      real*4, allocatable :: rout_pos(:), rout_fld(:)
      character*3 filein
      character*12 vtufile
      integer*4 icalld1
      save      icalld1
      data      icalld1 /0/
      integer*4 vtu,pth,prevs(2,ppiclf_np)
      integer*8 idisp_pos,idisp_fld,stride_len
      integer*4 iint, nnp, nxx, npt_total, jx, jy, jz, if_sz, isize,
     >          if_pos, if_fld, ic_pos, ic_fld, i, j, k, nps, nglob,
     >          nkey, ndum, icount_pos, icount_fld, iorank, ierr,
     >          ivtu_size
      integer*4 ppiclf_iglsum
      external ppiclf_iglsum
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!

      call ppiclf_printsi(' *Begin WriteParticleVTU$',ppiclf_cycle)

      ! the fields to write (ppiclf_user_OutputFields), built once
      if (ppiclf_nout .lt. 0) then
         ppiclf_nout = 0
         call ppiclf_user_OutputFields
      endif

      if (icalld1 .eq. 0) icalld1 = istartout

      icalld1 = icalld1+1

      nnp   = ppiclf_np
      nxx   = PPICLF_NPART
      allocate(rout_pos(3*max(nxx,1)),
     >         rout_fld(max(ppiclf_nout,1)*max(nxx,1)))

      npt_total = ppiclf_iglsum(nxx,1)

      jx    = 1
      jy    = 2
      jz    = 3

      if_sz = len(filein1)
      if (if_sz .lt. 3) then
         filein = 'par'
      else 
         write(filein,'(A3)') filein1
      endif

! --------------------------------------------------
! COPY PARTICLES TO OUTPUT ARRAY
! --------------------------------------------------

      isize = 4

      if_pos = 3*isize*npt_total
      if_fld = 1*isize*npt_total

      ic_pos = 0
      do i=1,nxx

         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jx,i))
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jy,i))
         ic_pos = ic_pos + 1
         rout_pos(ic_pos) = sngl(ppiclf_y(jz,i))
      enddo
      ! one block of nxx values per output field
      ic_fld = 0
      do k=1,ppiclf_nout
         j = ppiclf_out_slot(k)
         do i=1,nxx
            ic_fld = ic_fld + 1
            if (ppiclf_out_source(k) .eq. PPICLF_OUT_Y) then
               rout_fld(ic_fld) = sngl(ppiclf_y(j,i))
            elseif (ppiclf_out_source(k) .eq. PPICLF_OUT_RPROP) then
               rout_fld(ic_fld) = sngl(ppiclf_rprop(j,i))
            else
               rout_fld(ic_fld) = real(ppiclf_iprop(j,i))
            endif
         enddo
      enddo

! --------------------------------------------------
! FIRST GET HOW MANY PARTICLES WERE BEFORE THIS RANK
! --------------------------------------------------
      do i=1,nnp
         prevs(1,i) = i-1
         prevs(2,i) = nxx
      enddo

      nps   = 1 ! index of new proc for doing stuff
      nglob = 1 ! unique key to sort by
      nkey  = 1 ! number of keys (just 1 here)
      ndum = 2
      call pfgslib_crystal_ituple_transfer(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nnp,nps)
      call pfgslib_crystal_ituple_sort(ppiclf_cr_hndl,prevs,
     >                 ndum,nnp,nglob,nkey)

      stride_len = 0
      if (ppiclf_nid .ne. 0) then
      do i=1,ppiclf_nid
         stride_len = stride_len + prevs(2,i)
      enddo
      endif

! ----------------------------------------------------
! WRITE EACH INDIVIDUAL COMPONENT OF A BINARY VTU FILE
! ----------------------------------------------------
      write(vtufile,'(A3,I5.5,A4)') filein,icalld1,'.vtu'

      if (ppiclf_nid .eq. 0) then

      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='replace')

! ------------
! FRONT MATTER
! ------------
      write(vtu,'(A)',advance='no') '<VTKFile '
      write(vtu,'(A)',advance='no') 'type="UnstructuredGrid" '
      write(vtu,'(A)',advance='no') 'version="1.0" '
      if (ppiclf_iendian .eq. 0) then
         write(vtu,'(A)',advance='yes') 'byte_order="LittleEndian">'
      elseif (ppiclf_iendian .eq. 1) then
         write(vtu,'(A)',advance='yes') 'byte_order="BigEndian">'
      endif

      write(vtu,'(A)',advance='yes') ' <UnstructuredGrid>'

      write(vtu,'(A)',advance='yes') '  <FieldData>' 
      write(vtu,'(A)',advance='no')  '   <DataArray '  ! time
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="TIME" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(E14.7)',advance='no') ppiclf_time
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='no') '   <DataArray '  ! cycle
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="CYCLE" '
      write(vtu,'(A)',advance='no') 'NumberOfTuples="1" '
      write(vtu,'(A)',advance='no') 'format="ascii"> '
      write(vtu,'(I0)',advance='no') ppiclf_cycle
      write(vtu,'(A)',advance='yes') ' </DataArray> '

      write(vtu,'(A)',advance='yes') '  </FieldData>'
      write(vtu,'(A)',advance='no') '  <Piece '
      write(vtu,'(A)',advance='no') 'NumberOfPoints="'
      write(vtu,'(I0)',advance='no') npt_total
      write(vtu,'(A)',advance='yes') '" NumberOfCells="0"> '

! -----------
! COORDINATES 
! -----------
      iint = 0
      write(vtu,'(A)',advance='yes') '   <Points>'
      call ppiclf_io_WriteDataArrayVTU(vtu,"Position",3,iint)
      iint = iint + 3*isize*npt_total + isize
      write(vtu,'(A)',advance='yes') '   </Points>'

! ----
! DATA 
! ----
      write(vtu,'(A)',advance='yes') '   <PointData>'


      do k=1,ppiclf_nout
         call ppiclf_io_WriteDataArrayVTU(vtu,trim(ppiclf_out_name(k))
     >                                   ,1,iint)
         iint = iint + 1*isize*npt_total + isize
      enddo

      write(vtu,'(A)',advance='yes') '   </PointData> '

! ----------
! END MATTER
! ----------
      write(vtu,'(A)',advance='yes') '   <Cells> '
      write(vtu,'(A)',advance='no')  '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="connectivity" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="offsets" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Int32" '
      write(vtu,'(A)',advance='no') 'Name="types" '
      write(vtu,'(A)',advance='yes') 'format="ascii"/> '
      write(vtu,'(A)',advance='yes') '   </Cells> '
      write(vtu,'(A)',advance='yes') '  </Piece> '
      write(vtu,'(A)',advance='yes') ' </UnstructuredGrid> '

! -----------
! APPEND DATA  
! -----------
      write(vtu,'(A)',advance='no') ' <AppendedData encoding="raw">'
      close(vtu)

      open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >    ,position='append')
      write(vtu) '_'
      close(vtu)

      inquire(file=vtufile,size=ivtu_size)
      endif

      call ppiclf_bcast(ivtu_size, isize)

      ! byte-displacements
      idisp_pos = ivtu_size + isize*(3*stride_len + 1)

      ! how much to write
      icount_pos = 3*nxx
      icount_fld = 1*nxx

      iorank = -1

      ! integer write
      if (ppiclf_nid .eq. 0) then
        open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >      ,position='append')
        write(vtu) if_pos
        close(vtu)
      endif

      call mpi_barrier(ppiclf_comm,ierr)

      ! write
      call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
      call ppiclf_byte_set_view(idisp_pos,pth)
      call ppiclf_byte_write_mpi(rout_pos,icount_pos,iorank,pth,ierr)
      call ppiclf_byte_close_mpi(pth,ierr)

      call mpi_barrier(ppiclf_comm,ierr)

      do k=1,ppiclf_nout
         idisp_fld = ivtu_size + isize*(3*npt_total
     >                         + (k-1)*npt_total
     >                         + (1)*stride_len
     >                         + 1 + k)

         ! integer write
         if (ppiclf_nid .eq. 0) then
           open(unit=vtu,file=vtufile,access='stream',form="unformatted"
     >         ,position='append')
           write(vtu) if_fld
           close(vtu)
         endif

         call mpi_barrier(ppiclf_comm,ierr)

         j = (k-1)*nxx + 1

         ! write
         call ppiclf_byte_open_mpi(vtufile,pth,.false.,ierr)
         call ppiclf_byte_set_view(idisp_fld,pth)
         call ppiclf_byte_write_mpi(rout_fld(j),icount_fld,iorank,pth
     >                             ,ierr)
         call ppiclf_byte_close_mpi(pth,ierr)
      enddo

      if (ppiclf_nid .eq. 0) then
      vtu=867+ppiclf_nid
      open(unit=vtu,file=vtufile,status='old',position='append')

      write(vtu,'(A)',advance='yes') '</AppendedData>'
      write(vtu,'(A)',advance='yes') '</VTKFile>'

      close(vtu)
      endif

      call ppiclf_printsi(' *End WriteParticleVTU$',ppiclf_cycle)

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_AddOutputField(isource,islot,fname)
!
! Appends one field to the particle output (par*.vtu): slot islot of
! ppiclf_y (isource = PPICLF_OUT_Y), ppiclf_rprop (PPICLF_OUT_RPROP)
! or ppiclf_iprop (PPICLF_OUT_TAG), written under the name fname.
! Called from ppiclf_user_OutputFields. Slot 0 of y or rprop is a slot
! this build does not have (PPICLF_LAYOUT.h): the field is left out.
!
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Input:
!
      integer*4 isource, islot
      character (len = *) fname
!
! Internal:
!
      integer*4 k, nslot
!
      if (islot .eq. 0 .and. isource .ne. PPICLF_OUT_TAG) return
      if (ppiclf_nout .ge. ppiclf_max_out) call ppiclf_exittr(
     >   'AddOutputField: too many output fields$',0.0d0,ppiclf_nout)

      if (isource .eq. PPICLF_OUT_Y) then
         nslot = ppiclf_lrs
      elseif (isource .eq. PPICLF_OUT_RPROP) then
         nslot = ppiclf_lrp
      elseif (isource .eq. PPICLF_OUT_TAG) then
         nslot = PPICLF_LIP
      else
         call ppiclf_exittr('AddOutputField: unknown source$',0.0d0,
     >                      isource)
      endif
      if (islot .lt. 1 .or. islot .gt. nslot) call ppiclf_exittr(
     >   'AddOutputField: slot out of range$',0.0d0,islot)

      if (len_trim(fname) .lt. 1 .or. len_trim(fname) .gt.
     >    len(ppiclf_out_name(1))) call ppiclf_exittr(
     >   'AddOutputField: empty or too long name$',0.0d0,islot)
      do k=1,ppiclf_nout
         if (ppiclf_out_name(k) .eq. fname) call ppiclf_exittr(
     >      'AddOutputField: name used twice$',0.0d0,k)
      enddo

      ppiclf_nout = ppiclf_nout + 1
      ppiclf_out_source(ppiclf_nout) = isource
      ppiclf_out_slot(ppiclf_nout)   = islot
      ppiclf_out_name(ppiclf_nout)   = fname

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_WriteDataArrayVTU(vtu,dataname,ncomp,idist)
!
      implicit none
!
! Input:
!
      integer*4 vtu,ncomp
      integer*4 idist
      character (len = *) dataname
!
      write(vtu,'(A)',advance='no') '    <DataArray '
      write(vtu,'(A)',advance='no') 'type="Float32" '
      write(vtu,'(A)',advance='no') 'Name="'
      write(vtu,'(A)',advance='no') dataname
      write(vtu,'(A)',advance='no') '" NumberOfComponents="'
      write(vtu,'(I0)',advance='no') ncomp
      write(vtu,'(A)',advance='no') '" format="append" '
      write(vtu,'(A)',advance='no') 'offset="'
      write(vtu,'(I0)',advance='no') idist
      write(vtu,'(A)',advance='yes') '"/>'

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagAll
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
      call ppiclf_prints('*********** PPICLF OUTPUT *****************$')
      call ppiclf_io_OutputDiagGen
      call ppiclf_io_OutputDiagGhost
      if (ppiclf_overlap) call ppiclf_io_OutputDiagGrid

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGen
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot, npart_ide, nbin_total
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin General Info$')
         npart_max = ppiclf_iglmax(ppiclf_npart,1)
         npart_min = ppiclf_iglmin(ppiclf_npart,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart,1)
         npart_ide = npart_tot/ppiclf_np

         nbin_total = ppiclf_n_bins(1)*ppiclf_n_bins(2)*ppiclf_n_bins(3)

      call ppiclf_printsi('  -Cycle                  :$',ppiclf_cycle)
      call ppiclf_printsi('  -Output Freq.           :$',ppiclf_iostep)
      call ppiclf_printsr('  -Time                   :$',ppiclf_time)
      call ppiclf_printsr('  -dt                     :$',ppiclf_dt)
      call ppiclf_printsi('  -Global particles       :$',npart_tot)
      call ppiclf_printsi('  -Local particles (Max)  :$',npart_max)
      call ppiclf_printsi('  -Local particles (Min)  :$',npart_min)
      call ppiclf_printsi('  -Local particles (Ideal):$',npart_ide)
      call ppiclf_printsi('  -Total ranks            :$',ppiclf_np)
      call ppiclf_printsi('  -Problem dimensions     :$',ppiclf_ndim)
      call ppiclf_printsi('  -Integration method     :$',ppiclf_imethod)
      call ppiclf_printsi('  -Number of bins total   :$',nbin_total)
      call ppiclf_printsi('  -Number of bins (x)     :$',
     >                    ppiclf_n_bins(1))
      call ppiclf_printsi('  -Number of bins (y)     :$'
     >                    ,ppiclf_n_bins(2))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsi('  -Number of bins (z)     :$'
     >                    ,ppiclf_n_bins(3))
      call ppiclf_printsr('  -Bin xl coordinate      :$',ppiclf_binb(1))
      call ppiclf_printsr('  -Bin xr coordinate      :$',ppiclf_binb(2))
      call ppiclf_printsr('  -Bin yl coordinate      :$',ppiclf_binb(3))
      call ppiclf_printsr('  -Bin yr coordinate      :$',ppiclf_binb(4))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zl coordinate      :$',ppiclf_binb(5))
      if (ppiclf_ndim .gt. 2)
     >call ppiclf_printsr('  -Bin zr coordinate      :$',ppiclf_binb(6))

      call ppiclf_prints('  End General Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGrid
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 nel_max_orig, nel_min_orig, nel_total_orig, nel_max_map,
     >          nel_min_map, nel_total_map
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin Grid Info$')

         nel_max_orig   = ppiclf_iglmax(ppiclf_nFVCells,1)
         nel_min_orig   = ppiclf_iglmin(ppiclf_nFVCells,1)
         nel_total_orig = ppiclf_iglsum(ppiclf_nFVCells,1)

         nel_max_map   = ppiclf_iglmax(ppiclf_nCells_FV2PICL,1)
         nel_min_map   = ppiclf_iglmin(ppiclf_nCells_FV2PICL,1)
         nel_total_map = ppiclf_iglsum(ppiclf_nCells_FV2PICL,1)

      call ppiclf_printsi('  -Orig. Global cells     :$',nel_total_orig)
      call ppiclf_printsi('  -Orig. Local cells (Max):$',nel_max_orig)
      call ppiclf_printsi('  -Orig. Local cells (Min):$',nel_min_orig)
      call ppiclf_printsi('  -Map Global cells       :$',nel_total_map)
      call ppiclf_printsi('  -Map Local cells (Max)  :$',nel_max_map)
      call ppiclf_printsi('  -Map Local cells (Min)  :$',nel_min_map)

      call ppiclf_prints('  End Grid Info$')

      return
      end
!-----------------------------------------------------------------------
      subroutine ppiclf_io_OutputDiagGhost
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 npart_max, npart_min, npart_tot
      integer*4 ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
      external ppiclf_iglmax, ppiclf_iglmin, ppiclf_iglsum
!
      call ppiclf_prints(' *Begin Ghost Info$')

         npart_max = ppiclf_iglmax(ppiclf_npart_gp,1)
         npart_min = ppiclf_iglmin(ppiclf_npart_gp,1)
         npart_tot = ppiclf_iglsum(ppiclf_npart_gp,1)

      call ppiclf_printsi('  -Global ghosts          :$',npart_tot)
      call ppiclf_printsi('  -Local ghosts (Max)     :$',npart_max)
      call ppiclf_printsi('  -Local ghosts (Min)     :$',npart_min)

      call ppiclf_prints('  End Ghost Info$')

      return
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Particle dump (restart_data/particles_<step>.dat; layout in
! modflu/RFLU_ModRestartData.F90): one file for all ranks. Rank 0 writes
! the header with stream access; then every rank writes its own block
! with MPI-IO, after the blocks of the lower ranks: y, rprop, iprop, the
! viscous-unsteady history (rprop3) and the feedback to the gas (feedbk)
! of this rank's particles, then the fields projected to this rank's
! fluid cells. At the start of a stage the gas reads the projected fields
! and ppiclF then projects feedbk of the stage before: a restart needs both.
! A column with the same bits for every particle is in the header only;
! of the projected fields, only the cells with a value not 0 are written.
      subroutine ppiclf_io_WriteRestart(fname,istep,time)
      use ppiclf_data
      use ppiclf_unsteady
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
! Input:
!
      character (len = *) fname
      integer*8 istep
      real*8 time
!
! Internal:
!
      integer*4 iunit, fh, ierr, j, k, nhist, nsaved, nnz
      integer*4 saved(ppiclf_lrp)
      integer*4 status(MPI_STATUS_SIZE)
      integer*4 sty(ppiclf_lrs), str(ppiclf_lrp), sti(PPICLF_LIP),
     >          sth(max(ppiclf_lrp3,1)), stf(ppiclf_lrp_pro),
     >          shi(PPICLF_LIP)
      integer*8 nlocal, ntotal, npart_rank(ppiclf_np), header_size,
     >          ncells, ncell_rank(ppiclf_np), nnzl,
     >          nnz_rank(ppiclf_np), perpart
      integer(kind=MPI_OFFSET_KIND) offset
      real*8 shy(ppiclf_lrs), shr(ppiclf_lrp),
     >       shh(max(ppiclf_lrp3,1)), shf(ppiclf_lrp_pro)
      integer*4, allocatable :: cells(:)
      real*8, allocatable :: nonzero(:,:), props(:,:)
!
      ! the output fields (ppiclf_user_OutputFields), built once
      if (ppiclf_nout .lt. 0) then
         ppiclf_nout = 0
         call ppiclf_user_OutputFields
      endif

      ! rprop slots saved: those a restart needs and those of the output
      call ppiclf_io_SavedSlots(saved,nsaved)

      ! pack the viscous-unsteady history into rprop3
      if (ppiclf_lrp3 .gt. 0) call ppiclf_user_plag2prop
      nhist = max(ppiclf_nUnsteadyData,0)

      nlocal = ppiclf_npart
      call mpi_allgather(nlocal,1,MPI_INTEGER8,npart_rank,1,
     >                   MPI_INTEGER8,ppiclf_comm,ierr)
      ntotal = sum(npart_rank)
      ncells = ppiclf_nFVCells
      call mpi_allgather(ncells,1,MPI_INTEGER8,ncell_rank,1,
     >                   MPI_INTEGER8,ppiclf_comm,ierr)

      allocate(props(max(nsaved,1),max(nlocal,1_8)))
      if (nlocal .gt. 0 .and. nsaved .gt. 0)
     >   props(1:nsaved,1:nlocal) =
     >      ppiclf_rprop(saved(1:nsaved),1:nlocal)

      ! columns with the same bits for every particle of every rank
      sth = 1
      shh = 0.0d0
      call ppiclf_io_SharedColumns(ppiclf_y,ppiclf_lrs,nlocal,sty,shy)
      call ppiclf_io_SharedColumns(props,nsaved,nlocal,str,shr)
      call ppiclf_io_SharedTags(ppiclf_iprop,PPICLF_LIP,nlocal,sti,shi)
      if (ppiclf_lrp3 .gt. 0)
     >   call ppiclf_io_SharedColumns(ppiclf_rprop3,ppiclf_lrp3,nlocal,
     >                                sth,shh)
      call ppiclf_io_SharedColumns(ppiclf_feedbk,ppiclf_lrp_pro,nlocal,
     >                             stf,shf)
      perpart = 8*count(sty .eq. 1) + 8*count(str(1:nsaved) .eq. 1)
     >        + 4*count(sti .eq. 1)
     >        + 8*count(sth(1:ppiclf_lrp3) .eq. 1)
     >        + 8*count(stf .eq. 1)

      ! this rank's cells with a projected field not 0 (bit pattern)
      allocate(cells(max(ncells,1_8)))
      nnz = 0
      do j=1,int(ncells,4)
         if (any(transfer(ppiclf_pro_fld(j,1:ppiclf_lrp_pro),1_8,
     >                    ppiclf_lrp_pro) .ne. 0_8)) then
            nnz = nnz + 1
            cells(nnz) = j
         endif
      enddo
      nnzl = nnz
      call mpi_allgather(nnzl,1,MPI_INTEGER8,nnz_rank,1,
     >                   MPI_INTEGER8,ppiclf_comm,ierr)

      if (ppiclf_nid .eq. 0) then
         open(newunit=iunit,file=fname,form='unformatted',
     >        access='stream',status='replace',iostat=ierr)
         if (ierr .ne. 0)
     >      call ppiclf_exittr('Cannot open particle dump$',0.0d0,0)
         write(iunit) 'PPICLDMP', 4_4, istep, time,
     >      int(ppiclf_np,4), ntotal,
     >      int(ppiclf_lrs,4), int(ppiclf_lrp,4), int(PPICLF_LIP,4),
     >      int(ppiclf_lrp3,4), int(ppiclf_lrp_pro,4), int(nhist,4),
     >      int(ppiclf_nTimeBH,4)
         if (nhist .gt. 0) write(iunit) ppiclf_TimeBH(1:nhist)
         write(iunit) ppiclf_n_bins, ppiclf_binb, ppiclf_bins_dx
         write(iunit) int(ppiclf_nout,4)
         do k=1,ppiclf_nout
            write(iunit) int(ppiclf_out_source(k),4),
     >                   int(ppiclf_out_slot(k),4), ppiclf_out_name(k)
         enddo
         write(iunit) int(nsaved,4), saved(1:nsaved)
         write(iunit) sty, shy
         write(iunit) str(1:nsaved), shr(1:nsaved)
         write(iunit) sti, shi
         write(iunit) sth(1:ppiclf_lrp3), shh(1:ppiclf_lrp3)
         write(iunit) stf, shf
         write(iunit) npart_rank, ncell_rank, nnz_rank
         ! the next byte to write, 1-based: the header is one byte less
         inquire(unit=iunit,pos=header_size)
         header_size = header_size - 1
         close(iunit)
      endif
      call mpi_bcast(header_size,1,MPI_INTEGER8,0,ppiclf_comm,ierr)

      call mpi_file_open(ppiclf_comm,fname,MPI_MODE_WRONLY,
     >                   MPI_INFO_NULL,fh,ierr)
      if (ierr .ne. MPI_SUCCESS)
     >   call ppiclf_exittr('Cannot open particle dump$',0.0d0,0)

      offset = header_size
      do k=1,ppiclf_nid
         offset = offset + npart_rank(k)*perpart
     >                   + nnz_rank(k)*(4 + 8*ppiclf_lrp_pro)
      enddo

      call ppiclf_io_WriteColumns(fh,offset,ppiclf_y,ppiclf_lrs,
     >                            nlocal,sty)
      call ppiclf_io_WriteColumns(fh,offset,props,nsaved,nlocal,str)
      deallocate(props)
      call ppiclf_io_WriteTags(fh,offset,ppiclf_iprop,PPICLF_LIP,
     >                         nlocal,sti)
      if (ppiclf_lrp3 .gt. 0)
     >   call ppiclf_io_WriteColumns(fh,offset,ppiclf_rprop3,
     >                               ppiclf_lrp3,nlocal,sth)
      ! the feedback to the gas, which the next stage projects first
      call ppiclf_io_WriteColumns(fh,offset,ppiclf_feedbk,
     >                            ppiclf_lrp_pro,nlocal,stf)
      ! the fields projected to this rank's fluid cells, which the next
      ! stage reads: the cells not 0, then their values
      call mpi_file_write_at_all(fh,offset,cells,nnz,MPI_INTEGER4,
     >                           status,ierr)
      offset = offset + 4*nnzl
      allocate(nonzero(max(nnz,1),ppiclf_lrp_pro))
      if (nnz .gt. 0) nonzero(1:nnz,:) =
     >   ppiclf_pro_fld(cells(1:nnz),1:ppiclf_lrp_pro)
      call mpi_file_write_at_all(fh,offset,nonzero,
     >   nnz*ppiclf_lrp_pro,MPI_DOUBLE_PRECISION,status,ierr)
      deallocate(nonzero,cells)

      call mpi_file_close(fh,ierr)

      return
      end
!-----------------------------------------------------------------------
! Particle dump: which columns of a(ncol,n) hold the same bits for every
! particle of every rank (stored 0, the value in shared); the others
! (stored 1) are written per particle. Collective.
      subroutine ppiclf_io_SharedColumns(a,ncol,n,stored,shared)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
      integer*4 ncol, stored(*)
      integer*8 n
      real*8 a(ncol,*), shared(*)
!
      integer*4 ierr, j
      integer*8 i, bits, lo(max(ncol,1)), hi(max(ncol,1)),
     >          glo(max(ncol,1)), ghi(max(ncol,1))
!
      if (ncol .le. 0) return
      lo = huge(1_8)
      hi = -huge(1_8)
      do i=1,n
         do j=1,ncol
            bits = transfer(a(j,i),1_8)
            lo(j) = min(lo(j),bits)
            hi(j) = max(hi(j),bits)
         enddo
      enddo
      call mpi_allreduce(lo,glo,ncol,MPI_INTEGER8,MPI_MIN,
     >                   ppiclf_comm,ierr)
      call mpi_allreduce(hi,ghi,ncol,MPI_INTEGER8,MPI_MAX,
     >                   ppiclf_comm,ierr)
      do j=1,ncol
         stored(j) = 1
         shared(j) = 0.0d0
         if (glo(j) .eq. ghi(j)) then
            stored(j) = 0
            shared(j) = transfer(glo(j),1.0d0)
         endif
      enddo

      return
      end
!-----------------------------------------------------------------------
! ppiclf_io_SharedColumns for the integer columns (iprop)
      subroutine ppiclf_io_SharedTags(a,ncol,n,stored,shared)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
      integer*4 ncol, stored(*), a(ncol,*), shared(*)
      integer*8 n
!
      integer*4 ierr, j, lo(max(ncol,1)), hi(max(ncol,1)),
     >          glo(max(ncol,1)), ghi(max(ncol,1))
      integer*8 i
!
      if (ncol .le. 0) return
      lo = huge(1_4)
      hi = -huge(1_4)
      do i=1,n
         do j=1,ncol
            lo(j) = min(lo(j),a(j,i))
            hi(j) = max(hi(j),a(j,i))
         enddo
      enddo
      call mpi_allreduce(lo,glo,ncol,MPI_INTEGER4,MPI_MIN,
     >                   ppiclf_comm,ierr)
      call mpi_allreduce(hi,ghi,ncol,MPI_INTEGER4,MPI_MAX,
     >                   ppiclf_comm,ierr)
      do j=1,ncol
         stored(j) = 1
         shared(j) = 0
         if (glo(j) .eq. ghi(j)) then
            stored(j) = 0
            shared(j) = glo(j)
         endif
      enddo

      return
      end
!-----------------------------------------------------------------------
! Particle dump: write the stored columns of a(ncol,n) at offset, as
! (columns,n), and move offset past them. Collective.
      subroutine ppiclf_io_WriteColumns(fh,offset,a,ncol,n,stored)
!
      implicit none
!
      include 'mpif.h'
!
      integer*4 fh, ncol, stored(*)
      integer*8 n
      integer(kind=MPI_OFFSET_KIND) offset
      real*8 a(ncol,*)
!
      integer*4 ierr, j, k, nkept
      integer*4 status(MPI_STATUS_SIZE)
      integer*8 i
      real*8, allocatable :: packed(:,:)
!
      nkept = count(stored(1:ncol) .eq. 1)
      allocate(packed(max(nkept,1),max(n,1_8)))
      k = 0
      do j=1,ncol
         if (stored(j) .eq. 1) then
            k = k + 1
            do i=1,n
               packed(k,i) = a(j,i)
            enddo
         endif
      enddo
      call mpi_file_write_at_all(fh,offset,packed,int(nkept*n,4),
     >                           MPI_DOUBLE_PRECISION,status,ierr)
      offset = offset + 8*nkept*n
      deallocate(packed)

      return
      end
!-----------------------------------------------------------------------
! ppiclf_io_WriteColumns for the integer columns (iprop)
      subroutine ppiclf_io_WriteTags(fh,offset,a,ncol,n,stored)
!
      implicit none
!
      include 'mpif.h'
!
      integer*4 fh, ncol, stored(*), a(ncol,*)
      integer*8 n
      integer(kind=MPI_OFFSET_KIND) offset
!
      integer*4 ierr, j, k, nkept
      integer*4 status(MPI_STATUS_SIZE)
      integer*8 i
      integer*4, allocatable :: packed(:,:)
!
      nkept = count(stored(1:ncol) .eq. 1)
      allocate(packed(max(nkept,1),max(n,1_8)))
      k = 0
      do j=1,ncol
         if (stored(j) .eq. 1) then
            k = k + 1
            do i=1,n
               packed(k,i) = a(j,i)
            enddo
         endif
      enddo
      call mpi_file_write_at_all(fh,offset,packed,int(nkept*n,4),
     >                           MPI_INTEGER4,status,ierr)
      offset = offset + 4*nkept*n
      deallocate(packed)

      return
      end
!-----------------------------------------------------------------------
! Particle dump: read the stored columns of a(ncol,n) at offset (as
! ppiclf_io_WriteColumns wrote them), set the others to their shared
! value, and move offset past them. Collective.
      subroutine ppiclf_io_ReadColumns(fh,offset,a,ncol,n,stored,shared)
!
      implicit none
!
      include 'mpif.h'
!
      integer*4 fh, ncol, stored(*)
      integer*8 n
      integer(kind=MPI_OFFSET_KIND) offset
      real*8 a(ncol,*), shared(*)
!
      integer*4 ierr, j, k, nkept
      integer*4 status(MPI_STATUS_SIZE)
      integer*8 i
      real*8, allocatable :: packed(:,:)
!
      nkept = count(stored(1:ncol) .eq. 1)
      allocate(packed(max(nkept,1),max(n,1_8)))
      call mpi_file_read_at_all(fh,offset,packed,int(nkept*n,4),
     >                          MPI_DOUBLE_PRECISION,status,ierr)
      offset = offset + 8*nkept*n
      k = 0
      do j=1,ncol
         if (stored(j) .eq. 1) then
            k = k + 1
            do i=1,n
               a(j,i) = packed(k,i)
            enddo
         else
            do i=1,n
               a(j,i) = shared(j)
            enddo
         endif
      enddo
      deallocate(packed)

      return
      end
!-----------------------------------------------------------------------
! ppiclf_io_ReadColumns for the integer columns (iprop)
      subroutine ppiclf_io_ReadTags(fh,offset,a,ncol,n,stored,shared)
!
      implicit none
!
      include 'mpif.h'
!
      integer*4 fh, ncol, stored(*), a(ncol,*), shared(*)
      integer*8 n
      integer(kind=MPI_OFFSET_KIND) offset
!
      integer*4 ierr, j, k, nkept
      integer*4 status(MPI_STATUS_SIZE)
      integer*8 i
      integer*4, allocatable :: packed(:,:)
!
      nkept = count(stored(1:ncol) .eq. 1)
      allocate(packed(max(nkept,1),max(n,1_8)))
      call mpi_file_read_at_all(fh,offset,packed,int(nkept*n,4),
     >                          MPI_INTEGER4,status,ierr)
      offset = offset + 4*nkept*n
      k = 0
      do j=1,ncol
         if (stored(j) .eq. 1) then
            k = k + 1
            do i=1,n
               a(j,i) = packed(k,i)
            enddo
         else
            do i=1,n
               a(j,i) = shared(j)
            enddo
         endif
      enddo
      deallocate(packed)

      return
      end
!-----------------------------------------------------------------------
! Restart from a particle dump (ppiclf_io_WriteRestart): every rank reads
! its own block, so particles keep their rank, order, IDs, random streams
! and viscous-unsteady history. Call instead of reading points.dat, before
! ppiclf_solve_InitParticle (which then adds no particles). The dump must
! come from a run with the same number of ranks and the same particle
! physics (sizes of y, rprop, iprop and the history).
!
! Input:  fname, istartoutin (par*.vtu files already present: new ones
!         continue the numbering)
! Output: ntotal (particles, all ranks), dp_min, dp_max and pbox
!         (particle box x/y/z min,max with dp/2 margin, all ranks)
      subroutine ppiclf_io_ReadRestart(fname,istartoutin,ntotal,
     >                                 dp_min,dp_max,pbox)
      use ppiclf_data
      use ppiclf_unsteady
!
      implicit none
!
      include "PPICLF"
      include 'mpif.h'
!
      character (len = *) fname
      integer*4 istartoutin
      integer*8 ntotal
      real*8 dp_min, dp_max, pbox(6)
!
      character*8 magic
      character*24 fieldname
      integer*4 iunit, fh, ierr, k, nhist, nranks, nout, version,
     >          lrs, lrp, lip, lrp3, lrppro, ntimebh, ibad, idum(2),
     >          nbins(3), nsaved, needed(ppiclf_lrp), nneeded,
     >          saved(ppiclf_lrp)
      integer*4 status(MPI_STATUS_SIZE)
      integer*4 nnz
      integer*8 istep, nlocal, npart_rank(ppiclf_np), header_size,
     >          ncells, ncell_rank(ppiclf_np), nnz_rank(ppiclf_np),
     >          perpart
      integer(kind=MPI_OFFSET_KIND) offset
      real*8 time, rdum(9), dpl
      integer*4 j
      real*8, allocatable :: ytmp(:,:), rtmp(:,:), htmp(:,:),
     >                       ftmp(:,:), props(:,:), timebh(:),
     >                       nonzero(:,:)
      real*8, allocatable :: shy(:), shr(:), shh(:), shf(:)
      integer*4, allocatable :: itmp(:,:), cells(:)
      integer*4, allocatable :: sty(:), str(:), sti(:), sth(:),
     >                          stf(:), shi(:)
      integer*4 istartout
      common /ppiclf_io_restart/ istartout
!
      call ppiclf_solve_InitZero
      PPICLF_RESTART = .true.
      istartout = istartoutin

      ! column tables of the dump; ppiclf_lrp3 is set by InitZero
      allocate(sty(ppiclf_lrs), str(ppiclf_lrp), sti(PPICLF_LIP),
     >         sth(max(ppiclf_lrp3,1)), stf(ppiclf_lrp_pro),
     >         shi(PPICLF_LIP), shy(ppiclf_lrs), shr(ppiclf_lrp),
     >         shh(max(ppiclf_lrp3,1)), shf(ppiclf_lrp_pro))

      nhist = max(ppiclf_nUnsteadyData,0)
      allocate(timebh(max(nhist,1)))
      timebh = 0.0d0
      ! version 3: every column per particle, every cell's fields
      sty = 1
      str = 1
      sti = 1
      sth = 1
      stf = 1
      shy = 0.0d0
      shr = 0.0d0
      shi = 0
      shh = 0.0d0
      shf = 0.0d0
      nnz_rank = 0

      ibad = 0
      version = 0
      if (ppiclf_nid .eq. 0) then
         open(newunit=iunit,file=fname,form='unformatted',
     >        access='stream',status='old',iostat=ierr)
         if (ierr .ne. 0) then
            ibad = 1
         else
            read(iunit,iostat=ierr) magic, version, istep, time,
     >         nranks, ntotal, lrs, lrp, lip, lrp3, lrppro, k, ntimebh
            if (ierr .ne. 0 .or. magic .ne. 'PPICLDMP' .or.
     >          (version .ne. 3 .and. version .ne. 4)) then
               ibad = 2
            elseif (nranks .ne. ppiclf_np) then
               ibad = 3
            elseif (lrs .ne. ppiclf_lrs .or. lrp .ne. ppiclf_lrp .or.
     >              lip .ne. PPICLF_LIP .or. lrp3 .ne. ppiclf_lrp3 .or.
     >              lrppro .ne. ppiclf_lrp_pro .or. k .ne. nhist) then
               ibad = 4
            else
               if (nhist .gt. 0) read(iunit) timebh(1:nhist)
               read(iunit) nbins, rdum
               read(iunit) nout
               do k=1,nout
                  read(iunit) idum, fieldname
               enddo
               read(iunit) nsaved
               if (nsaved .gt. ppiclf_lrp) then
                  ibad = 4
               else
                  read(iunit) saved(1:nsaved)
                  ! every slot a restart needs must be there
                  call ppiclf_io_SavedSlots(needed,nneeded)
                  do k=1,ppiclf_lrp_restart
                     if (.not. any(saved(1:nsaved) .eq.
     >                   ppiclf_restart_slots(k))) ibad = 4
                  enddo
               endif
               if (ibad .eq. 0 .and. version .eq. 4) then
                  read(iunit) sty, shy
                  read(iunit) str(1:nsaved), shr(1:nsaved)
                  read(iunit) sti, shi
                  read(iunit) sth(1:ppiclf_lrp3), shh(1:ppiclf_lrp3)
                  read(iunit) stf, shf
                  read(iunit) npart_rank, ncell_rank, nnz_rank
               else
                  read(iunit) npart_rank, ncell_rank
               endif
               ! the next byte to read, 1-based
               inquire(unit=iunit,pos=header_size)
               header_size = header_size - 1
            endif
            close(iunit)
         endif
      endif

      call mpi_bcast(ibad,1,MPI_INTEGER4,0,ppiclf_comm,ierr)
      if (ibad .eq. 1) call ppiclf_exittr(
     >   'Cannot open the particle dump$',0.0d0,0)
      if (ibad .eq. 2) call ppiclf_exittr(
     >   'Not a version 3 or 4 particle dump$',0.0d0,0)
      if (ibad .eq. 3) call ppiclf_exittr(
     >   'Particle dump from a different number of ranks$',0.0d0,0)
      if (ibad .eq. 4) call ppiclf_exittr(
     >   'Particle dump from different particle physics$',0.0d0,0)

      call mpi_bcast(header_size,1,MPI_INTEGER8,0,ppiclf_comm,ierr)
      call mpi_bcast(npart_rank,ppiclf_np,MPI_INTEGER8,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(ncell_rank,ppiclf_np,MPI_INTEGER8,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(nsaved,1,MPI_INTEGER4,0,ppiclf_comm,ierr)
      call mpi_bcast(saved,ppiclf_lrp,MPI_INTEGER4,0,ppiclf_comm,ierr)
      call mpi_bcast(version,1,MPI_INTEGER4,0,ppiclf_comm,ierr)
      call mpi_bcast(nnz_rank,ppiclf_np,MPI_INTEGER8,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(sty,ppiclf_lrs,MPI_INTEGER4,0,ppiclf_comm,ierr)
      call mpi_bcast(str,ppiclf_lrp,MPI_INTEGER4,0,ppiclf_comm,ierr)
      call mpi_bcast(sti,PPICLF_LIP,MPI_INTEGER4,0,ppiclf_comm,ierr)
      call mpi_bcast(sth,max(ppiclf_lrp3,1),MPI_INTEGER4,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(stf,ppiclf_lrp_pro,MPI_INTEGER4,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(shy,ppiclf_lrs,MPI_DOUBLE_PRECISION,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(shr,ppiclf_lrp,MPI_DOUBLE_PRECISION,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(shi,PPICLF_LIP,MPI_INTEGER4,0,ppiclf_comm,ierr)
      call mpi_bcast(shh,max(ppiclf_lrp3,1),MPI_DOUBLE_PRECISION,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(shf,ppiclf_lrp_pro,MPI_DOUBLE_PRECISION,0,
     >               ppiclf_comm,ierr)
      call mpi_bcast(ntotal,1,MPI_INTEGER8,0,ppiclf_comm,ierr)
      ! Lines of points.dat per rank (restarts keep the rank count)
      ppiclf_idlines = int((ntotal + ppiclf_np - 1)/ppiclf_np,4)
      call mpi_bcast(istep,1,MPI_INTEGER8,0,ppiclf_comm,ierr)
      call mpi_bcast(time,1,MPI_DOUBLE_PRECISION,0,ppiclf_comm,ierr)
      call mpi_bcast(ntimebh,1,MPI_INTEGER4,0,ppiclf_comm,ierr)
      ! the bin box (rdum(1:6)) and bin counts of the dump
      call mpi_bcast(rdum,9,MPI_DOUBLE_PRECISION,0,ppiclf_comm,ierr)
      call mpi_bcast(nbins,3,MPI_INTEGER4,0,ppiclf_comm,ierr)
      ppiclf_restart_binb(1:6) = rdum(1:6)
      ppiclf_restart_bins = .true.
      call mpi_bcast(timebh,max(nhist,1),MPI_DOUBLE_PRECISION,0,
     >               ppiclf_comm,ierr)

      nlocal = npart_rank(ppiclf_nid+1)
      allocate(ytmp(ppiclf_lrs,max(nlocal,1_8)),
     >         rtmp(ppiclf_lrp,max(nlocal,1_8)),
     >         props(max(nsaved,1),max(nlocal,1_8)),
     >         itmp(PPICLF_LIP,max(nlocal,1_8)),
     >         htmp(max(ppiclf_lrp3,1),max(nlocal,1_8)),
     >         ftmp(ppiclf_lrp_pro,max(nlocal,1_8)))
      ncells = ncell_rank(ppiclf_nid+1)
      if (allocated(ppiclf_restart_pro_fld))
     >   deallocate(ppiclf_restart_pro_fld)
      allocate(ppiclf_restart_pro_fld(ncells,ppiclf_lrp_pro))

      ! the same layout as ppiclf_io_WriteRestart
      perpart = 8*count(sty .eq. 1) + 8*count(str(1:nsaved) .eq. 1)
     >        + 4*count(sti .eq. 1)
     >        + 8*count(sth(1:ppiclf_lrp3) .eq. 1)
     >        + 8*count(stf .eq. 1)
      offset = header_size
      do k=1,ppiclf_nid
         if (version .eq. 3) then
            offset = offset + npart_rank(k)*perpart
     >                      + ncell_rank(k)*8*ppiclf_lrp_pro
         else
            offset = offset + npart_rank(k)*perpart
     >                      + nnz_rank(k)*(4 + 8*ppiclf_lrp_pro)
         endif
      enddo

      call mpi_file_open(ppiclf_comm,fname,MPI_MODE_RDONLY,
     >                   MPI_INFO_NULL,fh,ierr)
      call ppiclf_io_ReadColumns(fh,offset,ytmp,ppiclf_lrs,nlocal,
     >                           sty,shy)
      call ppiclf_io_ReadColumns(fh,offset,props,nsaved,nlocal,
     >                           str,shr)
      ! slots not saved start at 0: interpolated before they are used
      rtmp = 0.0d0
      if (nlocal .gt. 0 .and. nsaved .gt. 0)
     >   rtmp(saved(1:nsaved),1:nlocal) = props(1:nsaved,1:nlocal)
      call ppiclf_io_ReadTags(fh,offset,itmp,PPICLF_LIP,nlocal,
     >                        sti,shi)
      if (ppiclf_lrp3 .gt. 0)
     >   call ppiclf_io_ReadColumns(fh,offset,htmp,ppiclf_lrp3,nlocal,
     >                              sth,shh)
      call ppiclf_io_ReadColumns(fh,offset,ftmp,ppiclf_lrp_pro,nlocal,
     >                           stf,shf)
      if (version .eq. 3) then
         call mpi_file_read_at_all(fh,offset,ppiclf_restart_pro_fld,
     >      int(ncells*ppiclf_lrp_pro,4),MPI_DOUBLE_PRECISION,status,
     >      ierr)
      else
         ! the cells not 0, then their values; the others are 0
         nnz = int(nnz_rank(ppiclf_nid+1),4)
         allocate(cells(max(nnz,1)),nonzero(max(nnz,1),ppiclf_lrp_pro))
         call mpi_file_read_at_all(fh,offset,cells,nnz,MPI_INTEGER4,
     >                             status,ierr)
         offset = offset + 4*nnz_rank(ppiclf_nid+1)
         call mpi_file_read_at_all(fh,offset,nonzero,
     >      nnz*ppiclf_lrp_pro,MPI_DOUBLE_PRECISION,status,ierr)
         if (nnz .gt. 0) then
            if (minval(cells(1:nnz)) .lt. 1 .or.
     >          maxval(cells(1:nnz)) .gt. ncells)
     >         call ppiclf_exittr('Particle dump: bad cell number$',
     >                            0.0d0,0)
         endif
         ppiclf_restart_pro_fld = 0.0d0
         if (nnz .gt. 0) ppiclf_restart_pro_fld(cells(1:nnz),:) =
     >      nonzero(1:nnz,1:ppiclf_lrp_pro)
         deallocate(cells,nonzero)
      endif
      call mpi_file_close(fh,ierr)

      ! state and properties (PPICLF_RESTART: the IDs are not reset)
      call ppiclf_solve_AddParticles(int(nlocal,4),ytmp,rtmp)
      if (nlocal .gt. 0) then
         ppiclf_iprop(:,1:nlocal) = itmp(:,1:nlocal)
         ppiclf_feedbk(:,1:nlocal) = ftmp(:,1:nlocal)
         if (ppiclf_lrp3 .gt. 0) then
            ppiclf_rprop3(:,1:nlocal) = htmp(:,1:nlocal)
         endif
      endif
      if (ppiclf_lrp3 .gt. 0) call ppiclf_user_prop2plag
      if (nhist .gt. 0) ppiclf_TimeBH(1:nhist) = timebh(1:nhist)
      ppiclf_nTimeBH = ntimebh

      ppiclf_time  = time
      ppiclf_cycle = int(istep,4)

      ! diameter range over all ranks
      dp_min = huge(1.0d0)
      dp_max = 0.0d0
      do k=1,int(nlocal,4)
         dp_min = min(dp_min,ppiclf_rprop(PPICLF_R_JDP,k))
         dp_max = max(dp_max,ppiclf_rprop(PPICLF_R_JDP,k))
      enddo
      dpl = dp_min
      call mpi_allreduce(dpl,dp_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     >                   ppiclf_comm,ierr)
      dpl = dp_max
      call mpi_allreduce(dpl,dp_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     >                   ppiclf_comm,ierr)

      ! particle box over all ranks, each particle with its radius
      do j=1,3
         pbox(2*j-1) =  huge(1.0d0)
         pbox(2*j)   = -huge(1.0d0)
         do k=1,int(nlocal,4)
            pbox(2*j-1) = min(pbox(2*j-1),
     >             ppiclf_y(j,k)-0.5d0*ppiclf_rprop(PPICLF_R_JDP,k))
            pbox(2*j)   = max(pbox(2*j),
     >             ppiclf_y(j,k)+0.5d0*ppiclf_rprop(PPICLF_R_JDP,k))
         enddo
         call mpi_allreduce(MPI_IN_PLACE,pbox(2*j-1),1,
     >        MPI_DOUBLE_PRECISION,MPI_MIN,ppiclf_comm,ierr)
         call mpi_allreduce(MPI_IN_PLACE,pbox(2*j),1,
     >        MPI_DOUBLE_PRECISION,MPI_MAX,ppiclf_comm,ierr)
      enddo

      deallocate(ytmp,rtmp,itmp,htmp,ftmp,props,timebh)
      deallocate(sty,str,sti,sth,stf,shi,shy,shr,shh,shf)

      return
      end
!-----------------------------------------------------------------------
! The rprop slots the particle dump holds, in increasing order: those a
! restart needs (ppiclf_restart_slots) and those of the particle output
! (ppiclf_user_OutputFields).
      subroutine ppiclf_io_SavedSlots(saved,nsaved)
      use ppiclf_data
!
      implicit none
!
      include "PPICLF"
!
      integer*4 saved(*), nsaved
!
      integer*4 k, j, slot
      logical keep(ppiclf_lrp)
!
      if (ppiclf_nout .lt. 0) then
         ppiclf_nout = 0
         call ppiclf_user_OutputFields
      endif
      keep = .false.
      do k=1,ppiclf_lrp_restart
         keep(ppiclf_restart_slots(k)) = .true.
      enddo
      do k=1,ppiclf_nout
         if (ppiclf_out_source(k) .eq. PPICLF_OUT_RPROP) then
            slot = ppiclf_out_slot(k)
            if (slot .ge. 1 .and. slot .le. ppiclf_lrp)
     >         keep(slot) = .true.
         endif
      enddo
      nsaved = 0
      do j=1,ppiclf_lrp
         if (keep(j)) then
            nsaved = nsaved + 1
            saved(nsaved) = j
         endif
      enddo

      return
      end
