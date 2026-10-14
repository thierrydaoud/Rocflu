!-----------------------------------------------------------------------
!
! Particle-wall contacts.
!
! Walls are triangles (ppiclf_wall_c: vertices, ppiclf_wall_n: unit
! normal and area). They come from filein.vtk (ppiclf_io_ReadWallVTK)
! or from Rocflu wall patches marked PARTWALL 1 in the .bc file.
! A wall with ppiclf_wall_side = 1 knows its fluid side (the normal
! points into the fluid) and only pushes particles towards it, so a
! particle whose centre has crossed the wall is pushed back. A wall
! with ppiclf_wall_side = 0 (filein.vtk, whose triangle winding is
! arbitrary) pushes particles away from whichever side they are on.
!
! Each RK stage the walls are filed into this rank's collision
! subbins (ppiclf_solve_WallSubbins), and ppiclf_solve_WallContacts
! finds the contacts of one particle from its own subbin's list.
!
!-----------------------------------------------------------------------
      MODULE ppiclf_wallbin
!
! Wall triangles of subbin s: wsb_list(wsb_start(s):wsb_start(s+1)-1)
! Subbin (l1,l2,l3) spans wsb_x0 + [l, l+1]*wsb_h in each direction,
! the same grid as ppiclf_user_subbinMap.
!
      IMPLICIT NONE
      INTEGER*4 :: wsb_n(3) = 0, wsb_tot = 0, wsb_nwall = -1
      REAL*8 :: wsb_x0(3) = 0.0d0, wsb_h = 0.0d0
      INTEGER*4, ALLOCATABLE :: wsb_start(:), wsb_list(:)
      END MODULE ppiclf_wallbin
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_AddWall(xp1,xp2,xp3,rnf,iside)
!
! Add triangle (xp1,xp2,xp3). iside = 1: rnf points into the fluid
! and the stored normal is oriented the same way; iside = 0: the
! normal follows the vertex order and rnf is not used.
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      REAL*8 xp1(3), xp2(3), xp3(3), rnf(3)
      INTEGER*4 iside
!
! Internal:
!
      REAL*8 e1(3), e2(3), rn(3), rmag
      INTEGER*4 l
!
      IF (ppiclf_ndim .NE. 3)
     >CALL ppiclf_exittr('Walls need a 3D ppiclF run$',0.d0,0)

      DO l = 1,3
         e1(l) = xp2(l) - xp1(l)
         e2(l) = xp3(l) - xp1(l)
      END DO
      rn(1) = e1(2)*e2(3) - e1(3)*e2(2)
      rn(2) = e1(3)*e2(1) - e1(1)*e2(3)
      rn(3) = e1(1)*e2(2) - e1(2)*e2(1)
      rmag  = sqrt(rn(1)**2 + rn(2)**2 + rn(3)**2)
      ! Degenerate triangle: no area, no contacts
      IF (rmag .LE. 0.0d0) RETURN

      ppiclf_nwall = ppiclf_nwall + 1
      IF (ppiclf_nwall .GT. PPICLF_LWALL)
     >CALL ppiclf_exittr('Increase PPICLF_LWALL in PPICLF_USER.h$'
     >                  ,0.d0,ppiclf_nwall)

      DO l = 1,3
         rn(l) = rn(l)/rmag
      END DO
      IF (iside .EQ. 1) THEN
         IF (rn(1)*rnf(1) + rn(2)*rnf(2) + rn(3)*rnf(3) .LT. 0.0d0)
     >   THEN
            DO l = 1,3
               rn(l) = -rn(l)
            END DO
         END IF
      END IF

      DO l = 1,3
         ppiclf_wall_c(l  ,ppiclf_nwall) = xp1(l)
         ppiclf_wall_c(l+3,ppiclf_nwall) = xp2(l)
         ppiclf_wall_c(l+6,ppiclf_nwall) = xp3(l)
         ppiclf_wall_n(l  ,ppiclf_nwall) = rn(l)
      END DO
      ppiclf_wall_n(4,ppiclf_nwall) = 0.5d0*rmag
      ppiclf_wall_side(ppiclf_nwall) = iside

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_WallSubbins(i_Bin,n_SBin,tot_SBin)
!
! File every wall triangle into the subbins of this rank's bin that
! lie within ppiclf_nndist of it (bounding box, then plane distance
! from the subbin centre), so that a particle finds every triangle
! within ppiclf_nndist of its centre in its own subbin's list. The
! lists are rebuilt only when the subbin grid moves.
!
      USE ppiclf_wallbin
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 i_Bin(3), n_SBin(3), tot_SBin
!
! Internal:
!
      REAL*8 x0(3), h, rlo(3), rhi(3), cc(3), rmax, dpl
      INTEGER*4 ilo(3), ihi(3), l, w, ipass, n, i1, i2, i3, sb, s
      LOGICAL same
!
      IF (ppiclf_nwall .LE. 0 .OR. ppiclf_ndim .NE. 3) RETURN

      h = ppiclf_nndist
      DO l = 1,3
         x0(l) = ppiclf_binb(2*l-1) + i_Bin(l)*ppiclf_bins_dx(l) - h
      END DO

      same = wsb_tot .EQ. tot_SBin .AND. wsb_h .EQ. h
     >       .AND. wsb_nwall .EQ. ppiclf_nwall
      DO l = 1,3
         same = same .AND. wsb_n(l) .EQ. n_SBin(l)
     >               .AND. wsb_x0(l) .EQ. x0(l)
      END DO
      IF (same) RETURN

      IF (ALLOCATED(wsb_start)) DEALLOCATE(wsb_start)
      ALLOCATE(wsb_start(0:tot_SBin))

      ! Distance from a subbin centre within which a triangle can be
      ! closer than h to some point of the subbin
      rmax = h + 0.5d0*sqrt(3.0d0)*h

      ! Pass 1 counts the entries of each subbin, pass 2 fills them
      DO ipass = 1,2
         IF (ipass .EQ. 1) THEN
            DO sb = 0,tot_SBin
               wsb_start(sb) = 0
            END DO
         END IF
         DO w = 1,ppiclf_nwall
            DO l = 1,3
               rlo(l) = min(ppiclf_wall_c(l,w),ppiclf_wall_c(l+3,w),
     >                      ppiclf_wall_c(l+6,w)) - h
               rhi(l) = max(ppiclf_wall_c(l,w),ppiclf_wall_c(l+3,w),
     >                      ppiclf_wall_c(l+6,w)) + h
               ilo(l) = max(0, FLOOR((rlo(l)-x0(l))/h))
               ihi(l) = min(n_SBin(l)-1, FLOOR((rhi(l)-x0(l))/h))
            END DO
            DO i3 = ilo(3),ihi(3)
            DO i2 = ilo(2),ihi(2)
            DO i1 = ilo(1),ihi(1)
               cc(1) = x0(1) + (i1+0.5d0)*h
               cc(2) = x0(2) + (i2+0.5d0)*h
               cc(3) = x0(3) + (i3+0.5d0)*h
               dpl = 0.0d0
               DO l = 1,3
                  dpl = dpl + ppiclf_wall_n(l,w)
     >                       *(cc(l) - ppiclf_wall_c(l,w))
               END DO
               IF (abs(dpl) .GT. rmax) CYCLE
               sb = i1 + n_SBin(1)*i2 + n_SBin(1)*n_SBin(2)*i3
               IF (ipass .EQ. 1) THEN
                  wsb_start(sb+1) = wsb_start(sb+1) + 1
               ELSE
                  wsb_list(wsb_start(sb)) = w
                  wsb_start(sb) = wsb_start(sb) + 1
               END IF
            END DO
            END DO
            END DO
         END DO
         IF (ipass .EQ. 1) THEN
            ! Running sum: wsb_start(sb) = first entry of subbin sb
            wsb_start(0) = 1
            DO sb = 1,tot_SBin
               wsb_start(sb) = wsb_start(sb) + wsb_start(sb-1)
            END DO
            n = wsb_start(tot_SBin) - 1
            IF (ALLOCATED(wsb_list)) DEALLOCATE(wsb_list)
            ALLOCATE(wsb_list(max(n,1)))
         ELSE
            ! Pass 2 advanced each start to the next subbin's; shift
            DO sb = tot_SBin,1,-1
               wsb_start(sb) = wsb_start(sb-1)
            END DO
            wsb_start(0) = 1
         END IF
      END DO

      wsb_tot   = tot_SBin
      wsb_h     = h
      wsb_nwall = ppiclf_nwall
      DO l = 1,3
         wsb_n(l)  = n_SBin(l)
         wsb_x0(l) = x0(l)
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_ClosestPtTri(p,a,b,c,q,ireg)
!
! Closest point q on triangle (a,b,c) to p (Ericson, Real-Time
! Collision Detection, 5.1.5). ireg = 0: q is inside the triangle,
! 1: on an edge, 2: a vertex.
!
      IMPLICIT NONE
!
      REAL*8 p(3), a(3), b(3), c(3), q(3)
      INTEGER*4 ireg
!
      REAL*8 ab(3), ac(3), ap(3), bp(3), cp(3), d1, d2, d3, d4, d5, d6,
     >       va, vb, vc, v, w, denom
      INTEGER*4 l
!
      DO l = 1,3
         ab(l) = b(l) - a(l)
         ac(l) = c(l) - a(l)
         ap(l) = p(l) - a(l)
         bp(l) = p(l) - b(l)
         cp(l) = p(l) - c(l)
      END DO
      d1 = ab(1)*ap(1) + ab(2)*ap(2) + ab(3)*ap(3)
      d2 = ac(1)*ap(1) + ac(2)*ap(2) + ac(3)*ap(3)
      IF (d1 .LE. 0.0d0 .AND. d2 .LE. 0.0d0) THEN
         q = a
         ireg = 2
         RETURN
      END IF
      d3 = ab(1)*bp(1) + ab(2)*bp(2) + ab(3)*bp(3)
      d4 = ac(1)*bp(1) + ac(2)*bp(2) + ac(3)*bp(3)
      IF (d3 .GE. 0.0d0 .AND. d4 .LE. d3) THEN
         q = b
         ireg = 2
         RETURN
      END IF
      vc = d1*d4 - d3*d2
      IF (vc .LE. 0.0d0 .AND. d1 .GE. 0.0d0 .AND. d3 .LE. 0.0d0) THEN
         v = d1/(d1 - d3)
         q = a + v*ab
         ireg = 1
         RETURN
      END IF
      d5 = ab(1)*cp(1) + ab(2)*cp(2) + ab(3)*cp(3)
      d6 = ac(1)*cp(1) + ac(2)*cp(2) + ac(3)*cp(3)
      IF (d6 .GE. 0.0d0 .AND. d5 .LE. d6) THEN
         q = c
         ireg = 2
         RETURN
      END IF
      vb = d5*d2 - d1*d6
      IF (vb .LE. 0.0d0 .AND. d2 .GE. 0.0d0 .AND. d6 .LE. 0.0d0) THEN
         w = d2/(d2 - d6)
         q = a + w*ac
         ireg = 1
         RETURN
      END IF
      va = d3*d6 - d5*d4
      IF (va .LE. 0.0d0 .AND. (d4-d3) .GE. 0.0d0
     >                  .AND. (d5-d6) .GE. 0.0d0) THEN
         w = (d4 - d3)/((d4 - d3) + (d5 - d6))
         q = b + w*(c - b)
         ireg = 1
         RETURN
      END IF
      denom = 1.0d0/(va + vb + vc)
      v = vb*denom
      w = vc*denom
      q = a + ab*v + ac*w
      ireg = 0

      RETURN
      END
!-----------------------------------------------------------------------
      SUBROUTINE ppiclf_solve_WallContacts(i)
!
! Contacts of particle i with the wall triangles and their forces.
!
! A contact with the inside of a triangle has the wall normal as its
! direction; one with an edge or vertex points from the centre to
! the closest point. A flat wall made of several triangles must act
! as one surface, so a contact is dropped when it repeats one already
! kept: an inside contact on the same plane, or an edge/vertex
! contact lying on the plane of a kept inside contact (the shared
! edge of two triangles) or at the same point as a kept one.
!
! A particle whose centre has crossed a one-sided wall stays in
! contact with it, and is pushed back, as long as the triangle is in
! its subbin's list (about ppiclf_nndist, at least 4 diameters,
! behind the wall). This happens when the impact is fast compared
! with the contact time ksp*dt: the overlap is about v*ksp*dt/pi.
! One-sided walls therefore assume that the solid behind them is at
! least ppiclf_nndist thick, which holds for the boundary of the
! fluid domain.
!
      USE ppiclf_wallbin
!
      IMPLICIT NONE
!
      INCLUDE "PPICLF"
!
! Input:
!
      INTEGER*4 i
!
! Internal:
!
      INTEGER*4 LCON
      PARAMETER (LCON = 64)
      REAL*8 p(3), a(3), b(3), c(3), q(3), rnf(3), rad, s, d, dlt, tol,
     >       cq(3,LCON), cnf(3,LCON), crn(3,LCON), cdl(LCON), dot
      INTEGER*4 creg(LCON), ncon, isb(3), sb, k, w, l, ireg, m, n,
     >          iord(LCON), nkeep, m1
      LOGICAL keep(LCON), ppiclf_solve_PointBefore
      EXTERNAL ppiclf_solve_PointBefore
!
      IF (wsb_tot .LE. 0) RETURN

      DO l = 1,3
         p(l) = ppiclf_cp_map(l,i)
         isb(l) = FLOOR((p(l) - wsb_x0(l))/wsb_h)
         IF (isb(l) .LT. 0 .OR. isb(l) .GE. wsb_n(l)) RETURN
      END DO
      sb  = isb(1) + wsb_n(1)*isb(2) + wsb_n(1)*wsb_n(2)*isb(3)
      rad = 0.5d0*ppiclf_cp_map(PPICLF_LRS+PPICLF_R_JDP,i)
      tol = 1.0d-8*rad

      ncon = 0
      DO k = wsb_start(sb),wsb_start(sb+1)-1
         w = wsb_list(k)
         DO l = 1,3
            a(l)   = ppiclf_wall_c(l  ,w)
            b(l)   = ppiclf_wall_c(l+3,w)
            c(l)   = ppiclf_wall_c(l+6,w)
            rnf(l) = ppiclf_wall_n(l,w)
         END DO
         CALL ppiclf_solve_ClosestPtTri(p,a,b,c,q,ireg)
         ! Signed distance of the centre from the triangle's plane
         s = (p(1)-a(1))*rnf(1) + (p(2)-a(2))*rnf(2)
     >     + (p(3)-a(3))*rnf(3)

         IF (ppiclf_wall_side(w) .EQ. 0 .AND. s .LT. 0.0d0) THEN
            ! Two-sided wall: the side the particle is on is the fluid
            DO l = 1,3
               rnf(l) = -rnf(l)
            END DO
            s = -s
         END IF

         ! A closest point on an edge or vertex that is the centre's
         ! own projection on the plane (centre right above or below
         ! the triangle's boundary) is an inside contact: behind a
         ! one-sided wall the edge/vertex branch would skip it
         IF (ireg .NE. 0) THEN
            d = 0.0d0
            DO l = 1,3
               d = d + (p(l) - q(l) - s*rnf(l))**2
            END DO
            IF (sqrt(d) .LT. tol) ireg = 0
         END IF

         IF (ireg .EQ. 0) THEN
            dlt = rad - s
            DO l = 1,3
               crn(l,ncon+1) = -rnf(l)
            END DO
         ELSE
            ! Behind the plane but outside the triangle: left to the
            ! triangle the centre projects into
            IF (s .LT. 0.0d0) CYCLE
            d = sqrt((q(1)-p(1))**2 + (q(2)-p(2))**2 + (q(3)-p(3))**2)
            IF (d .LE. 0.0d0) CYCLE
            dlt = rad - d
            DO l = 1,3
               crn(l,ncon+1) = (q(l) - p(l))/d
            END DO
         END IF
         IF (dlt .LE. 0.0d0) CYCLE

         ncon = ncon + 1
         IF (ncon .GT. LCON)
     >   CALL ppiclf_exittr('Too many wall contacts for a particle$'
     >                     ,0.d0,ncon)
         DO l = 1,3
            cq(l,ncon)  = q(l)
            cnf(l,ncon) = rnf(l)
         END DO
         cdl(ncon)  = dlt
         creg(ncon) = ireg
      END DO

      IF (ncon .EQ. 0) RETURN

      ! Inside contacts first: drop one on the plane of a kept one
      DO m = 1,ncon
         keep(m) = .FALSE.
         IF (creg(m) .NE. 0) CYCLE
         keep(m) = .TRUE.
         DO n = 1,m-1
            IF (.NOT. keep(n) .OR. creg(n) .NE. 0) CYCLE
            dot = cnf(1,m)*cnf(1,n) + cnf(2,m)*cnf(2,n)
     >          + cnf(3,m)*cnf(3,n)
            IF (dot .GT. 1.0d0 - 1.0d-10 .AND.
     >          abs((cq(1,m)-cq(1,n))*cnf(1,n)
     >            + (cq(2,m)-cq(2,n))*cnf(2,n)
     >            + (cq(3,m)-cq(3,n))*cnf(3,n)) .LT. tol) THEN
               keep(m) = .FALSE.
               EXIT
            END IF
         END DO
      END DO
      ! Edge and vertex contacts
      DO m = 1,ncon
         IF (creg(m) .EQ. 0) CYCLE
         keep(m) = .TRUE.
         DO n = 1,ncon
            IF (.NOT. keep(n) .OR. n .EQ. m) CYCLE
            IF (creg(n) .EQ. 0) THEN
               IF (abs((cq(1,m)-cq(1,n))*cnf(1,n)
     >               + (cq(2,m)-cq(2,n))*cnf(2,n)
     >               + (cq(3,m)-cq(3,n))*cnf(3,n)) .LT. tol) THEN
                  keep(m) = .FALSE.
                  EXIT
               END IF
            ELSE IF (n .LT. m) THEN
               IF (abs(cq(1,m)-cq(1,n)) .LT. tol .AND.
     >             abs(cq(2,m)-cq(2,n)) .LT. tol .AND.
     >             abs(cq(3,m)-cq(3,n)) .LT. tol) THEN
                  keep(m) = .FALSE.
                  EXIT
               END IF
            END IF
         END DO
      END DO

      ! Apply the kept contacts in order of their contact points, so
      ! the sum does not depend on how the walls were numbered (the
      ! order of the ranks that sent them)
      nkeep = 0
      DO m = 1,ncon
         IF (.NOT. keep(m)) CYCLE
         nkeep = nkeep + 1
         iord(nkeep) = m
         DO n = nkeep,2,-1
            IF (.NOT. ppiclf_solve_PointBefore(cq(1,iord(n)),
     >                                    cq(1,iord(n-1)))) EXIT
            m1 = iord(n)
            iord(n) = iord(n-1)
            iord(n-1) = m1
         END DO
      END DO
      DO n = 1,nkeep
         m = iord(n)
         CALL ppiclf_user_EvalWallContact(i,ppiclf_cp_map(1,i)
     >        ,ppiclf_cp_map(1+PPICLF_LRS,i),crn(1,m),cdl(m))
      END DO

      RETURN
      END
!-----------------------------------------------------------------------
      LOGICAL FUNCTION ppiclf_solve_PointBefore(a,b)
!
! Lexicographic order of points a and b (x, then y, then z)
!
      IMPLICIT NONE
      REAL*8 a(3), b(3)
!
      IF (a(1) .NE. b(1)) THEN
         ppiclf_solve_PointBefore = a(1) .LT. b(1)
      ELSE IF (a(2) .NE. b(2)) THEN
         ppiclf_solve_PointBefore = a(2) .LT. b(2)
      ELSE
         ppiclf_solve_PointBefore = a(3) .LT. b(3)
      END IF

      RETURN
      END
