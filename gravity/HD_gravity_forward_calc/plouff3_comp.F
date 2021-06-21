!==========================================================================
!
!   plouff_comp.F
!   Computational routine for plouff3.F (Fortran version of plouff1.m).
!   This version uses the input absolute positions of Pp and Pobs.
!
!     This is a MEX-file for MATLAB.
!     by Atsuhiro Muto, Dept. of Geosciences,
!     The Pennsylvania State University
!     aum34@psu.edu, atsumuto@gmail.com
!     Last updated 07/01/2011
!==========================================================================

!#include "fintrf.h"

!==========================================================================
      SUBROUTINE plouff3_comp(g,Pobs,Pp,zp,zp2,rho1,rho2,dx,dy,m,
     +                       n,l1,l2,dist)
!==========================================================================



!--------------------------------------------------------------------------
! initialization
      integer m, n, a, i, j, k, l1, b, id, si, sj, s, c, l2

! input variables
      real*8 Pobs(n,2), Pp(m,4), zp(m,l1+1), rho1(l1), dx, dy
      real*8 dist, zp2(m,l2+1), rho2(l2)

! output variables
      real*8 g(n)

! local variables
      real*8 R, U(l1+1,4), gl(l1), xp(2), yp(2), xpp(2), ypp(2)
      real*8 gg, zpt(l1+1,m), Q, P
      real*8 R2, U2(l2+1,4), gl2(l2), Q2, zpt2(l2+1,m)
!--------------------------------------------------------------------------


      ! gravitational constant
      gg = 6.673e-11
      
      zpt = transpose(zp)
      zpt2 = transpose(zp2)
      g = 0.0d8

!--------------------------------------------------------------------------
! Loop over data points
      do c=1,n
        !----------------------------------------------------------------
        ! Loop over prisms
        do a=1,m
          xp(1) = Pp(a,1)-Pobs(c,1)
          xp(2) = Pp(a,2)-Pobs(c,1)

          if(xp(1).LT.0 .AND. xp(1).LT.xp(2) .AND. xp(1)*xp(2).GE.0)then
!            xp(1:2:-1) = xp(1:2)
            xpp(1) = xp(2)
            xpp(2) = xp(1)
            xp = abs(xpp)
          endif

          yp(1) = Pp(a,3)-Pobs(c,2)
          yp(2) = Pp(a,4)-Pobs(c,2)

          if(yp(1).LT.0 .AND. yp(1).LT.yp(2) .AND. yp(1)*yp(2).GE.0)then
!            yp(1:2:-1) = yp(1:2)
            ypp(1) = yp(2)
            ypp(2) = yp(1)
            yp = abs(ypp)
          endif

          P = sqrt(xp(1)**2+yp(1)**2)
          if(P.GT.dist)then
            gl = 0
            gl2 = 0
            go to 10
          endif

          si = 1
          sj = 1

          !-----------------------------------------------------------
          ! Loop over corners of the prism
          id = 1
          do i=1,2
            si = si*(-1)
            do j=1,2
              sj = sj*(-1)
              s = si*sj
              do k=1,l1+1
                R = sqrt(xp(i)**2+yp(j)**2+zpt(k,a)**2)
                Q = atan(xp(i)*yp(j)/zpt(k,a)/R)
                U(k,id) = s * (zpt(k,a) * Q
     +                       - xp(i)*log(R+yp(j))
     +                       - yp(j)*log(R+xp(i)))
              enddo
              
              do k=1,l2+1
                R2 = sqrt(xp(i)**2+yp(j)**2+zpt2(k,a)**2)
                Q2 = atan(xp(i)*yp(j)/zpt2(k,a)/R2)
                U2(k,id) = s * (zpt2(k,a) * Q2
     +                       - xp(i)*log(R2+yp(j))
     +                       - yp(j)*log(R2+xp(i)))
              enddo

              id = id+1
            enddo
          enddo
          !-----------------------------------------------------------


          !-------------------------------------------------
          ! gravity effect of each layer in the prism
          do b=1,l1
            gl(b) = rho1(b) * (sum(U(b,:))*(-1) +
     +                sum(U(b+1,:)))
          enddo

          do b=1,l2
            gl2(b) = rho2(b) * (sum(U2(b,:))*(-1) +
     +                sum(U2(b+1,:)))
          enddo
          !-------------------------------------------------
        

          ! gravity effect of all layers in the prism
10          g(c) =  g(c) + gg*sum(gl) + gg*sum(gl2) 

        enddo
        !----------------------------------------------------------------

      enddo
!--------------------------------------------------------------------------

      return
      end