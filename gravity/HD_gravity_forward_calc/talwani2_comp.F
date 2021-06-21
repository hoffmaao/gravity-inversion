!==========================================================================
!
!   talwani_comp.F
!   Computational routine for talwani.F (Fortran version of twodpoly2.m).
!
!     This is a MEX-file for MATLAB.
!     by Atsuhiro Muto, Dept. of Geosciences,
!     The Pennsylvania State University
!     aum34@psu.edu, atsumuto@gmail.com
!     Last updated 05/06/2014
!==========================================================================

!#include "fintrf.h"

!==========================================================================
	  SUBROUTINE talwani2_comp(g,xobs,xp,zp,drho,mx,nx,nv)
!==========================================================================



!--------------------------------------------------------------------------
! initialization
      integer l, mx, nx, nv, ip, im, iv, i

! input variables
      real*8 xobs(mx), xp(nx), zp(nx), drho

! output variables
      real*8 g(mx)

! local variables
      real*8 U(nv), dx, dz
      real*8 zint, xint, theta, gg
!--------------------------------------------------------------------------


      ! gravitational constant
      gg = 6.673e-11


      ! -------------------------------------------------------------------------
      ! loop over observation points
      !do ip=1,l
    
        ! -------------------------------------------------------------
        ! loop over observation points  
        do im=1,mx   

          U = 0.0d8
		  
          ! --------------------------------------------------
          ! loop over each side of polygon
          do iv=1,nv   
            dx = (xp(iv+1)-xp(iv))/nx
            dz = (zp(iv+1)-zp(iv))/nx            

            ! -----------------------------------------
            ! Integrate over side    
            do i=1,nx
              zint = zp(iv)+i*dz
              xint = xp(iv)+i*dx
              theta = atan(abs(xint-xobs(im))/zint) 

              if(xobs(im) .GE. xint)then
                theta = theta*(-1.0)
              endif

              U(iv) = U(iv) + theta*dz
            enddo
            ! -----------------------------------------
            
          enddo
          ! --------------------------------------------------
        
          g(im) = 2*gg*drho*(sum(U))

          enddo
          ! -------------------------------------------------------------
    
      !enddo
      ! -------------------------------------------------------------------------

      !g = sum(dgg,dim=2)*1e5

!--------------------------------------------------------------------------
      return
      end