!==========================================================================
!
!   plouff.F
!   Gateway routine for plouff_comp.F (Fortran version of 
!   plouff8.m).
!   This version uses the input absolute positions of Pp and Pobs.
!
!     This is a MEX-file for MATLAB.
!     by Atsuhiro Muto, Dept. of Geosciences,
!     The Pennsylvania State University
!     aum34@psu.edu, atsumuto@gmail.com
!     Last updated 10/22/2014
!==========================================================================

#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

!-----------------------------------------------------------------------
      mwPointer plhs(*), prhs(*)
      mwPointer mxCreateDoubleMatrix, mxGetPr
      mwPointer Pobs_pr, Pp_pr, zp_pr, rho1_pr, dx_pr, dy_pr, g_pr
      mwPointer dist_pr, zp2_pr, rho2_pr
!-----------------------------------------------------------------------

      integer nlhs, nrhs
      mwsize mxGetM, mxGetN
      mwsize m, n, l1, l2
      real*8 Pobs(1000000), Pp(1000000), zp(750000), rho1(5), dx, dy
      real*8 g(500000), dist, zp2(750000), rho2(5)

! get the size of the input matrices
      n = mxGetM(prhs(1))
      m = mxGetM(prhs(2))
      l1 = mxGetM(prhs(5))
      l2 = mxGetM(prhs(6))

! create matrix for the return argument
      plhs(1) = mxCreateDoubleMatrix(n,1,0)
      g_pr = mxGetPr(plhs(1))

! create matrix for the input argument
      Pobs_pr = mxGetPr(prhs(1))
      Pp_pr = mxGetPr(prhs(2))
      zp_pr = mxGetPr(prhs(3))
      zp2_pr = mxGetPr(prhs(4))
      rho1_pr = mxGetPr(prhs(5))
      rho2_pr = mxGetPr(prhs(6))
      dx_pr = mxGetPr(prhs(7))
      dy_pr = mxGetPr(prhs(8))
      dist_pr = mxGetPr(prhs(9))

! load data into Fortran arrays
      call mxCopyPtrToReal8(Pobs_pr,Pobs,n*2)
      call mxCopyPtrToReal8(Pp_pr,Pp,m*4)
      call mxCopyPtrToReal8(zp_pr,zp,m*(l1+1))
      call mxCopyPtrToReal8(zp2_pr,zp2,m*(l2+1))
      call mxCopyPtrToReal8(rho1_pr,rho1,l1)
      call mxCopyPtrToReal8(rho2_pr,rho2,l2)
      call mxCopyPtrToReal8(dx_pr,dx,1)
      call mxCopyPtrToReal8(dy_pr,dy,1)
      call mxCopyPtrToReal8(dist_pr,dist,1)

! call computational subroutine
      call plouff3_comp(g,Pobs,Pp,zp,zp2,rho1,rho2,dx,dy,m,
     +                  n,l1,l2,dist)

! load output into Matlab array
      call mxCopyReal8ToPtr(g,g_pr,n)

      return
      end