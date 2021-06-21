!==========================================================================
!
!   talwani.F
!   Gateway routine for talwani_comp.F (Fortran version of 
!   twodpoly3.m).
!
!     This is a MEX-file for MATLAB.
!     by Atsuhiro Muto, Dept. of Geosciences,
!     The Pennsylvania State University
!     aum34@psu.edu, atsumuto@gmail.com
!     Last updated 05/06/2014
!==========================================================================

#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

!-----------------------------------------------------------------------
      mwPointer plhs(*), prhs(*)
      mwPointer mxCreateDoubleMatrix, mxGetPr
      mwPointer xobs_pr, xp_pr, zp_pr, drho_pr
      mwPointer g_pr
!-----------------------------------------------------------------------

      integer nlhs, nrhs, nv 
      mwsize mxGetM, mxGetN
      mwsize mx, nx
      real*8 xobs(10000), xp(10000), zp(10000), drho(1)
      real*8 g(10000)

! get the size of the input matrices
      mx = mxGetM(prhs(1))
      nx = mxGetM(prhs(2))
      !l = mxGetM(prhs(4))
	  nv = nx-1

! create matrix for the return argument
      plhs(1) = mxCreateDoubleMatrix(mx,1,0)
      g_pr = mxGetPr(plhs(1))

! create matrix for the input argument
      xobs_pr = mxGetPr(prhs(1))
      xp_pr = mxGetPr(prhs(2))
      zp_pr = mxGetPr(prhs(3))
      drho_pr = mxGetPr(prhs(4))

! load data into Fortran arrays
      call mxCopyPtrToReal8(xobs_pr,xobs,mx)
      call mxCopyPtrToReal8(xp_pr,xp,nx)
      call mxCopyPtrToReal8(zp_pr,zp,nx)
      call mxCopyPtrToReal8(drho_pr,drho,1)

! call computational subroutine
      call talwani2_comp(g,xobs,xp,zp,drho,mx,nx,nv)

! load output into Matlab array
      call mxCopyReal8ToPtr(g,g_pr,mx)

      return
      end