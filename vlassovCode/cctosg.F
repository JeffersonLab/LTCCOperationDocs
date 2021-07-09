      SUBROUTINE CCTOSG(ncc,x_in,x_fn)
c     ================================
c
c  Documentation for subroutine CCTOSG(ncc,x_in,x_fn)
c
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : To transfer point x_in in CC segment R.S.
C-                         to SGCC R.S. (x_fn)
C-
C-   Inputs  :  ncc     -  Segment #
C-              x_in(3) -  Point coordinates in the CC Segment R.S.
C-              /CCP/ is used to extract CC parameters
C-
C-   Outputs :  x_fn(3) -  Point coordinates in the SGCC volume  R.S.
C-   Comments:  x_in and x_fn could be the same
C-
C-   Created    17-JUN-1994   Alexander V. Vlassov
C-   Modified   21-AUG-1995   Alexander V. Vlassov
C-
C----------------------------------------------------------------------
c
c_end_doc
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
#include "ccp.inc"
c
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
      integer ncc
      real x_in(3),x_fn(3)
c
c  Local pre-defined variables:
c  ---------------------------
      integer i,j
      DOUBLE PRECISION T(3,3),DEGRAD
      DOUBLE PRECISION tet,a(3),un
c
      DATA un /1.0d000/
      PARAMETER (DEGRAD=0.0174532925199432958)
c_end_var
c
c  executable code for routine CCTOSG
c------------------------------------------------------------------------------
c
      if(ncc.le.0.or.ncc.gt.18) then
        print*, 'Wrong CC number! (1-18):',ncc
        stop
      end if
c
C----------------------------------------------------------------------
C
      tet = DEGRAD*th(ncc)
C
      T(1,1) = 0.
      T(2,1) = un
      T(3,1) = 0.
c
      T(1,2) = sin(tet)
      T(2,2) = 0.
      T(3,2) = cos(tet)
c
      T(1,3) = cos(tet)
      T(2,3) = 0.
      T(3,3) = -sin(tet)
c
      do j = 1,3
c
        a(j) = 0.
        do i = 1,3
          a(j) = a(j) + T(j,i)*x_in(i)
        end do
      end do
c
      do j = 1,3
        x_fn(j) = a(j)
      end do
c
      RETURN
      END
