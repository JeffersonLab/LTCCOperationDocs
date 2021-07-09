
CDECK  ID>, CCWCPAR.
      SUBROUTINE CCWCPAR(iwc,az,ar,S,P,wr)
c
c_begin_doc
c  RCS ID string
c  $Id: ccwcpar.F,v 1.2 1997/04/24 15:38:16 holtrop Exp $
c
c  Documentation for subroutine CCWCPAR(iwc,az,ar,S,P,wr)
c
C---------------------------------------------------------------------
C-
C-   Purpose and Methods : Calculates WC parameters
C-
C-   Inputs  :  iwc   - number of WC
C-              az,ar - radius of WC (ar) at the distance az from
C-                      the WC window
C-              /CCP/ - parameters of WC window and PMT
C-
C-   Outputs :  S     - coordinates of the WC center in the
C-                      Segment median plane R.S.
C-              P     - WC parameters in form of
C-                      P(1)*x**2 + P(2)*y**2 + P(3)*z**2 + 1 = 0
C-                      ( in it's coordinate system )
C-
C-   Controls:
C-
C-   Created     8-MAY-1994   Alexander V. Vlassov
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
c
#include "ccp.inc"
c
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
      integer iwc
      real az, ar, S(3), P(3), wr
c
c
c  Local pre-defined variables:
c  ---------------------------
C----------------------------------------------------------------------
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs rcs
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs rcs
c                                                                         rcs
c  RCS information:                                                       rcs
c                                                                         rcs
      CHARACTER*(132)  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile: ccwcpar.F,v $') 
      PARAMETER (CREVIS = '$Revision: 1.2 $')
      PARAMETER (CSTATE = '$State: Exp $')
      PARAMETER (CDATE = '$Date: 1997/04/24 15:38:16 $')
      PARAMETER (CAUTHO = '$Author: holtrop $')
      DATA CRCSID/ 
     1'$Id: ccwcpar.F,v 1.2 1997/04/24 15:38:16 holtrop Exp $'
     2/ 
c  Module information: 
      CHARACTER*(*)  CRNAME, CRAUTH
      PARAMETER (CRNAME='ccwcpar')                                   
      PARAMETER (CRAUTH='Alexander Vlassov')
c                                                                         rcs
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs rcs
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs rcs
c
c  Local User defined variables:
c  -----------------------------
c
      integer i
      double precision un,a1,a2,a3,a4,a5,a6,a7
      double precision q,b,c
c
      data un/1.0d000/
c_end_var
c
c  executable code for routine CCWCPAR
c------------------------------------------------------------------------------
c
      a1 = un *  az   *  az                                ! R1**2
      a2 = un *(pmt0(1,iwc)-wc0(1,iwc))*(pmt0(1,iwc)-wc0(1,iwc)) +
     &     un *(pmt0(2,iwc)-wc0(2,iwc))*(pmt0(2,iwc)-wc0(2,iwc)) +
     &     un *(pmt0(3,iwc)-wc0(3,iwc))*(pmt0(3,iwc)-wc0(3,iwc))
c                                                          ! R2**2
      a3 = un * wcr(iwc)* wcr(iwc)                         ! L1**2
      a4 = un *  ar   *  ar                                ! L2**2
      a5 = un *pmtr(iwc)*pmtr(iwc)                         ! L3**2
c
      a6 = a3 - a4
      a7 = a3 - a5
c
      q  = (a1*a7 - a2*a6) / ( (a6*sqrt(a2) - a7*sqrt(a1))*2.)
c
      b  = a3 + q*q*a6/(a2 + q*sqrt(a2)*2.)                ! a**2
      c  = b*(q + sqrt(a2))*(q + sqrt(a2))/(b - a5)        ! b**2
c
      a3 = sqrt(b)
      a4 = sqrt(c)
c
c  Ellipzoid parameters:
c
      P(1) = -un/b
      P(2) = -un/b
      P(3) = -un/c
c
c  Ellipzoid center coordinates:
c
      do i = 1,3
c
        S(i) = un*wc0(i,iwc) + q*(wc0(i,iwc) - pmt0(i,iwc))/sqrt(a2)
c
      end do
c
c  Winston cone enterance window
c
      a5 = b*(un - q*q/c)
      wr = sqrt(a5)
c
C----------------------------------------------------------------------
  999 RETURN
      END
