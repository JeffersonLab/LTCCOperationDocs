      SUBROUTINE CCDEIN
c     =================
c  Documentation for subroutine CCDEIN
c
C     CC GEOMETRICAL INITIALIZATION
C
C    Purpose and methods : To read geometry parameters for CC
C                          from the file, calculate some other parameters
C                          needed for CC and store common blocks /CCP/
C                          and /CCPW/ (see ccp.inc file)
C----------------------------------------------------------------------
C
C    Side effects : Common /CCP/ initialization.
C
C-
C-   Created    16-APR-1994   Alexander V. Vlassov
C-   Modified   18-AUG-1995   Alexander V. Vlassov
C-   Modified   06-OCT-2009   Alexander V. Vlassov
C
c_end_doc
c
C----------------------------------------------------------------------
      IMPLICIT NONE
      SAVE
C----------------------------------------------------------------------
c
c_begin_inc
c  include files :
c  ---------------------
#include "gsimpar.inc"
#include "ccp.inc"
#include "ccfrmp.inc"
#include "gconsp.inc"
#include "ccgeom12.inc"
c  ---------------------
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
c  Local pre-defined variables:
c  ---------------------------
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs
c                                                                        
c  RCS information:                                                      
c                                                                        
      CHARACTER*(132)  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile: ccdein.F,v $') 
      PARAMETER (CREVIS = '$Revision: 1.3 $')
      PARAMETER (CSTATE = '$State: Exp $')
      PARAMETER (CDATE = '$Date: 1999/03/24 20:30:07 $')
      PARAMETER (CAUTHO = '$Author: fklein $')
      DATA CRCSID/ 
     1'$Id: ccdein.F,v 1.3 1999/03/24 20:30:07 fklein Exp $'
     2/ 
c  Module information: 
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='ccdein')                                   
      PARAMETER (CRAUTH='Alexander Vlassov')
c                                                                        
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs
c rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rcs rsc rcs rcs
c
C  External references:
      real vdot, vdist
      external vdot, vdist
c
c  Local User defined variables:
c  -----------------------------
      INTEGER i,j,ncc,k
      INTEGER LUNG,IRES,ifirst
c
      character*80 filename
c
      REAL vlar,r(5),t(4),s,d,ve(3),xv(3),xu(3),ul(3),teta
      REAL tm(3),pm(3), xmi,xma
      real Z(3), al, wcz(18), wcxy(18), xms(3,6), pms(3,4,4)
c
      DATA vlar /1000001./
      DATA Z /0., 0., 0./
      data ifirst /0/
C
C  =================   END OF DECLARATIONS  ====================
C
c_end_var
c
c  executable code for routine CCDEIN
C
***************************************************8
C
C --- Input parameters for elliptical and hyperbolic mirrors
C
      LUNG = 21
      print *, ' CCGEOM12.dat lun', LUNG
      open (lung,file='ccgeom12.dat',status='old',err=2)
C
      DO i=1,18
        READ (LUNG,*,end=2) r
      END DO
c
      DO i=1,18
        READ (LUNG,*,end=2) r
      END DO
c
c     new way to define mirror parameters :
c
      do i = 1,18
        do j = 1,5
          pe0(j,i) = ppe0(j,i)
          ph0(j,i) = pph0(j,i)
        end do
      end do
c
      do i = 1,18
        read(LUNG,*,end=2) j,wc0(1,i),wc0(2,i),pmt0(1,i),pmt0(2,i),
     &  wcr(i),pmtr(i),dplwc(i),wcang(i)
        wc0 (3,i) = 0.
        pmt0(3,i) = 0.
      end do
c
C --- Input parameters for PMT wall
c
      DO i=1,18
        READ (LUNG,*) t
        x11(i)= t(1)
        y11(i)= t(2)
        x12(i)= t(3)
        y12(i)= t(4)
      END DO
c
C --- Input parameters for coil plane
c
      DO i=1,18
        READ (LUNG,*) t
        x21(i)= t(1)
        y21(i)= t(2)
        x22(i)= t(3)
        y22(i)= t(4)
      END DO
c
      DO i=1,18
        READ (LUNG,*) hwde(i), hwdh(i), wcz(i), wcxy(i),
     &  cmsx1(2,i), cmsx1(3,i), cmsx2(2,i), cmsx2(3,i),
     &  hmssz(1,i),  hmssz(2,i),  hmssz(3,i)
        cmsx1(1,i) = 0.
        cmsx2(1,i) = 0.
      END DO
c
      DO i=1,18
        READ (LUNG,*) j, thmin(i), th(i), scrnp1(1,i),
     &  scrnp1(2,i), scrnp2(1,i), scrnp2(2,i)
        scrnp1(3,i) = 0.  
        scrnp2(3,i) = 0.  
      END DO
c
      read (lung,*) j,thmin(19)
c
C---- end of reading -----
c
      CLOSE (LUNG)
C
c---   PARAMETER calculations   ----
c
      entry CCPARM
c     ============
c
C     Winstone Cone parameters
C
      DO i=1,18
C
C ---  Parameters for Winston Cone slope
C
        s = (wc0(1,i) - pmt0(1,i))*(wc0(1,i) - pmt0(1,i))
     &    + (wc0(2,i) - pmt0(2,i))*(wc0(2,i) - pmt0(2,i))
        d = ( pmt0(1,i) - wc0(1,i) )/sqrt(s)
        csa0(i)=asin(d)
c-
        call ccwcpar(i,wcz(i),wcxy(i),sw0(1,i),pw0(1,i),wcer(i))
c-
      END DO
C
c-    Let's redefine x22,y22
c     ======================
c-
      do i = 1,18
c-
        d = csa0(i)
        x22(i) = wc0(1,i) + wcr(i)*cos(d)
        y22(i) = wc0(2,i) - wcr(i)*sin(d)
c-
c-    Let's define xp1,yp1,xp2,yp2
c     ============================
c
        xp1(i) = pmt0(1,i) - pmtr(i)*cos(d)
        yp1(i) = pmt0(2,i) + pmtr(i)*sin(d)
        xp2(i) = pmt0(1,i) + pmtr(i)*cos(d)
        yp2(i) = pmt0(2,i) - pmtr(i)*sin(d)
c-
      end do
c
C --- Middle plane  ----
c
      p00(1)=vlar
      p00(2)=0.
      p00(3)=0.
c
C --- Up-Front line (wall of PMT)
c
      DO i=1,18
        d=x11(i)*y12(i)-x12(i)*y11(i)
        p10(1,i)=(y11(i)-y12(i))/d
        p10(2,i)=(x12(i)-x11(i))/d
        p10(3,i)=0.
      END DO
c
C --- Back-bottom (under PMT) line
c
      DO i=1,18
        d=x21(i)*y22(i)-x22(i)*y21(i)
        p20(1,i)=(y21(i)-y22(i))/d
        p20(2,i)=(x22(i)-x21(i))/d
        p20(3,i)=0.
      END DO
c
C --- Input in Winston Cone (WC window plane)
c
      DO i=1,18
        d = 0.
        hwcpmt(i) = 0.
        do j = 1,3
          ve(j) = pmt0(j,i) - wc0(j,i)
          d = d + wc0(j,i)*ve(j)
          hwcpmt(i) = hwcpmt(i) + ve(j)*ve(j)
        end do
        hwcpmt(i) = sqrt(hwcpmt(i))
        do j = 1,3
          p30(j,i) = -ve(j)/d
        end do
      END DO
c
C --- PhotoMultipier Tube window plane (PMT plane)
c
      DO i=1,18
        d = 0.
        do j = 1,3
          ve(j) = pmt0(j,i) - wc0(j,i)
          d = d + pmt0(j,i)*ve(j)
        end do
        do j = 1,3
          p40(j,i) = -ve(j)/d
        end do
      END DO
c
C --- bottom black wall around hyperbolic mirror
c
      DO i=1,18
        p50(1,i)=0.
        p50(2,i)=-1./y21(i)
        p50(3,i)=0.
      END DO
c
C --- close plane(mirror for #18),through it light comes to next CC sect.
c
      DO i=1,18
        d=(th(i)-thmin(i))*DEGRAD
        p60(1,i)=0.
        p60(2,i)=-vlar
        p60(3,i)=-vlar/tan(d)
      END DO
c
C --- fare plane(mirror for # 1),through it light comes to prev. CC sect.
c
      DO i=1,18
        d=(thmin(i+1)-th(i))*DEGRAD
        p70(1,i)=0.
        p70(2,i)=-vlar
        p70(3,i)=vlar/tan(d)
      END DO
c
C --- X-loc., where elliptical mirrow should be replaced by attached plane
C
      do i = 1,18
         x_pl(i) = 7.5
         if(i.ge.11) x_pl(i) = 10.0
      end do
c
C --- Plane parameters of this attached plane:
C     in form p(1)*x + p(2)*y + p(3)*z + 1 = 0.
C
      xv(2) = 0.
      xv(3) = 0.
c
      ve(1) = 0.
      ve(2) = 1.
      ve(3) = 0.
c
      do i = 1,18
c
c     crossing point :
c
        xv(1) = x_pl(i)
        call CCROSE(xv,ve,pe0(1,i),IRES,S,R)
        if(IRES.ne.1) then   
          print *,'WRONG parameters for CCROSE !!!'
          print *,'xv,ve',xv,ve
          print *,'i,pe0',i,pe0(1,i),pe0(2,i),pe0(3,i),pe0(4,i),pe0(5,i)
          print *,'ires,s,r',ires,s,r
          stop
        end if
C
        d = R(2)
C*******
        elt_pl(1,i) = 2*pe0(1,i)*x_pl(i) + pe0(3,i)*d + pe0(4,i)
        elt_pl(2,i) = 2*pe0(2,i)*d + pe0(3,i)*x_pl(i) + pe0(5,i)
        elt_pl(3,i) = 0.
C*******
        S = elt_pl(1,i)*x_pl(i) + elt_pl(2,i)*d
        elt_pl(1,i) = -elt_pl(1,i) / S
        elt_pl(2,i) = -elt_pl(2,i) / S
C*******
      end do
c
c    Plane parameters for WC bounding planes
c
      do i = 1,18
c
c *****  turn in contra-clockwise direction
c
        pwc10(1,i) =  2.*cos(degrad*wcang(i))/dplwc(i)
        pwc10(2,i) =  2.*sin(degrad*wcang(i))/dplwc(i)
        pwc10(3,i) =  0.
c
        pwc20(1,i) = -2.*cos(degrad*wcang(i))/dplwc(i)
        pwc20(2,i) = -2.*sin(degrad*wcang(i))/dplwc(i)
        pwc20(3,i) =  0.
c
c     Parameters of these planes in CC SEGment system
c
        ncc = i
        d     = degrad*wcang(i)   ! in radians
c
        ul(1) =  cos(d)
        ul(2) =  sin(d)
        ul(3) =  0.
c
c  ******  P80  ******
c
        xu(1) =  -dplwc(i)*cos(d)/2.
        xu(2) =  -dplwc(i)*sin(d)/2.
        xu(3) =   0.
c
c
c   From WC to CC Reference System :
c
        call CCFRWC(ncc,xu,ul,xv,ve)
c
c       Now xv and ve - point on the plane and plane unit vector
c
        s = vdot(xv,ve,3)
        do j = 1,3
          P80(j,i) = -ve(j)/s
        end do
c
c  ******  P90  ******
c
        xu(1) =  dplwc(i)*cos(d)/2.
        xu(2) =  dplwc(i)*sin(d)/2.
        xu(3) =  0.
c
c
c   From WC to CC Reference System :
c
        call CCFRWC(ncc,xu,ul,xv,ve)
c
c       Now xv and ve - point on the plane and plane unit vector
c
        d = vdot(xv,ve,3)
        do j = 1,3
          P90(j,i) = -ve(j)/d
        end do
c
      end do
c------------------------------------------------------------------
c
c    Parameters of the inner surface of the SideWall
c    ( 2.54 cm apart from the Sector Dividing plane )
c
c
c    Plane parameters
c
      xv(1) = sqrt(3.)/2.
      xv(2) = -0.5
      xv(3) =  0.
c
      ve(1) = -2.54 * sqrt(3.)/2.
      ve(2) =  2.54 / 2.
      ve(3) = 0.
c
c
c   TO NEW ( CC_i ) Reference System :
c
      do i = 1,18
        teta = th(i) - 90.
c
        tm(1) = 90.
        tm(2) = 90. - teta
        tm(3) = -teta
c
        pm(1) =  0.
        pm(2) = 90.
        pm(3) = 90.
c
C------------------------------------------------------------
        call CCREFS(ve,xv,tm,pm,Z,xu,ul)
c
c ==============================================================
c
        d = vdot(xu,ul,3)
        do j = 1,3
          PA0(j,i) = -ul(j)/d
        end do
c
      end do
c
C----------------------------------------------------------------
c    Magnetic shield planes definitions
C----------------------------------------------------------------
c
      do i = 1,18
c
c   ****   Point definitions in the R.S. of the "TURNED WC" ***
c  Center of MS :
c
        do j = 1,3
          xms(j,1) = 0.5*(cmsx1(j,i) + cmsx2(j,i))
        end do
        teta = atan((cmsx2(2,i)-cmsx1(2,i))/(cmsx2(3,i)-cmsx1(3,i)))
c
c  Centers of walls
c
        xms(1,2) = xms(1,1)
        xms(2,2) = xms(2,1) - hmssz(2,i)*cos(teta)
        xms(3,2) = xms(3,1) + hmssz(2,i)*sin(teta)
c
        xms(1,3) = xms(1,1) + hmssz(1,i)
        xms(2,3) = xms(2,1)
        xms(3,3) = xms(3,1)
c
        xms(1,4) = xms(1,1)
        xms(2,4) = xms(2,1) + hmssz(2,i)*cos(teta)
        xms(3,4) = xms(3,1) - hmssz(2,i)*sin(teta)
c
        xms(1,5) = xms(1,1) - hmssz(1,i)
        xms(2,5) = xms(2,1)
        xms(3,5) = xms(3,1)
c
c    Center of MS window = lower end of the MS axes
c
        call ucopy(cmsx1(1,i),xms(1,6),3)
c
c  Now transfer these coordinates to the CC R.S. :
c
        xu(1) = 0.
        xu(2) = 0.
        xu(3) = 1.
c
        do j = 1,6
          call CCWCFRS(i,xms(1,j),xms(1,j))
          call CCFRWC(i,xms(1,j),xu,xms(1,j),xv)
        end do
c
c   Plane parameters definitions :
c
c   Planes #1 - #4 :
c
        do j = 1,4
c
          k = j+1
          call vsub(xms(1,k),xms(1,1),ve,3)
          d = -1./vdot(ve,xms(1,k),3)
          call vscale(ve,d,plms0(1,j,i),3)
c
        end do
c
c       SAVE xms coordinates :
c       ====================
        call ucopy(xms,pcms0(1,1,i),18)   ! 18 = 3 * 6
c
      end do
c
C----------------------------------------------------------------
c               /CCFRMP/ parameters definition
C----------------------------------------------------------------
c
      do i = 1,18
c
        ncc = i
        call ccdefpar(i)
c
        xmi = 0.
        xma = x12(i)
        d   = (xma-xmi)/(mirrn - 1.)
c
        xu(2) = 0.
        xu(3) = 0.
 
        xv(1) = 0.
        xv(2) = 1.
        xv(3) = 0.
c ===  Ellips
        DO j = 1,mirrn
c       ==============
c
          xu(1) = xmi + d*(j-1)
          xell(j,i) = xu(1)
          call CCROSE(xu,xv,PE,IRES,S,R)
          yell(j,i) = R(2)
c
        END DO
c
        xu(2) = 700.
        xv(2) = -1.
        xma   = x21(i)
        d   = (xma-xmi)/(mirrn - 1.)
c ===  Hyperbola
        DO j = 1,mirrn
c       ==============
c
          xu(1) = xmi + d*(j-1)
          xhyp(j,i) = xu(1)
          call CCROSH(xu,xv,PH,IRES,S,R)
          yhyp(j,i) = R(2)
c
        END DO
c
c ===  Winston cone  ===
c
        teta  =  csa0(i)
        xv(1) = -cos(teta)
        xv(2) =  sin(teta)
        xv(3) =  0.
        xu(3) =  0.
c
        s = (wc0(1,i) - pmt0(1,i))*(wc0(1,i) - pmt0(1,i))
     &    + (wc0(2,i) - pmt0(2,i))*(wc0(2,i) - pmt0(2,i))
        d = sqrt(s)/(nwcon - 1.)
c
        do j = 1,nwcon
c       ==============
c
          xu(1) = wc0(1,i) + d*(j-1)*sin(teta)
          xu(2) = wc0(2,i) + d*(j-1)*cos(teta)
c
*          s  = d*(j-1)
          s  = vdist(xu,sw0(1,i),3)
          al = -(pw(3)*s*s + 1.)/pw(1)
          al = sqrt(al)
c
          xw1(j,i) = xu(1) + al*xv(1)
          yw1(j,i) = xu(2) + al*xv(2)
c
          xw2(j,i) = xu(1) - al*xv(1)
          yw2(j,i) = xu(2) - al*xv(2)
c
        end do
c
C----------------------------------------------------------------
c    Magnetic shield walls
C----------------------------------------------------------------
c
c       RESTORE xms coordinates :
c       =========================
        call ucopy(pcms0(1,1,i),xms,18)   ! 18 = 3 * 6
c
        call vsub(xms(1,6),xms(1,1),ve,3)
        call vsub(xms(1,2),xms(1,1),xu,3)
        call vsub(xms(1,3),xms(1,1),xv,3)
c
c   Wall #1
c
        call vadd(xms(1,2),ve,        ul,3)
        call vadd(      ul,xv,pms(1,1,1),3)
        call vsub(      ul,xv,pms(1,2,1),3)
        call vsub(xms(1,2),ve,        ul,3)
        call vsub(      ul,xv,pms(1,3,1),3)
        call vadd(      ul,xv,pms(1,4,1),3)
c
c   Wall #2
c
        call vadd(xms(1,3),ve,        ul,3)
        call vadd(      ul,xu,pms(1,1,2),3)
        call vsub(      ul,xu,pms(1,2,2),3)
        call vsub(xms(1,3),ve,        ul,3)
        call vsub(      ul,xu,pms(1,3,2),3)
        call vadd(      ul,xu,pms(1,4,2),3)
c
c   Wall #3
c
        call vadd(xms(1,4),ve,        ul,3)
        call vadd(      ul,xv,pms(1,1,3),3)
        call vsub(      ul,xv,pms(1,2,3),3)
        call vsub(xms(1,4),ve,        ul,3)
        call vsub(      ul,xv,pms(1,3,3),3)
        call vadd(      ul,xv,pms(1,4,3),3)
c
c   Wall #4
c
        call vadd(xms(1,5),ve,        ul,3)
        call vadd(      ul,xu,pms(1,1,4),3)
        call vsub(      ul,xu,pms(1,2,4),3)
        call vsub(xms(1,5),ve,        ul,3)
        call vsub(      ul,xu,pms(1,3,4),3)
        call vadd(      ul,xu,pms(1,4,4),3)
c
        do k = 1,4
          do j = 1,4
            xcms(j,k,i) = pms(1,j,k)
            ycms(j,k,i) = pms(2,j,k)
            zcms(j,k,i) = pms(3,j,k)
          end do
        end do
c
      end do
c
c =====================================================
c
c ==============   end of definitions   ===============
c
c =====================================================
c
      return
C
    2 continue
c
C  --- ERROR ---
c
      print *, 'Unable to retrieve the DETECTOR information'
      stop
C
      RETURN
      END
