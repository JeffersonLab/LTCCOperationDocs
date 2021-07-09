      PROGRAM AV
c     ==========
c
C    Purpose and methods : calculates new parameters for clas12 geometry
C
C    Created:   29-DEC-2009   Alexander V Vlassov
C
c=============================================================
c
      implicit none
      include 'ccp.inc'

c
      integer ib
      common /PAWC/ ib(100000)
c
************   +SEQ, CCP.    ********************************************
C
C --- Common /CCP/ - geometry ------------------------------------------
C
C --- Parameters of mirror definition p1*x2+p2*y2+p3*xy+p4*x+p5*y+1=0
C-   pe0 - elliptic mirror parameters
C-   ph0 - hiperbolic mirror parameters
C-   pwo - WC ellips parameters in its coordinate system:
C-     pw0(1)*x**2 + pw0(2)*y**2 + pw0(3)*z**2 + 1 = 0
C-     where z is along the biggest axis of the ellips.
C-   wc0 - coordinates of WC window center.
C-   pmt0 - coordinates of PMT center.
C-     PMT assume to be flat
C-   wcr  - WC max radius
C-   wcer - WC window radius
C-   pmtr - radius of PMT
C-   dplwc - distance between two planes in WC.
C-   wcang - angle between this planes and Segment median plane.
C-     PLANES ARE DEFINED as p(1)*x + p(2)*y + p(3)*z + 1.0 = 0
C-
C-   p00 - Sector middle plane
C-   p10 - PMT black wall
C-   p20 - Coil plane
C-   p30 - WC window plane
C-   p40 - PMT plane
C-   p50 - Bottom black wall
C-   p60 - edge plane between this CC segment and smaller number one
C-   p70 - edge plane between this CC segment and bigger number one
C-   p80 - bounding WC plane between this WC segment and smaller WC.
C-   p90 - bounding WC plane between this WC segment and bigger WC.
C-   pa0 - Plane, parallel to Sector Dividing plane,at dist. 1in.
C-         ( surface of the sidewall)
C-   pb0 - Additional mirror surface.
C-
c-   x11,y11 - x12,y12 : PMT side wall plane
c-   x11,y11 (in the WC window plane) is the point in OUTER
c-   WC surface (shield?), nearest to the middle plane
c-   x12,y12 is the distant edge of elliptic mirror
c-
c-   x21,y21 - x22,y22 : COIL plane
c-   x22,y22 is the distant point of inner surface of WC window.
c-   x21,y21 defines an edge point for all CC section when discribed
c-   as a trapezoid. If photon crosses this plane in positive
c-   direction, it suppose to be absorbed.
c-
c-   csa0 - angle between Y and WC axis (in radians).
c-   sw0  - WC elliptical surface center coordinates in SEGment
c-   Reference System. They are used for Reference transformation
c-   between CC SEGment R.S. and WC R.S.
c-
c-   th - polar angles (in degree) of the centers of cc segment
c-   thmin - polar angles of planes between CC segments.
c-   ****  should be redefined ! ****
c-
c-   xp1,yp1 - xp2,yp2  : PMT coordinates are calculated.
c-
c-   x_pl - X coordinate of the beginning of flat part of
c-   elliptic mirror: 7.0 or 10.0 cm
c-   elt_pl - this plane parameters 
c-
c-   hwcpmt - distance between WC window and PMT planes
c-
c-   pwc10  - parameters of WC bounding plane (near the smaller CC #)
c-   pwc20  - parameters of WC bounding plane (near the bigger  CC #)
c-
c-   hwde   - half of the elliptical width
c-   hwdh   - half of the hiperbolic width
c-
c-   cmsx1,cmsx2   - coordinates of the Magnetic Shield axis points
c- 
c-   hmssz  - half of the Magnetic Shield walls sizes
c-   
c-   plms0  - Magnetic shield planes parameters
c-      
c-   pcms0  - Magnetic shield points coordinates - center of MS (1),
c-          centers of MS walls (4) and center of MS enterance window (1)
c-   
c-   plmsw  - Magnetic shield window plane parameters
c- 
c-   scrnp1,scrnp2 - two points on the sidewall: first - center
c-                   of the screen, second - center of the hole
c- 
*******************   +SEQ, GCONSP.    ***************************** 
      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      DOUBLE PRECISION EMMU,PMASS,AVO
*
      PARAMETER (PI=3.14159265358979324)
      PARAMETER (TWOPI=6.28318530717958648)
      PARAMETER (PIBY2=1.57079632679489662)
      PARAMETER (DEGRAD=0.0174532925199432958)
      PARAMETER (RADDEG=57.2957795130823209)
      PARAMETER (CLIGHT=29979245800.)
      PARAMETER (BIG=10000000000.)
      PARAMETER (EMASS=0.0005109990615)
      PARAMETER (EMMU=0.105658387)
      PARAMETER (PMASS=0.9382723128)
      PARAMETER (AVO=0.60221367)
*******************   +SEQ, GCONSP.    ***************************** 
C
      integer i,j,k,ncc
      integer lunr,lunw
c
      real t(4),s,thmin12(19),th12(18),sz,zdiff,al
      real wc012(3,18),pmt012(3,18),Y12_11(18),y12_12(18)
      real y12_21(18),y12_22(18)
c
      double precision r(5)
c
      data zdiff /199.6/
c
      save
c
 900  format(5e25.15)
c
      call hlimit(99000)
c
      LUNR = 33
      OPEN( UNIT=LUNR, FILE='ccgeom.dat',
     &  STATUS='OLD',FORM='FORMATTED',ACCESS='SEQUENTIAL')
      write(*,*) ' File for CLAS geometry  is opened - Lun:',LUNR
c
      LUNW = 34
      OPEN( UNIT=LUNW, FILE='ccgeom_12.dat',
     &  STATUS='NEW',FORM='FORMATTED',ACCESS='SEQUENTIAL')
      write(*,*) ' File for CLAS12 geometry  is opened - Lun:',LUNW
C
      DO i=1,36
        READ (LUNR,*,end=2) r
        write(lunw,900,err=2) r
      END DO
c
c    WC and PMT parameters
c
      do i = 1,18
        read(LUNR,*,end=2) j,wc0(1,i),wc0(2,i),pmt0(1,i),pmt0(2,i),
     &  wcr(i),pmtr(i),dplwc(i),wcang(i)
        wc0 (3,i) = 0.
        pmt0(3,i) = 0.
      end do
c
C --- Input parameters for PMT wall
c
      DO i=1,18
        READ (LUNR,*) t
        x11(i)= t(1)
        y11(i)= t(2)
        x12(i)= t(3)
        y12(i)= t(4)
      END DO
c
C --- Input parameters for coil plane
c
      DO i=1,18
        READ (LUNR,*) t
        x21(i)= t(1)
        y21(i)= t(2)
        x22(i)= t(3)
        y22(i)= t(4)
      END DO
c
C --- Mirror parameters
c
      DO i=1,18
        READ (LUNR,*) hwde(i), hwdh(i), wcz(i), wcxy(i),
     &  cmsx1(2,i), cmsx1(3,i), cmsx2(2,i), cmsx2(3,i),
     &  hmssz(1,i),  hmssz(2,i),  hmssz(3,i)
        cmsx1(1,i) = 0.  
        cmsx2(1,i) = 0.  
      END DO
c
C --- Segment parameters
c
      DO i=1,18
        READ (LUNR,*) j, thmin(i), th(i), scrnp1(1,i),
     &  scrnp1(2,i), scrnp2(1,i), scrnp2(2,i)
        scrnp1(3,i) = 0.  
        scrnp2(3,i) = 0.  
      END DO
c
      read (lunr,*) j,thmin(19)
c
C---- end of reading -----
c
      CLOSE (LUNR)
c
      DO i=1,18
        s = pmt0(2,i)*cos(degrad*thmin(i))
        s = s + zdiff
        al = sqrt(s*s + (pmt0(2,i)*sin(degrad*thmin(i)))**2)
        thmin12(i) = raddeg*acos(s/al)
c        print *,i,thmin12(i)
        s = pmt0(2,i)*cos(degrad*th(i))
        s = s + zdiff
        al = sqrt(s*s + (pmt0(2,i)*sin(degrad*th(i)))**2)
        th12(i) = raddeg*acos(s/al)
c        print *,i,th12(i)
      END DO
c
      i = 19
      s = pmt0(2,i-1)*cos(degrad*thmin(i))
      s = s + zdiff
      al = sqrt(s*s + (pmt0(2,i-1)*sin(degrad*thmin(i)))**2)
      thmin12(i) = raddeg*acos(s/al)
c      print *,i,thmin12(i)
c
      do i = 1,18
        wc012(1,i)  =  wc0(1,i)
        wc012(2,i)  =  wc0(2,i) + zdiff*cos(degrad*th12(i))
        wc012(3,i)  =  wc0(3,i)
        pmt012(1,i) = pmt0(1,i)
        pmt012(2,i) = pmt0(2,i) + zdiff*cos(degrad*th12(i))
        pmt012(3,i) = pmt0(3,i)
      end do
c
      do i = 1,18
        y12_11(i) = y11(i) + zdiff*cos(degrad*th12(i))
        y12_12(i) = y12(i) + zdiff*cos(degrad*th12(i))
        y12_21(i) = y21(i) + zdiff*cos(degrad*th12(i))
        y12_22(i) = y22(i) + zdiff*cos(degrad*th12(i))
        scrnp1(2,i) = scrnp1(2,i) + zdiff*cos(degrad*th12(i))
        scrnp2(2,i) = scrnp2(2,i) + zdiff*cos(degrad*th12(i))
      end do
c
c    WC and PMT parameters
c
      do i = 1,18
        write(LUNW,901,err=2) i,wc012(1,i),wc012(2,i),pmt012(1,i),
     &  pmt012(2,i),wcr(i),pmtr(i),dplwc(i),wcang(i)
      end do
c
 901  format(i3,8f9.3)
c
      DO i=1,18
        write (LUNW,902) X11(i),Y12_11(i),X12(i),Y12_12(i)
      END DO
 902  format(4f12.4)
      DO i=1,18
        write (LUNW,902) X21(i),Y12_21(i),X22(i),Y12_22(i)
      END DO
c
      DO i=1,18
        write (LUNW,903) hwde(i), hwdh(i), wcz(i), wcxy(i),
     &  cmsx1(2,i), cmsx1(3,i), cmsx2(2,i), cmsx2(3,i),
     &  hmssz(1,i),  hmssz(2,i),  hmssz(3,i)
      END DO
 903  format(2f7.2,2f8.3,7f8.2)
c
      DO i=1,18
        write (LUNW,904) i, thmin12(i), th12(i), scrnp1(1,i),
     &  scrnp1(2,i), scrnp2(1,i), scrnp2(2,i)
      END DO
 904  format(i6,7f9.3)
c
      i = 19
      write (lunw,904) i,thmin12(i)
c
      close(lunw)
      stop
c
 2    continue
      print *,' Something is WRONG !!! '
      stop
c
      end
