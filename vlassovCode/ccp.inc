c  ccp.inc

c  $Id: ccp.inc,v 1.2 1997/04/24 15:38:10 holtrop Exp $

      character*132 ciccp
      data ciccp/
     1'$Id: ccp.inc,v 1.2 1997/04/24 15:38:10 holtrop Exp $'
     1/


c---  Alexander Vlassov  --------
C
C --- Common /CCP/ - geometry --------------------------------------------
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
c-            centers of MS walls (4) and center of MS enterance window (1)
c-
c-   plmsw  - Magnetic shield window plane parameters
c-
c-   scrnp1,scrnp2 - two points on the sidewall: first - center
c-                   of the screen, second - center of the hole
c-
      REAL pe0,ph0,pw0,p00,p10,p20,p30,p40,p50,p60,p70,p80,p90,pa0
      REAL x11,y11,x12,y12,x21,y21,x22,y22
      REAL csa0,wc0,sw0,th,thmin,xp1,yp1,xp2,yp2
      REAL x_pl,elt_pl,pmt0,wcr,wcer,pmtr,dplwc,wcang,hwcpmt
      REAL pwc10,pwc20,hwde,hwdh,cmsx1,cmsx2,hmssz,plms0,pcms0
      REAL plmsw,scrnp1,scrnp2
c
      COMMON/CCP/pe0(5,18),ph0(5,18),pw0(3,18),wc0(3,18),
     *  pmt0(3,18),wcr(18),wcer(18),pmtr(18),dplwc(18),wcang(18),
     *  p00(3),p10(3,18),p20(3,18),p30(3,18),p40(3,18),
     *  p50(3,18),p60(3,18),p70(3,18),p80(3,18),p90(3,18),pa0(3,18),
     *  x11(18),y11(18),x12(18),y12(18),x21(18),y21(18),x22(18),
     *  y22(18),xp1(18),yp1(18),xp2(18),yp2(18),csa0(18),sw0(3,18),
     *  th(18),thmin(19),
     *  x_pl(18),elt_pl(3,18),hwcpmt(18),pwc10(3,18),pwc20(3,18),
     *  hwde(18),hwdh(18),cmsx1(3,18),cmsx2(3,18),hmssz(3,18),
     *  plms0(3,4,18),pcms0(3,6,18),plmsw(3,18),scrnp1(3,18),
     *  scrnp2(3,18)
C
***************************************************8
C
C     Extracted parameters for ONE cc segment.
c-    Optionally filled by CCDEFPAR subroutine.
c-
c-    PE = pe0(i) ; PH = ph0(i) ; PW = pw0(i) ;
c-    P1,...P9,PA = p10(i),...p90(i),PA(i) ;
c-    PE_PL = elt_pl(i) ;
c-    WCC = wc0(i) ; SW = sw0(i) ;
c-    CSA = csa0(i) ;
c-    XE_PL = x_pl(i) ;
c-    pwc1,pwc2 = pwc10(i),pwc20(i)
c-    plms = plms0( , ,i)
c-    pcms = Magnetic Shield Point Coordinates - 4 points of 4 walls
c-    scrp1,scrp2 = scrnp1(i), scrnp2(i)
c-
      REAL PE,PH,PW,P1,P2,P3,P4,P5,P6,P7,P8,P9,PA,
     *  PE_PL,WCC,SW,CSA,XE_PL,pwc1,pwc2,plms,pcms,scrp1,scrp2
C-
      COMMON /CCPW/ PE(5),PH(5),PW(3),P1(3),P2(3),P3(3),P4(3),
     *  P5(3),P6(3),P7(3),P8(3),P9(3),PA(3),
     *  PE_PL(3),WCC(3),SW(3),CSA,XE_PL,
     *  pwc1(3),pwc2(3),plms(3,4),pcms(3,4,4),scrp1(3),scrp2(3)
C----
C
C *********************************************************************
C
