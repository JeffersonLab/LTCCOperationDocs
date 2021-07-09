      SUBROUTINE UGEOM
c     ================
*
************************************************************************
*                                                                      *
*             Routine to define the geometry of the set-up.            *
*                                                                      *
************************************************************************
*
c
      IMPLICIT NONE
      SAVE
#include "gttmed.inc"
#include "gctmed.inc"
#include "gtkine.inc"
#include "gckine.inc"
#include "gtvolu.inc"
#include "gcvolu.inc"
#include "gtcons.inc"
#include "gconst.inc"
#include "ccp.inc"
#include "ccph.inc"
C
      include 'gdinrmc.inc'
      include 'utarget.inc'
C
C=======================================================================
*
      integer mat_cc_gas,mat_sh1,mat_sh2,mat_mirr_cc
      integer med_cc_gas,med_sh1,med_sh2,med_mirr_cc
      integer ivolu,ivol,ICCRMS(36)
c
      REAL*4 PAR(10)
      REAL*4 ZCT(6),ACT(6),WCT(6),ZWA(2),AWA(2),WWA(2)
      REAL*4 ZLG(3),ALG(3),WLG(3),ZSC(2),ASC(2),WSC(2)
      REAL*4 ZGR(3),AGR(3),WGR(3),ZGL(2),AGL(2),WGL(2)
      REAL*4 ZCO2(2),ACO2(2),WCO2(2),DENS_CO2
      REAL*4 ZPOL(4),APOL(4),WPOL(4),DENS_POL
      REAL*4 ZRHC(4),ARHC(4),WRHC(4),DENS_RHC
      real*4 z1,z2,xx,yy,zz, ang0, ang1, ang2,r1,r2
      real*4 rmin1,rmin2,rmax1,rmax2, ztarget, tgl
      real x0ms(3),x1ms(3),x2ms(3),x3ms(3),x4ms(3),x5ms(3)
      real y0ms(3),y1ms(3),y2ms(3),y3ms(3),y4ms(3),y5ms(3)
      real tet1,tet2,tet3,phi1,phi2,phi3
c
      real cc_parm(19)
c
      real wsumGR,PHe3,he3tohe4,DMAXMS
c         
cc      data cc_parm/ -30.,  60.,  1., 5.,
cc     +    232., 276., 335.,
cc     +    339., 189., 335.,
cc     +    380., 130., 320.,
cc     +    420.,  40., 240.,
cc     +    510.,  40.,  45./
c
      data cc_parm/ -30.,  60.,  1., 5.,
     +    422., 276., 335.,
     +    529., 189., 335.,
     +    570., 130., 320.,
     +    610.,  40., 240.,
     +    700.,  40.,  45./
c
      Real*4 Xspace/2000./,Yspace/2000./,Zspace/2200./
      DATA ztarget / 0./   ! target position in cm
      DATA FIELDM /52./
C
C             Ground mixture parameters (Si02 + 5% H2O by weight)
C  Using ground density (RhoGrnd), and water content
C   - calculate WGR array in the program
C
      Real*4 RhoGrnd/1.8/,WinGrnd/0.05/
C                 H         O      Si
      DATA AGR/  1.00794, 16.0,   28.0855 /
      DATA ZGR/  1.0,      8.0,   14.0    /
      DATA WGR/  0.0,      0.0,    0.0    /
C
C             Concrete mixture parameters (from GCALOR writeup)
C                 H         O      Al     Si       Ca      Fe
      DATA ACT/  1.00794, 16.0,   27.0,  28.0855, 40.078, 55.847  /
      DATA ZCT/  1.0,      8.0,   13.0,  14.0,    20.0,   26.0    /
      DATA WCT/  0.004,    0.509,  0.034, 0.345,   0.070,  0.038  /
C
C             Lead glass mixture parameters

      DATA ZLG/  82.00,    14.00,   8.00  /
      DATA ALG/ 207.19,    28.086, 15.999 /
      DATA WLG/   1.0,      2.0,    9.0   /
C
C             Scintilator mixture parameters

      DATA ZSC/  6.00,    1.0   /
      DATA ASC/ 12.01,    1.008 /
      DATA WSC/ 19.00,   21.0   /
C
C             Water mixture parameters

      DATA ZWA/  8.00,    1.0   /
      DATA AWA/ 16.00,    1.008 /
      DATA WWA/  1.00,    2.0   /
C
C             SiO2 glass parameters

      DATA ZGL/ 14.00,    8.0   /
      DATA AGL/ 28.00,   16.0   /
      DATA WGL/  1.00,    2.0   /
c
C             CO2 gas parameters

      DATA ZCO2/ 6.00,    8.0   /
      DATA ACO2/ 12.01,   16.0   /
      DATA WCO2/  1.00,    2.0   /
      DATA DENS_CO2 /0.001977 /
c
c          Polyimide parameters
c
      DATA ZPOL /  6.0,  1.0,  7.0,   8.0 /     
      DATA APOL / 12.01, 1.0, 14.0,  16.0 /     
      DATA WPOL / 22.0, 10.0,  2.0,   5.0 /     
      DATA DENS_POL / 1.42 /     
c
c          ROHACELL parameters
c
      DATA ZRHC /  6.0,  1.0,  7.0,   8.0 /     
      DATA ARHC / 12.01, 1.0, 14.0,  16.0 /     
      DATA WRHC / 9.0, 13.0,  1.0,   2.0 /  ! C9H13O2N   
      DATA DENS_RHC / 0.1 /     
c      DATA DENS_RHC / 0.051 / 
c    
c  stainless steel
c
      integer mxstst
      parameter (mxstst=3)
      real astst(mxstst),zstst(mxstst),wstst(mxstst),dnstst
c                  Fe    Cr    Ni
c                  --    --    --
      data astst / 56. , 52. , 59. /
      data zstst / 26. , 24. , 28. /
      data wstst / 7.4 , 1.8 , 0.8 /
      data dnstst / 7.75 /
C
C--  Define Gas for CC ( C_4F_10 )
C
      INTEGER MXCCGS
      PARAMETER (MXCCGS = 2)
c
      REAL ACCGAS(MXCCGS), ZCCGAS(MXCCGS), WCCGAS(MXCCGS), DCCGAS
      DATA ACCGAS,         ZCCGAS,    WCCGAS,     DCCGAS
     &    /12.01,18.998,   6.0,9.0,  4.0,10.0,   0.009935/
C
C--  Define mirror materials (CH)
C
      REAL ASCIN(2), ZSCIN(2), WSCIN(2), DNSCIN
      DATA   ASCIN/ 12.01, 1.01/
      DATA   ZSCIN/  6.00, 1.00/
      DATA   WSCIN/  1.00, 1.00/
      DATA   DNSCIN/ 1.035/
C
C--  Define "Magnetic shield materials"
C
c       For CC segments 1 - 10 and 13 - 18
c
        REAL  ASHMAT1(2), ZSHMAT1(2), WSHMAT1(2), DSHMAT1
C                    Fe  +  Cu
	DATA ASHMAT1 / 55.85,  63.546 /
	DATA ZSHMAT1 / 26.00,  29.00  /
	DATA WSHMAT1 / 0.7725, 0.2275 /
        DATA DSHMAT1 /  8.79 /
C
C     For CC segments 11 - 12
c
        REAL  ASHMAT2(3), ZSHMAT2(3), WSHMAT2(3), DSHMAT2
c
C                    Fe  +  CH
c
	DATA ASHMAT2 / 55.85, 12.01, 1.01 /
	DATA ZSHMAT2 / 26.00,  6.00, 1.00 /
	DATA WSHMAT2 / 0.8374, 0.1490, 0.0136 /
        DATA DSHMAT2 /  7.49 /
C
C
      data ztarget / 0.0/
      data tgl /5.0/  ! 5.0 cm ( should be no more than 5.0 cm )
C
c
c--------------------------------------------------------------------
C                                                                   !
C  *           ------------------------------ - RMAX    *           !
C   \        /                                \        /|           !
C    \      *                                  *      /PHTRA        !
C     \      \                                /      /              !
C      \      \                              /      /               !
C       \      \                            /PHATRA/                !
C        \      \                          /    | /                 !
C         \*-----*                        *-----*/ - RATRA          !
C          \\                                  //                   !
C           \\                                //                    !
C            \\                        PHATRI//                     !
C             \\                           |//                      !
C              **--------------------------** - RMIN                !
C                                           |                       !
C                                        PHTRI                      !
C                                                                   !
c-------------------------------------------------------------------!
c
C
*     -----------------------------------------------------------------
C
      INTEGER NOFSEC,NOFSP1
      PARAMETER( NOFSEC = 18 )  ! Number of CC segments
      PARAMETER( NOFSP1 = NOFSEC + 1 )
      REAL SR3
      PARAMETER( SR3 = 1.732050807569 ) ! SQRT(3.)
c
      REAL RMIN(NOFSEC), RMAX(NOFSEC), PHTRI(NOFSEC), PHTRA(NOFSEC)
      INTEGER i,j,k,ICCROT(NOFSEC),iccroh(NOFSEC)
      real paredh(5),parede(5),Rh0,xshh,zshh,alenh,phih,
     &          xsh, ysh, zsh, Re0,xshe,zshe,alene,phie
      REAL THETAM,HTHETA,DR2,XSHIFT,ZSHIFT,TFI,TFA,PARCER(NOFSEC)
      real parsh1(5),parsh2(5)
c
      save
C -- Define Ground components
C -- H2O/SiO2 mol.ratio: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c
      WGR(1) = 2.*AGR(1)*WinGrnd*(AGR(3)+2.*AGR(2))/(2.*AGR(1)+AGR(2)) ! H
      WGR(2) = AGR(2)*(2.+WGR(1)/2./AGR(1))                            ! O
      WGR(3) = AGR(3)                                                  ! Si
      wsumGR = WGR(1)+WGR(2)+WGR(3)
      WGR(1) = WGR(1)/wsumGR
      WGR(2) = WGR(2)/wsumGR
      WGR(3) = WGR(3)/wsumGR
C
C*             Define materials
      CALL GSMATE( 1,'HYDROGEN$ ',  1.008, 1., 0.071,865., 0.0,0,0)
      CALL DRINIM( 1)
      CALL GSMATE( 2,'DEUTERIUM$',  2.01,  1., 0.162,757., 0.0,0,0)
      CALL DRINIM( 2)
      CALL GSMATE( 3,'HELIUM$   ',  4.00,  2., 0.125,755., 0.0,0,0)
      CALL DRINIM( 3)
      CALL GSMATE( 4,'LITHIUM$  ',  6.94,  3., 0.534,155., 0.0,0,0)
      CALL GSMATE( 5,'BERYLLIUM$',  9.01,  4., 1.848,35.3, 0.0,0,0)
      CALL DRINIM( 5)
      CALL GSMATE( 6,'CARBON$   ', 12.01,  6., 2.265,18.8, 0.0,0,0)
      CALL DRINIM( 6)
      CALL GSMATE( 9,'ALUMINIUM$', 26.98, 13., 2.7,   8.9, 0.0,0,0)
      CALL DRINIM( 9)
      CALL GSMATE(10,'IRON$     ', 55.85, 26., 7.87,  1.76,0.0,0,0)
      CALL DRINIM(10)
      CALL GSMATE(11,'COPPER$   ', 63.54, 29., 8.96,  1.43,0.0,0,0)
      CALL DRINIM(11)
      CALL GSMATE(12,'TUNGSTEN$ ',183.85, 74.,19.30,  0.35,0.0,0,0)
      CALL DRINIM(12)
      CALL GSMATE(13,'LEAD$     ',207.19, 82.,11.35,  0.56,0.0,0,0)
      CALL DRINIM(13)
      CALL GSMATE(14,'URANIUM$  ',238.03, 92.,18.95,  0.32,0.0,0,0)
      CALL DRINIM(14)
      CALL GSMATE(15,'AIR$      ', 14.610, 7.3, .001205, 30423.0,
     &                                                     0.0,0,0)
      CALL DRINIM(15)

      CALL GSMATE(16,'VACUUM$   ', 0.1E-16,0.1E-16,0.1E-16,0.1E+16,
     &                                                     0.0,0,0)
C###      CALL GSMATE(17,'HELIUM_GAS$ ',4.00,  2., 0.178e-3,529888.,0.,0,0)
      CALL GSMATE(17,'CADMIUM$  ',112.41,48., 8.65,1.0398,18.5,0,0)
      CALL DRINIM(17)
      CALL GSMIXT(18,'WATER$      ',AWA,  ZWA,  1.000,-2,WWA)
      CALL DRINIM(18)
      CALL GSMATE(19,'special ALUM$',26.98,13., 0.135, 178.0, 0.0,0,0)
      CALL DRINIM(19)
      CALL GSMIXT(20,'SCINTILATOR$',ASC,  ZSC,  1.032,-2,WSC)
      CALL DRINIM(20)
      CALL GSMIXT(22,'LEAD GLASS$' ,ALG,  ZLG,  4.7,  -3,WLG)
      CALL DRINIM(22)
c
      call gsmixt(23,' St. steel$',astst, zstst,dnstst,-mxstst,wstst)
      CALL DRINIM(23)
      CALL GSMIXT(24,'CONCRETE$   ',ACT,  ZCT,  2.35,      6,WCT)
      CALL DRINIM(24)
C
C -- He3 gas at 'PHe3' atm
      PHe3 = 10.
      he3tohe4 = PHe3*3.01603/4.0026
      CALL GSMATE(27,'HELIUM3_GAS$',3.00,  2., 0.178e-3*he3tohe4,
     &                                   529888./PHe3,0.,0,0)
      CALL DRINIM(27)
CC
CC -- He3 -- r.l._He3= r.l._He4*A_He3/A_He4*rho_He4/rho_He3
C      rho_He3 = 0.060
C      he3radl = 755.*3.01603/4.0026*0.125/rho_He3
C      CALL GSMATE(25,'HELIUM3$',3.00,  2., rho_He3,
C     &                                   he3radl,0.,0,0)
C      CALL DRINIM(25)
      CALL GSMIXT(30,'GROUND$     ',AGR,  ZGR,  RhoGrnd,   3,WGR)
      CALL DRINIM(30)
      CALL GSMATE(31,'OXYGEN$   ',15.9994,  8., 1.14, 28.7,0.0,0,0)
      CALL DRINIM(31)
      CALL GSMATE(32,'GLASS$', 83.0, 39.0, 2.76,   7.04, 0.0,0,0)
      CALL DRINIM(32)
      CALL GSMIXT(33,'GLASS_T$    ',AGL,  ZGL,  2.770,-2,WGL)
      CALL DRINIM(33)
      CALL GSMIXT(34,'CO2_GAS$    ',ACO2, ZCO2, DENS_CO2,-2,WCO2)
      CALL DRINIM(34)
      CALL GSMIXT(35,'Kapton$    ',APOL, ZPOL, DENS_POL,-4,WPOL)
      CALL DRINIM(35)
C
      CALL GSMATE(36,'Argon gas special$ ', 39.94, 18., .01784, 11133.0,
     &                                                     0.0,0,0)
      CALL DRINIM(36)
      CALL GSMATE(37,'Silicon$ ', 28.0855, 14.0, 2.33, 9.36,
     &                                                     0.0,0,0)
      CALL DRINIM(37)
c
C -- ROHACELL -- all parameters approximate
c
      CALL GSMIXT(38,'Rohacell$',ARHC, ZRHC, DENS_RHC,-4,WRHC)
      CALL DRINIM(38)
c
      CALL GSMATE(39,'Molibden$ ', 95.94, 42.0, 10.22, 0.96,
     &                                                     0.0,0,0)
      CALL DRINIM(39)
c
      CALL GSMATE(40,'Sn ( Tin)$ ', 118.71, 50.0, 7.31, 1.21,
     &                                                     0.0,0,0)
      CALL DRINIM(40)
C
C -- Boron 10 -- all parameters approximate
c
      CALL GSMATE(43,'BORON$', 10., 5., 2., 24.6, 0.0,0,0)
      CALL DRINIM(43)
c
      mat_cc_gas = 44
      CALL GSMIXT(mat_cc_gas,'Gas for CC$',
     & ACCGAS,ZCCGAS,DCCGAS,-MXCCGS,WCCGAS)
      CALL DRINIM(mat_cc_gas)
c
c --- Mirror material ( close to CH ) ---
c
      mat_mirr_cc = 45
      CALL GSMIXT(mat_mirr_cc,'CC mirror material $',
     &     ASCIN, ZSCIN,DNSCIN,-2, WSCIN)
      CALL DRINIM(mat_mirr_cc)
c
C --- Magnetic shield material (1) ---
c
      mat_sh1 = 46
      CALL GSMIXT(mat_sh1,'Magnet. Shield_1$',
     & ASHMAT1,ZSHMAT1,DSHMAT1,2,WSHMAT1)
      CALL DRINIM(mat_sh1)
C
C --- Magnetic shield material (2) ---
c
      mat_sh2 = 47
      CALL GSMIXT(mat_sh2,'Magnet. Shield_2$',
     & ASHMAT2,ZSHMAT2,DSHMAT2,3,WSHMAT2)
      CALL DRINIM(mat_sh2)

c
c ==============================================================
c
*             Defines tracking media parameters.
      DMAXMS =  0.50
      DEEMAX =  0.10
      EPSIL  =  0.0001
      STMIN  =  0.002
*
      CALL GSTMED( 1, 'HYDROGEN$     ', 1, 0, 2, 52., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED( 2, 'DEUTERIUM$    ', 2, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED( 3, 'HELIUM$       ', 3, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED( 4, 'LITHIUM$      ', 4, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED( 5, 'BERYLLIUM$    ', 5, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED( 6, 'CARBON$       ', 6, 0, 2, 52., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED( 9, 'ALUMINIUM$    ', 9,
     * 0, 2, 52., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(10, 'IRON$         ',10, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(11, 'COPPER$       ',11, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(12, 'TUNGSTEN$     ',12,
     * 0, 2, 52., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(13, 'LEAD$         ',13,
     * 0, 2, 52., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(14, 'URANIUM$      ',14, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(15, 'AIR$          ',15,
     * 0, 2, 52., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(16, 'VACUUM$       ',16,
     * 0, 2, 52., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(17, 'CADMIUM$      ',17, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(18, 'WATER$        ',18, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(19, 'special ALUMINIUM$    ', 19,
     * 0, 2, 52., 0., DMAXMS, DEEMAX, EPSIL, STMIN, 0, 0)

      CALL GSTMED(20, 'SCINTILLATOR$ ',20,
     * 1, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)

      CALL GSTMED(22, 'LEAD GLASS$   ',22, 1, 2, 50.0, 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      call gstmed(23,'St. Steel$',23,
     * 0, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(24, 'CONCRETE$     ',24, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
C      CALL GSTMED(25, 'HELIUM3$      ',25, 0, 0, 0., 0., DMAXMS,
C     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(27, 'HELIUM3-GAS$  ',27, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(30, 'GROUND$       ',30, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(31, 'OXYGEN$       ',31, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(32, 'GLASS$        ',32, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(33, 'GLASS_T$      ',33, 0, 0, 0., 0., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
c
c   CO2 definition - sensitive volume
c
      CALL GSTMED(34, 'CO2_GAS$      ',34,
     * 1, 2, 52., 3., DMAXMS,
     * DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(35, 'KAPTON$      ',35,
     * 0, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(36, ' Argon gas$  ',36,
     * 0, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(37, ' Silicon $  ',37,
     * 0, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(38, ' ROHACELL $  ',38,
     * 0, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(39, ' Molibden $  ',39,
     * 0, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
      CALL GSTMED(40, ' Sn (Tin) $  ',40,
     * 0, 2, 52., 3., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
c----------------------------------------------
      CALL GSTMED(43, 'BORON$    ',43, 0, 2, 52., 3., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
C
C -- Magnetic field
      CALL GSTMED(42, 'VACUUMB$      ',16, 0, 2, 50., 5., DMAXMS,
     *            DEEMAX, EPSIL, STMIN, 0, 0)
c
c --- Gas for CC ---
c
      med_cc_gas = 44
      IFIELD = 1
      CALL GSTMED(med_cc_gas,'GAS for CC$', mat_cc_gas,
     &     1,IFIELD,FIELDM,
     1     1., DMAXMS,DEEMAX, EPSIL, STMIN, 0, 0)
c ---
      med_mirr_cc = 45
      CALL GSTMED( med_mirr_cc,'CC mirror $', mat_mirr_cc,
     &     0, IFIELD,FIELDM,
     1     1.,    0.5,    1.,   0.01,    0.2, 0, 0)
c ---
      med_sh1 = 46
      CALL GSTMED(med_sh1,'Magnetic shield - 1$', mat_sh1,
     &     0,IFIELD,FIELDM,
     1     2.,    0.5,    1.,   0.01,    0.1, 0, 0)
c ---
      med_sh2 = 47
      CALL GSTMED(med_sh2,'Magnetic shield - 2$', mat_sh2,
     &     0,IFIELD,FIELDM,
     1     2.,    0.5,    1.,   0.01,    0.1, 0, 0)
c ---
C
C -- Birk's parameters for scintillator medium
      CALL GSTPAR(20,'BIRK1',1.)
      CALL GSTPAR(20,'BIRK2',0.013)
      CALL GSTPAR(20,'BIRK3',9.6E-6)
C
C --          Define geometry of the setup
C             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C -- Rotation matrices
      CALL GSROTM(1,90.,  0.,180.,0.,90., 90.)
      CALL GSROTM(2,90.,  0.,  0.,0.,90.,270.)
      CALL GSROTM(3,90., 90.,  0.,0.,90.,  0.)
      CALL GSROTM(4,90.,270.,  0.,0.,90.,180.)
      CALL GSROTM(5,90., 90., 90.,180.,0.,0.)
C
C -- Define the overall volume CLAS (air)
      PAR(1) = Xspace/2.
      PAR(2) = Yspace/2.
      PAR(3) = Zspace/2.
      CALL GSVOLU('CLAS', 'BOX ', 15, PAR, 3, IVOLU)
c ---> special      CALL GSVOLU('CLAS', 'BOX ', 16, PAR, 3, IVOLU)
c
c     some plate in imitating central silicon detector
c
      par(1) = 1.2
      par(2) = 25.0
      par(3) = 0.5*0.233  ! about 2.33 mm of Si
      xx = 0.
      yy = 0.
      zz = 12. + ztarget
c
c      CALL GSVOLU('CDET', 'TUBE', 37, PAR, 3, IVOLU)    ! Si
c      call GSPOS('CDET',1,'CLAS', xx, yy, zz, 0,'ONLY')
c
c ----   LH2 target ----- lenght tgl
c
c      include '/home/alex/.MCWORKS/JLAB/alex/ltcc/src/LH2_12.f'
c
C
C --  *** NEW SCATTERING CHAMBER ***
c
c      include '/home/alex/.MCWORKS/JLAB/alex/ltcc/src/scatt_chamb_a.f'

C
c
      ang0 = 25.0*degrad
      ang1 =  5.0*degrad
      ang2 = 45.0*degrad
      r1   = 400.
      r2   = 480.
c      z1   = r1*cos(ang0) + r1*sin(ang0)*sin(ang0)
      z1   = r1*cos(ang0)
      z2   = r2*cos(ang0)
      par(1) = 0.5*(z2 - z1)
      par(2) = z1*sin(ang1)
      par(3) = z1*sin(ang2)
      par(4) = z2*sin(ang1)
      par(5) = z2*sin(ang2)
      xx = 0.
      yy = 0.
      zz = 0.5 * ( z1 + z2 )
c
c      call gsvolu('CVOL','CONE',15,par,5,ivolu)
c      call gspos ('CVOL',1,'CLAS',xx,yy,zz,0,'MANY')
c
c
c  Define mother volume for Cerenkov (PGON)
c
      call gsvolu('CSEC','PGON',15,cc_parm,19,ivolu)
      call gspos ('CSEC',1,'CLAS',0.,0.,0.,0,'MANY')	
C
C====================================================================
c
      CALL CCDEIN
c     ===========
c
C====================================================================
c
C
C ==========================================================
C
C     Put CC-counts ('CCxx') into mother volume 'CC  '
C
      CALL GSVOLU ('CCVO', 'TRD2' , med_cc_gas , PARCER , 0 , IVOL )
      CALL GSVOLU ('CCHP', 'TUBS' , med_mirr_cc, PAREDH , 0 , IVOL )
      CALL GSVOLU ('CCEP', 'TUBS' , med_mirr_cc, PAREDE , 0 , IVOL )
      CALL GSVOLU ('CCHL', 'TUBS' , med_mirr_cc, PAREDH , 0 , IVOL )
      CALL GSVOLU ('CCEL', 'TUBS' , med_mirr_cc, PAREDE , 0 , IVOL )
c
c      call GSVOLU ('CCS1', 'BOX ' , med_sh1 , PARSH1 , 0 , IVOL )
c      call GSVOLU ('CCS2', 'BOX ' , med_sh2 , PARSH2 , 0 , IVOL )
c   *****  For direct electrons only   ******
      call GSVOLU ('CCS1', 'BOX ' , 15 , PARSH1 , 0 , IVOL )
      call GSVOLU ('CCS2', 'BOX ' , 15 , PARSH2 , 0 , IVOL )
c   *****  For direct electrons only   ******
c
c
c
      DO i=1,NOFSEC
C     =============
C
        RMIN(i) = Y21(i) - 0.5
        RMAX(i) = Y12(i) + 0.5
C
        THETAM  = th(i)
        TFI=RMIN(i)*SIN(THETAM*DEGRAD)/SR3
        TFA=RMAX(i)*SIN(THETAM*DEGRAD)/SR3
        HTHETA=TAN((THETAM-THMIN(i))*DEGRAD)
        DR2 = (RMAX(i)-RMIN(i)) / 2.
        XSHIFT=(DR2 + RMIN(i))*SIN(THETAM*DEGRAD)
        ZSHIFT=(DR2 + RMIN(i))*COS(THETAM*DEGRAD)
C
        PHTRI(i) = TFI - 4.5
        PHTRA(i) = TFA - 4.5
C
C --- Definition of cherenkov counter
C
        PARCER (1)=PHTRI(i)
        PARCER (2)=PHTRA(i)
        PARCER (3)=RMIN(i)*HTHETA
        PARCER (4)=RMAX(i)*HTHETA
        PARCER (5)=DR2
 
        Re0 = 280.
        alene = PHTRA(i)/2.
        phie  = raddeg*atan(alene/Re0)
        parede(1) = Re0 - 0.35
        parede(2) = Re0
        parede(3) = RMAX(i)*HTHETA
        parede(4) = -phie
        parede(5) =  phie
        zshe = DR2 - Re0
        xshe = alene
 
        Rh0 = 210.
        alenh = PHTRI(i)/2.
        phih  = raddeg*asin(alenh/Rh0)
        paredh(1) = Rh0
        paredh(2) = Rh0 + 0.3
        paredh(3) = RMIN(i)*HTHETA
        paredh(4) = -phih
        paredh(5) =  phih
        zshh      = -DR2 - sqrt(Rh0*Rh0 - alenh*alenh)
        xshh = alenh
C
C --- Define rotation to position CC section
c=====================================================
c
        CALL SGNROT(ICCROT(i))
        CALL GSROTM(ICCROT(i),90.,90.,THETAM-90.,0.,THETAM,0.)
        CALL GSPOSP( 'CCVO', i,'CSEC',XSHIFT,0.,ZSHIFT,ICCROT(i),
     &    'MANY', PARCER, 5)
c=====================================================
        XSH=(RMAX(i)-Re0)*SIN(THETAM*DEGRAD)
        YSH=alene
        ZSH=(RMAX(i)-Re0)*COS(THETAM*DEGRAD)
C
        CALL SGNROT(ICCROH(i))
C       ======================
        CALL GSROTM(ICCROH(i),THETAM,0.,90.,90.,THETAM-90.,0.)
        CALL GSPOSP( 'CCEP', i,'CSEC',xsh,ysh,zsh,
     &    ICCROH(i),'ONLY', parede, 5)
C
        YSH=-alene
        CALL GSPOSP( 'CCEL', i,'CSEC',xsh,ysh,zsh,
     &    ICCROH(i),'ONLY', parede, 5)
c=====================================================
        XSH=(DR2 + RMIN(i)+zshh)*SIN(THETAM*DEGRAD)
        YSH=alenh
        ZSH=(DR2 + RMIN(i)+zshh)*COS(THETAM*DEGRAD)
C
        CALL GSPOSP( 'CCHP', i,'CSEC',xsh,ysh,zsh,
     &    ICCROH(i),'ONLY', paredh, 5)
C
        YSH=-alenh
        CALL GSPOSP( 'CCHL', i,'CSEC',xsh,ysh,zsh,
     &    ICCROH(i),'ONLY', paredh, 5)
c
c=====================================================
c  Define the Magnetic shield walls
c
c  Center of MS box
c
        call ucopy(pcms0(1,1,i),x0ms,3)
c
c  Centers of walls
c
        call ucopy(pcms0(1,2,i),x1ms,3)
        call ucopy(pcms0(1,3,i),x2ms,3)
        call ucopy(pcms0(1,4,i),x3ms,3)
        call ucopy(pcms0(1,5,i),x4ms,3)
        call ucopy(pcms0(1,6,i),x5ms,3)
c
c  Remember now, that there are two sides of CC segment !
c
        call ucopy(x0ms,y0ms,3)
        call ucopy(x1ms,y1ms,3)
        call ucopy(x2ms,y2ms,3)
        call ucopy(x3ms,y3ms,3)
        call ucopy(x4ms,y4ms,3)
        call ucopy(x5ms,y5ms,3)
c
        y0ms(1) = -x0ms(1)
        y1ms(1) = -x1ms(1)
        y2ms(1) = -x2ms(1)
        y3ms(1) = -x3ms(1)
        y4ms(1) = -x4ms(1)
        y5ms(1) = -x5ms(1)
c
c  Now transfer that coordinates to the SG R.S.
c
c   From CC to SG Reference System :
c
        call CCTOSG(i,x0ms,x0ms)
        call CCTOSG(i,x1ms,x1ms)
        call CCTOSG(i,x2ms,x2ms)
        call CCTOSG(i,x3ms,x3ms)
        call CCTOSG(i,x4ms,x4ms)
        call CCTOSG(i,x5ms,x5ms)
        call CCTOSG(i,y0ms,y0ms)
        call CCTOSG(i,y1ms,y1ms)
        call CCTOSG(i,y2ms,y2ms)
        call CCTOSG(i,y3ms,y3ms)
        call CCTOSG(i,y4ms,y4ms)
        call CCTOSG(i,y5ms,y5ms)
c
c   Now we are to define magnetic shield walls in GEANT:
c
        call CCROTP(x0ms,x1ms,x2ms,tet1,tet2,tet3,phi1,phi2,phi3)
c
        k = 2*i - 1
        CALL SGNROT(ICCRMS(k))
        CALL GSROTM(ICCRMS(k),tet1,phi1,tet2,phi2,tet3,phi3)
c
c    The other side of CC segment:
c
        call CCROTP(y0ms,y1ms,y2ms,tet1,tet2,tet3,phi1,phi2,phi3)
c
        k = 2*i
*        ICCRMS(k) = get_rotm() + 20
        CALL SGNROT(ICCRMS(k))
        CALL GSROTM(ICCRMS(k),tet1,phi1,tet2,phi2,tet3,phi3)
c
C ================================================================
c
        if(i.ge.1.and.i.le.10) then
c
c ------   Variant 1 ------ * 9.93 * 7.76 * 21.08 *** l eff = 0.455
c
          parsh1(3) = 21.08
          k = 2*i - 1
c
          j = 8*(i-1) + 1
          parsh1(1) =  0.455
          parsh1(2) =  7.76
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x1ms(1),x1ms(2),x1ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-1) + 2
          parsh1(1) =  9.93
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x2ms(1),x2ms(2),x2ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-1) + 3
          parsh1(1) =  0.455
          parsh1(2) =  7.76
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x3ms(1),x3ms(2),x3ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-1) + 4
          parsh1(1) =  9.93
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x4ms(1),x4ms(2),x4ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          k = 2*i
c
          j = 8*(i-1) + 5
          parsh1(1) =  0.455
          parsh1(2) =  7.76
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y1ms(1),y1ms(2),y1ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-1) + 6
          parsh1(1) =  9.93
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y2ms(1),y2ms(2),y2ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-1) + 7
          parsh1(1) =  0.455
          parsh1(2) =  7.76
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y3ms(1),y3ms(2),y3ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-1) + 8
          parsh1(1) =  9.93
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y4ms(1),y4ms(2),y4ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
C ================================================================
c
        else if(i.ge.11.and.i.le.12) then
c
c ------   Variant 2 ------ * 11.46 * 8.46 * 24.13 *** l eff = 0.493
c
c
          parsh2(3) = 24.13
          k = 2*i - 1
c
          j = 8*(i-11) + 1
          parsh2(1) =  0.493
          parsh2(2) =  8.46
c
          CALL GSPOSP( 'CCS2', j,'CSEC',x1ms(1),x1ms(2),x1ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
          j = 8*(i-11) + 2
          parsh2(1) = 11.46
          parsh2(2) =  0.493
c
          CALL GSPOSP( 'CCS2', j,'CSEC',x2ms(1),x2ms(2),x2ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
          j = 8*(i-11) + 3
          parsh2(1) =  0.493
          parsh2(2) =  8.46
c
          CALL GSPOSP( 'CCS2', j,'CSEC',x3ms(1),x3ms(2),x3ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
          j = 8*(i-11) + 4
          parsh2(1) = 11.46
          parsh2(2) =  0.493
c
          CALL GSPOSP( 'CCS2', j,'CSEC',x4ms(1),x4ms(2),x4ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
          k = 2*i
c
          j = 8*(i-11) + 5
          parsh2(1) =  0.493
          parsh2(2) =  8.46
c
          CALL GSPOSP( 'CCS2', j,'CSEC',y1ms(1),y1ms(2),y1ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
          j = 8*(i-11) + 6
          parsh2(1) = 11.46
          parsh2(2) =  0.493
c
          CALL GSPOSP( 'CCS2', j,'CSEC',y2ms(1),y2ms(2),y2ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
          j = 8*(i-11) + 7
          parsh2(1) =  0.493
          parsh2(2) =  8.46
c
          CALL GSPOSP( 'CCS2', j,'CSEC',y3ms(1),y3ms(2),y3ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
          j = 8*(i-11) + 8
          parsh2(1) = 11.46
          parsh2(2) =  0.493
c
          CALL GSPOSP( 'CCS2', j,'CSEC',y4ms(1),y4ms(2),y4ms(3),
     &    ICCRMS(k),'ONLY', parsh2, 3)
c
C ================================================================
c
        else if(i.ge.13.and.i.le.18) then
c
c ------   Variant 3 ------ * 12.65 * 9.85 * 27.94 *** l eff = 0.455
c
          parsh1(3) = 27.94
          k = 2*i - 1
c
          j = 8*(i-13+10) + 1
          parsh1(1) =  0.455
          parsh1(2) =  9.85
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x1ms(1),x1ms(2),x1ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-13+10) + 2
          parsh1(1) = 12.65
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x2ms(1),x2ms(2),x2ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-13+10) + 3
          parsh1(1) =  0.455
          parsh1(2) =  9.85
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x3ms(1),x3ms(2),x3ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-13+10) + 4
          parsh1(1) = 12.65
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',x4ms(1),x4ms(2),x4ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          k = 2*i
c
          j = 8*(i-13+10) + 5
          parsh1(1) =  0.455
          parsh1(2) =  9.85
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y1ms(1),y1ms(2),y1ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-13+10) + 6
          parsh1(1) = 12.65
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y2ms(1),y2ms(2),y2ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-13+10) + 7
          parsh1(1) =  0.455
          parsh1(2) =  9.85
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y3ms(1),y3ms(2),y3ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
          j = 8*(i-13+10) + 8
          parsh1(1) = 12.65
          parsh1(2) =  0.455
c
          CALL GSPOSP( 'CCS1', j,'CSEC',y4ms(1),y4ms(2),y4ms(3),
     &    ICCRMS(k),'ONLY', parsh1, 3)
c
C ================================================================
c
        end if
c
c==========================================================
      ENDDO
c
c
c
c
      z1 = 701.0
      z2 = 703.0
      rmin1 = 0.
      rmax1 = 708.0
      PAR(1) = rmin1
      PAR(2) = rmax1
      PAR(3) = 0.5*(z2 - z1)
      xx = 0.
      yy = 0.
      zz = 0.5 * ( z1 + z2 )
      CALL GSVOLU('PSTP', 'TUBE', 15, PAR, 3, IVOLU)   ! Air
      call GSPOS('PSTP',1,'CLAS', xx, yy, zz, 0,'ONLY')

c
C -- Close geometry banks. Mandatory system routine.
C
      CALL GGCLOS
C
      END

