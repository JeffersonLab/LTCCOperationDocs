      PROGRAM TEST
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : TEST OF SOMETHING 
C-
C-
C-   Created  14-MAR-1994   Alexander V. Vlassov
C-
C----------------------------------------------------------------------
C
      double precision x1(2),x2(2),R,angled
C
      double precision p_in(5), p_fn(5), xshift(2), p_ff(5)
      double precision unit,d2,df2,ds 
      data unit /1.0d+000/
      real x_1(2,18),x_2(2,18),R_C(18),angle
C
C     Ellips first focus
CCC      data x_1 / 0.172,    9.926,
      data x_1 / 0.,    0.,
     &           0.263,    7.089,
     &           0.305,    5.346,
     &           0.288,    3.735,
     &           0.215,    2.225,
     &           0.088,    0.758,
*     &          -0.016,    0.116,
     &          -0.043,   -0.113,
*     &          -0.090,    0.584,
     &          -0.214,   -0.564,
     &           0.004,    0.024,
     *           0.146,    0.774,
*     &           0.375,    1.828, ! before the turn
*     &           0.343,    1.834, ! after the 1 degree turn
     &          -7.906,    3.602, ! after the 1 degree turn
*     &          -0.080,    0.362,
     &          -0.219,   -0.343,
     &           0.328,    1.364,
     &           0.434,    1.650,
     &           0.260,    0.906,
     &           0.082,    0.262,
     &           0.002,    0.005,
*     &          -0.258,    0.714/
     &          -0.627,    -0.696/
C
C     Ellips second focus
cccc      data x_2 / 7.151,  412.446,
      data x_2 / -91.06419,  213.29064,
     &          15.172,  409.393,
     &          23.272,  407.271,
     &          31.191,  405.128,
     &          38.930,  402.939,
     &          46.672,  400.634,
*     &          54.203,  398.796,
     &          54.207,  398.795,
*     &          60.942,  397.343,
     &          60.966,  397.343,
     &          67.997,  396.821,
     *          74.951,  396.344,
*     &          81.289,  396.193,  ! before turn
*     &          74.362,  397.551,   ! after 1 degree turn
     &          79.879,  396.495,   ! after 1 degree turn
*     &          85.387,  384.555,
     &          85.416,  384.560,
     &          92.399,  384.754,
     &          98.627,  375.156,
     &         106.813,  372.114,
     &         115.078,  368.942,
     &         123.882,  365.795,
*     &         131.037,  362.481/
     &         131.105,  362.564/
C
C     Semi- axes (large)
      data R_C / 259.9117,10*282.0, 2*289.5, 5*296.5 /
c      data R_C / 11*282.0, 2*289.5, 5*296.5 /
C
C ZEBRA and HBOOK initialization
C
      CALL MZEBRA(-1)
      CALL HLIMIT(10000)
      OPEN( UNIT=25, FILE='test.dat',STATUS='NEW',
     & ACCESS='SEQUENTIAL', FORM='FORMATTED' )
C
C=====   CYCLE on CC number
C
      do i = 1,18
C     ===========
C
C  Initialisation
c
      R     = R_C(i)
      x1(1) = x_1(1,i)
      x1(2) = x_1(2,i)
      x2(1) = x_2(1,i)
      x2(2) = x_2(2,i)
c
      df2 = (x1(1)-x2(1))*(x1(1)-x2(1)) + (x1(2)-x2(2))*(x1(2)-x2(2))
      write(*,*) df2
      ds  = (x2(2) - x1(2))/sqrt(df2) 
      write(*,*) ds
C      
      d2  = R*R - df2/4.
      p_in(1) = -unit/(R*R)
      p_in(2) = -unit/d2
      p_in(3) = 0.
      p_in(4) = 0.
      p_in(5) = 0.
C   =================  Equation in own system
C
C****   TURN  ****
C
      write(25,*) 'Parameters 1:',p_in
      angle = asin(ds)
      write(25,*) ' ANGLE: ',angle
      angled = angle
      call ccturn(angled,p_in,p_fn)
      write(25,*) 'Parameters 2:',p_fn
C
C****   SHIFT  ****
C
      xshift(1) = (x1(1) + x2(1))/2.
      xshift(2) = (x1(2) + x2(2))/2.
      write(25,*) ' SHIFT', xshift
      call ccshft(xshift,p_fn,p_ff)
      write(25,*) 'Parameters 3:'
      write(25,*) p_ff
C
C****   OUTPUT  ****
C
      write(*,*) 'New  mirror ****', i
      write(*,*) 'New parameters:', p_ff
C
      end do
C
      close (25)
      stop      
      end
      SUBROUTINE CCTURN(ANGLE,P_IN,P_FN)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : Calculation of new mirror parameters
C-                         after the turn on angle ANGLE
C-
C-   Inputs  : ANGLE, P_IN(5)
C-   Outputs : P_FN(5)
C-   Controls: 
C-
C-   Created  15-MAR-1994   Alexander V. Vlassov
C-
C----------------------------------------------------------------------
C
      double precision angle,p_in(5),p_fn(5)
      double precision c,s 
C
      c = cos(angle)
      s = sin(angle)
C
      p_fn(1) = p_in(1)*c*c + p_in(2)*s*s - p_in(3)*s*c
      p_fn(2) = p_in(1)*s*s + p_in(2)*c*c + p_in(3)*s*c
      p_fn(3) = 2.*(p_in(1) - p_in(2))*s*c + p_in(3)*(c*c - s*s)
      p_fn(4) = p_in(4)*c - p_in(5)*s
      p_fn(5) = p_in(5)*c + p_in(4)*s
C
  999 RETURN
      END
      SUBROUTINE CCSHFT(X,P_IN,P_FN)
C----------------------------------------------------------------------
C-
C-   Purpose and Methods : Calculation of new mirror parameters
C-                         after the SHIFT on vector X
C-
C-   Inputs  : X(2), P_IN(5)
C-   Outputs : P_FN(5)
C-   Controls: 
C-
C-   Created  22-MAR-1994   Alexander V. Vlassov
C-
C----------------------------------------------------------------------
C
      double precision x(2),p_in(5),p_fn(5)
      double precision d
C
      d = p_in(1)*x(1)*x(1) + p_in(2)*x(2)*x(2) + p_in(3)*x(1)*x(2)
     &  - p_in(4)*x(1) - p_in(5)*x(2) + 1.0
C     
      if(d.ne.0.) then
C 
        p_fn(1) = p_in(1)/d 
        p_fn(2) = p_in(2)/d
        p_fn(3) = p_in(3)/d
        p_fn(4) = (p_in(4) - 2.0*p_in(1)*x(1) - p_in(3)*x(2))/d
        p_fn(5) = (p_in(5) - 2.0*p_in(2)*x(2) - p_in(3)*x(1))/d
C
      else
        call ucopy(p_in,p_fn,5)
        write(*,*) '*** CCSHFT : *** INVALID parameters for shift'
        write(*,*) '*** CCSHFT : *** VECTOR ',x
        write(*,*) '*** CCSHFT : *** parameters',p_in
      end if
  999 RETURN
      END
