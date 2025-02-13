c This code computes zero temperature white
c dwarfs for the CD-ROM of Stellar Interiors, 2nd ed.
c It treats the surface radius as an eigenvalsue as
c was done for polytropes in Chapter 7.
c Carl Hansen, Sept 03
c *****************************************************
      PROGRAM WD
c$DEBUG
c******************************************************
c ******************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON X(400),YEIG(2,400),y0,eig,eigt,nsurf,verg,eps
      CHARACTER*20 OUTFILE
      DIMENSION amass(400),rout(400)
c
 3000 FORMAT (A20)
      WRITE (6,*)  ' YOU OUTPUT FILE IS '
      READ (5,3000) OUTFILE
      OPEN (11,FILE=OUTFILE,STATUS='UNKNOWN')
C
c    Set up x-grid  x between 0 and 1. eig is the outer boundary.
c
      nsurf=400
      dx=1.d0/399.d0
			 x(1)=0.d0
			 a=1.018466708d-2
			 do i=1,399
				x(i+1)=a*0.99d0**(i-1)+x(i)
			 end do
			 x(400)=1.00d0
   10 WRITE (6,2000)
 2000 FORMAT (' ENTER Central Density rho_0 (g/cm^3) and mue')
      WRITE (6,2001)
 2001 FORMAT (' (ENTER A rho_0 OF 0 TO STOP)')
      READ (5,*) rho0, amue
      IF (rho0.NE.0.0D0) GO TO 20
      WRITE (6,1000)
 1000 FORMAT (' CALCULATION COMPLETED')
      STOP
   20 WRITE (6,1010) rho0, amue
 1010 FORMAT (' rho0, mue =',2(1PE14.6))
      eig=4.3
c EIG is the first guess at the eigenvalue.
       x0=0.01009d0*(rho0/amue)**0.33333333d0
       y0=dsqrt(x0**2+1.0d0)
       alpha=7.77d8/(y0*amue*6.96d10)
       beta=5.741d33/(1.989d33*amue**2)
         NCONV = 0
         NSERCH=25
c
c Try to converge to a solution using Newton's method.
c
      VERG=1.00D-09
      EPS=1.00D-03
         DO 70 NTRY=1,NSERCH
            EIGT=EIG
               CALL GOOUT(DISC)
c
c Perturb eigenvalue and find change in boundary condition.
c DISC is the boundary condition. DISC should end up being
c 1/y0.
c
            EIGT=(1.D0+EPS)*EIG
               CALL GOOUT(DUM1)
     	    DEIGD=DUM1-DISC
          DEIG=(1/y0-DISC)*EPS*EIG/DEIGD
c            WRITE (6,1060) NTRY,EIG,DEIG
c           WRITE (6,1070) DISC,YEIG(1,NSURF)
c            WRITE (11,1060) NTRY,EIG,DEIG
c            WRITE (11,1070) DISC,YEIG(1,NSURF)
c 1070 FORMAT(6X,10H AND DISC=,1PE12.4,5H Y1= ,1PE12.4)
c 1060 FORMAT(' NTRY=',I3,' EIG=',1PE15.7,' DEIG=',
c     *         1PE13.5)
            ADEIG=DABS(DEIG/EIG)
c Check for convergence.
            IF (ADEIG.LT.VERG) NCONV=1
c Compute new guess for EIG.
            EPSEIG=DEIG/EIG
            IF (ADEIG.GT.0.05D0)
     *          EPSEIG=0.05D0*DEIG/DABS(DEIG)
            EIG=EIG*(1.D0+EPSEIG)
            IF (NCONV.EQ.1) GO TO 80
   70    CONTINUE
         WRITE (6,1080)
 1080 FORMAT (14H NOT CONVERGED)
         STOP
c
c Converged.
c
   80    EIGT=EIG
               CALL GOOUT(DISC)
c         WRITE ( 6,1090) EIG
c         WRITE (11,1090) EIG
c         WRITE (11,1100) DISC
c 1090 FORMAT (' FINAL VALUES:  EIG=',1PE17.9)
c 1100 FORMAT (28H BOUNDARY CONDITION IS DISC=,1PE11.3)
       totr=alpha*x(nsurf)*eig
       totmass=-beta*eig**2*x(nsurf)**2*yeig(2,nsurf)
 2099 FORMAT (' Total radius (Rsun)= ', 1PE12.4)
 2100 FORMAT (' Total mass (Msun)= ', 1PE12.4)
      write(11,2099) totr
      write(11,2100) totmass
	write(6,2099) totr
	write(6,2100) totmass
      write(11,1030)
      do i=1,nsurf
        amass(i)=-beta*eig**2*x(i)**2*yeig(2,i)
        rout(i)=x(i)*alpha*eig
         WRITE(11,1140) rout(i),amass(i)
      end do
 1030 FORMAT(' (r/Rsun ,  Mass/Msun)')
 1140 FORMAT(2(1PE12.4))
      GO TO 10
      END
c *****************************************************
      SUBROUTINE GOOUT(B2)
c This controls the integration.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON X(400),YEIG(2,400),y0,eig,eigt,nsurf,verg,eps
      DIMENSION Y(2),WORK(50),IWORK(5)
      EXTERNAL RKFCOW
      YEIG(1,1)=1.0D0
      YEIG(2,1)=0.0d0
c Expand around center in series to get 2nd & 3rd points.
c Expansion from Chandra, Chap XI, Eq 84.
      q=dsqrt((y0**2-1.0d0)/y0**2)
      t1=q**3*eigt**2/6.0d0
      t2=q**4*eigt**4/40.0d0
      t3=q**5*(5.0d0*q**2+14.0d0)*eigt**6/5040.0d0
      t4=q**6*(339.0d0*q**2+280.0d0)*eigt**8/1.08864d6
      t5=q**7*(1425.0d0*q**4+11436.0d0*q**2+4256.0d0)*eigt**10/1.9960d8
        x2=x(2)**2
      yeig(1,2)=1.0d0-x2*(t1-x2*(t2-x2*(t3-x2*(t4-t5*x2))))
      yeig(2,2)=-x(2)*(2.0d0*t1-x2*(4.0d0*t2-x2*(6.0d0*t3-
     *     x2*(8.0d0*t4-x2*10.0d0*t5))))/eigt
      x2=x(3)**2
      yeig(1,3)=1.0d0-x2*(t1-x2*(t2-x2*(t3-x2*(t4-t5*x2))))
      yeig(2,3)=-x(3)*(2.0d0*t1-x2*(4.0d0*t2-x2*(6.0d0*t3-
     *     x2*(8.0d0*t4-x2*10.0d0*t5))))/eigt
      Y(1)=YEIG(1,3)
      Y(2)=YEIG(2,3)
      IFLAG=1
      NEQN=2
      NSURF1=NSURF-1
   50 DO 20 J=3,NSURF1
         K=J+1
         XSTART=X(J)
         XFIN=X(K)
         RELERR=1.D-12
         ABSERR=DABS(Y(1))
      DO 30 I=2,NEQN
         TEMP=DABS(Y(I))
         ABSERR=DMIN1(ABSERR,TEMP)
   30 CONTINUE
         ABSERR=1.D-12*ABSERR
         ABSERR=DMAX1(ABSERR,1.D-20)
         CALL RKF(RKFCOW,NEQN,Y,XSTART,XFIN,RELERR,ABSERR,
     @        IFLAG,WORK,IWORK)
         IF (IFLAG.EQ.3.OR.IFLAG.EQ.4.OR.IFLAG.EQ.5)
     @              WRITE (11,1000) IFLAG,XFIN
         IFLAG=1
         DO 10 I=1,2
   10         YEIG(I,K)=Y(I)
   20 CONTINUE
c The following is the outer boundary condition. We
c iterate until B2 is close enough to zero that
c corrections to c the eigenvalue are within a small
c tolerance.
      B2=YEIG(1,NSURF)
 1000 FORMAT (17H WATCH OUT,IFLAG=,I4,6H XFIN=,1PE10.2)
      RETURN
      END
c ***************************************************
      SUBROUTINE RKFCOW(rt,Y,YP)
c
c This routine supplies the derivatives for use in RKF
c for the nonradial case.
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON X(400),YEIG(2,400),y0,eig,eigt,nsurf,verg,eps
      DIMENSION Y(2),YP(2)
      yp(1)=eigt*y(2)
	temp=y(1)**2-1/y0**2
      if (temp .lt. 1.0d-16) temp=0.d0
      yp(2)=-eigt*temp*dsqrt(temp)-2.d0*y(2)/rt
      RETURN
      END
c ****************************************************
C
      SUBROUTINE RKF(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK,IWORK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),WORK(50),IWORK(5)
      EXTERNAL F
      K1M=NEQN+1
      K1=K1M+1
      K2=K1+NEQN
      K3=K2+NEQN
      K4=K3+NEQN
      K5=K4+NEQN
      K6=K5+NEQN
      CALL RKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK(1),
     @    WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),
     @    WORK(K6),WORK(K6+1),IWORK(1),IWORK(2),IWORK(3),IWORK(4),
     @    IWORK(5))
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RKFS(F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,YP,H,F1,F2,F3,
     @         F4,F5,SAVRE,SAVAE,NFE,KOP,INIT,JFLAG,KFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL HFAILD,OUTPUT
      DIMENSION Y(NEQN),YP(NEQN),F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),
     @          F5(NEQN)
      EXTERNAL F
      DATA U26/2.D-13/ , REMIN/1.D-12/
      DATA MAXNFE/3000/
      IF (NEQN .LT. 1) GO TO 10
      IF ((RELERR .LT. 0.D0)  .OR.  (ABSERR .LT. 0.D0)) GO TO 10
      MFLAG=IABS(IFLAG)
      IF ((MFLAG .GE. 1) .AND. (MFLAG .LE. 7)) GO TO 20
   10 IFLAG=7
      RETURN
   20 IF (MFLAG .EQ. 1) GO TO 50
      IF (T .EQ. TOUT) GO TO 10
      IF(MFLAG .NE. 2) GO TO 25
      IF (INIT .EQ. 0) GO TO 45
      IF (KFLAG .EQ. 3) GO TO 40
      IF ((KFLAG .EQ. 4) .AND.  (ABSERR .EQ. 0.D0)) GO TO 30
      IF ((KFLAG .EQ. 5)  .AND. (RELERR .LE. SAVRE) .AND.
     @   (ABSERR .LE. SAVAE)) GO TO 30
      GO TO 50
   25 IF (IFLAG .EQ. 3) GO TO 40
      IF ((IFLAG .EQ. 4) .AND. (ABSERR .GT. 0.D0)) GO TO 45
   30 WRITE (6,1000) IFLAG,T
 1000 FORMAT (16H RKF SAYS IFLAG=,I5,3H X=,E12.4)
      STOP
   40 NFE=0
      IF (MFLAG .EQ. 2) GO TO 50
   45 IFLAG=JFLAG
   50 JFLAG=IFLAG
      KFLAG=0
      SAVRE=RELERR
      SAVAE=ABSERR
      RER=DMAX1(RELERR,REMIN)
      DT=TOUT-T
      IF (MFLAG .EQ. 1) GO TO 60
      IF (INIT .EQ. 0) GO TO 65
      GO TO 80
   60 INIT=0
      KOP=0
      A=T
      CALL F(A,Y,YP)
      NFE=1
      IF (T .NE. TOUT) GO TO 65
      IFLAG=2
      RETURN
   65 INIT=1
      YMAX=0.D0
      YPN=0.D0
      DO 70 K=1,NEQN
        YPN=DMAX1(ABS(YP(K)),YPN)
70      YMAX=DMAX1(ABS(Y(K)),YMAX)
      ETN=RER*YMAX+ABSERR
      H=ABS(DT)
      IF(ETN.GE.YPN*H**5) GO TO 80
      H=DMAX1((ETN/YPN)**0.2D0,U26*DMAX1(ABS(T),H))
80    H=SIGN(H,DT)
      IF (ABS(H).GE.ABS(DT)) KOP=KOP+1
      IF (KOP.NE.100) GO TO 85
      IFLAG=6
      RETURN
85    IF (ABS(DT).GT.U26*ABS(T))GO TO 95
      DO 90 K=1,NEQN
90      Y(K)=Y(K)+DT*YP(K)
      A=TOUT
      CALL F(A,Y,YP)
      NFE=NFE+1
      GO TO 300
95    OUTPUT=.FALSE.
      SCALE=2.D0/RER
      AE=SCALE*ABSERR
100   HFAILD=.FALSE.
      HMIN=U26*ABS(T)
      DT=TOUT-T
      IF (ABS(DT).GE.2.D0*ABS(H)) GO TO 200
      IF (ABS(DT).GT.ABS(H)/0.9D0) GO TO 150
      OUTPUT=.TRUE.
      H=DT
      GO TO 200
150   H=0.5D0*DT
200   IF (NFE.LE.MAXNFE) GO TO 220
      IFLAG=3
      KFLAG=3
      RETURN
220   CALL FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE=NFE+5
      EEOET=0.D0
      DO 250 K=1,NEQN
        ET=ABS(Y(K))+ABS(F1(K))+AE
        IF (ET.GT.0.D0) GO TO 240
        IFLAG=4
        KFLAG=4
        RETURN
240     EE=ABS((-2090.D0*YP(K)+(21970.D0*F3(K)-15048.D0*F4(K)))+
     1                       (22528.D0*F2(K)-27360.D0*F5(K)))
250     EEOET=DMAX1(EEOET,EE/ET)
      ESTTOL=ABS(H)*EEOET*SCALE/752400.D0
      IF (ESTTOL.LE.1.D0) GO TO 260
      HFAILD=.TRUE.
      OUTPUT=.FALSE.
      S=0.1D0
      IF (ESTTOL.LT.59049.D0)S=0.9D0/ESTTOL**0.2D0
      H=S*H
      IF (ABS(H).GT.HMIN) GO TO 200
      IFLAG=5
      KFLAG=5
      RETURN
260   T=T+H
      DO 270 K=1,NEQN
270     Y(K)=F1(K)
      A=T
      CALL F(A,Y,YP)
      NFE=NFE+1
      IF (HFAILD) GO TO 290
      S=5.D0
      IF (ESTTOL.GT.1.889568D-04) S=0.9D0/ESTTOL**0.2D0
      H=SIGN(DMAX1(S*ABS(H),HMIN),H)
290   IF (OUTPUT) GO TO 300
      IF (IFLAG.GT.0) GO TO 100
      IFLAG=-2
      RETURN
300   T=TOUT
      IFLAG=2
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE FEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(NEQN),YP(NEQN),F1(NEQN),F2(NEQN),F3(NEQN),F4(NEQN),
     @          F5(NEQN),S(NEQN)
      CH=0.25D0*H
      DO 10 K=1,NEQN
   10      F5(K)=Y(K)+CH*YP(K)
      CALL F(T+0.25D0*H,F5,F1)
      CH=0.09375D0*H
      DO 20 K=1,NEQN
   20      F5(K)=Y(K)+CH*(YP(K)+3.D0*F1(K))
      CALL F(T+0.375D0*H,F5,F2)
      CH=H/2197.D0
      DO 30 K=1,NEQN
   30      F5(K)=Y(K)+CH*(1932.D0*YP(K)+(7296.D0*F2(K)-7200.D0*F1(K)))
      CALL F(T+12.D0/13.D0*H,F5,F3)
      CH=H/4104.D0
      DO 40 K=1,NEQN
   40      F5(K)=Y(K)+CH*((8341.D0*YP(K)-845.D0*F3(K))+
     @                   (29440.D0*F2(K)-32832.D0*F1(K)))
      CALL F(T+H,F5,F4)
      CH=H/20520.D0
      DO 50 K=1,NEQN
   50    F1(K)=Y(K)+CH*((-6080.D0*YP(K)+(9295.D0*F3(K)-5643.D0*F4(K)))+
     @                             (41040.D0*F1(K)-28352.D0*F2(K)))
      CALL F(T+0.5D0*H,F1,F5)
      CH=H/7618050.D0
      DO 60 K=1,NEQN
   60     S(K)=Y(K)+CH*((902880.D0*YP(K)+(3855735.D0*F3(K)
     @         -1371249.D0*F4(K)))+(3953664.D0*F2(K)+277020.D0*F5(K)))
      RETURN
      END
C