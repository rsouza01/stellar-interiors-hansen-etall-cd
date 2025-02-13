          PROGRAM GETKAPPA
!$DEBUG
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL(8) KAPPA
      WRITE (6,1001)
          1001 Format(' INPUT T (K) THEN RHO (G/CM**3)')
      READ (5,*) TIN, RHO
      VIN=1.d0/RHO
      WRITE (6,1002)
	  1002 FORMAT(' INPUT X THEN Y')
      READ (5,*) X, Y
	  Z=1.d0-X-Y
      CALL EOS(X,Y,TIN,VIN,P,E,PE,PV,PT,ET,EV)
      CALL OPACITY(Y,Z,TIN,RHO,PE,KAPPA)
       WRITE (6,1003)
           1003 FORMAT(' Z=, KAPPA (CM**2/G)=:')
       WRITE(6,1000) Z,KAPPA
       1000 Format (1X,1P2D10.3)
               STOP ; END

      SUBROUTINE EOS (X,Y,TIN,VIN,P,E,PE,PV,PT,EV,ET)
! ***********************************************
! This routine and the other EOS routines are due
! to W. Dean Pesnell.
!
!     EQUATION OF STATE
!     ARGUMENTS...TIN (DEGREES K), VIN=1/RHO (CM**3/GM)
!     METALS...
!        NA,AL ALWAYS IONIZED
!        MG,SI,FE INCLUDED AS SINGLE ELEMENT
!        ALL OTHERS IGNORED
!
!          TABLE OF RETURNED QUANTITIES.
!
!      FUNCTION       NAME   DERIVATIVE WITH RESPECT TO
!                               TEMP.       SP. VOL.
!-------------------------------------------------------
!      PRESSURE     I    P I     PT      I     PV      I
!      INT. ENERGY  I    E I     ET      I     EV      I
!      ELEC. PRS.   I   PE I     PET     I     PEV     I
!      MOLAR DENSITYI   BP I     BPT     I     BPV     I
!-------------------------------------------------------
!
!
!          X = HYDROGEN MASS FRACTION
!          Y = HELIUM MASS FRACTION
!          Z = METALLIC MASS FRACTION
!          PARGRM = MEAN MOLECULAR MOLAR DENSITY WITHOUT
!                   ELECTRONS
!
!          R = GAS CONSTANT 8.31434E7
!          A = STEFAN-BOLTZMAN CONSTANT 7.56471E-15
!          BK = BOLTZMAN S CONSTANT 8.6170837E-5
!          AVAGD = AVAGADRO S NUMBER 6.02217E23
!          AD3 = A/3
!

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (ZERO = 0.D0,ONE = 1.D0, TWO=2.D0, &
       THRE=3.D0, FOR=4.D0, TEN=10.D0, AHF=0.5D0, &
       QRT=0.25D0 )
      PARAMETER ( R = 8.31434D7, A = 7.56471D-15, &
       BK = 8.6170837D-5, AVAGD = 6.02217D23, &
       AD3 = A/3.D0 )
      DATA T3OUT,T4OUT/1.665795163D-25,3.802592017D-28/
      DATA T2OUT,T5OUT/ 5.347896D-35,6.614536D-34/
!      IONIZATION POTENTIALS FOR HYDROGEN AND HELIUM
      DATA XH,XHE,XHE2/13.595D0,24.581D0,54.403D0/
      DATA C1,C2,C3/4.0092926D-9,1.00797D0,4.0026D0/
      DATA XM,CM,ZPZP/7.9D0,0.7D0,0.12014D0/
      DATA PREC /1.D-10/ ; DATA ONHLF/1.5D0/
      V = VIN ; T = TIN  ; Z=1-X-Y
      FRE = ZERO ;  ENT = ZERO
      PARGRM = X/C2 + Y/C3 ; RMUC = ONE/PARGRM
      RT=R*T ; TT4=T**4 ; TK=ONE/(T*BK); SQT=DSQRT(T)
             ! C1=ORIGINAL(C1(0.33334622))/R
      T1 = V*SQT**3*C1 ; T2 = T2OUT
      IF( T .GT. 2.D3 ) T2 = DEXP(-XH*TK)
          T3 = T3OUT
      IF( T .GT. 5.D3 ) T3 = DEXP(-XHE*TK)
          T4 = T4OUT
      IF( T .GT. 1.D4 ) T4 = DEXP(-XHE2*TK)
          T5 = T5OUT
      IF( T .GT. 1.2D3 ) T5 = DEXP(-XM*TK)
      D=T1*T2 ; B=FOR*T1*T3 ; C=B*T1*T4 ; DD=TWO*CM*T1*T5
      ZNA=Z*2.48D-3/24.969D0 ; ZMG=Z*ZPZP/45.807D0
!   CONVERGE ON ELECTRON DENSITY USING THE SAHA EQUATION.
!          GES IS THE MOLAR DENSITY OF ELECTRONS.
       GES = (X+Y*AHF)/(ONE+Y/(FOR*C))
      IF( GES .LT. X ) GES = AHF*(DSQRT(D*(D+FOR*X))-D)
      IF( GES .LT. 1.D-6*Z ) GES = 1.D-6*Z
      XC2 = X/C2 ; YC3 = Y/C3
!       NEWTON METHOD FOR ELECTRON DENSITY.
      DO 1 I=1,25
         T2 = C/GES+GES+B
         GEP = XC2*D/(GES+D)+YC3*(B+TWO*C/GES)/T2 &
          + ZMG*DD/(GES+DD) + ZNA
         T1 = ONE+XC2*D/(D+GES)**2+YC3/T2* &
          (TWO*C/GES**2+(B+TWO*C/GES)*(ONE-C/GES**2)/T2) &
          + ZMG*DD/(GES+DD)**2
         DGES = (GEP-GES)/T1 ; GES = GES+DGES
         IF( DABS(DGES)/GES .LT. PREC ) GOTO 3
   1  CONTINUE
      GOTO 12
   3  CONTINUE
!       ELECTRON PRESSURE
!      TOTLN = 1/MU = X/C2+Y/C4+Z/C3+GES
       PE=RT*GES/V ; TOTLN=PARGRM+GES ; XX=D/(GES+D)
      T2=GES+B+C/GES ; YY=B/T2 ; ZZ=C/(GES*T2) ; WW=DD/(GES+DD)
!          DERIVATIVES OF THE SAHA EQUATION
!          DERIVATIVES OF THE SAHA EQUATION FOR THE
!          PRESSURE AND INTERNAL ENERGY TEMPERATURE AND
!          DENSITY DERIVATIVES.

      T1 = YC3*(B+TWO*C/GES)
      QC0 = ONE+XC2*XX/(GES+D)+ZMG*WW/(GES+DD)+YC3/T2* &
        (TWO*C/GES**2+(B+TWO*C/GES)*(ONE-C/GES**2)/T2)
      QC1 = XC2*(ONE-XX)/(GES+D) ; QC4 = ZMG*(ONE-WW)/(GES+DD)
      QC2 = (YC3-T1/T2)/T2 ; QC3 = (YC3*TWO-T1/T2)/(GES*T2)
      QGV = (QC1*D+QC2*B+QC3*TWO*C+QC4*DD)/(QC0*V)
      QP1 = D*(ONHLF+XH*TK)/T ; QP2 = B*(ONHLF+XHE*TK)/T
      QP3 = C*(THRE+(XHE+XHE2)*TK)/T ; QP4 = DD*(ONHLF+XM*TK)/T
      QGT = (QC1*QP1+QC2*QP2+QC3*QP3+QC4* QP4)/QC0
!          ELECTRON PRESSURE DERIVATIVES.
      PET=ONE+QGT/GES ; PEV=-ONE+QGV/GES
!          PRESSURE DUE TO THE IDEAL GAS
      P=RT*TOTLN/V ; PT=P/T+RT*QGT/V ; PV=RT*QGV/V-P/V
      BP = R*TOTLN ; BPV = R*QGV ; BPT = R*QGT
!           ADD THE RADIATION PRESSURE
      P = P+AD3*TT4 ; PT = PT+FOR*AD3*TT4/T
!          IONIZATION ENERGY
      EI = (R/BK)*(XH*XX*XC2+YC3*(XHE*YY+(XHE+XHE2)*ZZ) &
        +ZMG*XM*WW + ZNA*5.524D0 )
!          TOTAL INTERNAL ENERGY
      E=ONHLF*RT*TOTLN+A*V*TT4+EI ; EV = T*PT-P
      QXT = ((ONE-XX)*QP1-XX*QGT)/(GES+D)
      DT2 = QGT*(ONE-C/GES**2)+QP2+QP3/GES
      QYT = (QP2-B*DT2/T2)/T2
      QZT = (QP3-C*QGT/GES-C*DT2/T2)/(T2*GES)
      QWT = ((ONE-WW)*QP4-WW*QGT)/(GES+DD)
      EIT = (R/BK)*(XH*QXT*XC2+YC3*(XHE*QYT+(XHE+XHE2)*QZT)+ &
        ZMG*XM*QWT)
      ET = ONHLF*R*(TOTLN+T*QGT)+FOR*A*V*TT4/T+EIT
      RETURN
  12  WRITE(6,102) T,V
 102  FORMAT(1H0,27HNO CONVERGENCE IN EOS    T=,1PE10.3,2X,&
         & 2HV=,E10.3)
      STOP ; END
! ***********************************************
      SUBROUTINE OPACITY (Y,Z,TIN,RHO,PE,KAPPA)
! ***********************************************
! NOTE: This opacity routine is the same as used in
! ZAMS.FOR except the input list has changed and
! no YGO(4) exists, etc.
! This routine computes the opacity (OP) given the
! density (RHO), temperature (T), Hydrogen mass fraction
! (X), Helium mass fraction (Y), and electron pressure (Pe
! from the equation of state routine called EOS).
! Taken from Stellingwerf 1975, ApJ, 195,
! 441., with corrections given in Stellingwerf 1975,
! ApJ, 199, 705.
! This routine is to be used for 0.6<X<0.8, 0.2<Y<0.4,
! and 0.001<Z<0.02.
! You may use it outside these bounds but, if you do, use
! with caution!
! Constructing these sorts of fits is not for the amateur.
! Try it!
! ***********************************************
! Temperature and volume in these fits are given in 10^4 K
! and 1/RHO.
! ***********************************************
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      REAL(8) KAPPA
! Set up variables.
      T4=TIN/1.0D4
      V=1.0D0/RHO
      V1=V**0.35D0
      V2=DSQRT(V1)
      Y1=6.294D-5-(6.0D-5*Y)
      Y2=3.53D6*Y-3.0447D5
      Z1=21.0D0*Z+0.979D0
      Z2=105.0D0*Z+0.895D0
! ***********************************************
!  Compute kappa  from equation D3 of Stellingwerf
!  1975 which involves some continued fractions
!  (temp1, temp2).
! ***********************************************
      TEMP1=Y1*V1*DSQRT(T4)*T4**3+1.0D0/(760.0D0*T4**5 &
        +316.0D0/V2)
      TEMP1=1.0D0/(10.0D0*T4**6+1.0D0/TEMP1)
      TEMP2=1780.0D0*DSQRT(T4)*T4*T4/Z1+1.0D0/(Z1*Y2/T4**10 &
        +2.13D-3*Z2*V2/DSQRT(T4)/T4**4)
      TEMP2=47.3D0/T4**8+1.0D0/TEMP2
      TEMP2=1.0D0/(4.0D3+1.0D0/TEMP2)
      KAPPA=PE*(4.819D-13*V/T4+TEMP1+TEMP2)
      RETURN
      END
! ***********************************************