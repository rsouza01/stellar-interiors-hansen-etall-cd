      program optest
c This program yields the total opacity (radiative combined
c with conductive) using Cox-Stewart tabulated opacities.
c The input is log1o(density), log1o(temperature), X and Y.

      implicit double precision (a-h,o-z)
      do 10 j=1,1000
        print *,' Enter log rho, log t'
        read *, dl,tl
        print *,' Enter x,y'
        read *,x,y
        z=1.d0-x-y
        c=z/10.d0
        o=z/8.d0
        call opacity(dl,tl,x,y,z,c,o,fkap,0,dum1,dum2)
        print *,fkap
10    continue
      stop
      end
      SUBROUTINE OPACITY(DL,TL,X,Y,Z,C,O,FKAP,IDERIV,DLKLT,DLKLD)
      implicit double precision (a-h,o-z)

C
C General suburoutine  to return the total opacity for a mixture with a
c composition specified by X,Y,C,O, and metallicity Z (=1-X-Y) at
c a given value of log10(rho) (=DL) and log10(T) (=TL).
c
c Returns opacity as fkap,
c   if IDERIV=1 then computes logarithmic derivatives DLKLD and DLKLT.
c
c First evaluate the radiative opacities by calling your favorite
c radiative opacity routine
      call csorad(dl,tl,x,y,z,c,o,fkrad,dlkrlt,dlkrld)
c and evaluate the conductive opacities with another favorite routine
c ***********************************
      call iocond(DL,TL,X,Y,FKCON)
c
      fkap=(fkcon*fkrad)/(fkcon+fkrad)
c
c logarithmic derivatives, if requested:
c
      if (ideriv.eq.1) then
        delta=1.d-4
c
        call iocond(dl,tl+delta,x,y,fkcon1)
        dlkclt=(dlog10(fkcon1)-dlog10(fkcon))/delta
        dlklt=fkap/fkrad*dlkrlt+fkap/fkcon*dlkclt
c
        call iocond(dl+delta,tl,x,y,fkcon1)
        dlkcld=(dlog10(fkcon1)-dlog10(fkcon))/delta
        dlkld=fkap/fkrad*dlkrld+fkap/fkcon*dlkcld
      endif
c
      return
      END
c
c
      SUBROUTINE CSORAD(DL,TL,X,Y,Z,C,O,FKRAD,DLKLTR,DLKLDR)
      implicit double precision (a-h,o-z)
c
c Interpolates within opacity tables of Cox and Stewart to obtain
c radiative opacity for a given log10(rho)=DL and log10(T)=TL and
c a mixture given by X,Y,C,O, and Z (where Z=1-X-Y)
c
      parameter(nz=5,nx=5,nt=29,nd=8,nzm=4,ntm=10,ndm=8)
      common/csop/csrad(nz,nx,nt,nd),tabd(nt,nd),tabt(nt),tabz(nz),
     1           tabx(nz,nx),
     2           csmet(nzm,ntm,ndm),tabdm(ntm,ndm),tabtm(ntm),
     3           tabzmy(nzm),tabzmc(nzm),tabzmo(nzm)
      save /csop/
      common/rdtab/iread
      dimension f(4,4),xf(4),yf(4),fk(4),fkx(4),fky(4),fkxy(4)
      dimension tempt(nd),temptm1(nd),temptp1(nd),temptp2(nd)
      dimension txzm1(nx),txzp1(nx),tx(nx)
      dimension fx1z1(4,4),fx2z1(4,4),fx2z2(4,4),fx1z2(4,4),zf(4)
      dimension dum(4)
c
c if table hasn't been read yet, then read in the opacity table
c read in opacity table header for normal stellar mixtures
c
      if (iread.ne.122) then
        call rdcsop(csrad,tabd,tabt,tabz,tabx,csmet,tabdm,tabtm,
     1              tabzmy,tabzmc,tabzmo)
        iread=122
      endif
c
c determine isotherm number
c
      call locate(tabt,nt,tl,it)
c
c determine Z table number
c
      call locate(tabz,nz,Z,iz)
c
c determine X table number for values of z surrounding the point
c
      do 10 i=1,nx
        tx(i)=tabx(iz,i)
        txzm1(i)=tabx(iz-1,i)
        txzp1(i)=tabx(iz+1,i)
10    continue
      call locate(tx,nx,X,ix)
      call locate(txzm1,nx,X,ixzm1)
      call locate(txzp1,nx,X,ixzp1)
c
c Determine density points: Note that for the isotherm above log T,
c the density point number is not the same as for the isotherm
c below log T.
c
      do 11 i=1,nd
        tempt(i)=tabd(it,i)
        temptm1(i)=tabd(it-1,i)
        temptp1(i)=tabd(it+1,i)
        temptp2(i)=tabd(it+2,i)
11    continue
      call locate(tempt,nd,DL,idt)
      call locate(temptm1,nd,DL,idtm1)
      call locate(temptp1,nd,DL,idtp1)
      call locate(temptp2,nd,DL,idtp2)
c
c
c      2-D interpolations to obtain radiative opacity.
c
c Labels for the corners of the interpolation rectangle are as follows:
c
c                 (1,2)---------(2,2)
c             ^     | 4         3 |
c           Y |     |             |
c                   | 1         2 |
c                 (1,1)---------(2,1)
c
c                       X ->
c
c Interpolations in log t and log rho (X->logT, Y->log rho)
c to determine opacity at (logT, log rho) for the values of Z and X
c that surround (Z,X) of the desired point.
c
      tlow=tabt(it)
      thi=tabt(it+1)
      dlow=tabd(it,idt)
      dhi=tabd(it,idt+1)
c
c Find the 4x4 matrices in the logT, log D plane that surround the
c input TL,DL point for the 4 compostions (X,Z) that surround the
c composition of the input point.
c
      call chunk(csrad,iz,ix,it,idt,idtp1,idtp2,idtm1,tabt,tabd,
     1            fx1z1,xf,yf)
      call chunk(csrad,iz,ix+1,it,idt,idtp1,idtp2,idtm1,tabt,tabd,
     1            fx2z1,dum,dum)
      call chunk(csrad,iz+1,ix+1,it,idt,idtp1,idtp2,idtm1,tabt,tabd,
     1            fx2z2,dum,dum)
      call chunk(csrad,iz+1,ix,it,idt,idtp1,idtp2,idtm1,tabt,tabd,
     1            fx1z2,dum,dum)
c Now interpolate in composition for the 16 (T,D) points
      xlow=tabx(iz,ix)
      xhi=tabx(iz,ix+1)
      zlow=tabz(iz)
      zhi=tabz(iz+1)
      call XZINTRP(fx1z1,fx2z1,fx2z2,fx1z2,x,z,xlow,xhi,zlow,zhi,f)
c	print *,((f(i,j),i=1,4),j=1,4)
c and use the composition-interpolated value FK of the opacity at the 16 (D,T)
c points to compute the desired opacity and its derivatives for the given point
c	print *,(fk(i),i=1,4)
      call SETINT(f,xf,yf,fk,fkx,fky,fkxy)
c	print *,(fk(i),i=1,4)
      call BCUINT(fk,fkx,fky,fkxy,tlow,thi,dlow,dhi,TL,DL,
     1            fkrad,dlkltr,dlkldr)
c
c convert opacity back into non-logarithmic form
c
      fkrad=10.0**fkrad
      return
      end
c
      SUBROUTINE LOCATE(xx,n,x,j)
      implicit double precision (a-h,o-z)
c
c given array XX of length N, and given a value X, returns a value J
c such that X is between XX(j) and XX(j+1).  XX must be monotonic,
c either increasing or decreasing.  J=0 or J=N is returned to indicate
c that X is out of range
      dimension xx(n)
      jl=0
      ju=n+1
 10   if (ju-jl.gt.1) then
        jm=(ju+jl)/2
        if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
          jl=jm
        else
          ju=jm
        endif
        go to 10
      endif
      j=jl
      return
      end
c
      subroutine chunk(c,i1,i2,iix,iiyx,iiyxp1,iiyxp2,iiyxm1,
     1                 tabx,taby,f,x,y)
      implicit double precision (a-h,o-z)
c
c Given the four dimensional array c, and the indices i1,i2,ix,iy
c that describe the interpolant, abstract a 4x4 square array (F) at the
c (fixed) values of i1 and i2, and with X and Y as the independent variables.
c The points of F surround the interpolating point, with the point
c lying within the box defined by indices (2,2),(3,2),(3,3),(2,3).
c F(2,2) is the value of c at ix,iy.
c
c In the process, if a subscript falls out of bounds, then use the last
c tabulated value for that point. This results in having zero
c derivatives at the boundary of the table.  Not the best solution, but at
c least whatever solution to the problem of skirting the edge of the table
c can be implemented in this subroutine alone.
c
c    Note that the subroutine LOCATE (called in the upper level routine
c to determine the values of i1,i2,ix,iyx for a given point) sets
c the subscript equal to zero if it falls off the low edge of the table.
c
c   Note also that we assume the temperature boundaries of the table should
c NEVER be crossed.  That is left to the upper level driver to check...
c need to use a Christy-like opacity in the case of low T.
c
      parameter(nz=5,nx=5,nt=29,nd=8)
      dimension c(nz,nx,nt,nd),f(4,4),x(4),y(4),tabx(nt),taby(nt,nd)
c
      iyxm1=iiyxm1
      iyx=iiyx
      iyxp1=iiyxp1
      iyxp2=iiyxp2
      ix=iix
c make sure the indices are in the table bounds.  If not, set them equal
c to the indices of the table top or bottom
      if (iyxm1.eq.0) iyxm1=1
      if (iyx.eq.0) iyx=1
      if (iyxp1.eq.0) iyxp1=1
      if (iyxp2.eq.0) iyxp2=1
      if (ix.eq.0) ix=1
      ixm1=max(ix-1,1)
      ixp1=min(ix+1,nt)
      ixp2=min(ix+2,nt)
c fill the F array, being careful at the edges of the table...
c
c Ben's fix: got to be careful at the low hydrogen edge; if i2 comes in
c at greater than 5, then that means Hydrogen was zero, and therefore it
c needs to be set back to 5 so that you don't go beyond that second index
c in the C array and really (potentially) muck things up...
c
      i2old=i2
      if (i2.eq.6) i2=5
c
      iy11=max(iyxm1-1,1)
      F(1,1)=c(i1,i2,ixm1,iy11)
      iy21=max(iyx-1,1)
      F(2,1)=c(i1,i2,ix,iy21)
      iy31=max(iyxp1-1,1)
      F(3,1)=c(i1,i2,ixp1,iy31)
      iy41=max(iyxp2-1,1)
      F(4,1)=c(i1,i2,ixp2,iy41)
c
      F(1,2)=c(i1,i2,ixm1,iyxm1)
      F(2,2)=c(i1,i2,ix,iyx)
      F(3,2)=c(i1,i2,ixp1,iyxp1)
      F(4,2)=c(i1,i2,ixp2,iyxp2)
c
      iy13=min(iyxm1+1,nd)
      F(1,3)=c(i1,i2,ixm1,iy13)
      iy23=min(iyx+1,nd)
      F(2,3)=c(i1,i2,ix,iy23)
      iy33=min(iyxp1+1,nd)
      F(3,3)=c(i1,i2,ixp1,iy33)
      iy43=min(iyxp2+1,nd)
      F(4,3)=c(i1,i2,ixp2,iy43)
c
      iy14=min(iyxm1+2,nd)
      F(1,4)=c(i1,i2,ixm1,iy14)
      iy24=min(iyx+2,nd)
      F(2,4)=c(i1,i2,ix,iy24)
      iy34=min(iyxp1+2,nd)
      F(3,4)=c(i1,i2,ixp1,iy34)
      iy44=min(iyxp2+2,nd)
      F(4,4)=c(i1,i2,ixp2,iy44)
c
c and now the values of the independent variables at the grid points
c being doubly careful about going over the edge, 'cuz you'll be
c dividing by the differences of these things soon... In that case
c you artificially extend the values of X and Y with the last table
c spacings...
c
      X(1)=tabx(ixm1)
      if (ix.eq.1) X(1)=tabx(1)-(tabx(2)-tabx(1))
      X(2)=tabx(ix)
      X(3)=tabx(ixp1)
      X(4)=tabx(ixp2)
      if (ix.eq.nt-1) then
        x(4)=tabx(nt)+(tabx(nt)-tabx(nt-1))
      elseif (ix.eq.nt) then
        x(4)=tabx(nt)+2.0*(tabx(nt)-tabx(nt-1))
        x(3)=tabx(nt)+(tabx(nt)-tabx(nt-1))
      endif
c
      if (iyx.eq.1) then
        Y(1)=taby(ix,1)-(taby(ix,2)-taby(ix,1))
      else
        Y(1)=taby(ix,iyx-1)
      endif
      Y(2)=taby(ix,iyx)
      if (iyx.eq.nd-1) then
        iy3=min(iyx+1,nd)
        Y(3)=taby(ix,nd)
        Y(4)=taby(ix,nd)+(taby(ix,nd)-taby(ix,nd-1))
      elseif (iyx.ge.nd) then
        Y(3)=taby(ix,nd)+(taby(ix,nd)-taby(ix,nd-1))
        Y(4)=taby(ix,nd)+2.0*(taby(ix,nd)-taby(ix,nd-1))
      else
        Y(3)=taby(ix,iyx+1)
        Y(4)=taby(ix,iyx+2)
      endif
c
      i2=i2old
      return
      end
c
      subroutine setint (f,xf,yf,fk,fkx,fky,fkxy)
      implicit double precision (a-h,o-z)
      dimension f(4,4),xf(4),yf(4),fk(4),fkx(4),fky(4),fkxy(4)
c
c takes a 4x4 array F of values of a function at mesh points xf(4),yf(4)
c and returns the value of the function at the four points surrounding
c the desired location in array FK, the X and Y derivatives at the same
c in FKX and FKY, and the cross derivative FKXY, at the same points
c
c values at grid points
      fk(1)=f(2,2)
      fk(2)=f(3,2)
      fk(3)=f(3,3)
      fk(4)=f(2,3)
c  derivatives w.r.t. x at grid points
      fkx(1)=(f(3,2)-f(1,2))/(xf(3)-xf(1))
      fkx(2)=(f(4,2)-f(2,2))/(xf(4)-xf(2))
      fkx(3)=(f(4,3)-f(2,3))/(xf(4)-xf(2))
      fkx(4)=(f(3,3)-f(1,3))/(xf(3)-xf(1))
c  derivatives w.r.t. y at grid points
      fky(1)=(f(2,3)-f(2,1))/(yf(3)-yf(1))
      fky(2)=(f(3,3)-f(3,1))/(yf(3)-yf(1))
      fky(3)=(f(3,4)-f(3,2))/(yf(4)-yf(2))
      fky(4)=(f(2,4)-f(2,2))/(yf(4)-yf(2))
c  cross derivatives at grid points
      fkxy(1)=(f(3,3)-f(3,1)-f(1,3)+f(1,1))/
     1           ((xf(3)-xf(1))*(yf(3)-yf(1)))
      fkxy(2)=(f(4,3)-f(4,1)-f(2,3)+f(2,1))/
     1           ((xf(4)-xf(2))*(yf(3)-yf(1)))
      fkxy(3)=(f(4,4)-f(4,2)-f(2,4)+f(2,2))/
     1           ((xf(4)-xf(2))*(yf(4)-yf(2)))
      fkxy(4)=(f(3,4)-f(3,2)-f(1,4)+f(1,2))/
     1           ((xf(3)-xf(1))*(yf(4)-yf(2)))
c
      return
      end
c
      subroutine XZINTRP(f11,f21,f22,f12,x,z,xlow,xhi,zlow,zhi,fout)
      implicit double precision (a-h,o-z)
      dimension f11(4,4),f21(4,4),f22(4,4),f12(4,4),fout(4,4),f(2,2)
      save f
c interpolation in (X,Z) to make an F array for bicubic interpolation
      do 10 it=1,4
        do 20 id=1,4
           f(1,1)=f11(it,id)
           f(2,1)=f21(it,id)
           f(2,2)=f22(it,id)
           f(1,2)=f12(it,id)
           call TWODLI(f,xlow,xhi,zlow,zhi,x,z,fout(it,id),dum,dum)
20      continue
10    continue
      return
      end
c
      SUBROUTINE TWODLI(FTAB,XLOW,XHI,YLOW,YHI,X,Y,F,DFDX,DFDY)
      implicit double precision (a-h,o-z)
c
c Performs two dimensional linear interpolation.
c Labels for the corners of the interpolation rectangle are as follows:
c
c                 (1,2)------------(2,2)  yhi
c             ^     | 4            3 |
c           Y |     |     +(x,y)     |
c                   |                |
c                   | 1            2 |
c                 (1,1)------------(2,1)  ylow
c                  xlow             xhi
c                       X ->
c FTAB(4) is the value of the function at the four grid points.
c Returns F, the function at (X,Y), DFDX,DFDY
c
      DIMENSION FTAB(4)
      T=(X-XLOW)/(XHI-XLOW)
      U=(Y-YLOW)/(YHI-YLOW)
      F=(1.0-T)*(1.0-U)*FTAB(1)+T*(1.0-U)*FTAB(2)+
     1              T*U*FTAB(3)+(1.0-T)*U*FTAB(4)
      DFDT=-(1.0-U)*FTAB(1)+(1.0-U)*FTAB(2)+
     1         U*FTAB(3)-U*FTAB(4)
      DFDX=DFDT/(XHI-XLOW)
      DFDU=-(1.0-T)*FTAB(1)-T*FTAB(2)+
     1         T*FTAB(3)+(1.0-T)*FTAB(4)
      DFDY=DFDU/(YHI-YLOW)
      RETURN
      END
c
      subroutine rdcsop(csrad,tabd,tabt,tabz,tabx,csmet,tabdm,tabtm,
     1                  tabzmy,tabzmc,tabzmo)
      implicit double precision (a-h,o-z)
      parameter(nz=5,nx=5,nt=29,nd=8,nzm=4,ntm=10,ndm=8)
      dimension csrad(nz,nx,nt,nd),tabd(nt,nd),tabt(nt),tabz(nz)
      dimension tabx(nz,nx)
      dimension csmet(nzm,ntm,ndm),tabdm(ntm,ndm),tabtm(ntm),
     1           tabzmy(nzm),tabzmc(nzm),tabzmo(nzm)
      dimension temp(nd)
c
c This subroutine reads in the Cox-Stewart opacity tables
c and converts the opacity values to log opacity for internal storage
c
c read in opacity table header
101   format (8F6.1)
102   format (8F8.5)
103   format (F8.5,5x,5f7.4)
      open (1,file='cs1970.dat',status='OLD')
      read (1,101) (tabd(i,1),i=1,nt)
      read (1,102) (tabt(i),i=1,nt)
c only the first density point is in the header; density spacing in this
c case is 1 in log10(d) so fill in the rest of the tabd(nt,nd) array
      do 3 it=1,nt
        do 4 id=2,nd
          tabd(it,id)=tabd(it,id-1)+1.0
4       continue
3     continue
c
      do 5 iz=1,nz
        read (1,103) tabz(iz),(tabx(iz,j),j=1,nx)
5     continue
c
c read in opacity data for normal stellar mixtures
c and convert opacities to log opacity
c
100   format(8E9.2)
      do 15 iz=1,nz
        do 25 ix=1,nx
          read (1,100)
          do 35 it=1,nt
            read (1,100) (temp(id),id=1,nd)
            do 45 id=1,nd
              csrad(iz,ix,it,id)=dlog10(temp(id))
45          continue
35        continue
25      continue
15    continue
c
c now read in opacity data for metallic stellar mixtures
c
      read (1,101) (tabdm(i,1),i=1,ntm)
      read (1,102) (tabtm(i),i=1,ntm)
      do 6 i=1,nzm
        read (1,103) tabzmy(i),tabzmc(i),tabzmo(i)
6     continue
c
c As above fill in the tabd(ntm,ndm) array since only first density point
c is given
      do 13 it=1,ntm
        do 14 id=2,ndm
          tabdm(it,id)=tabdm(it,id-1)+1.0
14      continue
13    continue
c
      do 16 izm=1,nzm
        read (1,100)
        do 36 itm=1,ntm
          read (1,100) (temp(idm),idm=1,ndm)
          do 46 idm=1,ndm
            csmet(izm,itm,idm)=dlog10(temp(idm))
46        continue
36      continue
16    continue
c
c Done reading table!
c
      return
      end
c
      SUBROUTINE BCUINT(Y,Y1,Y2,Y12,X1L,X1U,X2L,X2U,X1,X2,ANSY,ANSY1,
     *ANSY2)
      implicit double precision (a-h,o-z)
c
c Bicubic interpolation within a grid square.  Input quantities are Y,
c Y1, Y2, and Y12 (as described in BCUCOF); X1L and X1U, the lower and
c upper coordinates of the grid square in the abscissa direction; X2L and X2U
c likewise for the ordinate-direction; and X1,X2, the coordinates of the desired
c point for the interpolation.  The interpolated function value is
c returned as ANSY and the interpolated gradient values as ANSY1 and
c ANSY2.  This routine requires BCUCOF
c
c Numerical Recipes, p. 98-100
c
      DIMENSION Y(4),Y1(4),Y2(4),Y12(4),C(4,4)
c      print *,(y(i),i=1,4)
      CALL BCUCOF(Y,Y1,Y2,Y12,X1U-X1L,X2U-X2L,C)
      IF(X1U.EQ.X1L.OR.X2U.EQ.X2L)PAUSE 'bad input'
      T=(X1-X1L)/(X1U-X1L)
      U=(X2-X2L)/(X2U-X2L)
      ANSY=0.
      ANSY2=0.
      ANSY1=0.
      DO 11 I=4,1,-1
        ANSY=T*ANSY+((C(I,4)*U+C(I,3))*U+C(I,2))*U+C(I,1)
        ANSY2=T*ANSY2+(3.*C(I,4)*U+2.*C(I,3))*U+C(I,2)
        ANSY1=U*ANSY1+(3.*C(4,I)*T+2.*C(3,I))*T+C(2,I)
11    CONTINUE
      ANSY1=ANSY1/(X1U-X1L)
      ANSY2=ANSY2/(X2U-X2L)
      RETURN
      END
C
      SUBROUTINE BCUCOF(Y,Y1,Y2,Y12,D1,D2,C)
      implicit double precision (a-h,o-z)
c
c Given arrays Y, Y1, Y2, and Y12, each of length 4, containing the
c funcition, gradients, and cross derivative at the four grid points of a
c rectangular grid cell (numbered counterclockwise from the lower left),
c and given D1 and D2, the length of the grid cell in the X and Y
c directions, this routine returns the table C that is used by routine
c BCUINT for bicubic interpolation.
c
c Numerical Recipes, p. 98-100
c
      DIMENSION C(4,4),Y(4),Y1(4),Y2(4),Y12(4),CL(16),X(16),WT(16,16)
      SAVE wt
      DATA WT/1.,0.,-3.,2.,4*0.,-3.,0.,9.,-6.,2.,0.,-6.,
     *  4.,8*0.,3.,0.,-9.,6.,-2.,0.,6.,-4.,10*0.,9.,-6.,
     *  2*0.,-6.,4.,2*0.,3.,-2.,6*0.,-9.,6.,2*0.,6.,-4.,
     *  4*0.,1.,0.,-3.,2.,-2.,0.,6.,-4.,1.,0.,-3.,2.,8*0.,
     *  -1.,0.,3.,-2.,1.,0.,-3.,2.,10*0.,-3.,2.,2*0.,3.,
     *  -2.,6*0.,3.,-2.,2*0.,-6.,4.,2*0.,3.,-2.,0.,1.,-2.,
     *  1.,5*0.,-3.,6.,-3.,0.,2.,-4.,2.,9*0.,3.,-6.,3.,0.,
     *  -2.,4.,-2.,10*0.,-3.,3.,2*0.,2.,-2.,2*0.,-1.,1.,
     *  6*0.,3.,-3.,2*0.,-2.,2.,5*0.,1.,-2.,1.,0.,-2.,4.,
     *  -2.,0.,1.,-2.,1.,9*0.,-1.,2.,-1.,0.,1.,-2.,1.,10*0.,
     *  1.,-1.,2*0.,-1.,1.,6*0.,-1.,1.,2*0.,2.,-2.,2*0.,-1.,1./
      D1D2=D1*D2
c pack a temporary vector X
      DO 11 I=1,4
        X(I)=Y(I)
        X(I+4)=Y1(I)*D1
        X(I+8)=Y2(I)*D2
        X(I+12)=Y12(I)*D1D2
11    CONTINUE
c matrix multiply by the stored table
      DO 13 I=1,16
        XX=0.
        DO 12 K=1,16
          XX=XX+WT(I,K)*X(K)
12      CONTINUE
        CL(I)=XX
13    CONTINUE
      L=0
c unpack the result into the output table
      DO 15 I=1,4
        DO 14 J=1,4
          L=L+1
          C(I,J)=CL(L)
14      CONTINUE
15    CONTINUE
      RETURN
      END
      subroutine iocond(dlin,tlin,xin,yin,ochl)
      implicit double precision (a-h,o-z)
c
c Returns the conductive opacity for a mixture (X,Y,Z=1-X-Y) at
c a given value of log density and log temperature.  Uses the Iben
c fit to the Hubbard and Lampe opacity tables (Ap. J. 196, 545).
c Implementation by S.Kawaler, August 1990.  Checked against
c actual H&L tables for helium and carbon. Well within 20%
c agreement for log rho>0, can be way off at lower densities, so
c caveat emptor...
c
      dl=dlin
      tl=tlin
      x=xin
      y=yin
      cln=dlog(10.d0)
c
c   If log rho is greater than 6, use Lamb(1974)'s extrapolation
c   of the tables for the conductive opacity when relativistically
c   degenerate...  Though I suppose you could use the Itoh et al.
c   opacities in this regime...
c
      if (dlin.gt.6.d0) then
        ochl1=-1.d0
        dl=6.d0
      endif
c   coarse approximation to <Z^2/A> and ue:
      ue=1.d0/(0.5d0*(1.d0+X))
      avz2oa=x+y+(1.d0-x-y)*3.d0
c (the statement number here is for branching to when dlin > 6)
 10   fldelt=dl-dlog10(ue)-1.5d0*(tl-6.d0)
c      delta=10.d0**fldelt
      delta=dexp(cln*fldelt)
      fleta=-0.55255+(2.d0/3.d0)*fldelt
c      eta=10.d0**fleta
      eta=dexp(cln*fleta)
c Iben equations (A2)-(A5)
      if (fldelt.lt.0.645d0) then
        fltokt=-3.2862d0+dlog10(delta*(1.d0+0.024417*delta))
      else
        a1=-3.29243d0+dlog10(delta*(1.d0+0.02804*delta))
        if (fldelt.lt.2.0d0) then
          fltokt=a1
        else
          b1=-4.80946d0+dlog10(delta*delta*(1.d0+9.376d0/eta/eta))
          if (fldelt.gt.2.5d0) then
            fltokt=b1
          else
            fltokt=2.d0*a1*(2.5d0-fldelt)+2.d0*b1*(fldelt-2.d0)
          endif
        endif
      endif
c (A8) thru (A10) for Pe/NkT
      a2=dlog10(1.d0+0.021876*delta)
      if (fldelt.lt.1.5d0) then
        peonkt=a2
      else
        b2=dlog10(0.4d0*eta*(1.d0+4.1124d0/eta/eta))
        if (fldelt.gt.2.0d0) then
          peonkt=b2
        else
          peonkt=2.d0*a2*(2.d0-fldelt)+2.0d0*b2*(fldelt-1.5d0)
        endif
      endif
c (A11) and (A12) for nednedeta
      if (delta.lt.40.d0) then
        dlnede=1.d0-0.01d0*delta*(2.8966d0-0.034838d0*delta)
      else
        dlnede=(1.5d0/eta)*(1.d0-0.8225/eta/eta)
      endif
c and (A7) for alfa
      temp=dlog10(ue*avz2oa+dlnede)
      alfa=-2.033983d0+fldelt-0.5d0*(tl-6.d0)-peonkt+temp
c
c Now calculate the value of Theta
c
c for Hydrogen
      if (x.le.0.d0) then
        flthx=0.d0
      else
        if (alfa.le.-3.d0) then
          flthx=1.048d0-0.124d0*alfa
        elseif (alfa.gt.-1.d0) then
          flthx=0.185d0-0.585d0*alfa
c note that .585 in above is from Iben's code; paper says 0.558
        else
          flthx=0.13d0-alfa*(0.745d0+0.105d0*alfa)
        endif
      endif
c for Helium
      if (y.le.0.d0) then
        flthy=0.d0
      else
        if (alfa.le.-3.d0) then
          flthy=0.937d0-0.111d0*alfa
        elseif (alfa.gt.0.d0) then
          flthy=0.24d0-0.6d0*alfa
        else
          flthy=0.24d0-alfa*(0.55d0+0.0689d0*alfa)
        endif
      endif
c for Carbon
      if (alfa.lt.-2.5d0) then
        flthc=1.27d0-0.1d0*alfa
      elseif (alfa.gt.0.5d0) then
        flthc=0.843d0-0.785d0*alfa
      else
        flthc=0.727d0-alfa*(0.511d0+0.0778d0*alfa)
      endif
c
      Zhere=1.d0-x-y
c      ochl=x*10.d0**flthx+y*10.d0**flthy+zhere*10.d0**flthc
      ochl= x*dexp(cln*flthx)+y*dexp(cln*flthy)+zhere*dexp(cln*flthc)
c      ochl=ochl/10.d0**(tl-6.d0)/10.d0**fltokt
c      ochl=ochl/10.d0**(tl-6.d0+fltokt)
      ochl=ochl/dexp(cln*(tl-6.d0+fltokt))
c
c If input density greater than 10^6, then extrapolate H&L opacities
c linearly using the opacity at 10^5 and 10^6
c
      if (dlin.gt.6.d0.and.ochl1.le.0.d0) then
         ochl1=ochl
         dl=5.d0
         go to 10
      elseif (dlin.gt.6.d0.and.ochl1.gt.0.d0) then
        flochl=dlog10(ochl1)+(dlog10(ochl1)-dlog10(ochl))*(dlin-6.d0)
        ochl=10.d0**flochl
        return
      endif
c
      return
      end