      subroutine roplat(a,alat,alpha,b,beta,c,gamma,isym,platcv,platro)
C- Make platro from lattice parameters
C ----------------------------------------------------------------------
Ci Inputs:
Ci   isym  :numbers characterizing the symmetry of lattice
Ci   a     :lattice parameter, always required
Ci   b     :lattice parameter, required unless isym1.gt.3
Ci         :                   (triclinic,monoclinic,orthorhombic)
Ci   c     :lattice parameter, required unless isym1.gt.6 or
Ci         :                   isym1.eq.5.and.isym3.eq.1
Ci         :                   (cubic or ?triple hexagonal cell?)
Ci   alpha :angle(b,c), in degrees, required for isym1=1
Ci   beta  :angle(c,a), in degrees, required for isym1=1,2
Ci   gamma :angle(a,b), in degrees, required for isym1=1
Co Outputs:
Co   alat  :length scale
Co   platcv:lattice vectors of conventional unit cell
Co   platro:lattice vectors of standard unit cell.
Co          They can be rotated into platst using stplat.f
Cr Remarks:
Cr Bravais lattice:                          platro
Cr ---------------------------------------------------------------------
Cr                            ( ... ,      ...       ,     ...        )
Cr triclinic:                 (  0  ,       1        ,      0         )
Cr                            (  0  , c/b*cos(alpha) , c/b*sin(alpha) )
Cr
Cr                            ( a/b*sin(beta) , 0 ,  a/b*cos(beta)   )
Cr monoclinic primitive:      (      0        , 1 ,      0           )
Cr b-unique axis              (      0        , 0 ,     c/b          )
Cr
Cr                            ( a/2b*sin(beta) , -1/2 , a/2b*cos(beta) )
Cr monoclinic C-centred:      ( a/2b*sin(beta) ,  1/2 , a/2b*cos(beta) )
Cr b-unique axis              (        0       ,   0  ,      c/b       )
Cr
Cr                            ( a/b*sin(beta) ,  0  , a/b*cos(beta) )
Cr monoclinic A-centred:      (       0       , 1/2 ,    -c/2b      )
Cr b-unique axis              (       0       , 1/2 ,     c/2b      )
Cr
Cr                         (-a/2b*sin(beta) , 1/2 ,-a/2b*cos(beta)+c/2b)
Cr monoclinic I-centred:   ( a/2b*sin(beta) ,-1/2 , a/2b*cos(beta)+c/2b)
Cr b-unique axis           ( a/2b*sin(beta) , 1/2 , a/2b*cos(beta)-c/2b)
Cr
Cr                            ( a/b , 0 ,  0  )
Cr orthorhombic primitive:    (  0  , 1 ,  0  )
Cr                            (  0  , 0 , c/b )
Cr
Cr                            ( a/2b , -1/2 ,  0  )
Cr orthorhombic C-centred:    ( a/2b ,  1/2 ,  0  )
Cr                            (  0   ,   0  , c/b )
Cr
Cr                            ( a/2b ,  0  , -c/2b )
Cr orthorhombic B-centred:    (  0   ,  1  ,   0   )
Cr                            ( a/2b ,  0  ,  c/2b )
Cr
Cr                            ( a/b ,  0  ,   0   )
Cr orthorhombic A-centred:    (  0  , 1/2 , -c/2b )
Cr                            (  0  , 1/2 ,  c/2c )
Cr
Cr                            ( -a/2b ,  1/2 ,  c/2b)
Cr orthorhombic body-centred: (  a/2b , -1/2 ,  c/2b)
Cr                            (  a/2b ,  1/2 , -c/2b)
Cr
Cr                            (  0   , 1/2 , c/2b)
Cr orthorhombic face-centred: ( a/2b ,  0  , c/2b)
Cr                            ( a/2b , 1/2  , 0  )
Cr
Cr                            ( 1 , 0 , 0   )
Cr tetragonal primitive:      ( 0 , 1 , 0   )
Cr                            ( 0 , 0 , c/a )
Cr
Cr                            ( -1/2 ,  1/2 ,  c/a )
Cr tetragonal body-centred:   (  1/2 , -1/2 ,  c/a )
Cr                            (  1/2 ,  1/2 , -c/a )
Cr
Cr                            (   1  , 0 , c/a)
Cr rhombohedral:              ( -1/2 , C , c/a)     C=dsqrt(3)/2
Cr                            ( -1/2 , C , c/a)
Cr
Cr                            ( C , -1/2 , 0  )
Cr hexagonal:                 ( 0 ,   1 ,  0  )     C=dsqrt(3)/2
Cr                            ( 0 ,   0 , c/a )
Cr
Cr                            ( 1 , 0 , 0 )
Cr cubic primitive:           ( 0 , 1 , 0 )
Cr                            ( 0 , 0 , 1 )
Cr
Cr                            ( -1/2 ,  1/2 ,  1/2)
Cr cubic body-centred:        (  1/2 , -1/2 ,  1/2)
Cr                            (  1/2 ,  1/2 , -1/2)
Cr
Cr                            (  0  , 1/2 , 1/2)
Cr  cubic face-centred:       ( 1/2 ,  0  , 1/2)
Cr                            ( 1/2 , 1/2 ,  0 )
Cr
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer isym(*)
      double precision a,alat,b,c,alpha,beta,gamma,platcv(9),platro(3,*)
C Local variables:
      integer isym1,iprint,isym3,k,lgunit,m,stdo
      double precision a1,b1,c1,al,be,ga,facdeg,pi,qlatro(3,3),
     .tb(9,7),parm,vol,wk
      parameter(pi=3.14159265358979324d0,facdeg=180.d0/pi,parm=1.d-4)
C External calls:
Ckino 2011.05 del external,       external dmpy,iprint,lgunit,query
C Data statements:
      data tb/1.0d0, .0d0, .0d0,  .0d0, 1.0d0, .0d0,  .0d0,.0d0, 1.0d0,
     .0.5d0,-.5d0, .0d0,  .5d0, 0.5d0, .0d0,  .0d0,.0d0, 1.0d0,
     .0.5d0, .0d0,-.5d0,  .0d0, 1.0d0, .0d0,  .5d0,.0d0, 0.5d0,
     .1.0d0, .0d0, .0d0,  .0d0, 0.5d0,-.5d0,  .0d0,.5d0, 0.5d0,
     .-0.5d0, .5d0, .5d0,  .5d0,-0.5d0, .5d0,  .5d0,.5d0,-0.5d0,
     .0.0d0, .5d0, .5d0,  .5d0, 0.0d0, .5d0,  .5d0,.5d0, 0.0d0,
     .0.6666666666667d0,  0.3333333333333d0,  0.3333333333333d0,
     .-0.3333333333333d0,  0.3333333333333d0,  0.3333333333333d0,
     .-0.3333333333333d0, -0.6666666666667d0,  0.3333333333333d0/

      stdo = lgunit(1)
      isym1=isym(1)
      isym3=isym(3)
      a1=a
      b1=b
      c1=c
      al=alpha/facdeg
      be=beta /facdeg
      ga=gamma/facdeg

      if (iprint().ge.120) write(stdo,300) isym1,isym3
      if (iprint().ge.120) write(stdo,301) a,b,c,alpha,beta,gamma

      if (isym1.eq.5.and.isym3.eq.1) then
C ----- rhombohedral primitive cell will be obtained from
C ----- triple hexagonal cell
        c1=a1*dsqrt(3.d0+6.d0*dcos(al))
        a1=a1*dsqrt(2.d0-2.d0*dcos(al))
        isym3=7
      endif
      if (isym1.gt.3) b1=a1
      if (isym1.gt.6) c1=a1
      if (isym1.gt.1) al = 90.d0/facdeg
      if (isym1.gt.2) be = 90.d0/facdeg
      if (isym1.gt.1) ga = 90.d0/facdeg
      if (isym1.eq.5.or.isym1.eq.6) ga =120.d0/facdeg

C --- Sanity checks ---
      call rxx(a1.le.0,'roplat: unspecified  `a'' lattice parameter')
      call rxx(b1.le.0,'roplat: unspecified  `b'' lattice parameter')
      call rxx(c1.le.0,'roplat: unspecified  `c'' lattice parameter')
      call rxx(al.le.0,'roplat: unspecified angle `alpha''')
      call rxx(be.le.0,'roplat: unspecified angle `beta''')
      call rxx(ga.le.0,'roplat: unspecified angle `gamma''')

C --- Choose b as length unit
      alat=b1
      if (isym1.eq.5) alat=alat/dsqrt(3.d0)

      if (a1.lt.parm.or.b1.lt.parm.or.c1.lt.parm.or.
     .dabs(dsin(al)).lt.parm.or.dabs(dsin(be)).lt.parm.or.
     .dabs(dsin(ga)).lt.parm)
     .call rx('LATPAR: bad lattice parameters')

      wk=1.d0-dcos(al)**2-dcos(be)**2-dcos(ga)**2
     .+2.d0*dcos(al)*dcos(be)*dcos(ga)
      if (wk.le.0d0) call rx('ROPLAT: impossible angles')
      platcv(1)=a1/alat*dsqrt(wk)/dsin(al)
      platcv(2)=a1/alat*dcos(ga)
      platcv(3)=a1/alat*(dcos(be)-dcos(al)*dcos(ga))/dsin(al)
      platcv(4)=0.d0
      platcv(5)=b1/alat
      platcv(6)=0.d0
      platcv(7)=0.d0
      platcv(8)=c1/alat*dcos(al)
      platcv(9)=c1/alat*dsin(al)

      call dmpy(platcv,3,1,tb(1,isym3),3,1,platro,3,1,3,3,3)
C ----take rhombohedral primitive cell  not hexagonal triple cell
      if(isym(1).eq.5.and.isym(3).eq.1) call dcopy(9,platro,1,platcv,1)

      call dinv33(platro,1,qlatro,vol)
      if (iprint().ge.50) write(stdo,302)alat,
     .((platro(m,k),m=1,3),(qlatro(m,k),m=1,3),k=1,3)
      if (iprint().ge.60) write(stdo,304) (platcv(m),m=1,9)

  300 format(/' ROPLAT: isym1,isym3=',2i3)
  301 format(/' ROPLAT:      a=',f8.5,' ,   b=',f8.5,' ,    c=',f8.5,
     ./'          alpha=',f8.4,' ,beta=',f8.4,' ,gamma=',f8.4)
  302 format(/' ROPLAT: rotated standard unit cell: alat=',f10.6,
     ./17x,'Platro',32x,'Qlatro',3(/2x,3f11.7,5x,3f11.7))
  304 format(17x,'Platcv'/3(2x,3f11.7/))
      end
      subroutine cpplat(platin,platcp)
C- Finds the most compact unit cell
C ----------------------------------------------------------------------
Ci Inputs:
Ci   platin:any primitive lattice vectors
Ci Outputs:
Co   platcp:lattice vectors of most compact primitive unit cell
Cr Remarks:
Cr  The unit cell with the same volume as the original cell and
Cr  with minimal |p1|*|p2|*|p3|*|q1|*|q2|*|q3| is seeked.
Cr  where |pm| = dsqrt(plat(1,m)**2+plat(2,m)**2+plat(3,m)**2), m=1,2,3
Cr  and   |pn| = dsqrt(qlat(1,n)**2+qlat(2,n)**2+qlat(3,n)**2), n=1,2,3
Cr  platin and platcp can point to the same address.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision platin(3,*)
C Local variables:
      integer ikd,ikd1,ikd2,ikd3,iprint,k,ll1,ltmax,lgunit,m,mn,nkd,
     .nplat,stdo
      parameter(ltmax=1,nplat=3,ll1=ltmax*2+1)
      double precision d1cd2(3),d1mach,d2,d2min,danrm2,ddot,
     .dlat(3,ll1**nplat),plat(3,3),platcp(3,3),
     .pr1,pr2,prn1,prn2,pro1,pro2,qlat(3,3),
     .qlatcp(3,3),qlatin(3,3),tollat,vol,voln
      logical change,lmorec
      character*72 messg
      parameter(tollat=1.d-5)
C External calls:
Ckino 2011.05 del external,       external  cross,d1mach,danrm2,dcopy,ddot,dinv33,
Ckino 2011.05 del external,      .iprint,lgunit,prodln
C Intrinsic functions:
      intrinsic  dabs,dsqrt,idnint,mod
C Statement functions:
      mn(ikd,m)=(mod(ikd+(ll1/2)*ll1**(m-1),ll1**m)
     .-mod(ikd+(ll1/2)*ll1**(m-1),ll1**(m-1)))
     ./ll1**(m-1)-ltmax

      call dcopy(9,platin,1,platcp,1)
      call dinv33(platin,1,qlatin,vol)
      vol=dabs(vol)
      call prodln(platin,qlatin,vol,pro1,pro2)
      prn1=pro1
      prn2=pro2
      change=.true.
      lmorec=.false.
      nkd=ll1**nplat-1
      stdo = lgunit(1)

      if (iprint().ge.40) write(stdo,300)
     .pro1,pro2,((platin(m,k),m=1,3),(qlatin(m,k),m=1,3),k=1,3)
      if (iprint().gt.50) write(stdo,303)
     .(1.d0/dsqrt(ddot(3,qlatin(1,m),1,qlatin(1,m),1)),m=1,3)

C --- generate relevant lattice vectors
      do while (change)
        d2min=d1mach(2)
        change=.false.
        do ikd=1,nkd
          do k = 1, 3
            dlat(k,ikd)=0.d0
            do m = 1, nplat
              dlat(k,ikd)=dlat(k,ikd)+mn(ikd,m)*platcp(k,m)
            enddo
          enddo
          d2 = danrm2(dlat(1,ikd))
          if (d2.lt.d2min) then
            d2min=d2
            ikd1=ikd
          endif
        enddo

        do ikd2=1,nkd
          call cross(dlat(1,ikd1),dlat(1,ikd2),d1cd2)
          do ikd3=ikd2+1,nkd
            voln = dabs(ddot(3,d1cd2,1,dlat(1,ikd3),1))
            if (idnint(voln/vol).eq.1) then
              call dcopy(3,dlat(1,ikd1),1,plat(1,1),1)
              call dcopy(3,dlat(1,ikd2),1,plat(1,2),1)
              call dcopy(3,dlat(1,ikd3),1,plat(1,3),1)
              call dinv33(plat,1,qlat,voln)
              call prodln(plat,qlat,vol,pr1,pr2)
c             if ((pr1.lt.prn1-tollat).or.
c    .            (pr1*pr2.lt.prn1*prn2-d1mach(3))) then
              if(pr1*pr2.lt.prn1*prn2-d1mach(3)) then
                lmorec=.true.
                change=.true.
                prn1=pr1
                prn2=pr2
                call dcopy(9,plat,1,platcp,1)
              endif
            endif
          enddo
        enddo
      enddo

      call dinv33(platcp,1,qlatcp,voln)

C --- Printout
      if (lmorec) then
        write(messg,400)
        if (iprint().ge.40) write(stdo,301)
     .  prn1,prn2,((platcp(m,k),m=1,3),(qlatcp(m,k),m=1,3),k=1,3)
        if (iprint().gt.50) write(stdo,303)
     .  (1.d0/dsqrt(ddot(3,qlatcp(1,m),1,qlatcp(1,m),1)),m=1,3)

      else
        if (iprint().ge.40) write(stdo,302)
C        if (iprint().ge.50) write(stdo,301)
C     .    prn1,prn2,((platcp(m,k),m=1,3),(qlatcp(m,k),m=1,3),k=1,3)
      endif

  300 format(/' CPPLAT:   P1*P2*P3/vol=',f8.5,9x,'Q1*Q2*Q3*vol=',f8.5,
     ./17x,'Platin',32x,'Qlatin ',3(/2x,3f11.7,5x,3f11.7))
  301 format( '           P1*P2*P3/VOL=',f8.5,9x,'Q1*Q2*Q3*VOL=',f8.5,
     ./17x,'Platcp',32x,'Qlatcp',3(/2x,3f11.7,5x,3f11.7))
  302 format( '    The supplied unit cell is the most compact one.')
  303 format( '    Plane distances:',3f10.5)
  400 format(' CPPLAT: A more compact unit cell exists.$')
      end
      subroutine prodln(plat,qlat,vol,p1,p2)
C- Takes the product of the lenght of plat and qlat
C ----------------------------------------------------------------------
Ci Inputs:
Ci   plat  :primitive lattice vectors
Ci   qlat  :primitive reciprocal lattice vectors
Ci   vol   :volume of unit cell
Co Outputs:
Co   p1    :product of lenght of primitive lattice vectors
Co   p2    :product of lenght of primitive reciprocal lattice vectors
Cr Remarks:
Cr The optimal unit cell is chosen in such that it's volume is
Cr the same as the original unit cell, and that the product of
Cr all the lenght plat_1,plat_2,plat_3,qlat_1,qlat_2,qlat_3 is minmal.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision danrm2,plat(*),qlat(*),p1,p2,vol
C Local variables:
      double precision avol
C External calls:
Ckino 2011.05 del external,       external danrm2
C Intrinsic functions:
      intrinsic dabs,dsqrt

      avol=dabs(vol)
      p1=dsqrt(danrm2(plat(1))*danrm2(plat(4))*danrm2(plat(7)))/avol
      p2=dsqrt(danrm2(qlat(1))*danrm2(qlat(4))*danrm2(qlat(7)))*avol

      end
      double precision function danrm2(r)
C- Anisotropic 'norm'-square
C ----------------------------------------------------------------------
Ci Inputs:
Ci   r     :vector
Co Outputs:
Co   danrm2:anisotropic 'norm'-square
Cr Remarks:
Cr  A small anisotropy proportional to tiny is added to d2 to
Cr  distinguish between the three directions x,y,x and +x,-x ...
Cr  tiny should NOT be changed
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      double precision r(3)
C Local variables:
      double precision anis,d1mach,tiny,x,y,z,d2

      tiny=d1mach(3)*100d0
      x=r(1)
      y=r(2)
      z=r(3)
      d2=x*x+y*y+z*z
      anis=(30.34d0*x-4.23d0*dabs(x))*x+
     .(60.45d0*y-2.12d0*dabs(y))*y+
     .(90.56d0*z-1.01d0*dabs(z))*z
      danrm2=d2+anis*tiny

      end

