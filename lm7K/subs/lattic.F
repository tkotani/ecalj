      subroutine lattic(slat,sctrl,ssite,sarray)
C- Sets up the real and reciprocal space lattice vectors
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Cio  slat  :struct for lattice information; see routine ulat
Cio    Elts read: alat as tol nkdmx nkqmx gam plat platl platr ldist
Cio               dist opos
Cio    Stored:    vol plat0 plat qlat platl platr awald nkd nkq odlv
Cio               oqlv
Cio  sctrl :struct for program flow parameters; see routine uctrl
Cio    Elts read: nbas
Cio    Stored:    *
Cio  ssite :struct for site-specific information; see routine usite
Cio    Elts read: *
Cio    Stored:    pos
Cio  sarray:structure containing offsets to various arrays
Cio    Elts read: npadl npadr
Cio    Stored:    *
Cr Remarks
Cr    For historical reasons, lattice distortions may be EITHER
Cr    defined through gam (special-purpose volume conserving shear) OR
Cr    by one of the ldist modes:
Cr    ldist: 1: defgrd holds rot about spec'd angle
Cr           2, lattice deformed with a general linear transformation
Cr           3, lattice deformed by a shear.
Cu Updates
Cu   2 Mar 04 Pass rpad to lattc
Cu   5 Jun 01 (ATP) Now calls lattc after lattice transformation
Cu  19 Apr 00 Fixed rotations; new argument list
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision slat(1),sctrl(1),ssite(1),sarray(1)
C ... Local parameters
      integer ldist,lmxst,nkd,nkdmx,nkq,nkqmx,nbas,odlv,oqlv,opos
      integer nbaspp,npadl,npadr
      double precision alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol,
     .  xx1,xx2,dotprd,pi,rpad,
     .  platl(3,3),platr(3,3),plat0(3,3),plat(3,3),qlat(3,3),dist(3,3)
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
C ... Heap
      integer w(1)
      common /w/ w
C ... External calls
      external daxpy,dcopy,defps2,defrr,lattc,lattdf,pack1,pack5,rdistn,
     .         redfrr,spackv,upack,upack1,upack2

C     call info(30,1,0,' Real and recip space lattices:',0,0)
      call upack('lat alat as tol rpad',slat,alat,awald0,tol,rpad,0)
      call upack('lat nkdmx nkqmx gam alat plat',slat,nkdmx,nkqmx,gam,
     .  alat,plat0)
      call upack1('ctrl nbas',sctrl,nbas)
      call upack2('array npadl npadr',sarray,npadl,npadr)
      nbaspp = nbas + 2*(npadl+npadr)
      if (nbaspp .gt. nbas) then
        call dcopy(9,plat0,1,plat,1)
        call upack2('lat platl platr',slat,platl,platr)
        call dcopy(6,plat0,1,platl,1)
        call dcopy(6,plat0,1,platr,1)
        call daxpy(3,2d0,platl(1,3),1,plat0(1,3),1)
        call daxpy(3,2d0,platr(1,3),1,plat0(1,3),1)

        pi = 4*datan(1d0)
        xx1 = 180/pi*dotprd(3,platl(1,3),1,plat(1,3),1)
        xx2 = 180/pi*dotprd(3,platr(1,3),1,plat(1,3),1)

        call info8(30,1,0,
     .    ' LATTIC:  Padding Plat(3) with end principal layers: '//
     .    '%N%3;11,6D Plat(3) as input'//
     .    '%N%3;11,6D PlatL:  angle (deg) with Plat(3) = %d'//
     .    '%N%3;11,6D PlatR:  angle (deg) with Plat(3) = %d'//
     .    '%N%3;11,6D Plat(3) including padding',
     .    plat(1,3),platl(1,3),xx1,platr(1,3),xx2,plat0(1,3),0,0)

      endif

C ... Apply specified linear transformation of lattice and basis vectors
      call upack('lat ldist dist opos',slat,ldist,dist,opos,0,0)
C     call prmx('pos from lat',w(opos),3,3,nbaspp)

      if (abs(gt-1) .gt. 1d-10) then
        call rdistn(w(opos),w(opos),nbaspp,gx,gy,gz,gt)
      elseif (ldist .ne. 0) then
        call lattdf(ldist,dist,plat0,nbaspp,w(opos),0,0d0)
        if (nbaspp .gt. nbas) then
          call lattdf(ldist,dist,platl,0,w,0,0d0)
          call lattdf(ldist,dist,platr,0,w,0,0d0)
        endif
      else
        call dpzero(dist,9)
        dist(1,1) = 1
        dist(2,2) = 1
        dist(3,3) = 1
      endif

      call defrr(odlv,  3*nkdmx)
      call defrr(oqlv,  3*nkqmx)
      lmxst = 6
      call lattc(awald0,tol,rpad,alat,alat,plat0,gx,gy,gz,gt,plat,qlat,
     .   lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx)

      if (nbaspp .gt. nbas) then
        if (cos(pi*xx1/180) .lt. 0) call rx('lattic:  '//
     .    'angle betw/ PlatL and Plat(3) > 90 degrees ... fix PlatL')
        if (cos(pi*xx2/180) .lt. 0) call rx('lattic:  '//
     .    'angle betw/ PlatR and Plat(3) > 90 degrees ... fix PlatR')

      endif

      call pack1('lat vol',slat,vol)
      call pack5('lat plat0 plat qlat platl platr',slat,plat0,plat,qlat,
     .  platl,platr)
      call redfrr(oqlv,3*nkq)
      call defps2(odlv,oqlv)
      call redfrr(odlv,3*nkd)
      call pack5('lat awald nkd nkq odlv oqlv',slat,awald,nkd,nkq,odlv,
     .  oqlv)
      call pack5('lat dist',slat,dist,0,0,0,0)
      call spackv(11,'site pos',ssite,1,nbaspp,w(opos))
C     call prmx('repack pos',w(opos),3,3,nbaspp)

      end

