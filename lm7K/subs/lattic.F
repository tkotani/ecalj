Cdelw1 var ['opos', '-1', 'undef', '0', 'undef', 'x', 'lattic']
Cdelw1 var ['odlv', '0', 'odlv', '3*nkdmx', 'real(8)', 'rv_w_', 'lattic']
Cdelw1 var ['oqlv', '0', 'oqlv', '3*nkqmx', 'real(8)', 'rv_w_', 'lattic']
Cdelw1 undel [['undef', 'wref'], 'opos', 'lattic']
Cdelw1 undel [['wref'], 'odlv', 'lattic']
Cdelw1 undel [['wref'], 'oqlv', 'lattic']
Cgetarg...info...           structure ['lattic', 'lat', 'slat']
Cgetarg...info...           structure ['lattic', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['lattic', 'array', 'sarray']
Cgetarg...info...           structure ['lattic', 'site', 'ssite']
Cgetarg...info...           use_to_add ['lattic', 'm_struc_def']
      subroutine lattic(slat,sctrl,ssite,sarray)
       
       use m_struc_def  !Cgetarg

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
Cgetarg       double precision slat(1),sctrl(1),ssite(1),sarray(1)
       type(s_lat)::slat
       type(s_ctrl)::sctrl
       type(s_site)::ssite(*)
       type(s_array)::sarray

C ... Local parameters
Cdelw1 do not change  [['wref'], 'odlv', 'lattic']
Cdelw1 do not change  [['wref'], 'oqlv', 'lattic']
Cdelw1 do not change  [['undef', 'wref'], 'opos', 'lattic']
      integer ldist,lmxst,nkd,nkdmx,nkq,nkqmx,nbas,odlv,oqlv,opos
      integer nbaspp,npadl,npadr
      double precision alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol,
     .  xx1,xx2,dotprd,pi,rpad,
     .  plat0(3,3),plat(3,3),qlat(3,3),dist(3,3) !platl(3,3),platr(3,3)
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
C ... Heap
      integer w(1)
      common /w/ w
C ... External calls
      external daxpy,dcopy,defps2,defrr,lattc,lattdf,pack1,pack5,rdistn,
     .         redfrr,spackv,upack,upack1,upack2

      print *,'lattic:'
C     call info(30,1,0,' Real and recip space lattices:',0,0)
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: awald0 , real(8) :: slat%as 1
Cgetarg...info...           real(8) :: tol , real(8) :: slat%tol 1
Cgetarg...info...           real(8) :: rpad , real(8) :: slat%rpad 1
Cgetarg       call upack('lat alat as tol rpad',slat,alat,awald0,tol,rpad,0)
       
       alat=slat%alat
       awald0=slat%as
       tol=slat%tol
       rpad=slat%rpad

Cgetarg...info...           integer :: nkdmx , integer(8) :: slat%nkdmx 1
Cgetarg...info...           integer :: nkqmx , integer(8) :: slat%nkqmx 1
Cgetarg...info...           real(8) :: gam(4) , real(8) :: slat%gam 4
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat0(3,3) , real(8) :: slat%plat 9
Cgetarg       call upack('lat nkdmx nkqmx gam alat plat',slat,nkdmx,nkqmx,gam,
Cgetarg      .  alat,plat0)
       
       nkdmx=slat%nkdmx
       nkqmx=slat%nkqmx
       i_copy_size=size(slat%gam)
       call dcopy(i_copy_size,slat%gam,1,gam,1) 
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat0,1) 

Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg       call upack1('ctrl nbas',sctrl,nbas)
       
       nbas=sctrl%nbas

Cgetarg...info...           integer :: npadl , integer(8) :: sarray%npadl 1
Cgetarg...info...           integer :: npadr , integer(8) :: sarray%npadr 1
Cgetarg       call upack2('array npadl npadr',sarray,npadl,npadr)
       
c       npadl=sarray%npadl
c       npadr=sarray%npadr

       nbaspp = nbas  !+ 2*(npadl+npadr)
ctakao assumes npadl=npapr=0, thus nbaspp=nbas

C$$$      if (nbaspp .gt. nbas) then
C$$$        call dcopy(9,plat0,1,plat,1)
C$$$Cgetarg...info...           real(8) :: platl(3,3) , real(8) :: slat%platl 9
C$$$Cgetarg...info...           real(8) :: platr(3,3) , real(8) :: slat%platr 9
C$$$Cgetarg         call upack2('lat platl platr',slat,platl,platr)
         
C$$$         i_copy_size=size(slat%platl) 
C$$$         call dcopy(i_copy_size,slat%platl,1,platl,1) 
C$$$         i_copy_size=size(slat%platr) 
C$$$         call dcopy(i_copy_size,slat%platr,1,platr,1) 

C$$$        call dcopy(6,plat0,1,platl,1)
C$$$        call dcopy(6,plat0,1,platr,1)
C$$$        call daxpy(3,2d0,platl(1,3),1,plat0(1,3),1)
C$$$        call daxpy(3,2d0,platr(1,3),1,plat0(1,3),1)

C$$$        pi = 4*datan(1d0)
C$$$        xx1 = 180/pi*dotprd(3,platl(1,3),1,plat(1,3),1)
C$$$        xx2 = 180/pi*dotprd(3,platr(1,3),1,plat(1,3),1)

C$$$        call info8(30,1,0,
C$$$     .    ' LATTIC:  Padding Plat(3) with end principal layers: '//
C$$$     .    '%N%3;11,6D Plat(3) as input'//
C$$$     .    '%N%3;11,6D PlatL:  angle (deg) with Plat(3) = %d'//
C$$$     .    '%N%3;11,6D PlatR:  angle (deg) with Plat(3) = %d'//
C$$$     .    '%N%3;11,6D Plat(3) including padding',
C$$$     .    plat(1,3),platl(1,3),xx1,platr(1,3),xx2,plat0(1,3),0,0)

C$$$      endif


c      print *,'xxxxxxx1'
C ... Apply specified linear transformation of lattice and basis vectors
Cgetarg...info...           integer :: ldist , integer(8) :: slat%ldist 1
Cgetarg...info...           real(8) :: dist(3,3) , real(8) :: slat%dist 9
Cgetarg...info...           integer :: opos , integer(8) :: slat%opos 1
Cgetarg       call upack('lat ldist dist opos',slat,ldist,dist,opos,0,0)
       
       ldist=slat%ldist
       i_copy_size=size(slat%dist) 
       call dcopy(i_copy_size,slat%dist,1,dist,1) 
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 146
Cdelw1 do not change opos because of ['undef', 'wref']
       opos=slat%opos

c      print *,'xxxxxxx2'
C     call prmx('pos from lat',w(opos),3,3,nbaspp)

      if (abs(gt-1) .gt. 1d-10) then
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 152
Cdelw1 do not change opos because of ['undef', 'wref']
        call rdistn(w(opos),w(opos),nbaspp,gx,gy,gz,gt)
      elseif (ldist .ne. 0) then
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 154
Cdelw1 do not change opos because of ['undef', 'wref']
        call lattdf(ldist,dist,plat0,nbaspp,w(opos),0,0d0)
c        if (nbaspp .gt. nbas) then
c          call lattdf(ldist,dist,platl,0,w,0,0d0)
c          call lattdf(ldist,dist,platr,0,w,0,0d0)
c        endif
      else
        call dpzero(dist,9)
        dist(1,1) = 1
        dist(2,2) = 1
        dist(3,3) = 1
      endif

Cdelw1 do not change odlv because of ['wref']
      call defrr(odlv,  3*nkdmx)
Cdelw1 do not change oqlv because of ['wref']
      call defrr(oqlv,  3*nkqmx)
      lmxst = 6
Cdelw1 do not change odlv because of ['wref']
Cdelw1 do not change oqlv because of ['wref']
      call lattc(awald0,tol,rpad,alat,alat,plat0,gx,gy,gz,gt,plat,qlat,
     .   lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx)

c      if (nbaspp .gt. nbas) then
c        if (cos(pi*xx1/180) .lt. 0) call rx('lattic:  '//
c     .    'angle betw/ PlatL and Plat(3) > 90 degrees ... fix PlatL')
c        if (cos(pi*xx2/180) .lt. 0) call rx('lattic:  '//
c     .    'angle betw/ PlatR and Plat(3) > 90 degrees ... fix PlatR')
c
c      endif
c      print *,'xxxxxxx3'

Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call pack1('lat vol',slat,vol)
       
       slat%vol=vol 

Cgetarg...info...           real(8) :: plat0(3,3) , real(8) :: slat%plat0 9
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: qlat(3,3) , real(8) :: slat%qlat 9
Cgetarg...info...           real(8) :: platl(3,3) , real(8) :: slat%platl 9
Cgetarg...info...           real(8) :: platr(3,3) , real(8) :: slat%platr 9
Cgetarg       call pack5('lat plat0 plat qlat platl platr',slat,plat0,plat,qlat,
Cgetarg      .  platl,platr)
       
       i_copy_size=size(slat%plat0) 
       call dcopy(i_copy_size,plat0,1,slat%plat0,1) 
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,plat,1,slat%plat,1) 
       i_copy_size=size(slat%qlat) 
       call dcopy(i_copy_size,qlat,1,slat%qlat,1) 
c       i_copy_size=size(slat%platl) 
c       call dcopy(i_copy_size,platl,1,slat%platl,1) 
c       i_copy_size=size(slat%platr) 
c       call dcopy(i_copy_size,platr,1,slat%platr,1) 

Cdelw1 do not change oqlv because of ['wref']
      call redfrr(oqlv,3*nkq)
Cdelw1 do not change odlv because of ['wref']
Cdelw1 do not change oqlv because of ['wref']
      call defps2(odlv,oqlv)
Cdelw1 do not change odlv because of ['wref']
      call redfrr(odlv,3*nkd)
Cgetarg...info...           real(8) :: awald , real(8) :: slat%awald 1
Cgetarg...info...           integer :: nkd , integer(8) :: slat%nkd 1
Cgetarg...info...           integer :: nkq , integer(8) :: slat%nkq 1
Cgetarg...info...           integer :: odlv , integer(8) :: slat%odlv 1
Cgetarg...info...           integer :: oqlv , integer(8) :: slat%oqlv 1
Cgetarg       call pack5('lat awald nkd nkq odlv oqlv',slat,awald,nkd,nkq,odlv,
Cgetarg      .  oqlv)
       
       slat%awald=awald 
       slat%nkd=nkd 
       slat%nkq=nkq 
Cdelw1 do not change odlv because of ['wref']
       slat%odlv=odlv 
Cdelw1 do not change oqlv because of ['wref']
       slat%oqlv=oqlv 

Cgetarg...info...           real(8) :: dist(3,3) , real(8) :: slat%dist 9
Cgetarg       call pack5('lat dist',slat,dist,0,0,0,0)
       
       i_copy_size=size(slat%dist) 
       call dcopy(i_copy_size,dist,1,slat%dist,1) 

Cgetarg       call spackv(11,'site pos',ssite,1,nbaspp,w(opos))
       i_copy_size=size(ssite(1)%pos)
       do i_spackv=1,nbaspp 
Cdelw1 warning(1) , probably  opos  is not defined yet at linenumber= 231
Cdelw1 do not change opos because of ['undef', 'wref']
       call spackv_array_copy_r8_r8('p',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
       enddo

C     call prmx('repack pos',w(opos),3,3,nbaspp)

Cgetarg       end
Cdelw1 w_varlist remains: odlv oqlv 
Cdelw1 w_varlistundel: [opos odlv oqlv]
Cdelw1 w_varlist (undel), remains: [odlv oqlv]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate oqlv because of [wref]
Cdelw1 not deallocate odlv because of [wref]
       end subroutine lattic 


