Cchp1 cmdvar= real(8),rv_p_,odlv
Cchp1 cmdvar= real(8),rv_p_,oqlv
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
      type(s_lat)::slat
      type(s_ctrl)::sctrl
      type(s_site)::ssite(*)
      type(s_array)::sarray

C ... Local parameters
Cchp1       integer ldist,lmxst,nkd,nkdmx,nkq,nkqmx,nbas,odlv,oqlv,opos
Cchp1        integer:: ldist , lmxst , nkd , nkdmx , nkq , nkqmx , nbas , 
Cchp1      . odlv , opos 
        integer:: ldist , lmxst , nkd , nkdmx , nkq , nkqmx , nbas , 
     .  opos 
        real(8),pointer :: rv_p_odlv(:)

       real(8),pointer :: rv_p_oqlv(:)
       real(8),allocatable:: rv_a_tmp(:)

      integer nbaspp,npadl,npadr
      double precision alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol,
     .xx1,xx2,dotprd,pi,rpad,
     .plat0(3,3),plat(3,3),qlat(3,3),dist(3,3) !platl(3,3),platr(3,3)
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
C ... Heap
      integer w(1)
      common /w/ w
C ... External calls
      external daxpy,dcopy,defps2,defrr,lattc,lattdf,pack1,pack5,rdistn,
     .redfrr,spackv,upack,upack1,upack2

      print *,'lattic:'
C     call info(30,1,0,' Real and recip space lattices:',0,0)

      alat=slat%alat
      awald0=slat%as
      tol=slat%tol
      rpad=slat%rpad


      nkdmx=slat%nkdmx
      nkqmx=slat%nkqmx
      i_copy_size=size(slat%gam)
      call dcopy(i_copy_size,slat%gam,1,gam,1)
      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat0,1)


      nbas=sctrl%nbas


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

      ldist=slat%ldist
      i_copy_size=size(slat%dist)
      call dcopy(i_copy_size,slat%dist,1,dist,1)
      opos=slat%opos

c      print *,'xxxxxxx2'
C     call prmx('pos from lat',w(opos),3,3,nbaspp)

      if (abs(gt-1) .gt. 1d-10) then
        call rdistn(w(opos),w(opos),nbaspp,gx,gy,gz,gt)
      elseif (ldist .ne. 0) then
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

Cchp1       call defrr(odlv,  3*nkdmx)
       allocate(rv_p_odlv(abs(3*nkdmx)))
       if (3*nkdmx<0) rv_p_odlv(:)=0.0d0

Cchp1       call defrr(oqlv,  3*nkqmx)
       allocate(rv_p_oqlv(abs(3*nkqmx)))
       if (3*nkqmx<0) rv_p_oqlv(:)=0.0d0

      lmxst = 6
Cchp1       call lattc(awald0,tol,rpad,alat,alat,plat0,gx,gy,gz,gt,plat,qlat,
Cchp1      .lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx)
Cchp1        call lattc ( awald0 , tol , rpad , alat , alat , plat0 , gx , 
Cchp1      . gy , gz , gt , plat , qlat , lmxst , vol , awald , w ( odlv ) 
Cchp1      . , nkd , rv_p_oqlv , nkq , nkdmx , nkqmx ) 
        call lattc ( awald0 , tol , rpad , alat , alat , plat0 , gx , 
     .  gy , gz , gt , plat , qlat , lmxst , vol , awald , rv_p_odlv 
     .  , nkd , rv_p_oqlv , nkq , nkdmx , nkqmx ) 



c      if (nbaspp .gt. nbas) then
c        if (cos(pi*xx1/180) .lt. 0) call rx('lattic:  '//
c     .    'angle betw/ PlatL and Plat(3) > 90 degrees ... fix PlatL')
c        if (cos(pi*xx2/180) .lt. 0) call rx('lattic:  '//
c     .    'angle betw/ PlatR and Plat(3) > 90 degrees ... fix PlatR')
c
c      endif
c      print *,'xxxxxxx3'


      slat%vol=vol


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

Cchp1       call redfrr(oqlv,3*nkq)
       i_data_size=size(rv_p_oqlv); allocate(rv_a_tmp(i_data_size))
       rv_a_tmp=rv_p_oqlv; deallocate(rv_p_oqlv)
       i_data_size=min(i_data_size,3*nkq); allocate(rv_p_oqlv(3*nkq))
       rv_p_oqlv(:i_data_size)=rv_a_tmp(:i_data_size); deallocate(rv_a_tmp)

cki      call defps2(odlv,oqlv)
Cchp1       call redfrr(odlv,3*nkd)
       i_data_size=size(rv_p_odlv); allocate(rv_a_tmp(i_data_size))
       rv_a_tmp=rv_p_odlv; deallocate(rv_p_odlv)
       i_data_size=min(i_data_size,3*nkd); allocate(rv_p_odlv(3*nkd))
       rv_p_odlv(:i_data_size)=rv_a_tmp(:i_data_size); deallocate(rv_a_tmp)


      slat%awald=awald
      slat%nkd=nkd
      slat%nkq=nkq
Cchp1       slat%odlv=odlv
       slat%rv_p_odlv => rv_p_odlv 

Cchp1       slat%oqlv=oqlv
       slat%rv_p_oqlv => rv_p_oqlv 



      i_copy_size=size(slat%dist)
      call dcopy(i_copy_size,dist,1,slat%dist,1)

      i_copy_size=size(ssite(1)%pos)
      do i_spackv=1,nbaspp
        call spackv_array_copy_r8_r8('p',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
      enddo

C     call prmx('repack pos',w(opos),3,3,nbaspp)

      end subroutine lattic


