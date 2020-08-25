      subroutine lattic() !slat,nbas)
      use m_struc_def  !Cgetarg
      use m_lmfinit,only: slat=>v_slat,v_sctrl
      use m_work,only: rv_a_opos,rv_a_oqlv,rv_a_odlv
C     - Sets up the real and reciprocal space lattice vectors
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
c      type(s_lat)::slat
c      type(s_ctrl)::sctrl
c      type(s_site)::ssite(*)
c      type(s_array)::sarray
C ... Local parameters
      integer:: ldist , lmxst , nkd , nkdmx , nkq , nkqmx , nbas
      real(8),allocatable:: rv_a_tmp(:)
      integer nbaspp,npadl,npadr
      double precision alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol,
     .xx1,xx2,dotprd,pi,rpad,
     .plat0(3,3),plat(3,3),qlat(3,3),dist(3,3) !platl(3,3),platr(3,3)
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
      external daxpy,dcopy,lattc,lattdf,rdistn

      alat=slat%alat
      awald0=slat%as
      tol=slat%tol
      rpad=slat%rpad
      nkdmx=slat%nkdmx
      nkqmx=slat%nkqmx
      gam(:)=slat%gam(:)
      alat=slat%alat
c      i_copy_size=size(slat%plat)
c      call dcopy(i_copy_size,slat%plat,1,plat0,1)
      plat0=slat%plat
      nbas=v_sctrl%nbas
      nbaspp = nbas 

C ... Apply specified linear transformation of lattice and basis vectors
      ldist = slat%ldist
      dist(1:3,1)= slat%dist(1:3)
      dist(1:3,2)= slat%dist(4:6)
      dist(1:3,3)= slat%dist(7:9)

      if (abs(gt-1) .gt. 1d-10) then
        call rdistn ( rv_a_opos , rv_a_opos , nbaspp , gx , gy , gz , gt )
      elseif (ldist .ne. 0) then
        call lattdf ( ldist , dist , plat0 , nbaspp , rv_a_opos , 0 ,  0d0 )
      else
        dist=0d0 
        dist(1,1) = 1
        dist(2,2) = 1
        dist(3,3) = 1
      endif

      allocate(rv_a_odlv(abs(3*nkdmx)))
      allocate(rv_a_oqlv(abs(3*nkqmx)))
      lmxst = 6
      call lattc ( awald0 , tol , rpad , alat , alat , plat0 , gx ,
     .  gy , gz , gt , plat , qlat , lmxst , vol , awald , rv_a_odlv
     .  , nkd , rv_a_oqlv , nkq , nkdmx , nkqmx )

      slat%vol  =vol
      slat%plat0=plat0
      slat%plat =plat
      slat%qlat =qlat

!! reduce size. necessary?
      i_data_size=size(rv_a_oqlv)
      allocate(rv_a_tmp(i_data_size))
      rv_a_tmp=rv_a_oqlv
      deallocate(rv_a_oqlv)
      i_data_size=min(i_data_size,3*nkq)
      allocate(rv_a_oqlv(3*nkq))
      rv_a_oqlv(:i_data_size)=rv_a_tmp(:i_data_size)
      deallocate(rv_a_tmp)

!! reduce size. necessary?
      i_data_size=size(rv_a_odlv)
      allocate(rv_a_tmp(i_data_size))
      rv_a_tmp=rv_a_odlv
      deallocate(rv_a_odlv)
      i_data_size=min(i_data_size,3*nkd)
      allocate(rv_a_odlv(3*nkd))
      rv_a_odlv(:i_data_size)=rv_a_tmp(:i_data_size)
      deallocate(rv_a_tmp)

      slat%awald=awald
      slat%nkd=nkd
      slat%nkq=nkq
      slat%dist(1:3)=dist(1:3,1)
      slat%dist(4:6)=dist(1:3,2)
      slat%dist(7:9)=dist(1:3,3)
      end subroutine lattic


