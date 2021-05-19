      module m_lattic
      real(8) , allocatable,protected ::  rv_a_odlv (:)
      real(8) , allocatable,protected ::  rv_a_oqlv (:)
      real(8),protected:: lat_plat(3,3)
      real(8), allocatable :: rv_a_opos(:,:)
      contains
      subroutine m_lattic_init() !slat,nbas)
      use m_lmfinit,only: ctrl_nbas, !,v_sctrl,slat=>v_slat,
     &     lat_alat, lat_as, lat_tol,lat_rpad,lat_nkdmx,lat_nkqmx,lat_gam,lat_dist,v_ssite,
     &     lat_platin,lat_qlat,lat_awald,lat_ldist,lat_nkd,lat_nkq,lat_vol!,rv_a_opos
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
      implicit none
      integer:: ldist , lmxst , nkd , nkdmx , nkq , nkqmx , nbas,i_data_size,i_spackv
      real(8),allocatable:: rv_a_tmp(:)
      integer nbaspp,npadl,npadr
      double precision alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol,
     .xx1,xx2,dotprd,pi,rpad,
     .plat0(3,3),plat(3,3),qlat(3,3),dist(3,3) !platl(3,3),platr(3,3)
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
      external daxpy,dcopy,lattc,lattdf,rdistn
      
      alat=lat_alat
      awald0=lat_as
      tol=lat_tol
      rpad=lat_rpad
      nkdmx=lat_nkdmx
      nkqmx=lat_nkqmx
      gam = lat_gam
      alat = lat_alat
      plat0=lat_platin
      nbas=ctrl_nbas
      nbaspp = nbas 
C ... Apply specified linear transformation of lattice and basis vectors
      ldist = lat_ldist
      dist(1:3,1)= lat_dist(1:3)
      dist(1:3,2)= lat_dist(4:6)
      dist(1:3,3)= lat_dist(7:9)
      allocate(rv_a_opos(3,nbas))
      do i_spackv=1,nbas
         rv_a_opos(:,i_spackv)= v_ssite( i_spackv )%pos
      enddo
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
      lat_vol  =vol
c      lat_plat0=plat0
      lat_plat =plat
      lat_qlat =qlat
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
      lat_awald=awald
      lat_nkd=nkd
      lat_nkq=nkq
      lat_dist(1:3)=dist(1:3,1)
      lat_dist(4:6)=dist(1:3,2)
      lat_dist(7:9)=dist(1:3,3)
      end subroutine m_lattic_init
      end module m_lattic
 
