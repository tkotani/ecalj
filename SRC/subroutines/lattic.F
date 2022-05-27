      module m_lattic
      real(8) , allocatable,protected ::  rv_a_odlv (:)
      real(8) , allocatable,protected ::  rv_a_oqlv (:)
      real(8),protected:: lat_plat(3,3),lat_qlat(3,3),lat_awald,lat_vol
      real(8), allocatable,protected :: rv_a_opos(:,:)
c      real(8),protected:: lat_dist(3,3) !unused because no deformation of cell allowed currently.
      integer,protected:: lat_nkd,lat_nkq
      contains
      
      subroutine Setopos() !called from lmfp to revise atomic position by reading rst file.
      use m_lmfinit,only: ctrl_nbas,v_ssite
      integer:: i
      do i=1,ctrl_nbas
         rv_a_opos(:,i)= v_ssite(i)%pos
      enddo
      end
      
      subroutine m_lattic_init() !slat,nbas)
       use m_lmfinit,only: ctrl_nbas, !,v_sctrl,slat=>v_slat,
     &     lat_alat, lat_as, lat_tol,lat_rpad,lat_nkdmx,lat_nkqmx,lat_gam,v_ssite,
     &     lat_platin!, lat_dist0 ,lat_ldist!,rv_a_opos
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
      integer::  lmxst , nkd , nkdmx , nkq , nkqmx , nbas,i_data_size,i_spackv !ldist ,
      real(8),allocatable:: rv_a_tmp(:)
c      integer nbaspp,npadl,npadr
      real(8):: alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol,
     .xx1,xx2,dotprd,pi,rpad,
     .plat0(3,3),plat(3,3),qlat(3,3) !platl(3,3),platr(3,3),dist(3,3)
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
      call tcn('m_lattic_init')
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
c      nbaspp = nbas 
C ... Apply specified linear transformation of lattice and basis vectors
c      ldist = lat_ldist
c      dist= lat_dist0
      allocate(rv_a_opos(3,nbas))
      do i_spackv=1,nbas
         rv_a_opos(:,i_spackv)= v_ssite( i_spackv )%pos
      enddo
      if (abs(gt-1d0) .gt. 1d-10) then
        call rdistn ( rv_a_opos , rv_a_opos , nbas , gx , gy , gz , gt )
c        call rdistn ( rv_a_opos , rv_a_opos , nbaspp , gx , gy , gz , gt )
c      elseif (ldist .ne. 0) then
c        call lattdf ( ldist , dist , plat0 , nbaspp , rv_a_opos , 0 ,  0d0 )
c      else
c        dist=0d0 
c        dist(1,1) = 1
c        dist(2,2) = 1
c        dist(3,3) = 1
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
c      lat_dist=dist
      call tcx('m_lattic_init')
      end subroutine m_lattic_init
      
      subroutine lattdf(ldist,defgrd,plat,nbas,bas,ngen,gen)
C- Rotates or deforms lattice
C ----------------------------------------------------------------
Ci Inputs
Ci   ldist: 0, no deformations.  For abs(ldist):
Ci          1: defgrd holds rot about spec'd angle
Ci          2, lattice deformed with a general linear transformation
Ci          3, lattice deformed by a shear.
Ci          SIGN ldist <0 => suppress rotation of plat
Ci   defgrd:transformation matrix, whose form is specified by ldist
Ci   nbas  :size of basis
Ci   ngen  :number of generators of the symmetry group
Co Outputs
Cio  plat  :primitive lattice vectors, in units of alat
Cio        :On output, plat is transformed by defgrd
Cio  bas   :basis vectors, in units of alat
Cio        :On output, bas is transformed by defgrd
Cio  gen   :On input, generators of the symmetry group
Cio        :On output, generators are transformed by defgrd
Cio        :Note: gen can only be transformed by a rotation
Cio  defgrd:Output defgrd is actual transformation matrix.
Cu Updates
Cu   19 Mar 06 Blend Voigt strains into other types
C ----------------------------------------------------------------
      implicit none
      integer ldist,nbas,ngen
      double precision defgrd(3,3), plat(3,3), bas(3,1), gen(3,3,1)
      double precision work(3,3),rinv(3,3),det,gold(3,3)
      integer ipr,i,j,ib
      real(8) ,allocatable :: bwk_rv(:)
      if (ldist .eq. 0) return
      call getpr(ipr)
      allocate(bwk_rv(3*nbas))
      call dcopy(9,defgrd,1,work,1)
      if (iabs(ldist) .eq. 1) call makrot(work,defgrd)
      if (iabs(ldist) .eq. 3) then
        det = defgrd(1,3)
        if ( det .eq. 0 ) then
          if (allocated(bwk_rv)) deallocate(bwk_rv)
          return
        endif
        call shear(0,plat,bas,det,defgrd,defgrd)
      endif
      if (ipr .ge. 30) then
        print 333, ldist, ((defgrd(i,j), j=1,3), i=1,3)
  333   format(/' LATTDF:  deformation matrix for mode',i3,':'/
     .  (3f12.7))
      endif
C ... Rotate or shear plat
      if (ldist .ge. 1 .and. ldist .le. 3) then
        call dcopy(9,plat,1,work,1)
        call dinv33(defgrd,0,rinv,det)
        call dmpy(defgrd,3,1,work,3,1,plat,3,1,3,3,3)
        if (ipr .ge. 30) then
          print 334, ((work(i,j), i=1,3), (plat(i,j), i=1,3),j=1,3)
  334     format(10x,' Lattice vectors:',25x,'Transformed to:'/
     .    (3f12.7,2x,3f12.7))
        endif
      endif
      if (iabs(ldist) .ge. 1 .and. iabs(ldist) .le. 3) then
         if (ipr .ge. 30 .and. nbas .gt. 0)
     .        print '(10x,''  Basis vectors:'',25x,''Transformed to:'')'
         do  10  ib = 1, nbas
            call dcopy(3,bas(1,ib),1,work,1)
            call dmpy(defgrd,3,1,work,3,1,bas(1,ib),3,1,3,1,3)
            if (ipr .ge. 30) print '(3f12.7,2x,3f12.7)',
     .           (work(i,1), i=1,3), (bas(i,ib), i=1,3)
 10      continue
         if (ipr .ge. 30 .and. ngen .gt. 0)
     .        print '(15x,''Group ops:'',26x,''Rotated to:'')'
         call dinv33(defgrd,0,rinv,det)
         do  20  ib = 1, ngen
            call dcopy(9,gen(1,1,ib),1,gold,1)
            call dmpy(defgrd,3,1,gen(1,1,ib),3,1,work,3,1,3,3,3)
            call dmpy(work,3,1,rinv,3,1,gen(1,1,ib),3,1,3,3,3)
            if (ipr .ge. 30) print '(/(3f12.7,2x,3f12.7))',
     .           ((gold(j,i), i=1,3), (gen(j,i,ib), i=1,3),j=1,3)
 20      continue
      endif
      if (allocated(bwk_rv)) deallocate(bwk_rv)
      end
      
      end module m_lattic
 
