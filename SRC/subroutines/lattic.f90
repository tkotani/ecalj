module m_lattic
  real(8) , allocatable,protected ::  rv_a_odlv (:)
  real(8) , allocatable,protected ::  rv_a_oqlv (:)
  real(8),protected:: lat_plat(3,3),lat_qlat(3,3),lat_awald,lat_vol
  real(8), allocatable,protected :: rv_a_opos(:,:)
  !      real(8),protected:: lat_dist(3,3) !unused because no deformation of cell allowed currently.
  integer,protected:: lat_nkd,lat_nkq
contains

  subroutine Setopos() !called from lmfp to revise atomic position by reading rst file.
    use m_lmfinit,only: ctrl_nbas,v_ssite
    integer:: i
    do i=1,ctrl_nbas
       rv_a_opos(:,i)= v_ssite(i)%pos
    enddo
  end subroutine Setopos

  subroutine m_lattic_init() !slat,nbas)
    use m_lmfinit,only: ctrl_nbas,lat_alat, lat_as, lat_tol,lat_rpad,lat_nkdmx,lat_nkqmx,lat_gam,v_ssite, &
         lat_platin!, lat_dist0 ,lat_ldist!,rv_a_opos
    !     - Sets up the real and reciprocal space lattice vectors
    ! ----------------------------------------------------------------------
    ! o Inputs/Outputs
    ! o  slat  :struct for lattice information; see routine ulat
    ! o    Elts read: alat as tol nkdmx nkqmx gam plat platl platr ldist
    ! o               dist opos
    ! o    Stored:    vol plat0 plat qlat platl platr awald nkd nkq odlv
    ! o               oqlv
    ! o  sctrl :struct for program flow parameters; see routine uctrl
    ! o    Elts read: nbas
    ! o    Stored:    *
    ! o  ssite :struct for site-specific information; see routine usite
    ! o    Elts read: *
    ! o    Stored:    pos
    ! o  sarray:structure containing offsets to various arrays
    ! o    Elts read: npadl npadr
    ! o    Stored:    *
    !r Remarks
    !r    For historical reasons, lattice distortions may be EITHER
    !r    defined through gam (special-purpose volume conserving shear) OR
    !r    by one of the ldist modes:
    !r    ldist: 1: defgrd holds rot about spec'd angle
    !r           2, lattice deformed with a general linear transformation
    !r           3, lattice deformed by a shear.
    !u Updates
    !u   2 Mar 04 Pass rpad to lattc
    !u   5 Jun 01 (ATP) Now calls lattc after lattice transformation
    !u  19 Apr 00 Fixed rotations; new argument list
    ! ----------------------------------------------------------------------
    implicit none
    integer::  lmxst , nkd , nkdmx , nkq , nkqmx , nbas,i_data_size,i_spackv !ldist ,
    real(8),allocatable:: rv_a_tmp(:)
    !      integer nbaspp,npadl,npadr
    real(8):: alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol, &
         xx1,xx2,dotprd,pi,rpad, &
         plat0(3,3),plat(3,3),qlat(3,3) !platl(3,3),platr(3,3),dist(3,3)
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
    !      nbaspp = nbas
    ! ... Apply specified linear transformation of lattice and basis vectors
    !      ldist = lat_ldist
    !      dist= lat_dist0
    allocate(rv_a_opos(3,nbas))
    do i_spackv=1,nbas
       rv_a_opos(:,i_spackv)= v_ssite( i_spackv )%pos
    enddo
    if (abs(gt-1d0) > 1d-10) then
       call rdistn ( rv_a_opos , rv_a_opos , nbas , gx , gy , gz , gt )
       !        call rdistn ( rv_a_opos , rv_a_opos , nbaspp , gx , gy , gz , gt )
       !      elseif (ldist .ne. 0) then
       !        call lattdf ( ldist , dist , plat0 , nbaspp , rv_a_opos , 0 ,  0d0 )
       !      else
       !        dist=0d0
       !        dist(1,1) = 1
       !        dist(2,2) = 1
       !        dist(3,3) = 1
    endif

    allocate(rv_a_odlv(abs(3*nkdmx)))
    allocate(rv_a_oqlv(abs(3*nkqmx)))
    lmxst = 6
    call lattc ( awald0 , tol , rpad , alat , alat , plat0 , gx , &
         gy , gz , gt , plat , qlat , lmxst , vol , awald , rv_a_odlv &
         , nkd , rv_a_oqlv , nkq , nkdmx , nkqmx )
    lat_vol  =vol
    !      lat_plat0=plat0
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
    !      lat_dist=dist
    call tcx('m_lattic_init')
  end subroutine m_lattic_init

  subroutine lattdf(ldist,defgrd,plat,nbas,bas,ngen,gen)
    !- Rotates or deforms lattice
    ! ----------------------------------------------------------------
    !i Inputs
    !i   ldist: 0, no deformations.  For abs(ldist):
    !i          1: defgrd holds rot about spec'd angle
    !i          2, lattice deformed with a general linear transformation
    !i          3, lattice deformed by a shear.
    !i          SIGN ldist <0 => suppress rotation of plat
    !i   defgrd:transformation matrix, whose form is specified by ldist
    !i   nbas  :size of basis
    !i   ngen  :number of generators of the symmetry group
    !o Outputs
    ! o  plat  :primitive lattice vectors, in units of alat
    ! o        :On output, plat is transformed by defgrd
    ! o  bas   :basis vectors, in units of alat
    ! o        :On output, bas is transformed by defgrd
    ! o  gen   :On input, generators of the symmetry group
    ! o        :On output, generators are transformed by defgrd
    ! o        :Note: gen can only be transformed by a rotation
    ! o  defgrd:Output defgrd is actual transformation matrix.
    !u Updates
    !u   19 Mar 06 Blend Voigt strains into other types
    ! ----------------------------------------------------------------
    implicit none
    integer :: ldist,nbas,ngen
    double precision :: defgrd(3,3), plat(3,3), bas(3,1), gen(3,3,1)
    double precision :: work(3,3),rinv(3,3),det,gold(3,3)
    integer :: ipr,i,j,ib
    real(8) ,allocatable :: bwk_rv(:)
    if (ldist == 0) return
    call getpr(ipr)
    allocate(bwk_rv(3*nbas))
    call dcopy(9,defgrd,1,work,1)
    if (iabs(ldist) == 1) call makrot(work,defgrd)
    if (iabs(ldist) == 3) then
       det = defgrd(1,3)
       if ( det == 0 ) then
          if (allocated(bwk_rv)) deallocate(bwk_rv)
          return
       endif
       call shear(0,plat,bas,det,defgrd,defgrd)
    endif
    if (ipr >= 30) then
       print 333, ldist, ((defgrd(i,j), j=1,3), i=1,3)
333    format(/' LATTDF:  deformation matrix for mode',i3,':'/ &
            (3f12.7))
    endif
    ! ... Rotate or shear plat
    if (ldist >= 1 .AND. ldist <= 3) then
       call dcopy(9,plat,1,work,1)
       call dinv33(defgrd,0,rinv,det)
       call dmpy(defgrd,3,1,work,3,1,plat,3,1,3,3,3)
       if (ipr >= 30) then
          print 334, ((work(i,j), i=1,3), (plat(i,j), i=1,3),j=1,3)
334       format(10x,' Lattice vectors:',25x,'Transformed to:'/ &
               (3f12.7,2x,3f12.7))
       endif
    endif
    if (iabs(ldist) >= 1 .AND. iabs(ldist) <= 3) then
       if (ipr >= 30 .AND. nbas > 0) &
            print '(10x,''  Basis vectors:'',25x,''Transformed to:'')'
       do  10  ib = 1, nbas
          call dcopy(3,bas(1,ib),1,work,1)
          call dmpy(defgrd,3,1,work,3,1,bas(1,ib),3,1,3,1,3)
          if (ipr >= 30) print '(3f12.7,2x,3f12.7)', &
               (work(i,1), i=1,3), (bas(i,ib), i=1,3)
10     enddo
       if (ipr >= 30 .AND. ngen > 0) &
            print '(15x,''Group ops:'',26x,''Rotated to:'')'
       call dinv33(defgrd,0,rinv,det)
       do  20  ib = 1, ngen
          call dcopy(9,gen(1,1,ib),1,gold,1)
          call dmpy(defgrd,3,1,gen(1,1,ib),3,1,work,3,1,3,3,3)
          call dmpy(work,3,1,rinv,3,1,gen(1,1,ib),3,1,3,3,3)
          if (ipr >= 30) print '(/(3f12.7,2x,3f12.7))', &
               ((gold(j,i), i=1,3), (gen(j,i,ib), i=1,3),j=1,3)
20     enddo
    endif
    if (allocated(bwk_rv)) deallocate(bwk_rv)
  end subroutine lattdf

end module m_lattic

