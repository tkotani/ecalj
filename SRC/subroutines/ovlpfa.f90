subroutine ovlpfa(ssite,nbas,nxi,nxi0,exi,hfc,rsmfa,ng,ngmx,  gv,cv)
  use m_lmfinit,only:lat_alat,nsp
  use m_lattic,only: lat_vol
  use m_ftox
  use m_lgunit,only:stdo
  use m_struc_def
  !      use m_globalvariables
  !- Set up Fourier coeffs to overlap the smooth part of FA densities.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: *
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: alat vol
  !i     Stored:    *
  !i     Passed to: *
  !i   nbas  :size of basis
  !i   nxi   :number of Hankels
  !i   nxi0  :leading dimension of hfc
  !i   exi   :smoothed Hankel energies; see Remarks
  !i   hfc   :coefficients to smoothed Hankels
  !i   rsmfa :Hankel smoothing radius
  !i   ng    :number of G-vectors
  !i   ngmx  :leading dimension of gv
  !i   gv    :list of reciprocal lattice vectors G (glist.f)
  !o Outputs
  !o   cv    :Fourier coefficients
  !r Remarks
  !u Updates
  !u   12 May 07 parallelized (MPI)
  !u   01 Jul 05 Zero-radius empty spheres skip as having no local part
  !u   13 Jun 00 spin polarized
  !u   24 Apr 00 Adapted from nfp ovlpfa.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nbas,nxi(1),nxi0,ng,ngmx
  real(8):: gv(ngmx,3) , rsmfa(1) , exi(nxi0,1) , hfc(nxi0,2,1)
  !      type(s_lat)::slat
  type(s_site)::ssite(*)

  double complex cv(ng,*)
  ! ... Local parameters
  integer :: ipr,iprint,ib,is,nx,i,ixi,ig,isp
  double precision :: v(3),pi,y0,alat,vol,tpiba,sum(2),px,py,pz,pos(3), &
       sam(2),e,cof,rsm,gam,v2,aa,scalp
  double complex phase
  equivalence (px,pos(1)),(py,pos(2)),(pz,pos(3))
! #if MPI | MPIK
!   integer, dimension(:),allocatable :: kpproc
!   integer :: ierr
!   !     integer ifi,fopna
!   integer :: procid, master, mpipid, numprocs
!   logical :: mlog,cmdopt
!   character strn*120
! #endif
  integer:: ibini,ibend

  call tcn('ovlpfa')
  ipr  = iprint()
  !      stdo = lgunit(1)
  ! angenglob      nsp  = nglob('nsp')
  !      nsp  = globalvariables%nsp
  pi   = 4d0*datan(1d0)
  y0   = 1d0/dsqrt(4*pi)

  alat=lat_alat
  vol=lat_vol

  tpiba = 2*pi/alat
  call dpzero(cv,2*ng*nsp)
! #if MPI | MPIK
!   procid = mpipid(1)
!   numprocs = mpipid(0)
!   master = 0
!   mlog = cmdopt('--mlog',6,0,strn)
! #endif

  ! --- Loop over sites ---
  sum(1) = 0d0
  sum(2) = 0d0
  call info0(31,1,0,' ovlpfa: overlap smooth part of FA densities')
! #if MPI | MPIK
!   allocate (kpproc(0:numprocs), stat=ierr)
!   call pshpr(ipr-10)
!   call dstrbp(nbas,numprocs,1,kpproc(0))
!   call poppr
!   ipr = 0
!   !      do ib = kpproc(procid), kpproc(procid+1)-1
!   ibini = kpproc(procid)
!   ibend = kpproc(procid+1)-1
! #else
  !      do  ib = 1, nbas
  ibini=1
  ibend=nbas
!#endif

  do ib=ibini,ibend
     is=ssite(ib)%spec
     !        i_copy_size=size(ssite(ib)%pos)
     !        call dcopy(i_copy_size,ssite(ib)%pos,1,pos,1)
     pos=ssite(ib)%pos
     nx = nxi(is)
     !   ... Loop over Hankels at this site
     sam(1) = 0
     sam(2) = 0
     do  isp = 1, nsp
        do  ixi = 1, nx
           e = exi(ixi,is)
           cof = hfc(ixi,isp,is)
           rsm = rsmfa(is)
           gam = 0.25d0*rsm**2
           sam(isp) = sam(isp) - cof*y0*4d0*pi*dexp(gam*e)/e
           !       ... Loop over reciprocal lattice vectors
           do  ig = 1, ng
              v(1) = gv(ig,1)*tpiba
              v(2) = gv(ig,2)*tpiba
              v(3) = gv(ig,3)*tpiba
              v2 = v(1)**2+v(2)**2+v(3)**2
              aa = -4d0*pi*dexp(gam*(e-v2))/(e-v2)
              scalp = -alat*(px*v(1)+py*v(2)+pz*v(3))
              phase = dcmplx(dcos(scalp),dsin(scalp))
              cv(ig,isp) = cv(ig,isp) + cof*aa*phase*y0/vol
           enddo
        enddo
        sum(isp) = sum(isp) + sam(isp)
     enddo
     if (ipr > 30 .AND. nx > 0) then
        write(stdo,ftox)' site',ib,'spec',is,'pos',ftof(pos,4), &
             'Qsmooth',sam(1)+sam(2),'mom', sam(1)-sam(2)
        if (ipr >= 40) then
           write(stdo,700) 'energy:',(exi(i,is),i=1,nx)
           write(stdo,700) 'coeff:',(hfc(i,1,is),i=1,nx)
           if (nsp == 2) write(stdo,700) 'spin2:',(hfc(i,2,is),i=1,nx)
           write(stdo,'(1x)')
        endif
700     format(2x,a7,16f11.3)
     endif
  enddo

  ! ... Combine cv from separate threads
! #if MPI | MPIK
!   deallocate(kpproc, stat=ierr)
!   call mpibc2_real(cv,ng*nsp*2,'ovlpfa_cv')
!   call mpibc2_real(sum,2,'ovlpfa_sum')
!   !     Debugging
!   !     ifi = fopna('out',-1,0)
!   !     call ywrm(0,'cv',3,ifi,'(9f20.10)',cv,1,ng,ng,1)
!   !     call rx0('done')
! !#endif

  ipr  = iprint()
  call info5(31,0,0,' total smooth Q = %,6;6d'// &
       '%?#n==2#  moment = %,5;5d#%j#%?#n>40#  FT (0,0,0) = %,6;6d', &
       sum(1)+sum(2),nsp,sum(1)-sum(2),ipr, &
       (cv(1,1)+cv(1,nsp))*vol/(3-nsp))

  call tcx('ovlpfa')
end subroutine ovlpfa

