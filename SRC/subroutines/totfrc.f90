subroutine totfrc(leks,fes1,fes2,dfhf,fin,fout)! Add together and print contributions to force
  use m_lmfinit,only: nbas
  use m_mksym,only:  rv_a_osymgr , iv_a_oistab ,lat_nsgrp
  use m_struc_def
  use m_lgunit,only:stdo,stdl
  use m_ftox
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct containing site-specific information
  !i   leks: :<0 no forces calculated
  !i         : 0 Harris forces only
  !i         :>0 also print HKS forces
  !i         :>1 use HKS forces instead of HF forces
  !i   fes1  :contribution to HF forces from estat + xc potential   This is for input  density    !=3rd term in (B.5) in JPSJ.84.034705
  !i   fes2  :contribution to KS forces from estat + xc potential   This is for output density
  !i   dfhf  :2nd order corr. HF forces from ansatz density shift (dfrce) !=1st term  in (B.5)
  ! o Inputs/Outputs
  !i     fin  :On input, f is the contribution to force from eigval sum    !=2nd term in (B.5)
  ! o    fout :On output, f is the total force
  !c Remarks
  !u Updates
  !u   12 Apr 03 Prints out max correction to Harris force
  !u   30 May 00 Adapted from nfp totforce.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: leks,i_copy_size
  real(8):: fin(3,nbas) ,fout(3,nbas), f(3,nbas) , fes1(3,nbas) , fes2(3,nbas) , dfhf(3,nbas)
  integer:: ipr , ipl , ib , m , ng , ibmx , isw
  real(8) ,allocatable :: wk_rv(:)
  double precision :: fev1(3),fev2(3),fhar(3),fks(3),c,ddot,fmax,dfmax
  call tcn('totfrc')
  f=fin
  if (leks < 0) return
  c = 1000d0
  call getpr(ipr)
  ipl = 1
  fmax = -1
  dfmax = -1
  if (ipr >= 30) then
     if (ddot(3*nbas,dfhf,1,dfhf,1) /= 0) then
        write(stdo,'(/'' Forces, with eigenvalue correction'')')
     else
        write(stdo,'(/''Forces:'')')
     endif
     write(stdo,"('  ib',11x,'estatic',18x,'eigval',20x,'total')")
  endif
  do ib = 1, nbas
     fev1 = f(:,ib) + dfhf(:,ib)
     fev2 = f(:,ib)
     fhar = fev1(:) + fes1(:,ib)
     f(:,ib) = fhar(:)
     if(dsqrt(ddot(3,fhar,1,fhar,1)) > fmax) then
        ibmx = ib
        fmax = dsqrt(ddot(3,fhar,1,fhar,1))
     endif
     if (dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1)) > dfmax) then
        ibmx = ib
        dfmax = dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1))
     endif
     if(ipr >= 30) write(stdo,200) ib,(c*fes1(m,ib),m=1,3),(c*fev1(m),m=1,3),(c*fhar(m),m=1,3)
200  format(i4,3f8.2,1x,3f8.2,1x,3f8.2)
     if (leks == 0 .AND. ipl > 0 .AND. ipr >= 30) &
          write(stdl,"('fp ib',i4,'  fh ',3f8.2,2x,3f8.2)") ib,(c*fhar(m),m=1,3)
     if (leks >= 1) then
        fks  = fev2 + fes2(:,ib)
        if (leks > 1) f(:,ib) = fks
        if (ipr > 40) write(stdo,210) (c*fes2(m,ib),m=1,3),(c*fev2(m),m=1,3),(c*fks(m),m=1,3)
210     format('  KS',3f8.2,1x,3f8.2,1x,3f8.2)
        if (ipl > 0 .AND. ipr >= 30) &
             write (stdl,711) ib,(c*fhar(m),m=1,3),(c*fks(m),m=1,3)
711     format('fp ib',i4,'  fh ',3f8.2,'   fks ',3f8.2)
     endif
  enddo
  write(stdo,ftox)' Maximum Harris force = ',ftof(c*fmax),'mRy/au (site',ibmx,')'
  if(dfmax>0) write(stdo,ftox)' Max eval correction = ',ftof(c*dfmax)
  !     Symmetrize forces to machine precision
  ng=lat_nsgrp
  if (ng > 1) then
     if(ipr>=30) write(stdo,*)' Symmetrize forces ...'
     allocate(wk_rv(3*nbas))
     call symfor ( nbas , 1 , rv_a_osymgr , ng , iv_a_oistab , wk_rv , f )
     if (allocated(wk_rv)) deallocate(wk_rv)
  endif
  fout=f
  call tcx('totfrc')
end subroutine totfrc
subroutine symfor(nbas,mode,g,ng,istab,fwk,f)
  !- Symmetrize forces
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas:  number of atoms in basis
  !i   mode:  1 symmetrize by f(  ib ) = sum_g(ig) f(R(ib))
  !i          2 symmetrize by f(R(ib)) = sum_g(ig) f(  ib)
  !i   g,ng:  symmetry operations, and number
  !i   istab: site into which g,ag transforms site i
  !i   fwk:   work array of the same dimensions as f
  !o Outputs
  !o   f:  forces are symmetrized
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nbas,ng, istab(nbas,1),i,j,ib,ig,jb,mode
  double precision :: g(3,3,1),fwk(3,nbas),f(3,nbas)
  fwk=f/ng !!call dpcopy(f,fwk,1,3*nbas,1d0/ng)
  f=0d0
  if (mode == 1) then
     do  101  ig = 1, ng
        do  10  ib = 1, nbas
           jb = istab(ib,ig)
           do    i = 1, 3
              do    j = 1, 3
                 f(i,ib) = f(i,ib) + g(i,j,ig)*fwk(j,jb)
              enddo
           enddo
10      enddo
101  enddo
  else
     do  201  ig = 1, ng
        do  20  ib = 1, nbas
           jb = istab(ib,ig)
           do    i = 1, 3
              do    j = 1, 3
                 f(i,jb) = f(i,jb) + g(i,j,ig)*fwk(j,ib)
              enddo
           enddo
20      enddo
201  enddo
  endif
end subroutine symfor

