      subroutine totfrc(leks,fes1,fes2,dfhf,fin,fout)
      use m_lmfinit,only: nbas,ssite=>v_ssite
      use m_mksym,only:  rv_a_osymgr , iv_a_oistab ,lat_nsgrp
      use m_struc_def           
      use m_lgunit,only:stdo,stdl
C- Add together and print contributions to force
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci   leks: :<0 no forces calculated
Ci         : 0 Harris forces only
Ci         :>0 also print HKS forces   
Ci         :>1 use HKS forces instead of HF forces
Ci   fes1  :contribution to HF forces from estat + xc potential   This is for input  density    !=3rd term in (B.5) in JPSJ.84.034705
Ci   fes2  :contribution to KS forces from estat + xc potential   This is for output density
Ci   dfhf  :2nd order corr. HF forces from ansatz density shift (dfrce) !=1st term  in (B.5)
Cio Inputs/Outputs
Ci     fin  :On input, f is the contribution to force from eigval sum    !=2nd term in (B.5)
Cio    fout :On output, f is the total force
cc Remarks
Cu Updates
Cu   12 Apr 03 Prints out max correction to Harris force
Cu   30 May 00 Adapted from nfp totforce.f
C ----------------------------------------------------------------------
      implicit none
      integer leks,i_copy_size
      real(8):: fin(3,nbas) ,fout(3,nbas), f(3,nbas) , fes1(3,nbas) , fes2(3,nbas) , dfhf(3,nbas)
      integer:: ipr , ipl , ib , m , ng , ibmx , isw
      real(8) ,allocatable :: wk_rv(:)
      double precision fev1(3),fev2(3),fhar(3),fks(3),c,ddot,fmax,dfmax
      call tcn('totfrc')
      f=fin
c      stdo = lgunit(1)
c      stdl = lgunit(2)
      if (leks .lt. 0) return
      c = 1000d0
      call getpr(ipr)
      ipl = 1
      fmax = -1
      dfmax = -1

      if (ipr .ge. 30) then
        if (ddot(3*nbas,dfhf,1,dfhf,1) .ne. 0) then
          write(stdo,'(/'' Forces, with eigenvalue correction'')')
        else
          write(stdo,'(/''Forces:'')')
        endif
        write(stdo,201)
  201   format('  ib',11x,'estatic',18x,'eigval',20x,'total')
      endif
      
      do ib = 1, nbas
        fev1 = f(:,ib) + dfhf(:,ib)
        fev2 = f(:,ib)
        fhar = fev1(:) + fes1(:,ib)
        f(:,ib) = fhar(:)
        if(dsqrt(ddot(3,fhar,1,fhar,1)) .gt. fmax) then
          ibmx = ib
          fmax = dsqrt(ddot(3,fhar,1,fhar,1))
        endif
        if (dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1)) .gt. dfmax) then
          ibmx = ib
          dfmax = dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1))
        endif
        if(ipr .ge. 30) write(stdo,200) ib,(c*fes1(m,ib),m=1,3),(c*fev1(m),m=1,3),(c*fhar(m),m=1,3)
  200   format(i4,3f8.2,1x,3f8.2,1x,3f8.2)
        if (leks. eq. 0 .and. ipl .gt. 0 .and. ipr .ge. 30)
     .  write(stdl,"('fp ib',i4,'  fh ',3f8.2,2x,3f8.2)") ib,(c*fhar(m),m=1,3)
        if (leks >= 1) then
          fks  = fev2 + fes2(:,ib)
          if (leks .gt. 1) f(:,ib) = fks
          if (ipr .gt. 40) write(stdo,210) (c*fes2(m,ib),m=1,3),(c*fev2(m),m=1,3),(c*fks(m),m=1,3)
 210      format('  KS',3f8.2,1x,3f8.2,1x,3f8.2)
          if (ipl .gt. 0 .and. ipr .ge.30)
     .    write (stdl,711) ib,(c*fhar(m),m=1,3),(c*fks(m),m=1,3)
  711     format('fp ib',i4,'  fh ',3f8.2,'   fks ',3f8.2)
        endif
        ssite(ib)%force = f(:,ib)
      enddo
      call info5(10,0,0,' Maximum Harris force = %;3g mRy/au (site %i)'
     .//'%?#n#  Max eval correction = %;1d##',c*fmax,ibmx,isw(dfmax.gt.0),c*dfmax,0)

C     Symmetrize forces to machine precision
      ng=lat_nsgrp
      if (ng .gt. 1) then
        call info(30,1,0,' Symmetrize forces ...',0,0)
        allocate(wk_rv(3*nbas))
        call symfor ( nbas , 1 , rv_a_osymgr , ng , iv_a_oistab , wk_rv , f )
        if (allocated(wk_rv)) deallocate(wk_rv)
      endif
      fout=f
      call tcx('totfrc')
      end subroutine totfrc


