      module m_chkdmu
      real(8),protected:: eorb=0d0
      
      contains
      subroutine chkdmu(eks, dmatu,dmatuo,vorb) !,tolu,umix)
      use m_struc_def,only: s_site,s_spec
      use m_lmfinit,only: stdl,nbas,nsp,nlibu,lmaxu,ssite=>v_ssite,sspec=>v_sspec,lldau,
     &     tolu=>mix_tolu,umix=>mix_umix
      use m_MPItk,only: master_mpi
      use m_mksym,only: g=>rv_a_osymgr,istab=>iv_a_oistab, ng =>lat_nsgrp
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
C- LDA+U total energy and mixing of lda+U density matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu
Ci   lmaxu :dimensioning parameter for U matrix
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored: *
Ci     Passed to: symdmu rotycs
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idu uh jh
Ci     Stored: *
Ci     Passed to: symdmu rotycs
Ci   idvsh :0 dmatu, dmatuo, vorb input/output in real harmonics
Ci         :1 dmatu, dmatuo, vorb input/output in spherical harmonics
Ci   dmatu : dmatu produced in current iteration
Ci         : dmatu is passed in real harmonics
Ci   dmatuo: dmatu produced in prior iteration
Ci         : dmatuo is passed in real harmonics
Ci   tolu  :convergence tolerance density-matrix
Ci   umix  :linear mixing parameter for density matrix
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   istab :site istab(i,ig) is transformed into site i by grp op ig
Cio Inputs/Outputs
Cio  vorb  :orbital dependent potential matrices
Cio        :vorb is updated on output
Cl Local variables
Cl   eorb  : U contribution to LDA+U total energy
Cb Bugs
Cb   This routine should not update vorb
Cu Updates
Cu   12 Nov 07 Works with dmatu and vorb in either real or spher. harmonics
Cu   31 Jan 06 Printouts in spherical harmonics
Cu   09 Nov 05 Convert dmat to complex form
Cu   29 Oct 05 doesn't update vorb in tot. E eval; restores dmatu if conv.
Cu    2 Jun 05 Evaluates total energy contribution from output dmatu
Cu   27 Apr 05 Lambrecht first created
C ----------------------------------------------------------------------
      implicit none
      integer:: ierr
      include "mpif.h"
      integer:: idvsh=0,i_copy_size
      real(8):: eks, eorbxxx
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex dmatuo(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      integer l,idu(4),lmxa,ib,is,iblu,igetss,idmat,ivsiz,ifile_handle
      integer iprint,ipl,havesh
      double precision ddmat,uh(4),jh(4),eorbi,eterms(20),ddot,xx
      complex(8) ,allocatable :: dmwk_zv(:)
      logical:: fexist,mmtargetx,eee
      real(8),allocatable:: uhall(:,:)
      real(8):: mmsite(nbas),uhxx,mmhist(10000),uhhist(10000),mmtarget,uhdiff
      real(8),save::uhxnew,uhx,alpha,alphax
      integer:: nn,ibas,ifx,key,i,nit
      integer,save::ncount=0
      complex(8):: dmatuav(-lmaxu:lmaxu,-lmaxu:lmaxu)
      real(8):: sss(1),sigin
      real(8):: fac
C     print *, 'dmatu',dmatu(-lmaxu,-lmaxu,1,1)
C     call rx0('done')
C     call prmx('dmatu',dmatu,2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

      if (nlibu .eq. 0) return
      havesh = idvsh
      ipl = 1
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call info0(20,1,0,' chkdmu:  '//
     .'check LDA+U density-matrix for convergence and update ...')
C --- Symmetrize output dmatu (req. real harmonics); compare diff ---
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .' Unsymmetrized output dmats',dmatu)
      if (havesh .eq. 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      allocate(dmwk_zv(ivsiz))
      dmwk_zv(:)=0.0d0
      call symdmu ( dmatu , dmwk_zv , nbas , nsp , lmaxu , sspec 
     ., ssite , ng , g , istab , lldau , xx )
      deallocate(dmwk_zv)
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = idvsh
      endif
C     RMS change : dmatu-dmatuo; restore dmatuo
      call daxpy(2*ivsiz,-1d0,dmatu,1,dmatuo,1)
      ddmat = dsqrt(ddot(2*ivsiz,dmatuo,1,dmatuo,1)/(2*ivsiz))
      call daxpy(2*ivsiz,1d0,dmatu,1,dmatuo,1)
C --- Printout dmatu in real or spherical harmonics, fixed by idvsh ---
      call info2(30,0,0,' chkdmu:  RMS change in dmat'//
     .' from symmetrization = %,6d',xx,0)
C --- Compute U contribution to total energy; make vorb ---
C     This block requires dmatu to be in spherical harmonics
      if (havesh .ne. 1) then
        call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 1
      endif
      call info0(20,1,0,'%9pLDA+U total energy ...')
      eorb = 0
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)
          lmxa=sspec(is)%lmxa
          i_copy_size=size(sspec(is)%idu)
          call icopy(i_copy_size,sspec(is)%idu,1,idu,1)
          i_copy_size=size(sspec(is)%uh)
          call dcopy(i_copy_size,sspec(is)%uh,1,uh,1)
          i_copy_size=size(sspec(is)%jh)
          call dcopy(i_copy_size,sspec(is)%jh,1,jh,1)
          do  l = 0, min(lmxa,3)
            if (idu(l+1) .ne. 0) then
              iblu = iblu+1
              eorbi = 999
              call ldau(100+idu(l+1),l,iblu,uh(l+1),jh(l+1),dmatu,nsp,
     .        lmaxu,vorb,eorbi)
              eorb = eorb + eorbi
            endif
          enddo
        endif
      enddo
C --- LDA total energy terms ---
      call info5(20,0,0,'%9peks = %,6;6d  '//
     .'e[U] = %,6;6d  Etot(LDA+U) = %,6;6d',eks,eorb,eks+eorb,0,0)
      if (master_mpi .and. ipl.gt.0)
     .write (stdl,720) eks,eorb,eks+eorb
  720 format('ldau EHK ',f14.6,'  U',f12.6,'  ELDA+U ',f14.6)

C --- Restore dmatu, vorb to harmonics specified by idvsh ---
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs(2*idvsh-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = idvsh
      endif
C --- Case self-consistent within tolerance tolu ---
      if (ddmat .lt. tolu) then
        call info5(20,0,0,' LDA+U dmatu converged'//
     .  '  RMS diff (%;3g) < tolu (%;3g)',ddmat,tolu,0,0,0)
C       Restore dmatuo to dmatu
        call dcopy(2*ivsiz,dmatuo,1,dmatu,1)
        return
C --- Case not self-consistent ---
      else
        call info0(20,0,0,'%9pLDA+U update density matrix ...')
        call info5(20,0,0,'%9fRMS diff in dens mat(%;3g) > tolu (%;3g)'
     .  //' Linear mix with beta=%;3g',ddmat,tolu,umix,0,0)
C   ... Make new dmatu by mixing    new*umix + old*(1-umix)
        call dscal(2*ivsiz,umix,dmatu,1)
        call daxpy(2*ivsiz,1-umix,dmatuo,1,dmatu,1)
C   --- Make Vorb from mixed dmatu ---
C       This block requires dmatu to be in spherical harmonics
        if (havesh .ne. 1) then
          call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
          havesh = 1
       endif
       
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccc   
c$$$!! experimental block to keep magnetic moment for AF.
c$$$       inquire(file='mmtarget.aftest',exist=mmtargetx)
c$$$       if(mmtargetx) then
c$$$          if(master_mpi) then
c$$$             open(newunit=ifx,file='mmtarget.aftest')
c$$$             alpha  = .1d0
c$$$             read(ifx,*) mmtarget !,alpha
c$$$             close(ifx)
c$$$             open(newunit=ifx,file='uhval.aftest')
c$$$             nit=0
c$$$             mmsite=0d0
c$$$             uhx=0d0
c$$$             do
c$$$                read(ifx,*,end=1112,err=1112) uhx
c$$$                read(ifx,*,end=1112,err=1112) (mmsite(ibas),ibas=1,nbas)
c$$$                nit=nit+1
c$$$                mmhist(nit)=(mmsite(1)-mmsite(2))/2d0
c$$$                uhhist(nit)=uhx
c$$$             enddo
c$$$ 1112        continue
c$$$             close(ifx)
c$$$             write(6,"('uhval: ',i5,f10.6,2x,12f10.6)")nit,uhx,(mmsite(ibas),ibas=1,nbas)
c$$$!     ! Generate new uhx based on the  history of uhx mmsites for given mm
c$$$!     ! test uh          
c$$$             uhx= uhx + alpha*(mmhist(nit)-mmtarget)**2 - 2d0*(mmhist(nit)-mmtarget)
c$$$             write(6,"('mmhist0: UH mm',i5, 3d13.4)" ) nit,uhx,mmhist(nit)
c$$$             sss(1)=uhx
c$$$             call mixuh(sss)
c$$$             uhx=sss(1)
c$$$             write(6,"('mmhist:  UH mm ',i5, 3d13.4)") nit,uhx,mmhist(nit)
c$$$             if(ncount==0) then
c$$$                open(newunit=ifx,file="mixuh.aftest")
c$$$                close(ifx,status='delete')
c$$$                ncount=1
c$$$             endif
c$$$             open(newunit=ifx,file='uhval.aftest',position='append')
c$$$             write(ifx,"(d23.15,' !uhx')") uhx !, uhxnew
c$$$             close(ifx)
c$$$          endif
c$$$!! broadcast uhx       
c$$$          call MPI_Bcast(uhx, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
c$$$! Averaged dmatu up dn.
c$$$          dmatuav=0d0
c$$$          do i=-lmaxu,lmaxu
c$$$             dmatuav(i,i) =  1d0
c$$$          enddo
c$$$          do iblu =1,2
c$$$             fac=1d0
c$$$             if(iblu==2) fac=-1d0
c$$$c             dmatuav(:,:) =  .5d0*dmatu(:,:,1,iblu) +.5d0*dmatu(:,:,2,iblu)
c$$$             dmatu(:,:,1,iblu) =  fac*dmatuav(:,:)
c$$$             dmatu(:,:,2,iblu) = -fac*dmatuav(:,:)
c$$$          enddo
c$$$
c$$$       endif
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
        iblu = 0
        do  ib = 1, nbas
          if (lldau(ib) .ne. 0) then
            is = int(ssite(ib)%spec)
            lmxa=sspec(is)%lmxa
            i_copy_size=size(sspec(is)%idu)
            call icopy(i_copy_size,sspec(is)%idu,1,idu,1)
            i_copy_size=size(sspec(is)%uh)
            call dcopy(i_copy_size,sspec(is)%uh,1,uh,1)
            i_copy_size=size(sspec(is)%jh)
            call dcopy(i_copy_size,sspec(is)%jh,1,jh,1)
            do  l = 0, min(lmxa,3)
              if (idu(l+1) .ne. 0) then
                iblu = iblu+1
                call pshpr(iprint()-20)
c$$$                if(.not.mmtargetx) uhx = uh(l+1)
                uhx = uh(l+1)
                call ldau(idu(l+1),l,iblu,uhx,jh(l+1),dmatu,nsp, !uhxx,jh(l+1),dmatu,nsp, !
     .          lmaxu,vorb,eorbxxx) !eorbxxx is dummy?
                call poppr
              endif
            enddo
          endif
        enddo
       
C   ... Symmetrize vorb to check (symdmu requires real harmonics)
        allocate(dmwk_zv(ivsiz))
        dmwk_zv(:)=0.0d0
        call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        call symdmu ( vorb , dmwk_zv , nbas , nsp , lmaxu , sspec 
     .  , ssite , ng , g , istab , lldau , xx )
        call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        if (allocated(dmwk_zv)) deallocate(dmwk_zv)

!! ==>  At this point, dmatu and vorb are in spherical harmonics

        if(Iprint()>20) write(6,"(a,f9.4)")
     &    '  RMS change in vorb from symmetrization =',xx
        if(xx>.0001d0.and.iprint()>30)
     &   write(6,'(a)')' (warning) RMS change unexpectely large'
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite, ' Mixed dmats',dmatu)
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite, ' New vorb',vorb)
!! Write dmatu to file
        if (master_mpi) then
          idmat = ifile_handle()
          open(idmat,file='dmats.'//trim(sname))
          call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' mixed dmats',dmatu)
          close(idmat)
        endif
!! Exit with dmatu, vorb in real harmonics, depending on idvsh
        if (idvsh .eq. 0) then
          call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
          call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' Mixed dmats',dmatu)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' New vorb',vorb)
        endif
      endif
      end subroutine chkdmu



      
      end module
