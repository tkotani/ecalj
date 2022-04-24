      module m_chkdmu
      real(8),protected,public:: eorb=0d0
      public:: chkdmu,sudmtu
      private
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


      subroutine sudmtu(dmatu,vorb)
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
      use m_lmfinit,only: nbas,nsp,nlibu,lmaxu,lldau,ssite=>v_ssite,sspec=>v_sspec,stdo
      use m_mksym,only: g=>rv_a_osymgr,istab=>iv_a_oistab, ng =>lat_nsgrp
C- Initialize site density matrix and vorb  for LDA+U
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu : nlibu total number of U blocks
Ci   lmaxu :dimensioning parameter for U matrix
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:   *
Ci     Passed to: symdmu rotycs
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idu uh jh
Ci     Stored:   *
Ci     Passed to: symdmu rotycs
Ci   idvsh :0 dmatu and vorb returned in real harmonics
Ci         :1 dmatu and vorb returned in spherical harmonics
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   istab :site istab(i,ig) is transformed into site i by grp op ig
Co Outputs
Co   dmatu :density matrix for LDA+U orbitals
Co         :in real spherical harmonics basis
Co   vorb  :orbital dependent potential matrices
Co         :in real spherical harmonics basis
Cl Local variables
Cl   eorb  : U contribution to LDA+U total energy
Cr Remarks
Cr   Reads in diagonal occupation numbers from file occnum.ext or dmatu
Cr   given in order of m=-l,l, isp=1,2, and constructs initial vorb
Cu Updates
Cu   12 Nov 07 Generate dmatu and vorb in either real or spher. harmonics
Cu   07 May 07 Bug fix MPI mode, when reading occnum instead of dmats
Cu   31 Jan 06 Setup and printouts in spherical harmonics
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   27 Apr 05 Lambrecht first created
C-------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer:: idvsh=0,i_copy_size !nbas,nsp,nlibu,lmaxu,
c      integer lldau(nbas),istab(nbas,ng)
c      real(8):: g(9,ng)
c      type(s_site)::ssite(*)
c      type(s_spec)::sspec(*)

      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C ... Local parameters
      logical rdstrn,parstr,mmtargetx,eee
      integer i,isp,ib,l,lmxa,m,m2,foccn,rdm,havesh,ivsiz,ipr,ifx
C     integer nglob,stdo
      double precision nocc(-3:3,2),iv(7)
      integer idu(4),igetss,is,idmat,fxst,iblu,nlm,nlmu,a2vec,nn
      double precision uh(4),jh(4),eorb,tmp(2,7,7),xx
      character str*80,spid*8,aaa*24
      complex(8) ,allocatable :: dmwk_zv(:)
      real(8):: uhx,uhxx
C ... MPI
      include "mpif.h"
      integer procid,master,mpipid,ierr
      logical mlog,occe,dexist

!!
      real(8)::sss
      call rxx(nsp.ne.2,'LDA+U must be spin-polarized!')
      procid = mpipid(1)
      master = 0
      call getpr(ipr)
C     ipr = 55

!! When LDAU is dummy (usually just in order to print our dmats file).
      sss=0d0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)
          sss = sss + sum(abs(sspec(is)%uh))+sum(abs(sspec(is)%jh))
        endif
      enddo
      if(sss<1d-6) then
        dmatu = 0d0
        havesh = 1
        goto 1185
      endif

C --- Read in dmatu if file  dmats.ext  exists ---
      if(procid .eq. master) then
        inquire(file='dmats.'//trim(sname),exist=dexist)
        if(dexist) then
          open(newunit=idmat,file='dmats.'//trim(sname)) 
          rewind idmat
  825     continue
          if (rdstrn(idmat,str,len(str),.false.)) then
            if (str(1:1) .eq. '#') goto 825
            i = 0
            if (parstr(str,'sharm ',len(str)-6,5,' ',i,m)) then
              havesh = 1
            else
              havesh = 0
            endif
          endif
          call info2(20,0,0,' sudmtu:  reading density matrix from file'//
     .    ' dmats in %?#n#spherical#real# harmonics',havesh,0)
          backspace idmat
          iblu = 0
          do  ib = 1, nbas
            if (lldau(ib) .ne. 0) then
              is = int(ssite(ib)%spec)
              lmxa=sspec(is)%lmxa
              i_copy_size=size(sspec(is)%idu)
              call icopy(i_copy_size,sspec(is)%idu,1,idu,1)
              do l = 0, min(lmxa,3)
                if (idu(l+1) .ne. 0) then
                  iblu = iblu+1
                  nlm = 2*l+1
                  nlmu = 2*lmaxu+1
                  do  isp = 1, 2
                    if (rdm(idmat,40,2*nlm**2,' ',tmp,nlm,nlm) .ne. 2)
     .              call rxi('sudmtu failed to read dmats for site',ib)
                    call zmscop(0,nlm,nlm,nlm,nlmu,0,0,0,0,tmp,
     .              dmatu(-l,-l,isp,iblu))
                  enddo
                endif
              enddo
            endif
          enddo

C --- Otherwise read in occnumbers and construct diagonal dmatu ---
        else
          inquire(file='occnum.'//trim(sname),exist=occe)
          if(occe) then !if no occunum.* initial dmatu=0 
             call info0(20,1,0,' sudmtu:  initial (diagonal) density-matrix from occ numbers')
             open(newunit=foccn,file='occnum.'//trim(sname))
             havesh = 1
 12          if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
             if (str(1:1) .eq. '#') goto 12
             if (str(1:1) .eq. '%') then
                i = 0
                if (parstr(str,'real ',len(str)-5,4,' ',i,m)) havesh = 0
             else
                backspace foccn
             endif
             iblu = 0
             do  ib = 1, nbas
                if (lldau(ib) .ne. 0) then
                   is = int(ssite(ib)%spec)
                   lmxa=sspec(is)%lmxa
                   i_copy_size=size(sspec(is)%idu)
                   call icopy(i_copy_size,sspec(is)%idu,1,idu,1)
                   do l = 0,min(lmxa,3)
                      if (idu(l+1) .ne. 0) then
                         iblu = iblu+1
                         do  isp = 1, 2
 11                         continue
                            if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
C     Skip comment lines
                            if (str(1:1) .eq. '#') goto 11
                            i = 0
                            m = a2vec(str,len(str),i,4,', ',2,3,2*l+1,iv,nocc(-l,isp))
                            if (m .lt. 0) goto 99
                         enddo
                         call info8(20,0,0,' occ numbers, site %i l=%i:  '//
     .                        '%n:1d (spin 1)  %n:1d (spin 2)',
     .                        ib,l,2*l+1,nocc(-l,1),2*l+1,nocc(-l,2),0,0)
                         do isp = 1, 2
                         do m = -l, l
                         do m2 = -l, l
                            dmatu(m,m2,isp,iblu) = dcmplx(0d0,0d0)
                         enddo
                            dmatu(m,m,isp,iblu) = dcmplx(nocc(m,isp),0d0)
                         enddo
                         enddo
                      endif
                   enddo
                endif
             enddo
             close(foccn) !call fclose(foccn)
          else
             dmatu=0d0
          endif
        endif
      endif
!!
 1185 continue

C ... Initial printout
      call praldm(0,51,51,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .' dmats read from disk',dmatu)
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call mpibc1(dmatu,2*ivsiz,4,mlog,'sudmtu','dmatu')
      call mpibc1(havesh,1,2,.false.,' ',' ')

C ... Density matrix in real or spherical harmonics (fixed by idvsh)
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = idvsh
      endif

C ... Symmetrize dmatu (symdmu requires real harmonics)
      allocate(dmwk_zv(ivsiz))
      dmwk_zv=0d0
      if (havesh .eq. 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      call symdmu ( dmatu , dmwk_zv , nbas , nsp , lmaxu , sspec 
     ., ssite , ng , g , istab , lldau , xx )

      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs ( 2 * idvsh - 1 , dmwk_zv , nbas , nsp , lmaxu 
     .  , sspec , ssite , lldau )

        havesh = idvsh
      endif
      if (ng .ne. 0) then
        call info2(30,0,0,' sudmtu:  RMS change in dmats'//
     .  ' from symmetrization = %,6d',xx,0)
        if (xx .gt. .01d0) call info0(30,0,0,
     .  '          (warning) RMS change unexpectely large')
        call daxpy ( ivsiz * 2 , - 1d0 , dmatu , 1 , dmwk_zv , 1 )

        call info0(60,0,0,' change in dmat wrought by symmetrization')
        call praldm ( 0 , 60 , 60 , 0 , nbas , nsp , lmaxu , lldau , 
     .  sspec , ssite , ' ' , dmwk_zv )

        if (allocated(dmwk_zv)) deallocate(dmwk_zv)

      endif

C     Print dmats in specified harmonics
      allocate(dmwk_zv(abs(-ivsiz)))
      if (-ivsiz<0) dmwk_zv(:)=0.0d0

      call dcopy ( ivsiz * 2 , dmatu , 1 , dmwk_zv , 1 )

      if (havesh .ne. idvsh) then
        call rotycs ( 2 * idvsh - 1 , dmwk_zv , nbas , nsp , lmaxu 
     .  , sspec , ssite , lldau )

      endif
      call info0(30,0,0,' ')
      call praldm ( 0 , 30 , 30 , idvsh , nbas , nsp , lmaxu , lldau 
     ., sspec , ssite , ' Symmetrized dmats' , dmwk_zv )


C     Print dmats in complementary harmonics
c      if (ipr.ge.45 .or. ipr.ge.40 .and. idvsh.eq.0) then
        i = 1-idvsh
        call rotycs ( 2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , sspec 
     .  , ssite , lldau )

        call info0(30,0,0,' ')
        call praldm ( 0 , 30 , 30 , i , nbas , nsp , lmaxu , lldau , 
     .  sspec , ssite , ' Symmetrized dmats' , dmwk_zv )

c      endif
      if (allocated(dmwk_zv)) deallocate(dmwk_zv)


C ... Make Vorb (ldau requires spherical harmonics)
      if (havesh .ne. 1) then
        call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 1
      endif

c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c$$$!! experimental block to keep magnetic moment for AF.
c$$$       inquire(file='mmtarget.aftest',exist=mmtargetx)
c$$$       inquire(file='mmagfield.aftest',exist=eee)
c$$$       if(mmtargetx.and. (procid==master)) then
c$$$          uhx=0d0
c$$$          open(newunit=ifx,file='mmagfield.aftest')
c$$$          do
c$$$             read(ifx,*,end=1112,err=1112) uhxx,aaa
c$$$             uhx=uhxx
c$$$             if(trim(aaa)=='!Magfield') uhx=uhxx
c$$$          enddo
c$$$ 1112     continue
c$$$          close(ifx)
c$$$          write(stdo,"('sudmtu: mmtarget mode. Readin Magfield from mmagfield.aftest=',f10.6)")uhx
c$$$          open(newunit=ifx,file='mmagfield.aftest')
c$$$          write(ifx,"(d23.15,1x,'!Magfield is read from previous mmagfield.aftest')") uhx
c$$$          close(ifx)
c$$$       endif
c$$$       call mpibc1_real(uhx,1,'sudmtu_uhx')
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      iblu = 0
      do  20  ib = 1, nbas
        if (lldau(ib) .eq. 0) goto 20
        is = int(ssite(ib)%spec)
        lmxa=sspec(is)%lmxa
        i_copy_size=size(sspec(is)%idu)
        call icopy(i_copy_size,sspec(is)%idu,1,idu,1)
        i_copy_size=size(sspec(is)%uh)
        call dcopy(i_copy_size,sspec(is)%uh,1,uh,1)
        i_copy_size=size(sspec(is)%jh)
        call dcopy(i_copy_size,sspec(is)%jh,1,jh,1)

c        do i_spacks=is,is
c          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c        enddo
        spid=sspec(is)%name

        i = min(lmxa,3)
        call info8(30,1,0,' Species '//spid//
     .  '%a: mode=%n:1i    U=%n:1d    J=%n:1d',
     .  i+1,idu,i+1,uh,i+1,jh,0,0)
        do  22  l = 0, i
          if (idu(l+1) .ne. 0) then
            iblu = iblu+1
c            if(.not.mmtargetx) uhx=uh(l+1)
            uhx=uh(l+1)
            call ldau(idu(l+1),l,iblu,uhx,jh(l+1),dmatu,nsp,
     .      lmaxu,vorb,eorb)
C            call prdmts(0,30,30,0,' sudmtu:  Vorb in spherical harm',ib,
C     .        l,lmaxu,iblu,vorb,nsp,1)

          endif
   22   continue
   20 continue
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .' Unsymmetrized vorb',vorb)
C     At this point, dmatu and vorb are in spherical harmonics

C ... Symmetrize vorb to check (symdmu requires real harmonics)
      allocate(dmwk_zv(abs(-ivsiz)))
      if (-ivsiz<0) dmwk_zv(:)=0.0d0

      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
      call symdmu ( vorb , dmwk_zv , nbas , nsp , lmaxu , sspec 
     ., ssite , ng , g , istab , lldau , xx )

      if (allocated(dmwk_zv)) deallocate(dmwk_zv)


C ... Exit with vorb,dmatu in spherical or real harmonics
C     EITHER: vorb =>  spherical harmonics OR dmatu => real harmonics
      if (idvsh .eq. 1) then
        call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 1
      else
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      if (ng .ne. 0) then
        call info2(30,1,0,' sudmtu:  RMS change in vorb '//
     .  'from symmetrization = %,6d',xx,0)
        if (xx .gt. .01d0) call info0(30,0,0,
     .  '          (warning) RMS change unexpectely large')
      endif

C     Print vorb in specified harmonics
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .' Symmetrized vorb',vorb)

C     Print vorb in complementary harmonics
c      if (ipr.ge.45 .or. ipr.ge.40 .and. idvsh.eq.0) then
        i = 1-idvsh
        allocate(dmwk_zv(abs(-ivsiz)))
        if (-ivsiz<0) dmwk_zv(:)=0.0d0

C        call dcopy(ivsiz*2,dmatu,1,w(odmwk),1)
C        call rotycs(2*i-1,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,lldau)
C        call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,sspec,ssite,
C     .    ' dmats',w(odmwk))
        call dcopy ( ivsiz * 2 , vorb , 1 , dmwk_zv , 1 )

        call rotycs ( 2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , sspec 
     .  , ssite , lldau )

        call info0(30,0,0,' ')
        call praldm ( 0 , 30 , 30 , i , nbas , nsp , lmaxu , lldau , 
     .  sspec , ssite , ' Vorb' , dmwk_zv )

        if (allocated(dmwk_zv)) deallocate(dmwk_zv)

c      endif

      eorb = 0d0

C     call prmx('vorb',vorb(1,1,1,1),2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

C --- Error exit ---
      return
   99 continue
      call awrit2('bad occnum file, site %i, l=%i',str,len(str),0,ib,l)
      call rx(str)
      end subroutine sudmtu

      subroutine rotycs(mode,a,nbas,nsp,lmaxu,sspec,ssite,lldau)
      use m_struc_def  !Cgetarg
C- Rotate matrix a from real to spherical harmonics
C for LDA+U objects densmat and vorb
C-------------------------------------
Ci mode =1 from real to spherical
Ci      -1 from spherical to real
Ci a matrix to be transformed a(m,m,isp,iblu)  could be vorb or densmat
Ci nbas : number of sites
Ci nsp  : number of spins
Ci lmaxu: lmax for U
Ci sspec: species info
Ci ssite: sites info
Ci lldau  :lldau(ib)=0 => no U on this site otherwise
Ci        :U on site ib with dmat in dmats(*,lldau(ib))
Co a rotated in place
Cr Remarks
Cr order of cubic harmonics ls (l-1)s,ms...1s 0 1c mc... (l-1)c lc
Cr order of spherical harmonics -l:l
Cr Yl-m=(Ylmc-iYlms)/sqrt(2)  Ylm=(-1)**m*conjg(Yl-m)
Cu Updates
Cu   18 Jan 06 A. Chantis changed rotation matrices in accordance with
Cu             the definition of real harmonics used in the rest of
Cu             the code (Hund's rules satisfied as indicated by orb. moment)
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   30 Apr 05 Lambrecht first created
C----------------------------------------------------------------
      implicit none
C passed variables
      integer nbas,lldau(nbas),mode,lmaxu,nsp,i_copy_size
      double complex a(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)

C local variables
      integer ib,m,l,idu(4),lmxa,is,igetss,i,j,k,ll,isp,iblu
      double complex rot(2,2)
      double complex b(2,2),c(2,2),add
      double precision s2

      s2 = 1/dsqrt(2d0)
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)
          lmxa=sspec(is)%lmxa
          i_copy_size=size(sspec(is)%idu)
          call icopy(i_copy_size,sspec(is)%idu,1,idu,1)

          do  l = 0, min(lmxa,3)
            if (idu(l+1) .ne. 0) then
              iblu = iblu+1
              do  isp = 1, 2
                do  m = 1, l
                  b(1,1) = a(m,m,isp,iblu)
                  b(1,2) = a(m,-m,isp,iblu)
                  b(2,1) = a(-m,m,isp,iblu)
                  b(2,2) = a(-m,-m,isp,iblu)
                  if (mode .eq. -1) then
C rotation from spherical  to cubic basis
                    rot(1,1) = dcmplx(s2,0d0)
                    rot(1,2) = dcmplx(s2*(-1d0)**m,0d0)
C                   rot(2,1) = dcmplx(0d0,-s2)
C                   rot(2,2) = dcmplx(0d0,s2*(-1d0)**m)
                    rot(2,1) = dcmplx(0d0,s2)
                    rot(2,2) = dcmplx(0d0,-s2*(-1d0)**m)
                  elseif (mode .eq. 1) then
C rotation from cubic  to spherical  basis
                    rot(1,1) = dcmplx(s2,0d0)
C                   rot(1,2) = dcmplx(0d0,s2)
                    rot(1,2) = dcmplx(0d0,-s2)
                    rot(2,1) = dcmplx(s2*(-1d0)**m,0d0)
C                   rot(2,2) = dcmplx(0d0,-s2*(-1d0)**m)
                    rot(2,2) = dcmplx(0d0,s2*(-1d0)**m)
                  else
                    call rx('ROTYCS: mode must be 1 or -1')
                  endif
C calculate matrix product c=rot*b*rot^+
                  do  i = 1, 2
                    do  j = 1, 2
                      add = dcmplx(0d0,0d0)
                      do  k = 1, 2
                        do  ll = 1, 2
                          add = add + rot(i,k)*b(k,ll)*dconjg(rot(j,ll))
                        enddo
                      enddo
                      c(i,j) = add
                    enddo
                  enddo
C place c in appropriate place
                  a(m,m,isp,iblu) = c(1,1)
                  a(m,-m,isp,iblu) = c(1,2)
                  a(-m,m,isp,iblu) = c(2,1)
                  a(-m,-m,isp,iblu) = c(2,2)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      end subroutine rotycs
      end module m_chkdmu
