      module m_chkdmu
      use m_ftox
      real(8),protected,public:: eorb=0d0
      public:: chkdmu,sudmtu
      private
      contains
      subroutine chkdmu(eks, dmatu,dmatuo,vorb) !,tolu,umix)
      use m_struc_def,only: s_site,s_spec
      use m_lmfinit,only: stdl,nbas,nsp,nlibu,lmaxu,ssite=>v_ssite,sspec=>v_sspec,lldau,
     &     tolu=>mix_tolu,umix=>mix_umix,stdo
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
Ci         :xxxx1 dmatu, dmatuo, vorb input/output in spherical harmonics
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
Co  vorb  :orbital dependent potential matrices
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
c      complex(8) :: dmwk_zv(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      logical:: fexist,mmtargetx,eee
      real(8),allocatable:: uhall(:,:)
      real(8):: mmsite(nbas),uhxx,mmhist(10000),uhhist(10000),mmtarget,uhdiff
      real(8),save::uhxnew,uhx,alpha,alphax
      integer:: nn,ibas,ifx,key,i,nit
      integer,save::ncount=0
      complex(8):: dmatuav(-lmaxu:lmaxu,-lmaxu:lmaxu)
      real(8):: sss(1),sigin
      real(8):: fac
      if (nlibu .eq. 0) return
      havesh =0
      idvsh  =0  ! We assume real harmonics for i/o
      ipl = 1
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
C --- Symmetrize output dmatu (req. real harmonics); compare diff ---
      if(iprint()>=60) call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .' Unsymmetrized output dmats',dmatu)
      call symdmu(nlibu,dmatu, nbas,nsp, lmaxu, sspec, ssite, ng, g, istab, lldau, xx)
      ddmat = dsqrt(ddot(2*ivsiz,dmatu-dmatuo,1,dmatu-dmatuo,1)/(2*ivsiz))
      if(master_mpi)write(stdo,ftox)' chkdmu: LDA+U. RMSdiff of dmat from symmetrization =',ftod(xx,2)
      
C --- Compute U contribution to total energy; make vorb ---
      call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !mode=1, dmatu is from rh to sh 
      havesh = 1
      eorb = 0
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)
          lmxa=sspec(is)%lmxa
          idu=sspec(is)%idu
          uh=sspec(is)%uh
          jh=sspec(is)%jh
          do  l = 0, min(lmxa,3)
            if (idu(l+1) .ne. 0) then
              iblu = iblu+1
              eorbi = 999
              call ldau(100+idu(l+1),l,iblu,uh(l+1),jh(l+1),dmatu,nsp,lmaxu,vorb,eorbi)
              eorb = eorb + eorbi
            endif
          enddo
        endif
      enddo
C --- LDA total energy terms ---
      if(iprint()>=60) write(stdo,ftox)' eks =',ftof(eks),'e[U] =',ftof(eorb),'Etot(LDA+U) =',ftof(eks+eorb)
      if(master_mpi) write (stdl,ftox) 'ldau EHK',ftof(eks),'U',ftof(eorb),'ELDA+U',ftof(eks+eorb)
      
C --- Restore dmatu, vorb to real harmonics 
      call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !-1, from sh to rh idvsh=0
      havesh = 0
      if(master_mpi)write(stdo,ftox)' LDA+U update density matrix ...'
      if(master_mpi)write(stdo,ftox)' RMS diff in densmat',ftod(ddmat) 
      dmatu = dmatu+(1d0-umix)*dmatuo ! new*umix + old*(1-umix)
      call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !from rh to sh
      havesh = 1
      
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
      vorb=1d99
      iblu = 0
      do  ib = 1, nbas
         if (lldau(ib) .ne. 0) then
            is = int(ssite(ib)%spec)
            lmxa=sspec(is)%lmxa
            idu=sspec(is)%idu
            uh=sspec(is)%uh
            jh=sspec(is)%jh
            do  l = 0, min(lmxa,3)
               if (idu(l+1) .ne. 0) then
                  iblu = iblu+1
                  call pshpr(iprint()-20)
c$$$  if(.not.mmtargetx) uhx = uh(l+1)
                  uhx = uh(l+1)
                  call ldau(idu(l+1),l,iblu,uhx,jh(l+1),dmatu,nsp, !uhxx,jh(l+1),dmatu,nsp, !
     .                 lmaxu,vorb,eorbxxx) !eorbxxx is dummy?
                  call poppr
               endif
            enddo
         endif
      enddo
C  ... Symmetrize vorb to check (symdmu requires real harmonics)
      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau) !vorb from sh to rh
      call symdmu(nlibu,vorb,nbas , nsp , lmaxu , sspec, ssite , ng , g , istab , lldau , xx )
      call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau) !vorb from rh to sh
      
!!=>  At this point, dmatu and vorb are in spherical harmonics
      if(Iprint()>20) write(6,"(a,f9.4)")'  RMS change in vorb from symmetrization =',xx
      if(xx>.0001d0.and.iprint()>30) write(6,'(a)')' (warning) RMS change unexpectely large'
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite, ' Mixed dmats',dmatu)
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite, ' New vorb',vorb)
      if (master_mpi) then
         idmat = ifile_handle()
         open(idmat,file='dmats.'//trim(sname)) !havesh mush be 1
         call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' mixed dmats',dmatu)
         close(idmat)
      endif
!! write in real harmonics      
      call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau) !from sh to rh
      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)  !from sh to rh
      havesh=0                  !I recovered this 2022May8
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' Mixed dmats',dmatu)
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' New vorb',vorb)
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
      integer:: idvsh=0,i_copy_size !nbas,nsp,nlibu,lmaxu,
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      logical rdstrn,parstr,mmtargetx,eee
      integer i,isp,ib,l,lmxa,m,m2,foccn,rdm,havesh,ivsiz,ipr,ifx
      double precision nocc(-3:3,2),iv(7)
      integer idu(4),igetss,is,idmat,fxst,iblu,nlm,nlmu,a2vec,nn
      double precision uh(4),jh(4),eorb,xx !tmp(2,7,7)
      complex(8):: tmp(7,7)
      character str*80,spid*8,aaa*24
      complex(8) :: dmwk_zv(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      real(8):: uhx,uhxx
C ... MPI
      include "mpif.h"
      integer procid,master,mpipid,ierr
      logical mlog,occe,dexist
      real(8)::sss
      character(128):: bbb
      call rxx(nsp.ne.2,'LDA+U must be spin-polarized!')
      procid = mpipid(1)
      master = 0
      call getpr(ipr)
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
         if(havesh ==1) bbb='spherical harmonics'
         if(havesh ==0) bbb='real harmonics'
         write(stdo,*)' sudmtu: reading density matrix from file dmats in '//trim(bbb)
          backspace idmat
          iblu = 0
          do  ib = 1, nbas
            if (lldau(ib) .ne. 0) then
              is = int(ssite(ib)%spec)
              lmxa=sspec(is)%lmxa
              idu=sspec(is)%idu
              do l = 0, min(lmxa,3)
                if (idu(l+1) .ne. 0) then
                  iblu = iblu+1
                  nlm = 2*l+1
                  nlmu = 2*lmaxu+1
                  do  isp = 1, 2
                    if (rdm(idmat,40,2*nlm**2,' ',tmp(1:nlm,1:nlm),nlm,nlm) .ne. 2)
     .              call rxi('sudmtu failed to read dmats for site',ib)
                    dmatu(-l:l,-l:l,isp,iblu)=tmp(1:nlm,1:nlm)
                  enddo
                endif
              enddo
            endif
          enddo
        else
           dmatu=0d0 !0 at LDA 
        endif   
      endif   
c$$$        else
c$$$          inquire(file='occnum.'//trim(sname),exist=occe)
c$$$          if(occe) then !if no occunum.* initial dmatu=0 
c$$$             write(stdo,*) ' sudmtu:  initial (diagonal) density-matrix from occ numbers'
c$$$             open(newunit=foccn,file='occnum.'//trim(sname))
c$$$             havesh = 1
c$$$ 12          if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
c$$$             if (str(1:1) .eq. '#') goto 12
c$$$             if (str(1:1) .eq. '%') then
c$$$                i = 0
c$$$                if (parstr(str,'real ',len(str)-5,4,' ',i,m)) havesh = 0
c$$$             else
c$$$                backspace foccn
c$$$             endif
c$$$             iblu = 0
c$$$             do  ib = 1, nbas
c$$$                if (lldau(ib) .ne. 0) then
c$$$                   is = int(ssite(ib)%spec)
c$$$                   lmxa=sspec(is)%lmxa
c$$$                   idu=sspec(is)%idu
c$$$                   do l = 0,min(lmxa,3)
c$$$                      if (idu(l+1) .ne. 0) then
c$$$                         iblu = iblu+1
c$$$                         do  isp = 1, 2
c$$$ 11                         continue
c$$$                            if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
c$$$C     Skip comment lines
c$$$                            if (str(1:1) .eq. '#') goto 11
c$$$                            i = 0
c$$$                            m = a2vec(str,len(str),i,4,', ',2,3,2*l+1,iv,nocc(-l,isp))
c$$$                            if (m .lt. 0) goto 99
c$$$                         enddo
c$$$                         do isp=1,2
c$$$                            write(stdo,ftox)' occ num: site',ib,'l',l,'isp',isp,' ',nocc(-l:l,isp)
c$$$                         enddo   
c$$$                         do isp = 1, 2
c$$$                         do m = -l, l
c$$$                         do m2 = -l, l
c$$$                            dmatu(m,m2,isp,iblu) = dcmplx(0d0,0d0)
c$$$                         enddo
c$$$                            dmatu(m,m,isp,iblu) = dcmplx(nocc(m,isp),0d0)
c$$$                         enddo
c$$$                         enddo
c$$$                      endif
c$$$                   enddo
c$$$                endif
c$$$             enddo
c$$$             close(foccn) !call fclose(foccn)
c$$$          else
c$$$             dmatu=0d0
c$$$          endif
c$$$        endif
c$$$      endif
 1185 continue

C ... Initial printout
      call praldm(0,51,51,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' dmats read from disk',dmatu)
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call mpibc1(dmatu,2*ivsiz,4,mlog,'sudmtu','dmatu')
      call mpibc1(havesh,1,2,.false.,' ',' ')
C ... Density matrix in real or spherical harmonics (fixed by idvsh)
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = idvsh
      endif
C ... Symmetrize dmatu (symdmu requires real harmonics)
      dmwk_zv=dmatu
      if (havesh .eq. 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      call symdmu(nlibu,dmatu , nbas , nsp , lmaxu , sspec, ssite , ng , g , istab , lldau , xx )
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs(2*idvsh-1,dmwk_zv,nbas,nsp,lmaxu,sspec,ssite, lldau )
        havesh = idvsh
      endif
      if (ng .ne. 0) then
        write(stdo,ftox)' sudmtu:  RMS change in dmats'//
     .  ' from symmetrization',ftof(xx)
        if (xx .gt. .01d0) write(stdo,*)'(warning) RMS change unexpectely large'
        call daxpy ( ivsiz * 2 , - 1d0 , dmatu , 1 , dmwk_zv , 1 )
        if(ipr>=60) write(stdo,*)' change in dmat wrought by symmetrization'
        call praldm ( 0 , 60 , 60 , 0 , nbas , nsp , lmaxu , lldau , 
     .  sspec , ssite , ' ' , dmwk_zv )
      endif
C     Print dmats in specified harmonics
      dmwk_zv=dmatu
      if (havesh .ne. idvsh) then
        call rotycs ( 2 * idvsh - 1 , dmwk_zv , nbas , nsp , lmaxu, sspec , ssite , lldau )
      endif
      write(stdo,*)
      call praldm(0,30,30,idvsh,nbas,nsp,lmaxu,lldau,sspec , ssite , ' Symmetrized dmats' , dmwk_zv )
C     Print dmats in complementary harmonics
      i = 1-idvsh
      call rotycs(2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , sspec , ssite , lldau )
      write(stdo,*)' '
      call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,sspec,ssite, ' Symmetrized dmats' , dmwk_zv )
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
        idu=sspec(is)%idu
        uh=sspec(is)%uh
        jh=sspec(is)%jh
        spid=sspec(is)%name
        i = min(lmxa,3)
        write(stdo,ftox)'Species '//spid//'mode',idu(1:i+1),'U',ftof(uh(1:i+1),2),'J',ftof(jh(1:i+1),2)
        do  22  l = 0, i
          if (idu(l+1) .ne. 0) then
            iblu = iblu+1
c            if(.not.mmtargetx) uhx=uh(l+1)
            uhx=uh(l+1)
            call ldau(idu(l+1),l,iblu,uhx,jh(l+1),dmatu,nsp,lmaxu,vorb,eorb)
          endif
   22   continue
   20 continue
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' Unsymmetrized vorb',vorb)
!==>     At this point, dmatu and vorb are in spherical harmonics
C ... Symmetrize vorb to check (symdmu requires real harmonics)
      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
      call symdmu (nlibu, vorb, nbas , nsp , lmaxu , sspec , ssite , ng , g , istab , lldau , xx )
C     EITHER: vorb =>  spherical harmonics OR dmatu => real harmonics
      if (idvsh .eq. 1) then
        call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 1
      else
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      if (ng .ne. 0) then
        write(stdo,ftox)' sudmtu:  RMS change in vorb from symmetrization = ',ftof(xx)
        if (xx .gt. .01d0) write(stdo,*)'          (warning) RMS change unexpectely large'
      endif
C     Print vorb in specified harmonics
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,' Symmetrized vorb',vorb)
      i = 1-idvsh
      dmwk_zv=vorb
      call rotycs( 2 * i - 1 , dmwk_zv , nbas , nsp , lmaxu , sspec  , ssite , lldau )
      write(stdo,*) !vorb in complementary harmonics
      call praldm(0,30,30, i , nbas , nsp , lmaxu , lldau , sspec , ssite , ' Vorb' , dmwk_zv )
      eorb = 0d0
C --- Error exit ---
      return
   99 continue
      write(str,"('bad occnum file, site l= ',i5,i5)")ib,l
      call rx(str)
      end subroutine sudmtu
      
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
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
          idu=sspec(is)%idu
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
