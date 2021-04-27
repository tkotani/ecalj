      subroutine chkdmu(nbas,nsp,nlibu,lmaxu,ssite,sspec,idvsh, !sham,
     .dmatu,dmatuo,vorb,tolu,umix,lldau,ng,g,istab)
      use m_struc_def,only: s_site,s_spec  !Cgetarg
      use m_lmfinit,only: ham_ehk,ham_eterms!,procid,master
      use m_lmfinit,only: globalvariables
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
c      integer:: ierr
      include "mpif.h"
C ... Passed parameters
      integer nbas,nsp,nlibu,lmaxu,ng,idvsh,i_copy_size
      integer lldau(nbas),istab(nbas,ng)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
c      type(s_ham)::sham

      double precision tolu,umix,g(9,ng)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex dmatuo(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double precision eorb
C ... Local parameters
      integer l,idu(4),lmxa,ib,is,iblu,igetss,idmat,fopna,ivsiz
      integer nglob,stdl,iprint,ipl,havesh
      double precision ddmat,uh(4),jh(4),eorbi,eterms(20),eks,ddot,xx
      equivalence (eterms(2),eks)
      complex(8) ,allocatable :: dmwk_zv(:)
      logical:: newaftest=.true.,fexist,mmtargetx
      real(8),allocatable:: uhall(:,:)
      real(8):: mmsite(nbas),uhxx,mmhist(10000),uhhist(10000),mmtarget,uhdiff
      real(8),save::uhxnew,uhx,alpha
      integer:: nn,ibas,ifx,key,i,nit
      integer,save::ncount=0
C ... MPI
      integer:: procid,master,mpipid,mlog=-999,ierr
      complex(8):: dmatuav(-lmaxu:lmaxu,-lmaxu:lmaxu)
      real(8):: sss(1),sigin
      real(8):: fac
C     logical mlog,cmdopt
C ... Heap
c      integer w(1)
c      common /w/ w

C     print *, 'dmatu',dmatu(-lmaxu,-lmaxu,1,1)
C     call rx0('done')
C     call prmx('dmatu',dmatu,2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

      if (nlibu .eq. 0) return
C     stdo = nglob('stdo')
      havesh = idvsh
Changenglob      stdl = nglob('stdl')
      stdl = globalvariables%stdl
      ipl = 1
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call info0(20,1,0,' chkdmu:  '//
     .'check LDA+U density-matrix for convergence and update ...')
C ... MPI
      procid = mpipid(1)
      master = 0
C     mlog = cmdopt('--mlog',6,0,strn)

C --- Symmetrize output dmatu (req. real harmonics); compare diff ---
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .' Unsymmetrized output dmats',dmatu)
      if (havesh .eq. 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
      allocate(dmwk_zv(abs(-ivsiz)))
      if (-ivsiz<0) dmwk_zv(:)=0.0d0

      call symdmu ( dmatu , dmwk_zv , nbas , nsp , lmaxu , sspec 
     ., ssite , ng , g , istab , lldau , xx )

      if (allocated(dmwk_zv)) deallocate(dmwk_zv)

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

      i_copy_size=size(ham_eterms)
      call dcopy(i_copy_size,ham_eterms,1,eterms,1)

      eks = eks + eorb

      ham_ehk=eks
      i_copy_size=size(ham_eterms)
      call dcopy(i_copy_size,eterms,1,ham_eterms,1)

      call info5(20,0,0,'%9peks = %,6;6d  '//
     .'e[U] = %,6;6d  Etot(LDA+U) = %,6;6d',eks-eorb,eorb,eks,0,0)
      if (mpipid(1) .eq. 0 .and. ipl.gt.0)
     .write (stdl,720) eks-eorb,eorb,eks
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

       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
!     ! experimental block to keep magnetic moment for AF.
       inquire(file='mmtarget.aftest',exist=mmtargetx)
       if(mmtargetx) then
          if(procid==master) then
             open(newunit=ifx,file='mmtarget.aftest')
             alpha=1.0d0
             read(ifx,*)mmtarget !,alpha
             close(ifx)
             allocate(uhall(0:3,nbas))
             do ib = 1, nbas
                is = int(ssite(ib)%spec)
                nn = size(sspec(is)%uh)
                if(nn>4) call rx('chkdmu:nn>4')
                uhall(0:nn-1,ib) = sspec(is)%uh(1:nn)
             enddo
             mmhist=0d0
             nit=0
             open(newunit=ifx,file='uhval.aftest')
             do
                read(ifx,*,end=1112,err=1112) uhx
                read(ifx,*,end=1112,err=1112) (mmsite(ibas),ibas=1,nbas)
                nit=nit+1
                mmhist(nit)=(mmsite(1)-mmsite(2))/2d0
                uhhist(nit)=uhx
                write(6,"('uhval: ',i5,f10.6,2x,12f10.6)") nit,uhx, (mmsite(ibas),ibas=1,nbas)
             enddo
 1112        continue
             close(ifx)
!     ! Generate new uhx based on the  history of uhx mmsites for given mm
!     ! test uh          
             uhx= uhx + alpha*(mmhist(nit)-mmtarget)**2 - 2d0*(mmhist(nit)-mmtarget)
             write(6,"('mmhist0: uh mm uhx',i5, 3d13.4)") nit,uhhist(nit),mmhist(nit),uhx
             sss(1)=uhx
             call mixuh(sss)
             if(mmhist(nit)<0.or.ncount==0) then
                open(newunit=ifx,file="mixuh.aftest")
                close(ifx,status='delete')
                ncount=1
                uhx = uhx 
             endif
             write(6,"('mmhist: uh mm uhx',i5, 3d13.4)") nit,uhhist(nit),mmhist(nit),uhx
          endif
!! broadcast uhx       
          call MPI_Bcast(uhx, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
          do iblu =1,2
             fac=1d0
             if(iblu==2) fac=-1d0
             dmatuav(:,:) =  .5d0*dmatu(:,:,1,iblu) +.5d0*dmatu(:,:,2,iblu)
             dmatu(:,:,1,iblu) =  fac*uhx*dmatuav(:,:)
             dmatu(:,:,2,iblu) = -fac*uhx*dmatuav(:,:)
          enddo
       endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       
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
                call ldau(idu(l+1),l,iblu,uh(l+1),jh(l+1),dmatu,nsp, !uhxx,jh(l+1),dmatu,nsp, !
     .          lmaxu,vorb,eorb)
                call poppr
              endif
            enddo
          endif
        enddo
C       At this point, dmatu and vorb are in spherical harmonics
        
!!-------------------------------------------------------
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!! experimental block to keep magnetic moment for AF.       
       if(procid==master.and.newaftest) then
       open(newunit=ifx,file='uhval.aftest',position='append')
       write(ifx,"(d23.15,2x,i5,d23.15)",advance='no') uhx !, uhxnew
       close(ifx)
       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
       
C   ... Symmetrize vorb to check (symdmu requires real harmonics)
        allocate(dmwk_zv(abs(-ivsiz)))
        if (-ivsiz<0) dmwk_zv(:)=0.0d0

        call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        call symdmu ( vorb , dmwk_zv , nbas , nsp , lmaxu , sspec 
     .  , ssite , ng , g , istab , lldau , xx )

        call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
        if (allocated(dmwk_zv)) deallocate(dmwk_zv)

C       At this point, dmatu and vorb are in spherical harmonics

C   ... Printout
        call info2(20,0,0,'         RMS change in vorb '//
     .  'from symmetrization = %,6d',xx,0)
        if (xx .gt. .0001d0) call info0(30,0,0,
     .  '         (warning) RMS change unexpectely large')
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .  ' Mixed dmats',dmatu)
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .  ' New vorb',vorb)

C   ... Write dmatu to file
        if (procid .eq. master) then
          idmat = fopna('dmats',-1,0)
          rewind idmat
          call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .    ' mixed dmats',dmatu)
          call fclose(idmat)
        endif

C   ... Exit with dmatu, vorb in real harmonics, depending on idvsh
        if (idvsh .eq. 0) then
          call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
          call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .    ' Mixed dmats',dmatu)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .    ' New vorb',vorb)
        endif

      endif

      end subroutine chkdmu



C----------------------------------------------------------------------
      subroutine mixuh(sss)
c  subroutine pqmixa(nda,nmix,mmix,mxsav,beta,rms2,a,tj)
C- Mixing routine for sigma. Modified from pqmixa in subs/pqmix.f
C- Anderson mixing of a vector
Ci  mmix: number of iterates available to mix
Cio nmix: nmix > 0: number of iter to try and mix
Ci        nmix < 0: use mmix instead of nmix.
Co  nmix: (abs)  number of iter actually mixed.
Co        (sign) <0, intended that caller update nmix for next call.
C  MvS Feb 04 use sigin as input sigma if available (lsigin=T)
C             Add mixnit as parameter
      implicit none
      logical lsigin
      integer,parameter:: nda=1
      integer:: nmix,mmix
      integer(4),parameter:: mxsav=10
      double precision rms2,tj(mxsav),beta
      integer im,imix,jmix,amix,onorm,okpvt,oa
      integer iprintxx,ifi,nitr,ndaf
      real(8)::sss(nda),sigin(nda)
      real(8):: tjmax
      real(8),allocatable::norm(:),a(:,:,:)
      integer(4),allocatable:: kpvt(:)
      integer(4)::ret
      character*8 fff
      logical fexist
      real(8):: acc
      integer(4):: ido,ifile_handle
      iprintxx = 30
      beta=.3d0
      allocate ( a(nda,0:mxsav+1,2) )
      fff="mixuh.aftest"
      INQUIRE (FILE =fff, EXIST = fexist)
      if(fexist)      write(6,*)'... reading file mixsigma'
      if(.not.fexist) write(6,*)'... No file mixsigma'
      ifi=ifile_handle()
      open(ifi,file=fff,form='unformatted')
      if(fexist) then
        read(ifi,err=903,end=903) nitr,ndaf
        if (ndaf .ne. nda) goto 903
        read(ifi,err=903,end=903) a
        goto 902
      endif
      goto 901
  903 continue
      print 368
  368 format(5x,'(warning) file mismatch ... mixing file not read')
  901 continue
      nitr = 0
 902  continue
      a(:,0,1) = sss      !output
C     if input sigma available, use it instead of file a(:,0,2)
      if (lsigin) then
        write(6,*)'... using input sigma read from sigm file'
        a(:,0,2) = sigin  !input
      endif
      write(6,*)'sum sss=',sum(abs(sss))
      imix=9
      mmix = min(max(nitr-1,0),imix)
      if (mmix > mxsav) mmix = mxsav
C     this information already printed out by amix
C     write(6,*)'mixing parameters for amix are fixed in mixsigma'
C     write(6,*)'   beta       =', beta
C     write(6,*)'   tjmax      =', tjmax
C     write(6,*)'   mmix mxsav =', mmix,mxsav
c     call getkeyvalue("GWinput","mixtj",acc,default=0d0,status=ret)
      acc=0d0
      if(acc/=0d0) then
        write(6,*)' readin mixtj from GWinput: mixtj=',acc
        tjmax=abs(acc)+1d-3
        if(mmix==1) then
          tj(1)=1d0
        else
          tj(1)= acc
          tj(2)= 1-acc
          mmix=2
        endif
        ido=2
      else
        tjmax=5d0
        ido=0
      endif
      allocate(norm(mxsav**2),kpvt(mxsav))
      imix = amix(nda,mmix,mxsav,ido,dabs(beta),iprintxx,tjmax,
     .  norm,kpvt,a,tj,rms2)
      deallocate(norm, kpvt)
C ... Restore PQ array, updating new x
c      call dpscop(a,w(oa),nda,1+nda*(mxsav+2),1+nda*(mxsav+2),1d0)
c      call dcopy(nda*(mxsav+2)*2,w(oa),1,a,1)
c ...
      sss = a(:,0,2)
      rewind(ifi)
      write(ifi) nitr+1,nda
      write(ifi) a
      close(ifi)
      end
      
