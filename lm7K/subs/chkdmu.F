Cdelw1 var ['odmwk', '0', 'odmwk', '-ivsiz', 'complex(8)', 'zv_w_', 'chkdmu']
Cgetarg...info...           structure ['chkdmu', 'site', 'ssite']
Cgetarg...info...           structure ['chkdmu', 'spec', 'sspec']
Cgetarg...info...           structure ['chkdmu', 'ham', 'sham']
Cgetarg...info...           use_to_add ['chkdmu', 'm_struc_def']
      subroutine chkdmu(nbas,nsp,nlibu,lmaxu,ssite,sspec,sham,idvsh,
     .  dmatu,dmatuo,vorb,tolu,umix,lldau,ng,g,istab)
       
       use m_struc_def  !Cgetarg

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
C     implicit none
C ... Passed parameters
      integer nbas,nsp,nlibu,lmaxu,ng,idvsh
      integer lldau(nbas),istab(nbas,ng)
Cgetarg       double precision ssite(1),sspec(1),sham(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_ham)::sham

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
Cdelw1       integer odmwk
       complex(8) ,allocatable :: zv_w_odmwk(:)

C ... MPI
      integer procid,master,mpipid
C     logical mlog,cmdopt
C ... Heap
c      integer w(1)
c      common /w/ w

C     print *, 'dmatu',dmatu(-lmaxu,-lmaxu,1,1)
C     call rx0('done')
C     call prmx('dmatu',dmatu,2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (nlibu .eq. 0) return
C     stdo = nglob('stdo')
      havesh = idvsh
      stdl = nglob('stdl')
      ipl = 1
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call info0(20,1,0,' chkdmu:  '//
     .  'check LDA+U density-matrix for convergence and update ...')
C ... MPI
      procid = mpipid(1)
      master = 0
C     mlog = cmdopt('--mlog',6,0,strn)

C --- Symmetrize output dmatu (req. real harmonics); compare diff ---
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .  ' Unsymmetrized output dmats',dmatu)
      if (havesh .eq. 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
Cdelw1       call defcc(odmwk,-ivsiz)
       allocate(zv_w_odmwk(abs(-ivsiz)))
       if (-ivsiz<0) zv_w_odmwk(:)=0.0d0

Cdelw1       call symdmu(dmatu,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,ng,g,istab,
Cdelw1      .  lldau,xx)
       call symdmu ( dmatu , zv_w_odmwk , nbas , nsp , lmaxu , sspec 
     . , ssite , ng , g , istab , lldau , xx ) 

Cdelw1 rlse name= odmwk old_list= odmwk 
Cdelw1 rlse name= odmwk new_list= (None)
Cdelw1       call rlse(odmwk)
       if (allocated(zv_w_odmwk)) deallocate(zv_w_odmwk)

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
     .  ' from symmetrization = %,6d',xx,0)

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
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idu(4) , integer(8) :: sspec%idu 4
Cgetarg...info...           real(8) :: uh(4) , real(8) :: sspec%uh 4
Cgetarg...info...           real(8) :: jh(4) , real(8) :: sspec%jh 4
Cgetarg         call upack('spec lmxa idu uh jh',sspec,is,lmxa,idu,uh,jh)
         
         lmxa=sspec(is)%lmxa
         i_copy_size=size(sspec(is)%idu) 
         call i8icopy(i_copy_size,sspec(is)%idu,1,idu,1) 
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
Cgetarg...info...           real(8) :: eterms(20) , real(8) :: sham%eterms 20
Cgetarg       call upack1('ham eterms',sham,eterms)
       
       i_copy_size=size(sham%eterms) 
       call dcopy(i_copy_size,sham%eterms,1,eterms,1) 

      eks = eks + eorb
Cgetarg...info...           real(8) :: eks , real(8) :: sham%ehk 1
Cgetarg...info...           real(8) :: eterms(20) , real(8) :: sham%eterms 20
Cgetarg       call pack5('ham ehk eterms',sham,eks,eterms,0,0,0)
       
       sham%ehk=eks 
       i_copy_size=size(sham%eterms) 
       call dcopy(i_copy_size,eterms,1,sham%eterms,1) 

      call info5(20,0,0,'%9peks = %,6;6d  '//
     .  'e[U] = %,6;6d  Etot(LDA+U) = %,6;6d',eks-eorb,eorb,eks,0,0)
      if (mpipid(1) .eq. 0 .and. ipl.gt.0)
     .  write (stdl,720) eks-eorb,eorb,eks
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
     .    '  RMS diff (%;3g) < tolu (%;3g)',ddmat,tolu,0,0,0)
C       Restore dmatuo to dmatu
        call dcopy(2*ivsiz,dmatuo,1,dmatu,1)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        return

C --- Case not self-consistent ---
      else
        call info0(20,0,0,'%9pLDA+U update density matrix ...')
        call info5(20,0,0,'%9fRMS diff in dens mat(%;3g) > tolu (%;3g)'
     .    //' Linear mix with beta=%;3g',ddmat,tolu,umix,0,0)

C   ... Make new dmatu by mixing    new*umix + old*(1-umix)
        call dscal(2*ivsiz,umix,dmatu,1)
        call daxpy(2*ivsiz,1-umix,dmatuo,1,dmatu,1)

C   --- Make Vorb from mixed dmatu ---
C       This block requires dmatu to be in spherical harmonics
        if (havesh .ne. 1) then
          call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
          havesh = 1
        endif
        iblu = 0
        do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
Cgetarg           is = igetss('site spec',ib,ssite)
           is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idu(4) , integer(8) :: sspec%idu 4
Cgetarg...info...           real(8) :: uh(4) , real(8) :: sspec%uh 4
Cgetarg...info...           real(8) :: jh(4) , real(8) :: sspec%jh 4
Cgetarg           call upack('spec lmxa idu uh jh',sspec,is,lmxa,idu,uh,jh)
           
           lmxa=sspec(is)%lmxa
           i_copy_size=size(sspec(is)%idu) 
           call i8icopy(i_copy_size,sspec(is)%idu,1,idu,1) 
           i_copy_size=size(sspec(is)%uh) 
           call dcopy(i_copy_size,sspec(is)%uh,1,uh,1) 
           i_copy_size=size(sspec(is)%jh) 
           call dcopy(i_copy_size,sspec(is)%jh,1,jh,1) 

          do  l = 0, min(lmxa,3)
            if (idu(l+1) .ne. 0) then
              iblu = iblu+1
              call pshpr(iprint()-20)
              call ldau(idu(l+1),l,iblu,uh(l+1),jh(l+1),dmatu,nsp,
     .          lmaxu,vorb,eorb)
              call poppr
            endif
          enddo
        endif
        enddo
C       At this point, dmatu and vorb are in spherical harmonics

C   ... Symmetrize vorb to check (symdmu requires real harmonics)
Cdelw1         call defcc(odmwk,-ivsiz)
         allocate(zv_w_odmwk(abs(-ivsiz)))
         if (-ivsiz<0) zv_w_odmwk(:)=0.0d0

        call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
Cdelw1         call symdmu(vorb,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,ng,g,istab,
Cdelw1      .    lldau,xx)
         call symdmu ( vorb , zv_w_odmwk , nbas , nsp , lmaxu , sspec 
     .   , ssite , ng , g , istab , lldau , xx ) 

        call rotycs(1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
Cdelw1 rlse name= odmwk old_list= odmwk 
Cdelw1 rlse name= odmwk new_list= (None)
Cdelw1         call rlse(odmwk)
         if (allocated(zv_w_odmwk)) deallocate(zv_w_odmwk)

C       At this point, dmatu and vorb are in spherical harmonics

C   ... Printout
        call info2(20,0,0,'         RMS change in vorb '//
     .    'from symmetrization = %,6d',xx,0)
        if (xx .gt. .0001d0) call info0(30,0,0,
     .    '         (warning) RMS change unexpectely large')
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .    ' Mixed dmats',dmatu)
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .    ' New vorb',vorb)

C   ... Write dmatu to file
        if (procid .eq. master) then
          idmat = fopna('dmats',-1,0)
          rewind idmat
          call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .      ' mixed dmats',dmatu)
          call fclose(idmat)
        endif

C   ... Exit with dmatu, vorb in real harmonics, depending on idvsh
        if (idvsh .eq. 0) then
          call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
          call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .      ' Mixed dmats',dmatu)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .      ' New vorb',vorb)
        endif

      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine chkdmu 


