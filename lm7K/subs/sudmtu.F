Cdelw1 var ['odmwk', '0', 'odmwk', '-ivsiz', 'complex(8)', 'zv_w_', 'sudmtu']
Cgetarg...info...           structure ['sudmtu', 'site', 'ssite']
Cgetarg...info...           structure ['sudmtu', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['sudmtu', 'm_struc_def']
      subroutine sudmtu(nbas,nsp,nlibu,lmaxu,ssite,sspec,idvsh,lldau,
     .  ng,g,istab,dmatu,vorb)

      use m_struc_def  !Cgetarg

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
C     implicit none
C ... Passed parameters
      integer nbas,nsp,nlibu,lmaxu,ng,idvsh
      integer lldau(nbas),istab(nbas,ng)
Cgetarg       double precision ssite(1),sspec(1),g(9,ng)
      real(8):: g(9,ng)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C ... Local parameters
      logical rdstrn,parstr
      integer i,isp,ib,l,lmxa,m,m2,foccn,fopna,rdm,havesh,ivsiz,ipr
C     integer nglob,stdo
      double precision nocc(-3:3,2),iv(7)
      integer idu(4),igetss,is,idmat,fxst,iblu,nlm,nlmu,a2vec
      double precision uh(4),jh(4),eorb,tmp(2,7,7),xx
      character str*80,spid*8
Cdelw1       integer odmwk
      complex(8) ,allocatable :: zv_w_odmwk(:)

C ... MPI
      integer procid,master,mpipid
      logical mlog
C ... Heap
      integer w(1)
      common /w/ w

C     stdo = nglob('stdo')
      call rxx(nsp.ne.2,'LDA+U must be spin-polarized!')
      procid = mpipid(1)
      master = 0
      call getpr(ipr)
C     ipr = 55

C --- Read in dmatu if file  dmats.ext  exists ---
      if (procid .eq. master) then
        if (fxst('dmats') .eq. 1) then
          idmat = fopna('dmats',-1,1)
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
Cgetarg           is = igetss('site spec',ib,ssite)
              is = int(ssite(ib)%spec)

Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idu(4) , integer(8) :: sspec%idu 4
Cgetarg           call upack('spec lmxa idu',sspec,is,lmxa,idu,0,0)

              lmxa=sspec(is)%lmxa
              i_copy_size=size(sspec(is)%idu)
              call i8icopy(i_copy_size,sspec(is)%idu,1,idu,1)

              do l = 0, min(lmxa,3)
                if (idu(l+1) .ne. 0) then
                  iblu = iblu+1
                  nlm = 2*l+1
                  nlmu = 2*lmaxu+1
                  do  isp = 1, 2
                    if (rdm(idmat,40,2*nlm**2,' ',tmp,nlm,nlm) .ne. 2)
     .          call rxi('sudmtu failed to read dmats for site',ib)
C              print *, 'dmat for ib=',ib,' spin',isp
C              call zprm('sudmtu, dmatu',2,tmp,nlm,nlm,nlm)
C              call dmscop(dmatu(-l,-l,isp,iblu),
C     .          nlmu,tmp,nlm,1,nlm,1,nlm,1,1,1d0)
                    call zmscop(0,nlm,nlm,nlm,nlmu,0,0,0,0,tmp,
     .          dmatu(-l,-l,isp,iblu))
                  enddo
                endif
              enddo
            endif
          enddo

C --- Otherwise read in occnumbers and construct diagonal dmatu ---
        else
          call info0(20,1,0,
     .  ' sudmtu:  initial (diagonal) density-matrix from occ numbers')
          foccn = fopna('occnum',-1,1)
          rewind foccn
          havesh = 1
   12     if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
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
Cgetarg           is = igetss('site spec',ib,ssite)
              is = int(ssite(ib)%spec)

Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: idu(4) , integer(8) :: sspec%idu 4
Cgetarg           call upack('spec lmxa idu',sspec,is,lmxa,idu,0,0)

              lmxa=sspec(is)%lmxa
              i_copy_size=size(sspec(is)%idu)
              call i8icopy(i_copy_size,sspec(is)%idu,1,idu,1)

              do l = 0,min(lmxa,3)
                if (idu(l+1) .ne. 0) then
                  iblu = iblu+1
                  do  isp = 1, 2
   11               continue
                    if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
C             Skip comment lines
                    if (str(1:1) .eq. '#') goto 11
                    i = 0
                    m = a2vec(str,len(str),i,4,', ',2,3,2*l+1,iv,nocc(-l,isp))
                    if (m .lt. 0) goto 99
                  enddo
                  call info8(20,0,0,' occ numbers, site %i l=%i:  '//
     .        '%n:1d (spin 1)  %n:1d (spin 2)',
     .        ib,l,2*l+1,nocc(-l,1),2*l+1,nocc(-l,2),0,0)
                  do  isp = 1, 2
                    do  m = -l, l
                      do  m2 = -l, l
                        dmatu(m,m2,isp,iblu) = dcmplx(0d0,0d0)
                      enddo
                      dmatu(m,m,isp,iblu) = dcmplx(nocc(m,isp),0d0)
                    enddo
                  enddo
                endif
              enddo
            endif
          enddo
          call fclose(foccn)
        endif
      endif
C ... Initial printout
      call praldm(0,51,51,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .  ' dmats read from disk',dmatu)

      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call mpibc1(dmatu,2*ivsiz,4,mlog,'sudmtu','dmatu')
      call mpibc1(havesh,1,2,.false.,' ',' ')

C ... Density matrix in real or spherical harmonics (fixed by idvsh)
      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = idvsh
      endif

C ... Symmetrize dmatu (symdmu requires real harmonics)
Cdelw1       call defcc(odmwk,-ivsiz)
      allocate(zv_w_odmwk(abs(-ivsiz)))
      if (-ivsiz<0) zv_w_odmwk(:)=0.0d0

      if (havesh .eq. 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 0
      endif
Cdelw1       call symdmu(dmatu,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,ng,g,istab,
Cdelw1      .  lldau,xx)
      call symdmu ( dmatu , zv_w_odmwk , nbas , nsp , lmaxu , sspec 
     . , ssite , ng , g , istab , lldau , xx )

      if (havesh .ne. idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
Cdelw1         call rotycs(2*idvsh-1,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs ( 2 * idvsh - 1 , zv_w_odmwk , nbas , nsp , lmaxu 
     .   , sspec , ssite , lldau )

        havesh = idvsh
      endif
      if (ng .ne. 0) then
        call info2(30,0,0,' sudmtu:  RMS change in dmats'//
     .  ' from symmetrization = %,6d',xx,0)
        if (xx .gt. .01d0) call info0(30,0,0,
     .  '          (warning) RMS change unexpectely large')
Cdelw1       call daxpy(ivsiz*2,-1d0,dmatu,1,w(odmwk),1)
        call daxpy ( ivsiz * 2 , - 1d0 , dmatu , 1 , zv_w_odmwk , 1 )

        call info0(60,0,0,' change in dmat wrought by symmetrization')
Cdelw1       call praldm(0,60,60,0,nbas,nsp,lmaxu,lldau,sspec,ssite,' ',
Cdelw1      .  w(odmwk))
        call praldm ( 0 , 60 , 60 , 0 , nbas , nsp , lmaxu , lldau , 
     . sspec , ssite , ' ' , zv_w_odmwk )

Cdelw1 rlse name= odmwk old_list= odmwk
Cdelw1 rlse name= odmwk new_list= (None)
Cdelw1       call rlse(odmwk)
        if (allocated(zv_w_odmwk)) deallocate(zv_w_odmwk)

      endif

C     Print dmats in specified harmonics
Cdelw1       call defcc(odmwk,-ivsiz)
      allocate(zv_w_odmwk(abs(-ivsiz)))
      if (-ivsiz<0) zv_w_odmwk(:)=0.0d0

Cdelw1       call dcopy(ivsiz*2,dmatu,1,w(odmwk),1)
      call dcopy ( ivsiz * 2 , dmatu , 1 , zv_w_odmwk , 1 )

      if (havesh .ne. idvsh) then
Cdelw1         call rotycs(2*idvsh-1,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs ( 2 * idvsh - 1 , zv_w_odmwk , nbas , nsp , lmaxu 
     .   , sspec , ssite , lldau )

      endif
      call info0(30,0,0,' ')
Cdelw1       call praldm(0,30,30,idvsh,nbas,nsp,lmaxu,lldau,sspec,ssite,
Cdelw1      .  ' Symmetrized dmats',w(odmwk))
      call praldm ( 0 , 30 , 30 , idvsh , nbas , nsp , lmaxu , lldau 
     . , sspec , ssite , ' Symmetrized dmats' , zv_w_odmwk )


C     Print dmats in complementary harmonics
      if (ipr.ge.45 .or. ipr.ge.40 .and. idvsh.eq.0) then
        i = 1-idvsh
Cdelw1         call rotycs(2*i-1,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs ( 2 * i - 1 , zv_w_odmwk , nbas , nsp , lmaxu , sspec 
     .   , ssite , lldau )

        call info0(30,0,0,' ')
Cdelw1         call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,sspec,ssite,
Cdelw1      .    ' Symmetrized dmats',w(odmwk))
        call praldm ( 0 , 30 , 30 , i , nbas , nsp , lmaxu , lldau , 
     .   sspec , ssite , ' Symmetrized dmats' , zv_w_odmwk )

      endif
Cdelw1 rlse name= odmwk old_list= odmwk
Cdelw1 rlse name= odmwk new_list= (None)
Cdelw1       call rlse(odmwk)
      if (allocated(zv_w_odmwk)) deallocate(zv_w_odmwk)


C ... Make Vorb (ldau requires spherical harmonics)
      if (havesh .ne. 1) then
        call rotycs(1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
        havesh = 1
      endif
      iblu = 0
      do  20  ib = 1, nbas
        if (lldau(ib) .eq. 0) goto 20
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

Cgetarg         call spacks(0,'spec name',sspec,spid,is,is)
        do i_spacks=is,is
          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
        enddo

        i = min(lmxa,3)
        call info8(30,1,0,' Species '//spid//
     .    '%a: mode=%n:1i    U=%n:1d    J=%n:1d',
     .    i+1,idu,i+1,uh,i+1,jh,0,0)
        do  22  l = 0, i
          if (idu(l+1) .ne. 0) then
            iblu = iblu+1
            call ldau(idu(l+1),l,iblu,uh(l+1),jh(l+1),dmatu,nsp,
     .                lmaxu,vorb,eorb)
C            call prdmts(0,30,30,0,' sudmtu:  Vorb in spherical harm',ib,
C     .        l,lmaxu,iblu,vorb,nsp,1)

          endif
   22   continue
   20 continue
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,sspec,ssite,
     .  ' Unsymmetrized vorb',vorb)
C     At this point, dmatu and vorb are in spherical harmonics

C ... Symmetrize vorb to check (symdmu requires real harmonics)
Cdelw1       call defcc(odmwk,-ivsiz)
      allocate(zv_w_odmwk(abs(-ivsiz)))
      if (-ivsiz<0) zv_w_odmwk(:)=0.0d0

      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
Cdelw1       call symdmu(vorb,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,ng,g,istab,
Cdelw1      .  lldau,xx)
      call symdmu ( vorb , zv_w_odmwk , nbas , nsp , lmaxu , sspec 
     . , ssite , ng , g , istab , lldau , xx )

Cdelw1 rlse name= odmwk old_list= odmwk
Cdelw1 rlse name= odmwk new_list= (None)
Cdelw1       call rlse(odmwk)
      if (allocated(zv_w_odmwk)) deallocate(zv_w_odmwk)


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
     .  ' Symmetrized vorb',vorb)

C     Print vorb in complementary harmonics
      if (ipr.ge.45 .or. ipr.ge.40 .and. idvsh.eq.0) then
        i = 1-idvsh
Cdelw1         call defcc(odmwk,-ivsiz)
        allocate(zv_w_odmwk(abs(-ivsiz)))
        if (-ivsiz<0) zv_w_odmwk(:)=0.0d0

C        call dcopy(ivsiz*2,dmatu,1,w(odmwk),1)
C        call rotycs(2*i-1,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,lldau)
C        call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,sspec,ssite,
C     .    ' dmats',w(odmwk))
Cdelw1         call dcopy(ivsiz*2,vorb,1,w(odmwk),1)
        call dcopy ( ivsiz * 2 , vorb , 1 , zv_w_odmwk , 1 )

Cdelw1         call rotycs(2*i-1,w(odmwk),nbas,nsp,lmaxu,sspec,ssite,lldau)
        call rotycs ( 2 * i - 1 , zv_w_odmwk , nbas , nsp , lmaxu , sspec 
     .   , ssite , lldau )

        call info0(30,0,0,' ')
Cdelw1         call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,sspec,ssite,
Cdelw1      .    ' Vorb',w(odmwk))
        call praldm ( 0 , 30 , 30 , i , nbas , nsp , lmaxu , lldau , 
     .   sspec , ssite , ' Vorb' , zv_w_odmwk )

Cdelw1 rlse name= odmwk old_list= odmwk
Cdelw1 rlse name= odmwk new_list= (None)
Cdelw1         call rlse(odmwk)
        if (allocated(zv_w_odmwk)) deallocate(zv_w_odmwk)

      endif

      eorb = 0

C     call prmx('vorb',vorb(1,1,1,1),2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

C --- Error exit ---
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return
   99 continue
      call awrit2('bad occnum file, site %i, l=%i',str,len(str),0,ib,l)
      call rx(str)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine sudmtu


