      module m_subzi ! Obtain weight wtkb(ib,isp,iq) for brillowine zone integation
Co Outputs
Co   lwtkb :0 weights are neither required nor available a priori
Co         :1 weights are available
Co         :-1 weights are not yet
Co   lswtk :Flags whether to make 'spin weights' swtk
Co         :-2 do not make spin weights
Co         : 1 spin weights array allocated; make them
      integer,protected:: lwtkb=-1
      real(8),allocatable,protected :: rv_a_owtkb(:,:,:) ! wtkb : integration weights. it might be from wkp.*
      integer,protected:: nevmx=0 !for band mode plbnd=1. --->nev=min(nevmx,ndimhx) =0
      integer,protected:: lswtk=0

      contains
      
      subroutine m_subzi_setlwtkb(lwtkbin )
      integer:: lwtkbin
      lwtkb = lwtkbin
      end
      
      subroutine m_subzi_init(lwt) 
      use m_ext,only: sname
      use m_lmfinit, only: nsp,nspc, nevmxin=>bz_nevmx, lmet=>bz_lmet, qbg=>zbak,stdo
      use m_mkqp,only: ntet=> bz_ntet ,bz_nkp
      use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx
      use m_mkpot,only:  qval
      use m_MPItk,only: master_mpi
C- Brillouin-integration setup
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmet  :See Remarks
Ci         :0 assume insulator
Ci         :1 save eigenvectors to disk
Ci         :2 read weights from file wkp.*
Ci         :3 always make two band passes; weights never needed a priori
Ci         :xxx 4 BZ integration with 3-point scheme
Ci   ltet  :T allocate space for tetrahedron weights
Ci   lwt   :F weights are not needed until all bands are obtained
Ci         :T weights are needed a priori (eg output density generated)
Ci   ndham :leading dimension of owtkb
Ci         :Hamiltonian should not exceed this dimension
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Cio Inputs/Outputs
Cio  nevmx :On input, maximum number of eigenvectors to find
Cio         input nevmx<0 => do not generate eigenvectors
Cio         input nevmx=0 => choose a default value
Co   swtk  :memory is allocated for spin weights, nspc=2
!!      
Cr   To integrate the output density in the Brillouin Zone, integration
Cr   weights are needed for each band and qp, but they are not known
Cr   until all the bands are obtained.  The problem is solved in one of
Cr   the following ways:
Cr -------------------------------------
Cr   zval  :total valence charge
Cr     lmet=0 system assumed to be an insulator; weights known a priori
Cr
Cr     lmet=1 eigenvectors are written to disk, in which case the
Cr            integration for the charge density can be deferred until
Cr            all the bands are obtained
Cr
Cr     lmet=2 integration weights are assumed from a prior band pass
Cr
Cr     lmet=3 After m_bandcal_init, we calcutlate wtkb, which is passed to m_bandcal_2nd to accumulate quantities of sum in BZ.
Cr
c lmet=4 is removed now ------------.      
Cu Updates
Cu   09 Jun 07 Setup for spin weights (noncollinear case)
Cu   11 Oct 02 (ATP) MPI
Cu   21 Mar 01 Added printout; argument list changed
C ----------------------------------------------------------------------
      implicit none
      logical ltet,lwt,tdos,cmdopt0,PROCARon,fullmesh
      integer nkp,mpsord,ifile_handle
      double precision zval,ef0,def,esmear
      integer ifi,lerr,iprint,isw,n
      character*11 strni(2)
      integer procid,master,mpipid,nevx0,nq0,nsp0
      data strni /'sampling','tetrahedron'/
      call tcn('m_subzi_init')
!! set nevmx for plbnd=1 Note: initial setting nevmx=0 is returned for band mode.
      if(cmdopt0('--cls').or.cmdopt0('--tdos').or.cmdopt0('--mkprocar'))  nevmx = ndhamx
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') 
      if(cmdopt0('--band').or.fullmesh) return
      nkp=bz_nkp
      if(allocated(rv_a_owtkb)) deallocate(rv_a_owtkb)
      zval=qval-qbg
      ltet = ntet>0
      nevmx=nevmxin !initial contdition except  cmdopt0('--band').or.fullmesh 
      procid = mpipid(1)
      master = 0
      lwtkb = 0
      if ( lmet .gt. 0 ) then
        allocate(rv_a_owtkb(ndham,nsp,nkp))
      endif
      if (nevmx .ge. 0) then
        if (lmet .eq. 2 .or. lmet .eq. 3) then  
          lwtkb = -1
C     ... Attempt to use existing weights
          if (lmet .eq. 2 .and. lwt) then
            if(master_mpi) then
              open(newunit=ifi,file='wkp.'//trim(sname),form='unformatted')
              read(ifi,end=8080,err=8080) nevx0,nq0,nsp0
              if(.not.(ndham*nspc==nevx0 .and. nkp==nq0 .and. nsp==nsp0)) goto 8080 
              read(ifi) rv_a_owtkb
              lwtkb=1
 8080         continue
              close(ifi)
            endif
            call mpibc1(lwtkb,1,2,.false.,'subzi','lwtkb')
            if (lwtkb==1) call mpibc1(rv_a_owtkb, ndham * nsp * nkp , 4 , .false. , 'subzi' , 'wtkb' )
          endif
        endif
      endif
      lswtk = -2
      if (nspc ==2.and.lwtkb .eq. 1) lswtk = 1
      if (nevmx .eq. 0) then
        nevmx = (int(zval) + 1)/2
        if (lmet .ne. 0) nevmx = max(nevmx+nevmx/2,9)
        nevmx = min(nevmx,ndham)
        if (nspc .eq. 2) nevmx = 2*nevmx
        nevmx=nevmx+5 !+5 is for safer setting. At least +1 is required...
      endif
cC ... Printout
c      if (nevmx .ge. 0 .and. iprint() .gt. 30) then
c        if (lmet .gt. 0) then
c          call awrit0('%N subzi: '//strni(isw(ltet)+1)//'%a integration of bands; '//
c     .    strni(isw(lmet.ne.4.and.ltet)+1)//'%a integration of density',' ',80,stdo)
c        else
c          call info(20,0,0,' subzi : nonmetal',0,0)
c        endif
c        write(stdo,'(1x)')
c      endif
      call tcx('m_subzi_init')
      end subroutine m_subzi_init
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine m_subzi_bzintegration(evlall,swtk, eferm,sev,sumqv,vnow)
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq
      use m_lmfinit, only: lso,nsp,ham_scaledsigma,nlibu,lmaxu,bz_w,
     & lmet=>bz_lmet,stdo,nbas,epsovl=>ham_oveps,nspc,bz_n,bz_fsmommethod,qbg=>zbak,fsmom=>bz_fsmom,ndos=>bz_ndos
      use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
      use m_lmfinit,only: ham_ldham , nchan=>pot_nlma, nvl=>pot_nlml
      use m_qplist,only: m_qplist_init, nkp,xdatt,labeli,labele,dqsyml,etolc,etolv,
     &     nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml
      use m_mkpot,only: qval
      use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx,nspx=>ham_nspx 
      use m_ext,only: sname
      implicit none
      intent(in)::                     evlall,swtk
      intent(out)::                                 eferm,sev,sumqv,vnow
      logical:: lfill=.false.,ltet
      logical:: debug,cmdopt0 !goto99,
      integer:: ierr,ifimag,i,ifi,ifile_handle,unlink !,iobzwt
      real(8):: dosrng,evlall(*),sev,sumqv(3,*),eferm,vnow,ef0,bz_ef
      real(8),parameter::    NULLR =-99999
      real(8):: swtk(*)
      call tcn('m_subzi_bzintegration')
      sev=0d0
      ltet = ntet>0
c      qbg = ctrl_zbak(1) !homogenious background charge
      debug    = cmdopt0('--debugbndfp')
      if (master_mpi) ierr=unlink('MagField') !delete
!! --- BZ integration for fermi level, band sum and qp weights ---
      dosrng = 8
      if (bz_n .lt. 0) dosrng = 16
      if(bz_fsmommethod == 1) then !takao dec2010
!     ! vnow june22013 !vnow !(in Ry) contains magnetic field
!     ! For eigenvalus, add  -vnow/2 for isp=1, and +vnow/2 for isp=2.
         call bzwtsf2 ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .        ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval-qbg, !note qval is output
     .        fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .        , dosrng , rv_a_owtkp , evlall ,  lswtk , swtk 
     .        , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill,vnow)
      else
!     ! vnow june22013
         call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .        ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval-qbg , 
     .        fsmom , lmet.ne.0 , ltet , bz_n , ndos , bz_w 
     .        , dosrng , rv_a_owtkp , evlall , lswtk , swtk 
     .        , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill, vnow)
      endif
!     ! june2013 magfield is added
      if(fsmom/=NULLR .and. master_mpi) then
         open(newunit=ifimag,file='MagField',status='unknown')
         write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow
         close(ifimag)
      endif
!     !         Store val charge & magnetic moment in sumqv(1..2)
      sumqv(1,1) = sumqv(1,2)
      sumqv(2,1) = sumqv(2,2)
      if (lmet > 0) then
         if (master_mpi) then
            open(newunit=ifi,file='wkp.'//trim(sname),form='unformatted')
            write(ifi) ndhamx,nkp,nspx
            write(ifi) rv_a_owtkb
            close(ifi)    
         endif
      endif
      call tcx('m_subzi_bzintegration')
      end subroutine
      end module m_subzi
