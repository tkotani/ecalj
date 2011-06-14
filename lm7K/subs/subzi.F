      module m_subzi
      contains
c      subroutine subzi(lmet,ltet,lwt,ndham,nsp,nspc,nkp,zval,nevmx,
c     .  lwtkb,efermi,oevl,owtkb,numq,lswtk,oswtk,ef0,def,mpsord,esmear)
      subroutine subzi(lmet,ltet,lwt,ndham,nsp,nspc,nkp,zval,nevmx,
     .lwtkb,efermi,rv_p_owtkb,numq,lswtk,rv_p_oswtk,ef0,def,mpsord,esmear)
C- Brillouin-integration setup
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmet  :See Remarks
Ci         :0 assume insulator
Ci         :1 save eigenvectors to disk
Ci         :2 read weights from file, if they are needed
Ci         :3 always make two band passes; weights never needed a priori
Ci         :4 BZ integration with 3-point scheme
Ci   ltet  :T allocate space for tetrahedron weights
Ci   lwt   :F weights are not needed until all bands are obtained
Ci         :T weights are needed a priori (eg output density generated)
Ci   ndham :leading dimension of owtkb, oevl, oswtk
Ci         :Hamiltonian should not exceed this dimension
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   zval  :total valence charge
Ci   ef0   :(for printout only) trial fermi level (sampling)
Ci   def   :(for printout only) Fermi level window
Ci   mpsord:(for printout only) Methfessel-Paxton order, sampling integ.
Ci   esmear:(for printout only) gaussian broadening, sampling integ.
Cio Inputs/Outputs
Cio  nevmx :On input, maximum number of eigenvectors to find
Cio         input nevmx<0 => do not generate eigenvectors
Cio         input nevmx=0 => choose a default value
Co Outputs
Co   lwtkb :0 weights are neither required nor available a priori
Co         :1 weights are required a priori, and are read from disk
Co         :-1 weights are required a priori, but were not read
Co   efermi:(lwtkb=1) Fermi level corresonding weights; otherwise
Co         :efermi is set to -99.
Co   oevl  :memory is allocated for energy bands
Co   owtkb :memory is allocated for integration weights, and they
Co         :are possibly read from disk.
Co         :Not allocated if lwtkb=0.
Co   numq  :number of Fermi levels. Usu. 1, but numq=3 if lmet=4
Co   lswtk :Flags whether to make 'spin weights' swtk
Co         :-2 do not make spin weights
Co         : 1 spin weights array allocated; make them
Co   oswtk :memory is allocated for spin weights, nspc=2
Cr   To integrate the output density in the Brillouin Zone, integration
Cr   weights are needed for each band and qp, but they are not known
Cr   until all the bands are obtained.  The problem is solved in one of
Cr   the following ways:
Cr
Cr     lmet=0 system assumed to be an insulator; weights known a priori
Cr
Cr     lmet=1 eigenvectors are written to disk, in which case the
Cr            integration for the charge density can be deferred until
Cr            all the bands are obtained
Cr
Cr     lmet=2 integration weights are assumed from a prior band pass
Cr
Cr     lmet=3 two band passes are made; the first generates only evals
Cr
Cr     lmet=4 information is retained for three distinct Fermi levels.
Cr            After the Fermi level is determined, the density is
Cr            obtained by interpolation of the three points.  (This
Cr            scheme is suitable for sampling only, since in that case
Cr            just the Fermi level is needed to set integration weights.
Cr            When this scheme is used in conjunction with the
Cr            tetrahedron method, the charge density is calculated with
Cr            sampling.
Cu Updates
Cu   09 Jun 07 Setup for spin weights (noncollinear case)
Cu   25 Apr 04 subzi returns efermi=-99, or file value if lwtkb=1.
Cu             Altered argument list
Cu   11 Oct 02 (ATP) MPI
Cu   21 Mar 01 Added printout; argument list changed
Cu   23 Jan 01 set numq=3 for lmet=4 always
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical ltet,lwt
      integer lmet,ndham,nsp,nspc,nkp,lwtkb,lswtk,numq,nevmx,mpsord
c        real(8),pointer :: rv_p_owtkb(:)
ckino change rv_p_owtkb allocatable        real(8),pointer :: rv_p_owtkb(:,:,:)
        real(8),allocatable :: rv_p_owtkb(:,:,:)

ckino change rv_p_oswtk allocatable       real(8),pointer :: rv_p_oswtk(:) 
       real(8),allocatable :: rv_p_oswtk(:) 

      double precision zval,ef0,def,esmear,efermi
C ... Local parameters
      integer ifi,fopna,lerr,iobzwt,iprint,stdo,lgunit,isw,n
      character*11 strni(2)
      integer procid,master,mpipid
C ... Heap
      data strni /'sampling','tetrahedron'/

      procid = mpipid(1)
      master = 0

      n = isign(1,mpsord) * mod(iabs(mpsord),100)
      efermi = -99
      numq = 1
      if (lmet .eq. 4) numq = 3
      lwtkb = 0
       if ( lmet .gt. 0 .and. ( lmet .ne. 4 .or. ltet ) ) then 
       allocate(rv_p_owtkb(ndham,nsp,nkp))
c       allocate(rv_p_owtkb(ndham*nsp*nkp))
       endif 

      if (nevmx .ge. 0) then
        if (lmet .eq. 2 .or. lmet .eq. 3) then
          lwtkb = -1
C     ... Attempt to use existing weights
          if (lmet .eq. 2 .and. lwt) then
            if (procid .eq. master) then
              ifi = fopna('wkp',-1,4)
               if ( nspc .eq. 1 ) lerr = iobzwt ( 0 , ndham , nkp , nsp , efermi 
     .         , rv_p_owtkb , ifi ) 

               if ( nspc .eq. 2 ) lerr = iobzwt ( 0 , ndham * 2 , nkp , 1 , 
     .         efermi , rv_p_owtkb , ifi ) 

              call fclose(ifi)
              if (lerr .eq. 0) lwtkb = 1
              if (lerr .ne. 0) efermi = -99
            endif
C           Broadcast lwtkb,efermi,wtkb
            call mpibc1(lwtkb,1,2,.false.,'subzi','lwtkb')
            if (lwtkb .eq. 1) then
              call mpibc1(efermi,1,4,.false.,'subzi','efermi')
              call mpibc1 ( rv_p_owtkb , ndham * nsp * nkp , 4 , .false. , 
     .         'subzi' , 'wtkb' ) 

            endif
          endif
        endif
      endif
      lswtk = -2
cki      oswtk = 1
ckino change rv_p_oswtk allocatable      nullify(rv_p_oswtk)
      if (nspc .eq. 2) then
        if (lwtkb .eq. 1) lswtk = 1
         if (allocated(rv_p_oswtk)) deallocate( rv_p_oswtk ) 
         allocate(rv_p_oswtk(abs(ndham*nsp*nkp)))
         if (ndham*nsp*nkp<0) rv_p_oswtk(:)=0.0d0

      endif

      if (nevmx .eq. 0) then
        nevmx = (int(zval) + 1)/2
        if (lmet .ne. 0) nevmx = max(nevmx+nevmx/2,9)
        nevmx = min(nevmx,ndham)
        if (nspc .eq. 2) nevmx = 2*nevmx
      endif

C ... Printout
      if (nevmx .ge. 0 .and. iprint() .gt. 30) then
        stdo = lgunit(1)
        if (lmet .gt. 0) then
          call awrit0('%N subzi: '//strni(isw(ltet)+1)//
     .    '%a integration of bands; '//
     .    strni(isw(lmet.ne.4.and.ltet)+1)//
     .    '%a integration of density',' ',80,stdo)
          if (lmet .eq. 4 .or. .not. ltet) then
            call info8(30,0,0,
     .      '%7p sampling integration uses:  '//
     .      '%?#(n<0)#Fermi distribution (T=%;4d)#%-1jN=%i  W=%;4d#'//
     .      '  ef0=%;6;6d%?#(n==4)#  def=%;6d##',
     .      n,esmear,ef0,lmet,def,0,0,0)
          endif
        else
          call info(20,0,0,' subzi : nonmetal',0,0)
        endif
        write(stdo,'(1x)')
      endif

      end subroutine subzi
      end module m_subzi

      integer function iobzwt(mode,nevx,nq,nsp,efermi,wtkb,ifi)
C- File I/O of dos weights
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 File I/O header and band weights
Ci         :  Header consists of parms nevx,nq,nsp,efermi
Ci         :  For file read, nevx,nq,nsp must match passed values
Ci         :1 File I/O header information only
Ci         :  For file read, nevx,nq,nsp must match passed values
Ci   nevx  :leading dimension of wtkb
Ci   nq    :number of irreducible k-points (bzmesh.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  efermi:ifi>0: Fermi level written to disk
Cio        :ifi<0: Fermi level read from disk
Cio  wtkb  :ifi>0: weights for k-point integration written to disk
Cio        :ifi<0: weights for k-point integration read from disk
Co Outputs
Co   iobzwt:0  File I/O was successful
Co   iobzwt:-1 File I/O was not successful
Cr Remarks
Cu Updates
Cu   16 May 01 when writing header info only, write nq=0
Cu   18 Feb 01 Added Fermi level to file (and argument list)
Cu    5 May 00 Adapted from nfp pvwts
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nevx,nsp,nq,ifi
      double precision efermi,wtkb(nevx,nsp,nq)
C ... Local parameters
      integer jfi,nevx0,nq0,nsp0,iprint,lgunit
      logical:: isanrg, l_dummy_isanrg

Ckino isanrg is logical function,       call isanrg(mode,0,1,' iobzwt','mode',.false.)
      l_dummy_isanrg=isanrg(mode,0,1,' iobzwt','mode',.false.)

      iobzwt = 0
C ... File read
      if (ifi .gt. 0) then
        rewind ifi
        read(ifi,end=80,err=80) nevx0,nq0,nsp0,efermi
        if (mode .eq. 1) then
          call info(20,1,0,
     .    ' Read efermi from weights file : ef = %,6;6d',efermi,0)
          return
        endif
        if (nevx .ne. nevx0 .or. nq .ne. nq0 .or. nsp .ne. nsp0) goto 80
        read(ifi) wtkb
        call info(20,1,0,' Read qp weights ...  ef=%;6d',efermi,0)
        return
   80   continue
        iobzwt = -1
        call info(20,1,0,' Incompatible or missing qp weights file ...',
     .  0,0)
        return

C ... File write
      else
        jfi = -ifi
        rewind jfi
        if (mode .ne. 1) write(jfi) nevx,nq,nsp,efermi
        if (mode .eq. 1) write(jfi) nevx,0,nsp,efermi
        if (mode .eq. 1) then
          if (iprint() .ge. 20) then
            call awrit1('%N Saved Fermi level to weights file ... '//
     .      'ef = %,6;6d',' ',80,lgunit(1),efermi)
          endif
          return
        endif
        write(jfi) wtkb
        call info(20,1,0,' Saved qp weights ...',0,0)
      endif
      end



