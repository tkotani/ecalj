!> Get the weights and index for tetrahedron method for the Lindhard function.
!!    - nbnb = total number of weight.
!!    - n1b  = band index for occ.   1\ge n1b \ge nband+nctot.
!!       "Valence index->core index" ordering(Core index follows valence index).
!!    - n2b  = band index for unocc. 1\ge n2b \ge nband
!!    - wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
!!
!! - NOTE: 'call getbzdata1' generates nteti,ntetf,... See mkqg.F about how to call it.
!!   
      module m_tetwt
      real(8),allocatable :: whw(:)
      integer,allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      integer:: nbnbx,nhwtot
      integer,allocatable :: n1b(:,:,:),n2b(:,:,:),nbnb(:,:)
!!
      contains !! ----------------------------------------------------------------------
!! routine
      subroutine gettetwt(q,iq,is,isf,nwgt) !this routine set output data above.
!! input data; read only
      use m_readeigen,only:   readeval !we assume init_readeval is called already
      use m_genallcf_v3,only: ecore,nctot    !we assume genallcf_v3 called already.
      use m_read_bzdata,only: nqbz,qbas,ginv,nqbzw,nteti,ntetf,idtetf,qbzw,ib1bz,nqibz,qbz
        ! we assume read_bzdata is called already
      use m_freq,only:                   !we assume getfreq is called already.
     &   frhis, nwhis,npm !output of getfreq
      use m_zmel,only: nband

      implicit none
      real(8),intent(in):: q(3)
      integer,intent(in):: is,isf,iq,nwgt(*)

      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      logical,allocatable :: iwgt(:,:,:,:)
      integer,allocatable:: nbnbtt(:,:),noccxvv(:) !  &         idtetf(:,:),ib1bz(:)
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:) !qbzw(:,: )
      logical :: eibzmode,tetra,tmpwwk=.false.,debug,eibz4x0
      integer:: kx,ncc,job,jpm,noccxvx(2)=-9999,ik,jhwtot,ib1,ib2,ibib,noccx,noccxv,verbose,ifief,ifile_handle
      real(8)::ef

      tetra=.true.
      eibzmode = eibz4x0()
      debug=.false.
      if(verbose()>=100) debug=.true.

      if(.not.allocated(nbnb)) allocate( nbnb(nqbz,npm))
      allocate( nbnbtt(nqbz,npm),ekxx1(nband,nqbz),ekxx2(nband,nqbz))

!!===========tetraini block tetra==.true.===============================1ini
c      if(tetra) then
      write(6,"(' tetra mode nqbz nband ispin q=',2i7,i2,3f13.6)") nqbz,nband,is,q
!!     ekxx1 for rk
!!     ekxx2 for q+rk See tetwt4
      do kx = 1, nqbz
        call readeval(qbz(:,kx),   is,  ekxx1(1:nband, kx) ) 
        call readeval(q+qbz(:,kx), isf, ekxx2(1:nband, kx) )
      enddo
c     takao-feb/2002 I replaced tetwt4 (1d30) with tetwt5(job=0) -----
C     ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
c     the pairs are not dependent on the energy omega
c     in the denominator of the dielectric function.
      write(6,"(' -- First tetwt5 is to get size of array --')")
      job = 0
      if(npm==1) then
        ncc=0
      else
        ncc=nctot
      endif
      allocate( demin(nband+nctot,nband+ncc,nqbz,npm),
     &          demax(nband+nctot,nband+ncc,nqbz,npm) )
      allocate( iwgt (nband+nctot,nband+ncc,nqbz,npm) )
!     wgt, demin, demax may require too much memory in epsilon mode.
!     We will have to remove these memory allocations in future.
!     tetwt5x_dtet2 can be very slow because of these poor memory allocation.
      if(nctot==0) then
        deallocate(ecore)
        allocate(ecore(1,2))    !this is dummry
      endif
      allocate(ibjb(1,1,1,1),ihw(1,1,1),jhw(1,1,1),nhw(1,1,1),whw(1)) !dummy
c--- EFERMI
      ifief=ifile_handle()
      open(ifief,file='EFERMI')
      read(ifief,*) ef
      close(ifief)
      call tetwt5x_dtet4(npm,ncc,  
     i q, ekxx1, ekxx2, qbas,ginv,ef, 
     d ntetf,nqbzw, nband,nqbz,
     i nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i job,
     o iwgt,nbnb,               !job=0
     o demin,demax,             !job=0
     i frhis, nwhis,            ! job=1    not-used
     i nbnbx,ibjb,nhwtot,       ! job=1    not-used
     i ihw,nhw,jhw,             ! job=1    not-used
     o whw,                     ! job=1    not-used
     i iq,is,isf,nqibz, eibzmode,nwgt) 
      deallocate(ibjb,ihw,jhw,nhw,whw) !dummy
      nbnbx = maxval(nbnb(1:nqbz,1:npm)) !nbnbx = nbnbxx
      if(debug) write(6,*)' nbnbx=',nbnbx
      allocate(  n1b(nbnbx,nqbz,npm)
     &          ,n2b(nbnbx,nqbz,npm))
      n1b=0; n2b=0 
      do jpm=1,npm
        call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx,
     o   n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
      enddo
      if(debug) then
        do kx  = 1, nqbz
          do jpm = 1, npm
            write(6,"('jpm kx  minval n1b n2b=',4i5)")jpm,kx,
     &       minval(n1b(1:nbnb(kx,jpm),kx,jpm)),
     &       minval(n2b(1:nbnb(kx,jpm),kx,jpm))
          enddo
        enddo
      endif
      if(sum(abs(nbnb-nbnbtt))/=0)then
        do ik=1,nqbz
          write(6,*)
          write(6,*)"nbnb  =",nbnb(ik,:)
          write(6,*)"nbnbtt=",nbnbtt(ik,:)
        enddo
        call rx( 'hx0fp0:sum(nbnb-nbnbtt)/=0')
      endif
      noccxv = maxval(noccxvx)
      noccx  = nctot + noccxv
      write(6,*)' Tetra mode: nctot noccxv= ',nctot,noccxv
      deallocate(iwgt)
c      endif
c=========end of tetraini block==========================================1end

!! TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
c     if(ixc==11) then !sf 21May02
C     --- METHOD (tetwt5) for the tetrahedron weight
!     Histogram secstions are specified by frhis(1:nwp)
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].
!     nwp=nw+1; frhis(1)=0
!     takao-feb/2002
      if(frhis(1)/=0d0) call rx( ' hx0fp0: we assume frhis(1)=0d0')
      write(6,*)' ----------------nbnbx nqbz= ',nbnbx,nqbz
!!     ... make index sets
      allocate(ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm))
      ihw=0; nhw=0; jhw=0
      jhwtot = 1
      do jpm =1,npm
        do ik   = 1,nqbz
          do ibib = 1,nbnb(ik,jpm)
            call hisrange( frhis, nwhis,  
     i       demin(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     i       demax(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     o       ihw(ibib,ik,jpm),nhw(ibib,ik,jpm))
            jhw(ibib,ik,jpm)= jhwtot
            jhwtot = jhwtot + nhw(ibib,ik,jpm)
          enddo
        enddo
      enddo
      nhwtot = jhwtot-1
      write(6,*)' nhwtot=',nhwtot
      deallocate(demin,demax)
      allocate( whw(nhwtot),    ! histo-weight
     & ibjb(nctot+nband,nband+ncc,nqbz,npm) )
      whw=0d0
      ibjb = 0
      do jpm=1,npm
        do ik   = 1,nqbz
          do ibib = 1,nbnb(ik,jpm)
            ib1  = n1b(ibib,ik,jpm)
            ib2  = n2b(ibib,ik,jpm)
            ibjb(ib1,ib2,ik,jpm) = ibib
          enddo
        enddo
      enddo
!!     ... Generate the histogram weights whw
      job=1
      write(6,*) 'goto tetwt5x_dtet4 job=',job
      allocate(demin(1,1,1,1),demax(1,1,1,1),iwgt(1,1,1,1)) !dummy
      call tetwt5x_dtet4(  npm,ncc, 
     i q, ekxx1, ekxx2, qbas,ginv,ef, 
     d ntetf,nqbzw, nband,nqbz,
     i nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i job,
     o iwgt,nbnb,               ! job=0
     o demin,demax,             ! job=0
     i frhis,nwhis,             ! job=1
     i nbnbx,ibjb,nhwtot,       ! job=1
     i ihw,nhw,jhw,             ! job=1
     o whw,                     ! job=1
     i iq,is,isf,nqibz, eibzmode,nwgt) 
      deallocate(demin,demax,iwgt) !duumy
!! ======TetrahedronWeight_5 block end =========
      end subroutine
      end module
