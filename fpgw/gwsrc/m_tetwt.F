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
      real(8),allocatable,protected :: whw(:)
      integer,allocatable,protected:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      integer,protected:: nbnbx,nhwtot
      integer,allocatable,protected :: n1b(:,:,:),n2b(:,:,:),nbnb(:,:)
!!
      contains !! ----------------------------------------------------------------------
      subroutine tetdeallocate()
      deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b,nbnb)
      end subroutine

!! routine --------------------------------------------------------
      subroutine gettetwt(q,iq,is,isf,nwgt,frhis,nwhis,npm, 
     i   qbas,ginv, ef, nqibz_mtet, nband,ekxx1,ekxx2, nctot,ecore,
     i   nqbz,qbz,nqbzw,qbzw, ntetf,idtetf,ib1bz,
     i   nbmx,ebmx,mtet,eibzmode, !nov2016
     i   wan ) !Jan2019 okumura's option
!! INPUT DATA; read only
!! nqibz_mtet: is only for mtet/=(/1,1,1/) --->(we usually use only this case)
!! 
!! output data in returened in the module variables above.

!! we assume read_bzdata is called already
c      use m_read_bzdata,only: qbas,ginv, ntetf,idtetf,ib1bz !, qbzw,nqbzw,qbz, nqibz

c      use m_readeigen,only:   readeval !we assume init_readeval is called already
c      use m_genallcf_v3,only: ecore,nctot    !we assume genallcf_v3 called already.
c      use m_read_bzdata,only: nqbz,qbas,ginv,nqbzw,nteti,ntetf,idtetf,qbzw,ib1bz,nqibz,qbz
c      use m_freq,only:                   !we assume getfreq is called already.
c     &   frhis, nwhis,npm !output of getfreq
c      use m_zmel,only: nband
c      use m_ReadEfermi,only: readefermi,ef

      implicit none
      integer,intent(in):: is,isf,iq,nwgt(:),nqibz_mtet,nqbz,nqbzw,nband,npm,nwhis,nctot,nbmx,mtet(3)
      integer,intent(in):: ntetf,idtetf(0:3,ntetf),ib1bz(nqbzw)
      real(8),intent(in):: q(3),qbas(3,3),ginv(3,3),ef,qbz(3,nqbz),qbzw(3,nqbzw),ebmx
      real(8),intent(in):: ekxx1(nband,nqbz),ekxx2(nband,nqbz) !qbzw(:,: )
      real(8),intent(in):: frhis(1:nwhis+1),ecore(nctot,2)

      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      logical,allocatable :: iwgt(:,:,:,:)
      integer,allocatable:: nbnbtt(:,:),noccxvv(:) !  &         idtetf(:,:),ib1bz(:)
      logical :: eibzmode,tetra,tmpwwk=.false.,debug,eibz4x0
      integer::kx,ncc,job,jpm,noccxvx(2)=-9999,ik,jhwtot,ib1,ib2,ibib,noccx,noccxv,verbose,ifief,ifile_handle
      real(8),allocatable:: ecore_(:,:)
      integer:: ix
      logical,optional:: wan
      logical:: wan1

      if(nctot==0) then
        allocate(ecore_(1,2))    !this is dummy
      else
        allocate(ecore_(nctot,2))
        ecore_=ecore
      endif

      tetra=.true.
c      eibzmode = eibz4x0()
      debug=.false.
      if(verbose()>=100) debug=.true.
c      if(.not.allocated(nbnb)) 
      allocate( nbnb(nqbz,npm)   )
      allocate( nbnbtt(nqbz,npm) ) !,ekxx1(nband,nqbz),ekxx2(nband,nqbz))
!!===========tetraini block tetra==.true.===============================1ini
      write(6,"(' tetra mode nqbz nband ispin q=',2i7,i2,3f13.6)") nqbz,nband,is,q

c     takao-feb/2002 I replaced tetwt4 (1d30) with tetwt5(job=0) -----
C     ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
c     the pairs are not dependent on the energy otemega
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
c      if(nctot==0) then
c        deallocate(ecore)
c        allocate(ecore(1,2))    !this is dummry
c      endif
      allocate(ibjb(1,1,1,1),ihw(1,1,1),jhw(1,1,1),nhw(1,1,1),whw(1)) !dummy
c--- EFERMI
c      ifief=ifile_handle()
c      open(ifief,file='EFERMI')
c      read(ifief,*) ef
c      close(ifief)
c      call readefermi() !comment out,since ef is passed nov2016 
ccccccccccccccccc
c      print *,'nqbz,nqbzw,nteti,ntetf,nqibz_mtet=',nqbz,nqbzw,nteti,ntetf,nqibz_mtet
      wan1=.false.
      if (present(wan) .and. wan) wan1=.true.
      call tetwt5x_dtet4(npm,ncc,  
     i q, ekxx1, ekxx2, qbas,ginv,ef, 
     d ntetf,nqbzw, nband,nqbz,
     i nctot,ecore_(1,is),idtetf,qbzw,ib1bz, 
     i job,
     o iwgt,nbnb,               !job=0
     o demin,demax,             !job=0
     i frhis, nwhis,            ! job=1    not-used
     i nbnbx,ibjb,nhwtot,       ! job=1    not-used
     i ihw,nhw,jhw,             ! job=1    not-used
     o whw,                     ! job=1    not-used
     i iq,is,isf,nqibz_mtet, eibzmode,nwgt,
     i   nbmx,ebmx,mtet, !nov2016
     i     wan1 ) !Jan2019 for Wannier
 
      deallocate(ibjb,ihw,jhw,nhw,whw) !dummy
      nbnbx = max(maxval(nbnb(1:nqbz,1:npm)),1) !nbnbx = nbnbxx
      if(debug) write(6,*)' nbnbx=',nbnbx
      allocate(  n1b(nbnbx,nqbz,npm)
     &          ,n2b(nbnbx,nqbz,npm))
      n1b=0
      n2b=0 
      do jpm=1,npm
        call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx,
     o   n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
      enddo
!!
      if(debug) then
        do kx  = 1, nqbz
        do jpm = 1, npm
           if( nbnb(kx,jpm) >0) then
              write(6,"('jpm kx minval n1b,n2b=',5i5)")jpm,kx,nbnb(kx,jpm),
     &             minval(n1b(1:nbnb(kx,jpm),kx,jpm)),
     &             minval(n2b(1:nbnb(kx,jpm),kx,jpm))
           else
              write(6,"('jpm kx minval nbnb=0',5i5)")jpm,kx
           endif
        enddo
        enddo
c$$$cccccccccccccccccccccc
c$$$          do kx  = 1, nqbz
c$$$          do jpm = 1, npm
c$$$          do ix=1,nbnb(kx,jpm)
c$$$            write(6,"('jpm kx  minval n1b n2b=',3i5,2x,2i10)")jpm,kx,ix,
c$$$     &       n1b(ix,kx,jpm),n2b(ix,kx,jpm)
c$$$          enddo
c$$$          enddo
c$$$          enddo
c          call MPI__Finalize
c          call rx0( ' OK! hx0fp0_sc ixc=11 Sergey F. mode')
c          stop 'vvvvvvvvvvvvvv111'
cccccccccccccccccccccc
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
      if(abs(frhis(1))>1d-12) call rx( ' hx0fp0: we assume frhis(1)=0d0')
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

      wan1=.false.
      if (present(wan) .and. wan) wan1=.true.
      call tetwt5x_dtet4(  npm,ncc, 
     i q, ekxx1, ekxx2, qbas,ginv,ef, 
     d ntetf,nqbzw, nband,nqbz,
     i nctot,ecore_(1,is),idtetf,qbzw,ib1bz, 
     i job,
     o iwgt,nbnb,               ! job=0
     o demin,demax,             ! job=0
     i frhis,nwhis,             ! job=1
     i nbnbx,ibjb,nhwtot,       ! job=1
     i ihw,nhw,jhw,             ! job=1
     o whw,                     ! job=1
     i iq,is,isf,nqibz_mtet, eibzmode,nwgt,
     i nbmx,ebmx,mtet, !nov2016
     i wan1)           !Jan2019
      deallocate(demin,demax,iwgt,nbnbtt)

!! ======TetrahedronWeight_5 block end =========
      end subroutine gettetwt
      end module
