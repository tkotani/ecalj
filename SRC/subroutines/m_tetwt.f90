!> Tetrahedron weights are stored in this module
module m_tetwt 
!! output of gettetwt is passed to x0kf_v4hz
!!     nbnbx
!!     ihw(ibjb,kx): omega index, to specify the section of the histogram.
!!     nhw(ibjb,kx): the number of histogram sections
!!     jhw(ibjb,kx): pointer to whw
!!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
!!     : histogram weights for given ib,jb,kx for histogram sections
!!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.

!> Get the weights and index for tetrahedron method for the Lindhard function.
!!    - nbnb = total number of weight.
!!    - n1b  = band index for occ.   1\ge n1b \ge nband+nctot.
!!       "Valence index->core index" ordering(Core index follows valence index).
!!    - n2b  = band index for unocc. 1\ge n2b \ge nband
!!    - wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
!!
!! - NOTE: 'call getbzdata1' generates nteti,ntetf,... See mkqg.F about how to call it.
!!
  implicit none
  !! output ------------------------
  real(8),allocatable,protected,public :: whw(:)
  integer,allocatable,protected,public :: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
  integer,allocatable,protected,public :: n1b(:,:,:),n2b(:,:,:),nbnb(:,:)
  integer,protected,public:: nbnbx,nhwtot
  public:: Gettetwt, Tetdeallocate
  private
  !!----------------------------------------------
contains
  subroutine Tetdeallocate()
    deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b,nbnb)
  end subroutine Tetdeallocate
  subroutine Gettetwt(q,iq,is,isf,ekxx1,ekxx2,nband,wan )  
    use m_genallcf_v3,only: niw_in=>niw,ecore,nctot,nspin
    use m_freq,only: Getfreq2, frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,niw !output of getfreq
    use m_read_bzdata,only: qlat,ginv, ntetf,idtetf,ib1bz,nqibz_mtet=>nqibz,nqbz,qbz,nqbzw,qbzw, idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
    use m_ReadEfermi,only: ef
    use m_readgwinput,only: ebmx,nbmx,mtet
    use m_tetwt5,only:tetwt5x_dtet4,rsvwwk00_4,hisrange
    intent(in)::      q,iq,is,isf,ekxx1,ekxx2,nband,wan
    !! nqibz_mtet: is only for mtet/=(/1,1,1/) --->(we usually use only this case)
    !!
    !! output data in returened in the module variables above.

    !! we assume read_bzdata is called already
    !      use m_read_bzdata,only: qlat,ginv, ntetf,idtetf,ib1bz !, qbzw,nqbzw,qbz, nqibz

    !      use m_readeigen,only:   readeval !we assume init_readeval is called already
    !      use m_genallcf_v3,only: ecore,nctot    !we assume genallcf_v3 called already.
    !      use m_read_bzdata,only: nqbz,qlat,ginv,nqbzw,nteti,ntetf,idtetf,qbzw,ib1bz,nqibz,qbz
    !      use m_freq,only:                   !we assume getfreq is called already.
    !     &   frhis, nwhis,npm !output of getfreq
    !      use m_zmel,only: nband
    !      use m_ReadEfermi,only: readefermi,ef

    integer:: is,isf,iq,nband !,nwgt(:)
    !      integer:: ntetf,idtetf(0:3,ntetf),ib1bz(nqbzw)
    real(8):: q(3) !,qlat(3,3),ginv(3,3),ef,qbz(3,nqbz),qbzw(3,nqbzw),ebmx
    real(8):: ekxx1(nband,nqbz),ekxx2(nband,nqbz) !qbzw(:,: )
    !      real(8):: frhis(1:nwhis+1),ecore(nctot,2)
    real(8),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
    logical,allocatable :: iwgt(:,:,:,:)
    integer,allocatable:: nbnbtt(:,:),noccxvv(:) !  &         idtetf(:,:),ib1bz(:)
    logical :: eibzmode,tetra,tmpwwk=.false.,debug
    integer::kx,ncc,job,jpm,noccxvx(2)=-9999,ik,jhwtot,ib1,ib2,ibib,noccx,noccxv,verbose,ifief
    real(8),allocatable:: ecore_(:,:)
    integer:: ix,iqx
    logical,optional:: wan
    logical:: wan1,cmdopt0

    if(nctot==0) then
       allocate(ecore_(1,2))    !this is dummy
    else
       allocate(ecore_(nctot,2))
       ecore_=ecore
    endif

    tetra=.true.
    !      eibzmode = eibz4x0()
    debug=cmdopt0('--debug')
    !      if(.not.allocated(nbnb))
    allocate( nbnb(nqbz,npm)   )
    allocate( nbnbtt(nqbz,npm) )
!    do iqx=1,nqbz
!      write(6,*)'eee1=',ekxx1(:,iqx)-ef
!      write(6,*)'eee2=',ekxx2(:,iqx)-ef
!    enddo  
    !!===========tetraini block tetra==.true.===============================1ini
    write(6,"(' tetra mode: nqbz nband=',2i7,'nctot ispin q=',i2,i2,3f13.6)") nqbz,nband,nctot,is,q

    !     takao-feb/2002 I replaced tetwt4 (1d30) with tetwt5(job=0) -----
    !     ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
    !     the pairs are not dependent on the energy otemega
    !     in the denominator of the dielectric function.
    write(6,"(' -- First tetwt5 is to get size of array --')")
    job = 0
    if(npm==1) then
       ncc=0
    else
       ncc=nctot
    endif
    allocate( demin(nband+nctot,nband+ncc,nqbz,npm), &
         demax(nband+nctot,nband+ncc,nqbz,npm) )
    allocate( iwgt (nband+nctot,nband+ncc,nqbz,npm) )
    !     wgt, demin, demax may require too much memory in epsilon mode.
    !     We will have to remove these memory allocations in future.
    !     tetwt5x_dtet2 can be very slow because of these poor memory allocation.
    !      if(nctot==0) then
    !        deallocate(ecore)
    !        allocate(ecore(1,2))    !this is dummry
    !      endif
    allocate(ibjb(1,1,1,1),ihw(1,1,1),jhw(1,1,1),nhw(1,1,1),whw(1)) !dummy
    !--- EFERMI
    !      open(ifief,file='EFERMI')
    !      read(ifief,*) ef
    !      close(ifief)
    !      call readefermi() !comment out,since ef is passed nov2016
    ! ccccccccccccccc
    !      print *,'nqbz,nqbzw,nteti,ntetf,nqibz_mtet=',nqbz,nqbzw,nteti,ntetf,nqibz_mtet
    wan1=.false.
    if (present(wan)) then
       if(wan) wan1= .TRUE. 
    endif
    call tetwt5x_dtet4(npm,ncc, q, ekxx1, ekxx2, qlat,ginv,ef, ntetf,nqbzw, nband,nqbz, nctot,ecore_(1,is),idtetf,qbzw,ib1bz,job, &
         iwgt,nbnb,   demin,demax,                          & !job=0
         frhis, nwhis,nbnbx,ibjb,nhwtot,  ihw,nhw,jhw, whw, & ! job=1    not-used
         iq,is,isf,nqibz_mtet,&
         nbmx,ebmx,mtet, wan1 ) !Jan2019 for Wannier
    deallocate(ibjb,ihw,jhw,nhw,whw) !dummy
    nbnbx = max(maxval(nbnb(1:nqbz,1:npm)),1) !nbnbx = nbnbxx
    write(6,*)' nnnnnnnnn nbnbx=',maxval(nbnb(1:nqbz,1:npm))

    
    allocate(  n1b(nbnbx,nqbz,npm) ,n2b(nbnbx,nqbz,npm))
    n1b=0
    n2b=0
    do jpm=1,npm
       call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx, &
            n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
    enddo
    !!
    if(debug) then
       do kx  = 1, nqbz
          do jpm = 1, npm
             if( nbnb(kx,jpm) >0) then
                write(6,"('jpm kx minval n1b,n2b=',5i5)")jpm,kx,nbnb(kx,jpm), &
                     minval(n1b(1:nbnb(kx,jpm),kx,jpm)), &
                     minval(n2b(1:nbnb(kx,jpm),kx,jpm))
             else
                write(6,"('jpm kx minval nbnb=0',5i5)")jpm,kx
             endif
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
    !      endif
    !=========end of tetraini block==========================================1end

    !! TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
    !     if(ixc==11) then !sf 21May02
    !     --- METHOD (tetwt5) for the tetrahedron weight
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
             call hisrange( frhis, nwhis, &
                  demin(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm), &
                  demax(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm), &
                  ihw(ibib,ik,jpm),nhw(ibib,ik,jpm))
             jhw(ibib,ik,jpm)= jhwtot
             jhwtot = jhwtot + nhw(ibib,ik,jpm)
          enddo
       enddo
    enddo
    nhwtot = jhwtot-1
    write(6,*)' nhwtot=',nhwtot
    deallocate(demin,demax)
    allocate( whw(nhwtot), ibjb(nctot+nband,nband+ncc,nqbz,npm) )
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
    if (present(wan)) then
       if(wan) wan1= .TRUE. 
    endif
    call tetwt5x_dtet4(  npm,ncc, &
         q, ekxx1, ekxx2, qlat,ginv,ef, &
         ntetf,nqbzw, nband,nqbz, &
         nctot,ecore_(1,is),idtetf,qbzw,ib1bz, &
         job, &
         iwgt,nbnb,           &  !job=0
    demin,demax,              &  !job=0
    frhis,nwhis,              &  !job=1
    nbnbx,ibjb,nhwtot,        &  !job=1
    ihw,nhw,jhw,              &  !job=1
    whw,                      &  !job=1
    iq,is,isf,nqibz_mtet, &!eibzmode,nwgt, &
         nbmx,ebmx,mtet, wan1)           !Jan2019
    deallocate(demin,demax,iwgt,nbnbtt)

    !! ======TetrahedronWeight_5 block end =========
  end subroutine gettetwt
end module m_tetwt
