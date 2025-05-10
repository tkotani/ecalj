module m_genallcf_v3 ! Readin starting data dat in GWinput 
 use m_lgunit,only:stdo
  use m_mpi,only: ipr
  implicit none
  public:: setesmr, genallcf_v3
  integer,protected,public:: nrx,lcutmx
  real(8),allocatable,public::cutbase(:)
  integer,allocatable,protected,public::  & 
       nindx(:,:),konf(:,:),icore(:,:), ncore(:), nlnm(:),nlnmv(:), nlnmc(:), il(:,:), in(:,:), im(:,:),&
       nocc(:,:,:),nunocc(:,:,:),nindxc(:,:),lcutmxa(:),lmxa(:)
  integer,protected,public:: natom,nspin,nl,nn,nnv,nnc, nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot, niw,ndimanspc,nlmto
  real(8),protected,public::  plat(3,3),alat,deltaw,esmr,delta,tpioa,qval
  real(8), allocatable,protected,public:: pos(:,:),z(:),ecore(:,:) !,symgg(:,:,:)
  character(8),allocatable,protected,public :: spid(:)
  integer,protected,public:: nprecb,mrecb,mrece,ndima,nqbzt,nband,mrecg,nspc,nspx 
  logical,protected,public:: laf !: antiferro switch
  integer,allocatable,protected,public:: ibasf(:) ! specify AF pair atom.
  private
  logical,protected,private:: done_genallcf_v3=.false.
contains
  subroutine setesmr(esmr_in)
    intent(in)::     esmr_in
    real(8):: esmr_in
    esmr=esmr_in
  end subroutine setesmr
  subroutine genallcf_v3(incwfx)
    use m_keyvalue,only: getkeyvalue
    implicit none
    intent(in)::           incwfx
    ! incwfx: product basis setting switch
    !! input: efin,incwfx,
    !!        GWinput, MLOindex, ECORE
    !! output: All the output are given in the declear section above.
    ! Original idea of product basis is from F.Aryasetiawan. Some subroutines are written by him.
    ! We may need to clean them up in modern fortran.
    !! --------------------------------------------------------
    integer:: ifhbe
    integer::incwfx,ifec,i,j, ifi,ig,is,ix,ixoff,lx, infwfx,ret, n1,n2,n3,imagw,n,iatom,iatomt
    logical :: nocore,readon
    character(1000) :: tolchar
    real(8),   allocatable:: ecoret(:,:,:,:)
    integer,allocatable::ncwf2(:,:,:), nindxv(:,:),occv(:,:,:),unoccv(:,:,:), occc(:,:,:),unoccc(:,:,:),ncwf(:,:,:)!,iclass(:)
    integer:: ia,l,m,ic1,isp,lt,nt,nr,ncorex,ifix,nwdummy,ind,l2,lm,lmxax
    real(8):: a,b,zz, efdummy,dw,diw,efin
    real(8),parameter:: pi=4d0*datan(1d0)
    if(done_genallcf_v3) call rx('genallcf_v3 is already called')
    done_genallcf_v3=.true.
    open(newunit=ifi,file='__MTOindex',form='unformatted')
    read(ifi) natom,alat,plat,nspin,lmxax,nnv,nnc,nrx,qval,nspc,nlmto !,n1,n2,n3
    allocate(pos(3,natom),z(natom),spid(1:natom),ibasf(natom),lmxa(natom))
    read(ifi) pos,z(1:natom),spid(1:natom),lmxa(1:natom)
    read(ifi) nprecb,mrecb,mrece,ndima,nqbzt,nband,mrecg
    read(ifi) laf,ibasf
    close(ifi)
    nl   = lmxax+1
    tpioa= 2d0*pi/alat
    call getkeyvalue("GWinput","niw",   niw ) ! FREQUENCIES
    call getkeyvalue("GWinput","delta", delta )
    call getkeyvalue("GWinput","deltaw",deltaw )
    call getkeyvalue("GWinput","esmr",  esmr )
    if(ipr) write(stdo,*)' --- Freq ---'
    if(ipr) write(stdo,"(a,i6)")   '    niw  =',niw
    if(ipr) write(stdo,"(a,f12.6)")'    delta=',delta
    if(ipr) write(stdo,"(a,f12.6)")'    esmr =',esmr
    ReadProductBasis: block
      allocate(nindxv(nl,natom), nindxc(nl,natom), &
           occv(nl,nnv,natom),unoccv(nl,nnv,natom), &
           occc(nl,nnc,natom),unoccc(nl,nnc,natom),source=0)
      allocate(ncwf2(nl,nnc,natom),ncwf(nl,nnc,natom),source=0)
      allocate(cutbase(0:2*(nl-1)),source=0d0)
      ncwf  =99 !This is for counting the number of nctot in gencor.
      ncwf2 =99
      if(ipr) write(stdo,*)' reading <PRODUCT_BASIS> section'
      call getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifi,status=ret)!open GWinput and locate file position
      read(ifi,*)
      read(ifi,"(a)") tolchar !tolerance in percentage for optimal product basis
      readon=.false.
      lx=0
      do ix=1,1000
         if( .NOT. readon .AND. tolchar(ix:ix)/=' ') then
            readon=.true.
            ixoff=ix
         endif
         if(readon .AND. tolchar(ix:ix)==' ') then
            read(tolchar(ixoff:ix),*,err=1097) cutbase(lx)
            if(lx==2*(nl-1)) goto 1098
            readon=.false.
            lx=lx+1
         endif
      enddo
1097  continue
      cutbase(lx:)=cutbase(lx-1)
1098  continue
      do lx=0,2*(nl-1)
         if(ipr) write(stdo,"(' lx=',i3,' readin tolerance=',d11.3)") lx, cutbase(lx)
      enddo
      read(ifi,*)
      allocate(lcutmxa(1:natom))
      read(ifi,*) lcutmxa(1:natom)
      lcutmx=lcutmxa(1)
      if(ipr) write(stdo,'(20i3)') lcutmxa(1:natom)
      if(ipr) write(stdo,"(' --- prod section: lcutmx cutbase='i3,100d11.3)") lcutmx,cutbase
      read(ifi,*)
      do    iatom = 1,natom
         do l  = 0,lmxa(iatom) !nl-1
            read(ifi,*) iatomt,lt,nindxv(l+1,iatom),nindxc(l+1,iatom)
            if(ipr) write(stdo,*)iatomt,lt,nindxv(l+1,iatom),nindxc(l+1,iatom)
         enddo
      enddo
      if(ipr) write(stdo,*)' --- valence product basis section'
      occv=0
      unoccv=0
      read(ifi,*)
      do       iatom = 1,natom
         do     l = 0,lmxa(iatom)
            do  n = 1,nindxv(l+1,iatom)
               read(ifi,*)                  iatomt,lt,nt,occv(l+1,n,iatom),unoccv(l+1,n,iatom)
               if(ipr) write(stdo,"(100i3)") iatomt,lt,nt,occv(l+1,n,iatom),unoccv(l+1,n,iatom)
            enddo
         enddo
      enddo
      if(ipr) write(stdo,*)' --- core product basis section'
      read(ifi,*)
      do       iatom = 1,natom
         do    l  = 0,lmxa(iatom)
            do n  = 1,nindxc(l+1,iatom)
               read(ifi,*)                  iatomt,lt,nt,occc(l+1,n,iatom),unoccc(l+1,n,iatom),ncwf(l+1,n,iatom),ncwf2(l+1,n,iatom)
               if(ipr)write(stdo,"(100i3)") iatomt,lt,nt,occc(l+1,n,iatom),unoccc(l+1,n,iatom),ncwf(l+1,n,iatom),ncwf2(l+1,n,iatom)
            enddo
         enddo
      enddo
      close(ifi)
      if( incwfx==-1 ) then ;  if(ipr) write(stdo,*)' ### incwf=-1 Use ForSxc for core' !----- product basis setting
         ncwf = ncwf2
      elseif( incwfx==-2 ) then; if(ipr) write(stdo,*)' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
         ncwf = merge(1-ncwf2,ncwf2, ncwf2==0.or.ncwf2==1) 
         occc = ncwf
         unoccc= 0
         unoccv= merge(1,unoccv, occv==1.or.unoccv==1)
      elseif( incwfx==-3 ) then ; if(ipr) write(stdo,*)' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
         occc=1
         ncwf=1
         unoccc= 0
      elseif( incwfx==-4 ) then; if(ipr) write(stdo,*)' ### incwf=-4  occ=0 and unocc=0 for all core '
         occc=0
         unoccc=0
         ncwf=0
      elseif(incwfx==0) then;    if(ipr) write(stdo,*)' ### Use unocc occ ForX0 for core'
      else;    call rx( ' ### proper incwf is not given for genallcf_v3:rgwinf ')
      endif
      deallocate(ncwf2)
    endblock ReadProductBasis
    
    indexcoremto: block ! new nindx !-------------------------------------------- lmxa(iatom) instead of nl-1
      integer:: nnn1(natom),nnn2(natom),nnn3(natom),nnn4(natom),nnn5(natom),nnn6(natom),nlx
      ndima=0
      do iatom=1,natom
         ndima=ndima+sum([((2*l+1)*nindxv(l+1,iatom),l=0,lmxa(iatom))]) !l=0,nl-1)])
      enddo
      nn  =  maxval(nindxv(1:nl,1:natom)+nindxc(1:nl,1:natom))
      allocate(nindx(nl,natom),nocc(nl,nn,natom),nunocc(nl,nn,natom),source=0)
      reindxblock: block
        integer:: nval,ncore
        do    iatom = 1,natom
           do  l = 0,lmxa(iatom) !nl-1
              ncore  = nindxc(l+1,iatom)
              nval   = nindxv(l+1,iatom)
              nindx(l+1,iatom)= ncore + nval
              nocc(l+1,1:,iatom)   = [  (occc(l+1,n,iatom),n=1,ncore),  (occv(l+1,n,iatom),n=1,nval)]
              nunocc(l+1,1:,iatom) = [(unoccc(l+1,n,iatom),n=1,ncore),(unoccv(l+1,n,iatom),n=1,nval)]
           enddo
        enddo
      endblock reindxblock
      block
        do iatom=1,natom
           nlx=lmxa(iatom)+1
           nnn1(iatom)=sum(nindx(1:nlx,iatom))
           nnn2(iatom)=sum([(nindx(l+1,iatom)*(2*l+1),l=0,nlx-1)])
           nnn3(iatom)=sum(  nindxv(1:nlx,iatom))
           nnn4(iatom)=sum([(nindxv(l+1,iatom)*(2*l+1),l=0,nlx-1)])
           nnn5(iatom)=sum(  nindxc(1:nlx,iatom))
           nnn6(iatom)=sum([(nindxc(l+1,iatom)*(2*l+1),l=0,nlx-1)])
        enddo
        nlnx    = maxval(nnn1) ;  nlnmx   = maxval(nnn2);  nlnxv   = maxval(nnn3)
        nlnmxv  = maxval(nnn4) ;  nlnxc   = maxval(nnn5);  nlnmxc  = maxval(nnn6)
      endblock
      allocate(il(nlnmx,natom), in(nlnmx,natom), im(nlnmx,natom)) ! index for core and MTO basis =====================
      do iatom = 1,natom
        ind  = 0
        do l = 0,lmxa(iatom) !nl-1 ! core
          do  n = 1,nindxc(l+1,iatom)
            do  m = 1,2*l+1
              ind       = ind + 1
              lm        = l**2 + m
              il(ind,iatom)= l;   in(ind,iatom) = n;  im(ind,iatom)= m-l-1 !; ilnm(n,lm,iatom) = ind
            enddo
          enddo
        enddo
        do  l = 0,lmxa(iatom) !nl-1 ! valence
          ncorex  = nindxc(l+1,iatom)
          do    n = 1,nindxv(l+1,iatom)
            do      m = 1,2*l+1
              ind = ind + 1
              lm = l**2 + m
              il(ind,iatom)  =l;  in(ind,iatom) = ncorex + n; im(ind,iatom)  = m-l-1 !; ilnm(ncorex+n,lm,iatom) = ind
            enddo
          enddo
        enddo
      enddo
      allocate(nlnmv(natom),nlnmc(natom),nlnm(natom))
      do iatom=1,natom
         nlx=lmxa(iatom)+1
         nlnmv(iatom) = sum([(nindxv(l+1,iatom)*(2*l+1),l=0,nlx-1)])
         nlnmc(iatom) = sum([(nindxc(l+1,iatom)*(2*l+1),l=0,nlx-1)])
         nlnm(iatom)  = sum([(nindx(l+1,iatom)*(2*l+1),l=0,nlx-1)])
      enddo
    endblock indexcoremto
    coreblock: block
      real(8),external:: rydberg
      allocate(icore(nl**2*nnc,natom),ncore(natom),source=99999)
      do iatom = 1,natom ! index for allowed core states
         i  = 0
         j  = 0
         do       l = 0,lmxa(iatom) !nl-1
            do    n = 1,nindxc(l+1,iatom)
               do m = -l,l
                  j = j + 1
                  if (ncwf(l+1,n,iatom) == 1) then
                     i = i + 1
                     icore(i,iatom)= j
                  endif
               enddo
            enddo
         enddo
         ncore(iatom)  = i
      end do
      nctot = sum(ncore(1:natom))
      open(newunit=ifec,file='ECORE')         ! core energies ==========================
      read(ifec,*)
      allocate(konf(nl,natom),source=0)
      allocate(ecoret(0:nl-1,nnc,2,natom),ecore(nctot,2))
      ecoret=0d0
      do iatom = 1,natom
         if(ipr) write(stdo,*) ' read ECORE : iatom lmxa =',iatom,lmxa(iatom)
         read (ifec,*)
         read (ifec,*) (konf(l+1,iatom),l=0,lmxa(iatom))                ! number of cores for l 
         konf(1:lmxa(iatom)+1,iatom)=[(konf(l+1,iatom)+l+1,l=0,lmxa(iatom))]  ! minor change of ECORE but konf unchanged. 2025-5-10
         do  l = 0,lmxa(iatom) !nl-1
            ncorex = konf(l+1,iatom) -l-1 
            do n = 1,ncorex
               read (ifec,*) lt,nt,(ecoret(l,n,isp,iatom),isp=1,nspin) !takao
               if(nspin==1) ecoret(l,n,2,iatom) = ecoret(l,n,1,iatom)
               if(ipr) write(stdo,"(' read ecore=',3i4,2d13.5)")l,n,iatom,ecoret(l,n,1:nspin,iatom)
            end do
         end do
      end do
      ecoret=ecoret/rydberg() !readin ECORE is in eV from 2025-5-9
      close(ifec)
      i = 0
      do ia = 1,natom
         iatom  = ia
         do l = 0,lmxa(iatom) !nl-1
            do n = 1,nnc
               do m = -l,l
                  if (ncwf(l+1,n,iatom) == 1) then
                     i = i + 1
                     ecore(i,1:nspin) = ecoret(l,n,1:nspin,iatom)
                     if(ipr) write(stdo,"(' ecore=',4i4,2d13.5)")i, l,n,iatom,ecore(i,1:nspin)
                  endif
               enddo
            enddo
         enddo
      enddo
      if(size(ecore)==0) then !jun2020 moved from sxcf
         deallocate(ecore)
         allocate(ecore(1,2))
      endif
      deallocate(ecoret)
    endblock coreblock
    ndimanspc=ndima*nspc
    nspx=nspin/nspc
    call cputid(0); if(ipr) write(stdo,*) 'genallcf_v3'
  end subroutine genallcf_v3
end module m_genallcf_v3

module m_ReadEfermi
  use m_lgunit,only:stdo
  use m_mpi,only: ipr
  real(8),protected:: bandgap, ef, ef_kbt
  public:: readefermi,readefermi_kbt,setefermi
contains
  subroutine setefermi(efin)
    real(8)::efin
    ef=efin
  endsubroutine setefermi
  subroutine readefermi()
    implicit none
    integer:: ifief
    open(newunit=ifief,file='EFERMI')
    read(ifief,*) ef,bandgap
    close(ifief)
    if(ipr) write(stdo,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
  end subroutine readefermi
  subroutine readefermi_kbt()
    implicit none
    integer:: ifief_kbt
    open(newunit=ifief_kbt,file='EFERMI_kbt')
    read(ifief_kbt,*) ef_kbt,bandgap
    close(ifief_kbt)
    if(ipr) write(stdo,"(a,f12.6)")' --- READIN ef from EFERMI_kbt. ef=',ef_kbt
  end subroutine readefermi_kbt
end module m_ReadEfermi

