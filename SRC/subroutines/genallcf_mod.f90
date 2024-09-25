module m_genallcf_v3 ! Readin starting data dat in GWinput
  use m_lgunit,only:stdo
  implicit none
  public:: setesmr, genallcf_v3
  integer,protected,public:: nrx,lcutmx
  real(8),allocatable,public::cutbase(:)
  integer,allocatable,protected,public:: iclass(:), &
       nindx(:,:),konf(:,:),icore(:,:), ncore(:), &
       nlnm(:),nlnmv(:), nlnmc(:), il(:,:), in(:,:), im(:,:),&
       nocc(:,:,:),nunocc(:,:,:),nindxc(:,:),lcutmxa(:),lmxa(:)
  integer,protected,public:: nclass,natom,nspin,nl,nn,nnv,nnc,&
       nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot, niw,ndimanspc !ndima,
  real(8),protected,public::  plat(3,3),alat,deltaw,esmr,delta,tpioa,qval
  real(8), allocatable,protected,public:: pos(:,:),z(:),ecore(:,:) !,symgg(:,:,:)
  character(8),allocatable,protected,public:: spid(:)
  character(8),allocatable,protected,public :: clabl(:)
  integer,protected,public:: nprecb,mrecb,mrece,ndima,nqbzt,nband,mrecg,nspc,nspx !nspc=2 for so=1, zero otherwize.
  logical,protected,public:: laf !! - laf: antiferro switch
  integer,allocatable,protected,public:: ibasf(:) !! - ibasf(ibas) specify AF pair atom.
  private
  logical,protected,private:: done_genallcf_v3=.false.
!  integer,allocatable,protected,private:: &
!       ilv(:,:),inv(:,:),imv(:,:),  ilnmv(:,:,:),  &
!       ilc(:,:),inc(:,:),imc(:,:),  ilnmc(:,:,:) ,  ilnm(:,:,:)
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
    !! Return iclass=ibas.
    !    incwfx: product basis setting switch
    !! input: efin,incwfx,
    !!        GWinput, MLOindex, ECORE
    !! output: All the output are given in the declear section above.
    ! Original idea of product basis is from F.Aryasetiawan. Some subroutines are written by him.
    ! We may need to clean them up in modern fortran.
    !! --------------------------------------------------------
    integer:: ifhbe
    integer::incwfx,ifec,i,j, ifi,ig,is,ix,ixoff,lx
    integer:: infwfx,ret, n1,n2,n3,imagw,n,ic
    logical :: nocore,readon
    real(8)::efin
    character(1000) :: tolchar
    real(8),   allocatable:: ecoret(:,:,:,:)
    integer,allocatable::ncwf2(:,:,:), nindxv(:,:),occv(:,:,:),unoccv(:,:,:), occc(:,:,:),unoccc(:,:,:),ncwf(:,:,:)
    integer:: ia,l,m,ic1,isp,lt,nt,nr,ncorex,ifix
    real(8)::a,b,zz, efdummy,dw,diw
    integer:: nwdummy,ict,ind,l2,lm,lmxax1
    real(8),parameter:: pi=4d0*datan(1d0)
    if(done_genallcf_v3) call rx('genallcf_v3 is already called')
    done_genallcf_v3=.true.
    open(newunit=ifi,file='MTOindex',form='unformatted')
    read(ifi) natom,alat,plat,nspin,lmxax1,nnv,nnc,nrx,qval,nspc !,n1,n2,n3
    allocate(pos(3,natom),clabl(natom),z(natom),spid(1:natom),ibasf(natom),lmxa(natom))
    read(ifi) pos,z(1:natom),spid(1:natom),lmxa(1:natom)
    read(ifi) nprecb,mrecb,mrece,ndima,nqbzt,nband,mrecg
    read(ifi) laf,ibasf
    close(ifi)
    nclass = natom  !We set nclass = natom through the GW calculations
    nl=lmxax1
    clabl=spid
    tpioa=2d0*pi/alat
    call getkeyvalue("GWinput","niw",   niw ) ! FREQUENCIES
    call getkeyvalue("GWinput","delta", delta )
    call getkeyvalue("GWinput","deltaw",deltaw )
    call getkeyvalue("GWinput","esmr",  esmr )
    write(stdo,*)' --- Freq ---'
    write(stdo,"(a,i6)")   '    niw  =',niw
    write(stdo,"(a,f12.6)")'    delta=',delta
    write(stdo,"(a,f12.6)")'    esmr =',esmr
    allocate(iclass(natom),source=[(n,n=1,natom)]) !We set nclass = natom through the GW calculations
    ReadProductBasis: block
      allocate(nindxv(nl,nclass), nindxc(nl,nclass), &
           occv(nl,nnv,nclass),unoccv(nl,nnv,nclass), &
           occc(nl,nnc,nclass),unoccc(nl,nnc,nclass))
      allocate(ncwf2(nl,nnc,nclass),ncwf(nl,nnc,nclass))
      allocate(cutbase(0:2*(nl-1)))
      ncwf  =99 !This is for counting the number of nctot in gencor.
      ncwf2 =99
      !===================================================
      write(stdo,*)' reading <PRODUCT_BASIS> section'
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
         write(stdo,"(' lx=',i3,' readin tolerance=',d11.3)") lx, cutbase(lx)
      enddo
      read(ifi,*)
      allocate(lcutmxa(1:natom))
      read(ifi,*) lcutmxa(1:natom)
!      read(ifi,*)lcutmx
!  call Getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifinin,status=ret)
!  do
!     read(ifinin,*,err=980) aaa
!     if(aaa=='lcutmx(atom)') then
!        read(ifinin,*) lcutmxa(1:natom)
!        goto 990
!     endif
!  enddo
!980 continue
!  lcutmxa=lcutmx
!990 continue
      !  close(ifinin)
      lcutmx=lcutmxa(1)
      write(stdo,'(20i3)') lcutmxa(1:natom)
      write(stdo,"(' --- prod section: lcutmx cutbase='i3,100d11.3)") lcutmx,cutbase
      read(ifi,*)
      do    ic = 1,nclass
         do l  = 0,nl-1
            read(ifi,*) ict,lt,nindxv(l+1,ic),nindxc(l+1,ic)
            if(lt  /= l ) call rx( 'genallcf_mod /=l ')
         enddo
      enddo
      write(stdo,*)' --- valence product basis section'
      occv=0
      unoccv=0
      read(ifi,*)
      do       ic = 1,nclass
         do     l = 0,nl-1
            do  n = 1,nindxv(l+1,ic)
               read(ifi,*)           ict,lt,nt,occv(l+1,n,ic),unoccv(l+1,n,ic)
               write(stdo,"(100i3)") ict,lt,nt,occv(l+1,n,ic),unoccv(l+1,n,ic)
               if(lt  /= l )call rx( 'genallcf: wrong l valence')
               if(nt  /= n )call rx( 'genallcf: wrong n valence')
            enddo
         enddo
      enddo
      write(stdo,*)' --- core product basis section'
      read(ifi,*)
      do       ic = 1,nclass
         do    l  = 0,nl-1
            do n  = 1,nindxc(l+1,ic)
               read(ifi,*)           ict,lt,nt,occc(l+1,n,ic),unoccc(l+1,n,ic),ncwf(l+1,n,ic),ncwf2(l+1,n,ic)
               write(stdo,"(100i3)") ict,lt,nt,occc(l+1,n,ic),unoccc(l+1,n,ic),ncwf(l+1,n,ic),ncwf2(l+1,n,ic)  !ncwf2 is for Sigma calcuation
               if(lt /= l )call rx( 'rgwina: 2nd wrong l core')
               if(nt /= n )call rx( 'rgwina: wrong n core')
            enddo
         enddo
      enddo
      close(ifi)
      if( incwfx==-1 ) then !!----- product basis setting
         write(stdo,*)' ### incwf=-1 Use ForSxc for core'
         ncwf = ncwf2
      elseif( incwfx==-2 ) then
         write(stdo,*)' ### incwf=-2 Use NOT(ForSxc) for core and Pro-basis '
         ncwf = merge(1-ncwf2,ncwf2,ncwf2==0.or.ncwf2==1) 
         occc = ncwf
         unoccc= 0
         unoccv= merge(1,unoccv,occv==1.or.unoccv==1)
      elseif( incwfx==-3 ) then 
         occc=0
         ncwf=0
         do ic = 1,nclass
            do l = 0,nl-1
               occc(l+1,n,:)=[(1,i=1,nindxc(l,ic))]
               ncwf(l+1,n,:)=[(1,i=1,nindxc(l,ic))]
            end do
         end do
         unoccc= 0
         write(stdo,*)' ### incwf=-3  occ=1 unocc=0 incwf=1 for all core '
      elseif( incwfx==-4 ) then
         write(stdo,*)' ### incwf=-4  occ=0 and unocc=0 for all core '
         occc=0
         unoccc=0
         ncwf=0
      elseif(incwfx==0) then
         write(stdo,*)' ### Use unocc occ ForX0 for core'
      else
         call rx( ' ### proper incwf is not given for genallcf_v3:rgwinf ')
      endif
      deallocate(ncwf2)
    endblock ReadProductBasis
    
    indexcoremto: block ! new nindx !-------------------------------------------- lmxa(ic) instead of nl-1
      integer:: nnn1(nclass),nnn2(nclass),nnn3(nclass),nnn4(nclass),nnn5(nclass),nnn6(nclass),nlx
      ndima=0
      do ic=1,natom
         ndima=ndima+sum([((2*l+1)*nindxv(l+1,iclass(ic)),l=0,lmxa(ic))]) !l=0,nl-1)])
      enddo
      nn  =  maxval(nindxv(1:nl,1:nclass)+nindxc(1:nl,1:nclass))
      allocate(nindx(nl,nclass),nocc(nl,nn,nclass),nunocc(nl,nn,nclass),source=0)
      reindxblock: block
        integer:: nval,ncore
        do    ic = 1,nclass
           do  l = 0,nl-1
              ncore  = nindxc(l+1,ic)
              nval   = nindxv(l+1,ic)
              nindx(l+1,ic)= ncore + nval
              if (ncore+nval > nn) call rx( 'reindx: ncore+nval > nn')
              nocc(l+1,1:,ic)   = [  (occc(l+1,n,ic),n=1,ncore),  (occv(l+1,n,ic),n=1,nval)]
              nunocc(l+1,1:,ic) = [(unoccc(l+1,n,ic),n=1,ncore),(unoccv(l+1,n,ic),n=1,nval)]
           enddo
        enddo
      endblock reindxblock
      block
        do ic=1,nclass
           nlx=lmxa(ic)+1
           nnn1(ic)=sum(nindx(1:nlx,ic))
           nnn2(ic)=sum([(nindx(l+1,ic)*(2*l+1),l=0,nlx-1)])
           nnn3(ic)=sum(  nindxv(1:nlx,ic))
           nnn4(ic)=sum([(nindxv(l+1,ic)*(2*l+1),l=0,nlx-1)])
           nnn5(ic)=sum(  nindxc(1:nlx,ic))
           nnn6(ic)=sum([(nindxc(l+1,ic)*(2*l+1),l=0,nlx-1)])
        enddo
        nlnx    = maxval(nnn1) 
        nlnmx   = maxval(nnn2)
        nlnxv   = maxval(nnn3)
        nlnmxv  = maxval(nnn4)
        nlnxc   = maxval(nnn5)
        nlnmxc  = maxval(nnn6)
      endblock
      allocate(il(nlnmx,nclass), in(nlnmx,nclass), im(nlnmx,nclass)) ! index for core and MTO basis =====================
      do ic = 1,nclass
         ind  = 0
         do l = 0,lmxa(ic) !nl-1 ! core
            do  n = 1,nindxc(l+1,ic)
               do  m = 1,2*l+1
                  ind       = ind + 1
                  if(ind > nlnmx) call rx( 'idxlnmc: ind > nlnmx')
                  lm        = l**2 + m
                  il(ind,ic)= l;   in(ind,ic) = n;  im(ind,ic)= m-l-1 !; ilnm(n,lm,ic) = ind
               enddo
            enddo
         enddo
         do  l = 0,lmxa(ic) !nl-1 ! valence
            ncorex  = nindxc(l+1,ic)
            do    n = 1,nindxv(l+1,ic)
               if (ncorex+n > nn) call rx( 'idxlnmc: ncore+n > nn')
               do      m = 1,2*l+1
                  ind = ind + 1
                  if (ind > nlnmx) call rx( 'idxlnmc: ind > nlnmx')
                  lm = l**2 + m
                  il(ind,ic)  =l;  in(ind,ic) = ncorex + n; im(ind,ic)  = m-l-1 !; ilnm(ncorex+n,lm,ic) = ind
               enddo
            enddo
         enddo
      enddo
      allocate(nlnmv(nclass),nlnmc(nclass),nlnm(nclass))
      do ic=1,nclass
         nlx=lmxa(ic)+1
         nlnmv(ic) = sum([(nindxv(l+1,ic)*(2*l+1),l=0,nlx-1)])
         nlnmc(ic) = sum([(nindxc(l+1,ic)*(2*l+1),l=0,nlx-1)])
         nlnm(ic)  = sum([(nindx(l+1,ic)*(2*l+1),l=0,nlx-1)])
      enddo
    endblock indexcoremto
    coreblock: block
      allocate(icore(nl**2*nnc,nclass),ncore(nclass),source=99999)
      do ic = 1,nclass ! index for allowed core states
         i  = 0
         j  = 0
         do       l = 0,nl-1
            do    n = 1,nindxc(l+1,ic)
               do m = -l,l
                  j = j + 1
                  if (ncwf(l+1,n,ic) == 1) then
                     i = i + 1
                     icore(i,ic)= j
                  endif
               enddo
            enddo
         enddo
         ncore(ic)  = i
      end do
      nctot = sum(ncore(iclass(1:natom)))
      open(newunit=ifec,file='ECORE')         ! core energies ==========================
      allocate(konf(nl,nclass),ecore(nctot,2))
      konf=0
      allocate(ecoret(0:nl-1,nnc,2,nclass))
      ecoret=0d0
      do ic = 1,nclass
         write(stdo,*) ' read ECORE : ic=',ic
         read (ifec,*)
         read (ifec,*)
         read (ifec,*)
         read (ifec,*) !zz,ic1,nr ,a,b,nsp
         read (ifec,*)
         read (ifec,*) (konf(l+1,ic),l=0,lmxa(ic)) !nl-1)
         read (ifec,*)
         do  l = 0,nl-1
            ncorex = konf(l+1,ic)-l-1
            if (ncorex > nnc) call rx( 'ECORE: wrong nnc')
            do n = 1,ncorex
               read (ifec,*) lt,nt,(ecoret(l,n,isp,ic),isp=1,nspin) !takao
               if(nspin==1) ecoret(l,n,2,ic) = ecoret(l,n,1,ic)   !     write(stdo,"(' read ecore=',3i4,2d13.5)")l,n,ic,ecoret(l,n,1:nspin,ic)
               if (lt /= l) call rx( 'rcore: wrong l')
               if (nt /= n) call rx( 'rcore: wrong n')
            end do
         end do
      end do
      close(ifec)
      i = 0
      do ia = 1,nclass
         ic  = iclass(ia)
         do l = 0,nl-1
            do n = 1,nnc
               do m = -l,l
                  if (ncwf(l+1,n,ic) == 1) then
                     i = i + 1
                     if (i > nctot) call rx( 'genalloc_mod: wrong nctot')
                     ecore(i,1:nspin) = ecoret(l,n,1:nspin,ic)
                     write(stdo,"(' ecore=',4i4,2d13.5)")i, l,n,ic,ecore(i,1:nspin)
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
!hbe
!    open(newunit=ifhbe, file='hbe.d', action='read')
!    read (ifhbe,*) nprecb,mrecb,mrece,ndima,nqbzt,nband,mrecg,nspc
!    close(ifhbe)
    ndimanspc=ndima*nspc
    nspx=nspin/nspc
    call cputid(0); write(stdo,*) 'genallcf_v3'
  end subroutine genallcf_v3
end module m_genallcf_v3

module m_ReadEfermi
  use m_lgunit,only:stdo
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
    write(stdo,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
  end subroutine readefermi
  subroutine readefermi_kbt()
    implicit none
    integer:: ifief_kbt
    open(newunit=ifief_kbt,file='EFERMI_kbt')
    read(ifief_kbt,*) ef_kbt,bandgap
    close(ifief_kbt)
    write(stdo,"(a,f12.6)")' --- READIN ef from EFERMI_kbt. ef=',ef_kbt
  end subroutine readefermi_kbt
end module m_ReadEfermi

