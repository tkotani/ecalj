!> -- Read HamiltionanMTOinfo and HamiltonianMTO. Then convert HamMTO to HamRsMTO  ------
module m_hamMTO
!  public :: ReadHamMTOInfo,resetMTOInfo,geninitmloinfo
!  private
  real(8),allocatable,protected:: plat(:,:),pos(:,:),qplist(:,:),qlat(:,:),qprs(:,:),weight(:)
  integer,allocatable,protected:: nlat(:,:,:,:),npair(:,:),ipivot(:) &
       ,ib_table(:),l_table(:),k_table(:),m_table(:),nqwgt(:,:,:),itor(:),ir(:,:),ilo(:,:),ib(:)
  integer,protected:: nkk1,nkk2,nkk3,nbas,nkp,npairmx,NMTO,jsp,lso,nsp,NS,nqs,niq,niqs,nlo,nlo2,nr,NM
  real(8),protected:: epsovl,alat,eferm
  complex(8),dimension(:,:,:),allocatable,protected:: ovlr,hamr
  character*20,protected,allocatable:: orbch(:,:),orbsym(:),orbsymlo(:)
contains

  subroutine ReadHamMTOInfo(iflg)!(infoname)
    !!  read information for crystal strucre, k points, neighbor pairs.
    implicit none
    !      character*200 ,intent(in) :: infoname
    integer:: ififft,i,lold,kold,m=999999,aold,access,ifnm,ifqpl,iflg
    character*4:: cccx
    logical*4 ::debug=.false.

    open(newunit=ifqpl,file='qplist.dat')
    open(ifqpl,file='qplist.dat')
    read(ifqpl,*) eferm
    close(ifqpl)

    !     open(newunit=ififft,file=trim(adjustl(infoname))//".Info",form='unformatted')
    if(iflg==0)open(newunit=ififft,file="HamiltonianPMTInfo",form='unformatted')
    if(iflg==1)open(newunit=ififft,file="initmlo.Info",form='unformatted')

    allocate(plat(3,3),qlat(3,3)) !plat primitive vectors, qlat:primitive in reciprocal space
    read(ififft) plat,nkk1,nkk2,nkk3,nbas,qlat
    nkp = nkk1*nkk2*nkk3
    allocate(qplist(3,nkp))
    allocate(pos(3,nbas))
    read(ififft) pos,alat  !atomic positions, unit of the scale.
    read(ififft) qplist    !qpoint list. all mesh points in the BZ mesh
    allocate(npair(nbas,nbas)) ! pair table of atoms corresponding to the mesh points
    read(ififft) npair,npairmx
    allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
    read(ififft) nlat,nqwgt
    read(ififft) NMTO,lso,nsp   ! size of Hamiltonian: MTO part

    allocate(ib_table(NMTO),l_table(NMTO),k_table(NMTO),m_table(NMTO))
    read(ififft)ib_table,l_table,k_table
    close(ififft)
    allocate(orbch(0:3,-3:3))
    orbch(0,0)="s" ; orbch(1,-1)="py" ; orbch(1,0)="pz" ; orbch(1,1)="px"
    orbch(2,-2)="dxy" ; orbch(2,-1)="dyz" ; orbch(2,0)="dz2" ; orbch(2,1)="dxz" ; orbch(2,2)="dx2-y2"
    orbch(3,-3)="fy(3x2-y2)" ; orbch(3,-2)="fxyz"      ; orbch(3,-1)="fy(5z2-r2)"
    orbch(3,0 )="fz(5z2-3r2)"  ; orbch(3,1 )="fx(5z2-r2)" ; orbch(3, 2)="fz(x2-y2)"; orbch(3,3)="fx(x2-3y2)"

    lold=-999
    kold=-999
    aold=0
    do i = 1,NMTO
       if(l_table(i)/= lold .OR. k_table(i)/=kold .OR. aold/=ib_table(i)) then !reset m of lm
          m=-l_table(i)
          lold=l_table(i)
          kold=k_table(i)
       else
          m=m+1
       endif
       m_table(i)=m
       if(aold /= ib_table(i))then
          aold=ib_table(i)
          if(debug)write(6,"(a,i4,a)") "[input info]======",aold,"atom =========="
       end if
       if(k_table(i)==3)then
          if(debug)write(6,"(i4,3x,a12,a,i2,a)") i,orbch(l_table(i),m_table(i)),"of EH=",-k_table(i),"LO"
       else
          if(debug)write(6,"(i4,3x,a12,a,i2)") i,orbch(l_table(i),m_table(i)),"of EH=",-k_table(i)
       end if
    enddo
    if(debug)write(6,*) " "
    return
  end subroutine ReadHamMTOInfo

  subroutine resetMTOInfo()
    if(allocated(plat)) deallocate(plat)
    if(allocated(qlat)) deallocate(qlat)
    if(allocated(nlat)) deallocate(nlat)
    if(allocated(qplist)) deallocate(qplist)
    if(allocated(pos)) deallocate(pos)
    if(allocated(npair)) deallocate(npair)
    if(allocated(nqwgt)) deallocate(nqwgt)
    if(allocated(ib_table)) deallocate(ib_table)
    if(allocated(l_table)) deallocate(l_table)
    if(allocated(k_table)) deallocate(k_table)
    if(allocated(m_table)) deallocate(m_table)
    if(allocated(orbch)) deallocate(orbch)
    return
  end subroutine resetMTOInfo

  subroutine geninitmloinfo() !need ReadHamMTOInfo in advance!
    implicit none
    integer:: ififft,i,access,if2,ifind,lold,kold,aold,m,ifis
    character*4:: cccx

    open(newunit=if2,file="initmlo.Info",form='unformatted')
    write(if2) plat,nkk1,nkk2,nkk3,nbas,qlat
    nkp = nkk1*nkk2*nkk3
    write(if2) pos,alat  !atomic positions, unit of the scale.
    write(if2) qplist    !qpoint list. all mesh points in the BZ mesh
    write(if2) npair,npairmx
    write(if2) nlat,nqwgt
    write(if2) NMTO,lso,nsp   ! size of Hamiltonian: MTO part
    write(if2) ib_table,l_table,k_table ! changed via ib_new etc.
    close(if2)

    lold=-999
    kold=-999
    aold=0
    do i = 1,NMTO
       if(aold /= ib_table(i))then
          aold=ib_table(i)
          write(6,"(a,i4,a)") "[output info]======",aold,"atom========"
       end if
       if(k_table(i)==3)then
          write(6,"(i4,3x,a12,a,i2,a)") i,orbch(l_table(i),m_table(i)),"of EH=",-k_table(i),"LO"
       else
          write(6,"(i4,3x,a12,a,i2)") i,orbch(l_table(i),m_table(i)),"of EH=",-k_table(i)
       end if
    enddo
    write(6,*) " "
    return
  end subroutine geninitmloinfo


  subroutine ibzweight() !need ReadHamMTOInfo in advance!
    implicit none
    integer :: ifqpl,ifqpl2,ifibz,iflog,iq,i,j,jsp,ldum
    integer,allocatable:: itori(:)
    real(8) :: wdum,lend
    real(8),allocatable:: qpi(:,:)

    nqs=nkp*nsp
    allocate(qprs(3,nqs))
    open(newunit=ifqpl,file="QPLIST.RBZ")
    do i=1,nkp
       do jsp=1,nsp
          write(ifqpl,*) i,qplist(1:3,i),1.0d0/dble(nkp),jsp
          qprs(1:3,nsp*(i-1)+jsp)=qplist(1:3,i)
       end do
    end do
    close(ifqpl)

    open(newunit=ifqpl2,file="QPLIST.RBZ.spin")
    do i=1,nqs
       write(ifqpl2,*) i,qprs(1:3,i)
    end do
    close(ifqpl2)

    open(newunit=ifibz,file="QPLIST.IBZ")
    niq=0
    do
       read(ifibz,*,end=1999)
       niq=niq+1
    end do
1999 close(ifibz)
    niqs=niq*nsp
    allocate(qpi(3,niq),weight(niqs),itor(niqs),itori(niqs))
    write(*,*) "niq(irred. k-points) and niq*spin=",niq,niqs
    open(newunit=ifibz,file="QPLIST.IBZ")
    do i=1,niq
       read(ifibz,*) ldum,qpi(1:3,i),wdum
       !         write(*,*) "iq=",i,qpi(1:3,i)
       do j=1,nkp
          lend=(qplist(1,j)-qpi(1,i))**2d0+(qplist(2,j)-qpi(2,i))**2d0+(qplist(3,j)-qpi(3,i))**2d0
          !            write(*,*) "nq=",j,qplist(1:3,j),lend
          if(lend <= 1d-10)then
             do jsp=1,nsp
                itor(nsp*(i-1)+jsp)=nsp*(j-1)+jsp
                itori(nsp*(i-1)+jsp)=j
                weight(nsp*(i-1)+jsp)=wdum
             end do
          end if
       end do
    end do
    close(ifibz)

    open(newunit=iflog,file="IBZvsRBZandWEIGHT.dat")
    write(iflog,*) "ikp,ikpx,ikpy,ikpz,rkp,rkpx,rkpy,rkpz,weight"
    wdum=0d0
    do iq=1,niqs
       wdum=wdum+weight(iq)
       write(iflog,"(i4,3f10.5,a,2i4,3f10.5,a,f10.5)") &
            (iq+1)/2,qprs(1:3,itor(iq)),":",itor(iq),itori(iq),qpi(1:3,(iq+1)/2),":",weight(iq)
    end do
    close(iflog)
    if(abs(wdum/dble(nsp)-1.0d0) > 1d-10) stop "error! sum of weight is NOT 1"
  end subroutine ibzweight

  subroutine mkpivot()   !for initmlo.F
    implicit none
    integer::lold,i,j,n,if1,idum,access,ifi,ii,jj,iflg,ib1,l1,k1,m1,ib2,l2,k2,m2,ip,jp,jpf &
         ,di,dib,dk,dm,look(1000,3,16),ifim,len,ifo,lk,lk2,ioir,nmodel,outmodel,nlo,nk2,nlosub=0
    integer(4) ,allocatable :: ib_new(:),l_new(:),k_new(:),m_new(:) &
         ,ipnk2(:),irold(:,:),iro(:,:),ipLO(:),irl(:)!,ilo(:)
    character(2) add,a(23)
    character(10) tgs
    character(20) ,allocatable :: orbold(:)
    character(200) cdum

    open(newunit=ifim,file="meta_mto") !generated in mkinput.F
    look=0
    i=0
    do
       read(ifim,fmt="(4i4)",end=666) di,dib,dk,dm
       look(dib,dk,dm)=di
       i=i+1
    end do
666 close(ifim)
    len=i  !total number of (first-principles) MTO orbitals
    if(len /= NMTO) stop "error! inconsistent length of meta_mto read in m_hamMTO.F!"


    if(access("mtoctrl","") /= 0) stop "error! no mtoctrl"
    open(newunit=ifi,file="mtoctrl")

    do
       read(ifi,"(a10)") tgs
       !         write(*,*) tgs(1:1)
       if(tgs=="<mkhamloc>") exit
    end do
    do i=1,5
       read(ifi,*)
    end do

    open(newunit=ifo,file="status-model") !new version

    allocate(ipivot(1:NMTO))
    i=0
    do
       read(ifi,fmt="(a)",end=999) cdum
       if(cdum(1:1)==" " .OR. cdum(1:1)=="<" .OR. cdum(1:2)==" <") goto 999
       a="vo" !vo means "void"
       read(cdum,*,end=1919) (a(j),j=1,23)
1919   continue

       if(a(2) /= ":" .OR. a(4) /= ":" .OR. a(6) /= ":") stop "error! division by ':' is needed!"
       read(a(3),*) ib1
       read(a(5),*) k1

       do n=7,23
          if(a(n)/="vo")then
             read(a(n),*) m1 !transformation from character to integer!
             i=i+1
             ipivot(i)=look(ib1,k1,m1)
             lk2=look(ib1,1,m1) !sanity
             if(k1==2 .AND. lk2==0) stop "error! EH2 need corresponding EH1!" !sanity
             write(ifo,fmt="(i4)") ipivot(i)
          end if
          if(a(n)=="vo" .AND. k1==3)then
             nlosub=nlosub+(n-7)
             exit
          end if
       end do
    end do
999 close(ifi)
    close(ifo)
    nmodel=i
    write(*,*) "number of model orbital", nmodel
    !     stop "aaaaa"

    ! emproary ipivot
    outmodel=0
    do i=1,NMTO
       jpf=0
       do j=1,nmodel
          if(ipivot(j)==i) jpf=j
       end do
       if(jpf==0)then
          outmodel=outmodel+1
          ipivot(nmodel+outmodel)=i
       end if
    end do
    write(*,*) outmodel,NMTO-nmodel
    if(outmodel /= NMTO-nmodel)stop "error! wrong count for outmodel!(maybe duplicated orbitals in mtoctrl)"


    allocate(ipnk2(outmodel))
    nk2=0  !the number of EH2 basis, out of the model
    do i=nmodel+1,NMTO  !search nk2
       ip=ipivot(i)
       ib1=ib_table(ip)
       k1=k_table(ip)
       l1=l_table(ip)
       m1=m_table(ip)
       do j=1,nmodel
          jp=ipivot(j)
          ib2=ib_table(jp)    !in the model
          k2=k_table(jp)      !model
          l2=l_table(jp)   !model
          m2=m_table(jp)   !model
          if(ib1==ib2 .AND. l1==l2 .AND. m1==m2 .AND. k1 /= k2 .AND. k1 /= 3 .AND. k2 /= 3)then
             !               write(*,*) "unk",i,ip,j,jp
             nk2=nk2+1
             ipnk2(nk2)=ip
          end if
       end do
    end do
    do i=1,nk2
       ipivot(nmodel+i)=ipnk2(i)
    end do
    deallocate(ipnk2)

    jj=0 !the other orbitals(out of model, just trunctated)
    do j=1,NMTO
       iflg=0
       do i=1,nmodel+nk2!nr
          if(ipivot(i)==j) iflg=1
       end do
       if(iflg==0)then
          jj=jj+1
          !            ipivot(nmodel+nr+jj)=j
          ipivot(nmodel+nk2+jj)=j
       end if
    end do

    open(newunit=if1,file="IPVT") !just a record
    do i=1,NMTO
       write(if1,*) ipivot(i)
    end do

    allocate(ib_new(1:NMTO),l_new(1:NMTO),k_new(1:NMTO),m_new(1:NMTO))
    do i=1,NMTO
       ib_new(i)= ib_table(ipivot(i))
       l_new(i) =  l_table(ipivot(i))
       k_new(i) =  k_table(ipivot(i))
       m_new(i) =  m_table(ipivot(i))
    end do

    ib_table(:)=ib_new(:)
    l_table(:)=l_new(:)
    k_table(:)=k_new(:)
    m_table(:)=m_new(:)
    deallocate(ib_new,l_new,k_new,m_new)

    ! cccccccccccccc after the pivot turn ccccccccccccccccccccccc
    open(newunit=ifim,file="meta_pivot") !k,l,m index after pivoting
    do i=1,NMTO
       write(ifim,"(i4,a,3i4)") i,":",ib_table(i),k_table(i),1+(l_table(i))**2+l_table(i)+m_table(i)
    end do
    close(ifim)

    if(nk2==0)then
       write(*,*) "no orbital is truncated(nk2==0)."
       open(newunit=ioir,file="ir_list") !preparation for mkhamloc.F
       write(ioir,*) nmodel,nk2,0
       close(ioir)
       return
    end if

    allocate(orbsym(nk2),ir(nk2,3),irl(nk2))
    nlo=0
    ii=0
    ir=0
    do i=1,nmodel  !in the model space
       ib1=ib_table(i)
       k1=k_table(i)
       l1=l_table(i)
       m1=m_table(i)
       do j=nmodel+1,NMTO !out of the model space
          ib2=ib_table(j)
          k2=k_table(j)
          l2=l_table(j)
          m2=m_table(j)
          if(ib1==ib2 .AND. l1==l2 .AND. m1==m2 .AND. k1==1 .AND. k2==2)then
             ii=ii+1
             !               write(*,*) "b",ib2,k2,l2,m2,ib1,k1,l1,m1
             ir(ii,1)=i
             ir(ii,2)=j
             irl(ii)=l1
             orbsym(ii)=orbch(l2,m2)
          end if
          !            if(ib1==ib2.and.l1==l2.and.m1==m2.and.k1.ne.k2.and.k1.eq.3.and.k2.ne.3)then
          !               nlo=nlo+1
          !               ir(nlo,3)=i
          !            end if !under construction!
       end do
    end do
    if(ii /= nk2) stop "error! pivoting failed!" !sanity check

    open(newunit=ioir,file="ir_list.org") !preparation for mkhamloc.F
    write(ioir,"(3i6,a)") nmodel,nk2,nlo,"   size of matrix in mkhamloc and model space"
    write(*,*) "===========ir_list(original)==========="
    write(*,*) "i,ir(i,1),ir(i,2)"
    do i=1,nk2
       write(*,"(3i6,3x,a10,3i4)") ir(i,1),ir(i,2),ir(i,3),orbsym(i)
       !         if(k_table(ipivot(ir(i,2))).ne.3)then
       write(ioir,"(3i6,3x,a10)") ir(i,1),ir(i,2),ir(i,3),trim(adjustl(orbsym(i)))
       !         else !under construction
       !            write(ioir,"(3i6,3x,a10,a4)") ir(i,1),ir(i,2),ir(i,3),trim(adjustl(orbsym(i))),"(LO)"
       !         end if
    end do
    write(*,*) "===========ir_list(original)==========="
    close(ioir)

    allocate(irold(nk2,3),orbold(nk2))
    irold=ir
    orbold=orbsym

    !      ii=0
    !      do j=1,nr
    !         if(k_table(iro(j,2))==3)then
    !            ii=ii+1
    !            ir(ii,1)=irold(j,1)
    !            ir(ii,2)=irold(j,2)
    !            orbsym(ii)=orbold(j)
    !         end if
    !      end do

    ir=0d0
    ii=0
    do i=3,0,-1
       do j=1,nk2
          if(irl(j)==i)then
             ii=ii+1
             ir(ii,1)=irold(j,1)
             ir(ii,2)=irold(j,2)
             if(irold(j,3) /= 0) ir(ii,3)=irold(j,3)
             !               ir(ii,3)=irold(j,3)
             orbsym(ii)=orbold(j)
          end if
       end do
    end do


    open(newunit=ioir,file="ir_list") !preparation for mkhamloc.F
    write(*,*) "===========ir_list(practical)==========="
    write(ioir,"(3i6,a)") nmodel,nk2,nlo,"   size of matrix in mkhamloc, model space, lo"
    write(*,*) "i,ir(i,1),ir(i,2),ir(i,3)"
    do i=1,nk2
       write(*,"(3i6,3x,a10)") ir(i,1),ir(i,2),ir(i,3),orbsym(i)
       !         if(k_table(ipivot(ir(i,2))).ne.3)then
       write(ioir,"(3i6,3x,a10,a4)") ir(i,1),ir(i,2),ir(i,3),trim(adjustl(orbsym(i))),"    "
       !         else
       !            write(ioir,"(3i6,3x,a10,a4)") ir(i,1),ir(i,2),ir(i,3),trim(adjustl(orbsym(i))),"(LO)"
       !         end if
    end do
    write(*,*) "===========ir_list(practical)==========="
    close(ioir)

    return
  end subroutine mkpivot

  subroutine orbavg(rvec,jsp)
    implicit none
    integer,intent(in)::jsp
    real(8),intent(inout)::rvec(1:nr+nlo2,nsp)
    integer :: i,natm,at,ns,np,nd
    real(8) sv,pv,dv

    natm=maxval(ib_table)

    do at=1,natm
       ns=0   ; np=0   ; nd=0
       sv=0d0 ; pv=0d0 ; dv=0d0
       do i=1,nr
          if(at==ib_table(i))then
             if(l_table(ir(i,1))==0)then
                ns=ns+1
                sv=sv+rvec(i,jsp)
             else if(l_table(ir(i,1))==1)then
                np=np+1
                pv=pv+rvec(i,jsp)
             else if(l_table(ir(i,1))==2)then
                nd=nd+1
                dv=dv+rvec(i,jsp)
             end if
          end if
       end do
       do i=1,nr
          if(at==ib_table(i))then
             if(l_table(ir(i,1))==0)then
                rvec(i,jsp)=sv/ns
             else if(l_table(ir(i,1))==1)then
                rvec(i,jsp)=pv/np
             else if(l_table(ir(i,1))==2)then
                rvec(i,jsp)=dv/nd
             end if
          end if
       end do
    end do

  end subroutine orbavg

  !! read RealSpace MTO Hamiltonian
  subroutine mkRsMTO(fname,jsp)
    implicit none
    character*200,intent(in):: fname
    integer,intent(in)   ::jsp
    integer:: ifih,i,j,n,it,lso,ikp,jspin
    ! cccccccc
    logical:: debug=.false.,lprint=.true.,savez=.false.,getz=.false.
    integer(4) :: nev,ifig=-999
    real(8), allocatable:: evl(:)
    complex(8), allocatable:: c_C(:,:)
    ! cccccccccccc
    real(8)::qp(3),pi=4d0*atan(1d0)
    complex(8) :: phase,img=(0d0,1d0)
    complex(8),allocatable:: ovlM(:,:),hamM(:,:)

    open(newunit=ifih,file =trim(adjustl(fname)),form='unformatted')

    if(allocated(ovlr)) deallocate(ovlr)
    if(allocated(hamr)) deallocate(hamr)
    if( .NOT. allocated(ib)) allocate(ib(1:NMTO))

    do i=1,NMTO
       ib(i) = mod(ib_table(i),NMTO)
    end do
    write(6,*) " "
    write(6,*) "********** Making Real Space Hamiltonian(FT from k-space) ************"

    ikp=0
    do
       read(ifih,end=2019) qp,NS,lso,epsovl,jspin !NM and epsovl are module variable
       if( .NOT. allocated(ovlr) .AND. .NOT. allocated(hamr))then
          allocate(ovlr(1:NS,1:NS,1:npairmx),hamr(1:NS,1:NS,1:npairmx))
          hamr=0d0 ; ovlr=0d0
       end if

       allocate(ovlM(1:NS,1:NS),hamM(1:NS,1:NS))
       read(ifih) ovlM(1:NS,1:NS)
       read(ifih) hamM(1:NS,1:NS)
       if(jspin /= jsp)then
          deallocate(ovlM,hamM)
          cycle
       end if
       ikp=ikp+1
       write(6,"(i4)",advance="no") ikp
       if(mod(ikp,10)==0) write(6,"(a,i4,a)") "  : ",idnint((dble(ikp)/dble(nkp)*1d2))," %"

       !$$$         if(debug)then
       !$$$            allocate(evl(1:NS),c_C(1:NS,1:NS))
       !$$$            if(epsovl<1.000001d-14) then !epsovl is the trancation to remove poor linear-dependency basis
       !$$$               call zhev_tk2(NS,hamM,ovlM,NS,nev,evl,c_C,lprint,savez,getz,ifig)
       !$$$            else
       !$$$               call zhev_tk3(NS,hamM,ovlM,NS,nev,evl,c_C,lprint,savez,getz,ifig,epsovl)
       !$$$            endif
       !$$$            do i=1,nev
       !$$$               write(7777,"(2i4,4f13.7)") ikp,jsp,evl(i)
       !$$$            end do
       !$$$            write(7777,*) " "
       !$$$            deallocate(evl,c_C)
       !$$$         end if

       !         write(*,*) "jspin,jsp=",jspin,jsp
       do i=1,NS
          do j=1,NS
             do it =1,npair(ib(i),ib(j)) !hammr_ij (T)= \sum_k hamm(k) exp(ikT).
                phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib(i),ib(j)))))
                hamr(i,j,it)= hamr(i,j,it) + hamM(i,j)*phase
                ovlr(i,j,it)= ovlr(i,j,it) + ovlM(i,j)*phase
                !                  write(*,*) i,ib(i),j,ib(j),it,"b"
             enddo
          enddo
       enddo
       deallocate(ovlM,hamM)
    enddo
2019 close(ifih)
    return
  end subroutine mkRsMTO

  subroutine startprint(inp)
    integer(4),intent(in) :: inp
    character(200) :: chr
    if(inp==0)then
       chr="lmfham"
    else if(inp==1)then
       chr="mkhamloc"
    else if(inp==2)then
       chr="initmlo"
    else if(inp==3)then
       chr="plotmodel"
    else if(inp==4)then
       chr="mkmodel"
    else if(inp==-1)then
       chr="mkLO"
    end if

    write(6,*) "============================================================"
    write(6,*) "                        START ",trim(adjustl(chr))
    write(6,*) "============================================================"
    return
  end subroutine startprint

end module m_HamMTO

