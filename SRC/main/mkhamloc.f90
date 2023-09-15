! !     include "mkhamloc_inc.F"

! subroutine mkdiff(jp,ilo,neh2,nlo,ir,vec,M,L,Omat,Hmat,etr,dth)
!   implicit none
!   intent(in)  :: jp,ilo,neh2,nlo,ir,vec,M,L,Omat,Hmat
!   intent(out) :: dth,etr
!   logical     :: dbg,debug=.false.
!   integer ::  jp,ilo,neh2,nlo,ir(neh2,3),M,L,j1,j2,j3,jj,i,i1,i2,i3,ME
!   real(8)     :: vec(neh2,2),etr
!   complex(8)  :: Omat(L,L),Hmat(L,L),dth(1:7,1:2)
!   complex(8),allocatable,dimension(:,:) ::Ot,Ht,Oinv,rmat,OdL,OL,HL,emat
!   ME=M
!   allocate(Ot(ME,ME),Ht(ME,ME),Oinv(ME,ME),rmat(L,ME),OdL(ME,L),OL(ME,L),HL(ME,L),emat(ME,ME))
!   rmat=0d0
!   do i=1,ME
!      rmat(i,i)=1d0
!   end do

!   do i=1,neh2
!      i1=ir(i,1)
!      i2=ir(i,2)
!      i3=ir(i,3)
!      !         write(*,*) jp,i,i1,i2,i3
!      !         write(*,*) dcos(vec(i,1)),dcos(vec(i,2))
!      !         if(i3==0) write(*,*)  dsin(vec(i,2))
!      if(i3==0)then
!         if(dcos(vec(i,2)) /= 1d0) stop "error! csn"
!      end if
!      rmat(i1,i1)=dcos(vec(i,1))*dcos(vec(i,2))
!      rmat(i2,i1)=dsin(vec(i,1))*dcos(vec(i,2))
!      if(i3 /= 0) rmat(i3,i1)=  -dsin(vec(i,2))
!   end do

!   OL=matmul(transpose(rmat),Omat)
!   HL=matmul(transpose(rmat),Hmat)

!   Ot=matmul(OL,rmat)
!   Ht=matmul(HL,rmat)

!   call invs(M,Ot,Oinv) !Oinv is matrix inversion of Ot

!   emat=matmul(Oinv,Ht)
!   etr=0d0
!   do i=1,ME
!      etr=etr+dble(emat(i,i))
!   end do

!   OdL=matmul(Oinv,HL)-matmul(matmul(Oinv,Ht),matmul(Oinv,OL))

!   if(ilo==-1)then
!      j1=ir(jp,1)
!      j2=ir(jp,2)
!      dth(1,1)=OdL(j1,j1)*(-dsin(vec(jp,1)))*(dcos(vec(jp,2))) &
!           +OdL(j1,j2)*  dcos(vec(jp,1)) *(dcos(vec(jp,2)))
!   else if(ilo==-2)then
!      j1=ir(jp,1)
!      j2=ir(jp,2)
!      j3=ir(jp,3)
!      dth(1,1)=OdL(j1,j1)*dcos(vec(jp,1))*(-dsin(vec(jp,2))) &
!           +OdL(j1,j2)*dsin(vec(jp,1))*(-dsin(vec(jp,2))) &
!           +OdL(j1,j3)*                (-dcos(vec(jp,2)))
!   else if(ilo > 0)then
!      !         write(*,*) jp,ilo
!      do i=jp,ilo            !=jskip
!         j1=ir(i,1)
!         j2=ir(i,2)
!         j3=ir(i,3)
!         if(j3 <= 0) stop "error! 59"
!         dth(i,1)=OdL(j1,j1)*(-dsin(vec(i,1)))*dcos(vec(i,2)) &
!              +OdL(j1,j2)*  dcos(vec(i,1)) *dcos(vec(i,2))

!         dth(i,2)=OdL(j1,j1)*  dcos(vec(i,1)) *(-dsin(vec(i,2))) &
!              +OdL(j1,j2)*  dsin(vec(i,1)) *(-dsin(vec(i,2))) &
!              +OdL(j1,j3)*                  (-dcos(vec(i,2)))

!      end do
!   end if

!   return
! end subroutine mkdiff

! subroutine invs(len,mat,inv)
!   implicit none
!   intent(in)  :: len,mat
!   intent(out) :: inv
!   integer ::  len,ipv(len),info
!   complex(8)  :: mat(len,len),inv(len,len),work(len*3)
!   inv=mat
!   call zgetrf(len, len, inv, len, ipv(1:len), info)
!   call zgetri(len, inv, len, ipv(1:len), work, len*3, info)
!   return
! end subroutine invs

! program mkhamloc
!   use m_hamMTO,only:ReadHamMTOInfo,qplist,npair,ib_table,l_table,k_table,m_table,nsp,nkp,NMTO,lso,orbch, &
!        ibzweight,nqs,qprs,weight,itor,niqs,eferm,orbavg,startprint
!   implicit none
!   logical:: debug=.false.,lprint=.true.,savez=.false.,getz=.false.,dbg=.false.,LLmode=.false.
!   logical,allocatable::eflg(:),loflg(:)

!   integer(4):: ifmto,ifqpl,ifinput,ifout,i,j,jj,n,m,t,l,k,access,loop,nq,Lmax,i1,i2,i3,LD,imax,itr &
!        ,ikpd,ikp,ifih,it,iq,rq,lold,kold,ibold,ifig=-999,ii,id,jd,fflg,LLdum,idum,ifo2 &
!        ,jsp,nn,LL,LE,lso1,flg_osym,flg_ssym,lenMTO,nmax,nmin,nevd,j1,j2,j3,ioncut,jorb,lenlo,jskip &
!        ,ifir,NR,NS,NM,nlo,neh2,NL,jp,ilo,ifjorb,njorb,iffix,ifmax,lcut,rcut,nbsec,maxitr,jdum
!   integer(4),allocatable:: nev(:),jspl(:),isptab(:),llist(:),ir(:,:),ilorb(:),jlist(:,:),llen(:)

!   real(8)::qp(3),convc,h,ral,rimg,ryd,emax,emin,pi,epsdum,iden,fac,prdN,prdF,prdM,norm1,norm2 &
!        ,dum,r,r2,etr,etr0,et,tmp,adum,fD,dth,dthL,dthR,thmin &
!        , etrL,etrR,thL,thR,thM,thmax,cbsec,absd,theta,thLB(2),thUB(2),LB,UB
!   real(8),allocatable,dimension(:)    :: epsovl,evl_D,elist,dlist,thlist
!   real(8),allocatable,dimension(:,:)  :: vM,vL,vR,vec,diff
!   real(8),allocatable,dimension(:,:):: totvec,oldvec

!   complex(8),allocatable,dimension(:,:)   :: hamF,ovlF,hamD,ovlD,c_D,lmat1,lmat2,dmat,bdb,rmatout
!   complex(8),allocatable,dimension(:,:,:) :: ovlM,hamM,hamE,ovlE,damM,damE
!   complex(8) :: dt(1:7,1:2)

!   character*40,allocatable:: orbsym(:)
!   character*10:: tgs
!   character(200) :: iname,indexname,ifc

!   pi=4d0*atan(1d0)
!   ryd=13.605693d0
!   !      thLB(1)=-pi*0.5d0
!   !      thUB(1)=pi*0.5d0
!   !      thLB(2)=0d0       !lo
!   !     thUB(2)=pi !lo
!   LB=-pi*0.5d0
!   UB=pi*0.5d0


!   open(newunit=ifmto,file="mtoctrl")
!   do
!      read(ifmto,"(a10)") tgs
!      if(tgs=="<mkhamloc>") exit
!   end do
!   !      read(ifmto,*) convc,emin,emax,flg_osym,flg_ssym,maxitr,nbsec,cbsec,lcut,rcut
!   !     write(*,*) convc,emin,emax,flg_osym,flg_ssym,maxitr,nbsec,cbsec,lcut,rcut
!   read(ifmto,*) convc,cbsec
!   read(ifmto,*) emin,emax
!   read(ifmto,*) flg_osym,flg_ssym
!   read(ifmto,*) maxitr,nbsec
!   read(ifmto,*) lcut,rcut

!   close(ifmto)

!   emin=emin/ryd+eferm
!   emax=emax/ryd+eferm

!   call startprint(1)
!   !      call ReadHamMTOInfo(iname)
!   call ReadHamMTOInfo(1) !initmlo.Info
!   allocate(lmat1(NMTO,NMTO),lmat2(NMTO,NMTO))
!   LL=NMTO

!   open(newunit=ifir,file="ir_list")
!   read(ifir,*) NM,neh2,nlo

!   if(neh2==0 .AND. nlo==0)then
!      open(newunit=ifout,file="rotmat")
!      do i=1,NMTO
!         do j=1,NMTO
!            if(i==j .AND. i <= NM)then
!               write(ifout,"(2i4,2x,2(f12.5,x))") i,j,1d0,0d0
!            else
!               write(ifout,"(2i4,2x,2(f12.5,x))") i,j,0d0,0d0
!            end if
!         end do
!      end do
!      close(ifout)
!      stop "return 0(no trunctation)"
!   end if
!   NR=neh2+nlo
!   NL=NR+NM
!   write(*,*) "N(model),NL(model+EH2+lo)=",NM,NL
!   write(*,*) "NR(EH2+lo),N(EH2),N(lo)=",NR,neh2,nlo

!   if(nr <= 0) stop "error! wrong size of nr in ir_list!"
!   allocate(ir(neh2,3),orbsym(nr))
!   write(*,*) "=========ir list==========="
!   do i=1,neh2
!      read(ifir,*) ir(i,1),ir(i,2),ir(i,3),orbsym(i)
!      write(*,*) ir(i,1),ir(i,2),ir(i,3),orbsym(i)
!   end do
!   write(*,*) "=========ir list==========="
!   close(ifir)

!   !      do i=1,neh2
!   !         write(*,*) i,ir(i,1),ir(i,2),ir(i,3)
!   !      end do

!   call ibzweight()

!   ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   write(6,*) "====================open HamI file! ===================="
!   allocate(jspl(nqs),llist(nqs))
!   open(newunit=ifih,file="HamI",form='unformatted')
!   i=0 ; j=0; rq=0           !sanity check
!   do
!      read(ifih,end=2019) qp,LD,lso1,epsdum,jsp
!      if(lso /= lso1) stop "error! inconsistent lso!"
!      if(jsp==1) i=i+1       !for check
!      rq=rq+1
!      llist(rq)=LD
!      jspl(rq)=jsp
!      if((qp(1)-qprs(1,rq))**2d0+(qp(2)-qprs(2,rq))**2d0+(qp(3)-qprs(3,rq))**2d0 > 1d-10) stop "error! inconsistent qp!"
!      allocate(ovlF(1:LD,1:LD),hamF(1:LD,1:LD))
!      read(ifih) ovlF(1:LD,1:LD)
!      read(ifih) hamF(1:LD,1:LD)
!      deallocate(ovlF,hamF)
!   end do
! 2019 close(ifih)
!   if(debug) write(*,*) "**************sanity check is ok************************"
!   if(i /= nkp .OR. rq /= nqs)then
!      write(*,*) i,nq,rq,nqs
!      stop "inconsistent nq between initmlo.Info and HamI"
!   end if

!   write(6,*)'Read: total # of q for Ham=',nkp,i,"(from .Info and .K)"
!   write(6,*) "====================close HamI file! ===================="
!   open(newunit=ifih,file="HamI",form='unformatted')

!   Lmax=maxval(llist)
!   if(LLmode)then
!      ! list(:)=NMTO
!      llist(:)=NL
!      Lmax=NL
!   end if
!   allocate(ovlM(1:Lmax,1:Lmax,1:nqs),hamM(1:Lmax,1:Lmax,1:nqs),damM(1:Lmax,1:Lmax,1:nqs))
!   ! cccccccccccccccccccc diagonalization process! ccccccccccccccccccccccccccc
!   allocate(nev(nqs),isptab(niqs))

!   damM=0d0
!   open(newunit=ioncut,file="ncut.kpoint")
!   do rq=1,nqs
!      read(ifih) qp,LD,lso1,epsdum,jsp

!      allocate(ovlF(1:LD,1:LD),hamF(1:LD,1:LD))
!      read(ifih) ovlF(1:LD,1:LD)
!      read(ifih) hamF(1:LD,1:LD)

!      if(LLmode) LD=NL

!      ovlM(1:LD,1:LD,rq)=ovlF(1:LD,1:LD)
!      hamM(1:LD,1:LD,rq)=hamF(1:LD,1:LD)

!      allocate(hamD(LD,LD),ovlD(LD,LD),c_D(LD,LD),evl_D(LD))
!      ovlD(1:LD,1:LD)=ovlF(1:LD,1:LD)
!      hamD(1:LD,1:LD)=hamF(1:LD,1:LD)
!      call zhev_tk4(LD,hamD,ovlD,LD,nev(rq),evl_D,c_D,epsdum)

!      nmax=1
!      do i=1,nev(rq)
!         if(evl_D(i) > emax) exit
!         nmax=nmax+1
!      end do
!      nmin=1
!      do i=1,nev(rq)
!         if(evl_D(i) > emin) exit
!         nmin=nmin+1
!      end do

!      write(ioncut,"(a,2i5,3f15.5)") "rq,nmin,emin,maxeigen=",rq,nmin,emin*ryd,(evl_D(nmin)-eferm)*ryd
!      write(ioncut,"(a,2i5,3f15.5)") "rq,nmax,emax,maxeigen=",rq,nmax,emax*ryd,(evl_D(nmax)-eferm)*ryd

!      allocate(bdb(1:LD,1:LD))
!      bdb=matmul(c_D(1:LD,nmin:nmax),transpose(dconjg(c_D(1:LD,nmin:nmax))))
!      damM(1:LD,1:LD,rq)=matmul(ovlD,matmul(bdb,ovlD))
!      !         damM(1:LD,1:LD,rq)=-hamM(1:LD,1:LD,rq)
!      deallocate(ovlF,hamF,ovlD,hamD,c_D,evl_D,bdb)
!   end do
!   close(ioncut)
!   ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!   allocate(vM(neh2,2),vR(neh2,2),vL(neh2,2),totvec(neh2,2),oldvec(neh2,2))
!   write(6,*) "====================start main loop of mkhamloc!================"
!   open(newunit=idum,file="mkhamloc.log")
!   !      totvec=-asin(1.0d0)
!   !      totvec=-0.8d0
!   totvec=0d0

!   totvec(nlo+1:neh2,2)=0d0  !reset for local orbital
!   oldvec=totvec
!   if(access("jorblist","")==0)then
!      open(newunit=ifjorb,file="jorblist")
!      read(ifjorb,*) njorb
!      if(njorb <= 0) goto 2626
!   allocate(vec(neh2,2))
!   allocate(llen(njorb))
!   do i=1,njorb
!      read(ifjorb,*) llen(i)
!      write(*,*)  "llen(i)=",i,llen(i)
!   end do

! !  c         if(llen==1) 
! !  c         read(ifjorb,*) jlist(i,1),jlist(i,2),jlist(i,3),jlist(i,4),jlist(i,5) !(5,4)=local orbital(1) or not(0)
! !  c        if(jlist(i,5).gt.2.or.jlist(i,5).le.0) stop "error in jlist(i,5)"
! !  c      end do
! !  c      close(ifjorb)       

!   write(*,*) "[map.F]njorb,NMTO,NM,neh2,nlo,",njorb,NMTO,NM,neh2,nlo
!   open(newunit=ifmax,file="emax")
!   lcut=200
!   allocate(elist(0:lcut),dlist(0:lcut),thlist(0:lcut))
!   do n=1,njorb
!      if(allocated(jlist)) deallocate(jlist)
!      allocate(jlist(llen(n),2))
!      do j=1,llen(n)
!         read(ifjorb,*)  jlist(j,1), jlist(j,2)
!         write(*,*)"n,j,jlist(1),jlist(2=lo)",n,j,jlist(j,1), jlist(j,2)
!      end do
!      jp =jlist(1,1) !first orbital component
!      ilo=jlist(1,2)
!      j1=ir(jp,1)
!      j2=ir(jp,2)
!      j3=ir(jp,3)
!      write(*,*) "jpj1j2j3",jp,j1,j2,j3
!      emax=-1d0
!      emin=1000d0         
!      do loop=0,lcut
!         open(newunit=iffix,file="totvec1")
! !        c            write(*,*) "totvec1"
!         do i=1,neh2         
!            read(iffix,*) vec(i,1),vec(i,2)
! !           c              write(*,*) i,vec(i,1),vec(i,2)
!         end do
!         close(iffix)       
! !        c            theta=(dble(loop)/dble(lcut)*2d0-1d0)*pi*0.5 ![-pi/2,pi/2]
!         theta=dble(loop)/dble(lcut)*pi ![0,pi]
! !        c     theta=dble(loop)/dble(lcut)*pi*0.5-pi*0.5d0
!         do i=1,llen(n)
!            vec(jlist(i,1),jlist(i,2))=theta
! !           c               write(*,"(f6.2)",advance="no") vec(jlist(i,1),jlist(i,2))
!         end do
!         do i=1,neh2
!            write(*,"(f7.2)",advance="no") dsin(vec(i,1))
!         enddo
!         write(*,*)
!         do i=1,neh2
!            write(*,"(f7.2)",advance="no") dsin(vec(i,2))
!         enddo
!         write(*,*)
!         etr=0d0
!         dth=0d0
!         do rq=1,nqs         !begin rq loop
!            LD=llist(rq)
!            jsp=jspl(rq)            
!            call mkdiff(jp,-ilo,neh2,nlo,ir,vec,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
!            etr=etr+et/dble(nqs)               
!            dth=dth+dt(1,1)/dble(nqs)
!         end do !end rq loop
!  !       c     write(*,"(2f6.2)") vec(jorb,1),vec(jorb,2)
!         if(etr.gt.emax)then
!            emax=etr
!            thmax=theta
!         end if
!         if(etr.le.emin) emin=etr
!         elist(loop)=etr
!         dlist(loop)=dth
!         thlist(loop)=theta
!      end do
!      write(ifmax,*) "emax,theta-max,sin(theta-max)",emax,thmax,dsin(thmax)
!      do loop=0,lcut
!         etr=(elist(loop)-emin)/(emax-emin)-5d-1
!         write(8000+n,"(4f16.10)") dsin(thlist(loop)),etr,dlist(loop),thlist(loop) !,n,jlist(n,:)
!         write(7000+n,"(4f16.10)") dcos(thlist(loop)),etr,dlist(loop),thlist(loop) !,n,jlist(n,:)
!         write(9000+n,"(4f16.10)") thlist(loop),etr,dlist(loop),thlist(loop) !,n,jlist(n,:)
!      end do
!   end do
!   stop "checkprogram stoped. For return 0, change drawmap .false."
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !     include "map.f90"
!      !         call mkmap(neh2,nlo,NM,llist,NMTO,nqs,ir,ovlM,damM)
!   end if
! 2626 continue

!   do itr=1,maxitr
!      write(*,"(a,i6,i6)") "itr/max=",itr,maxitr
!      lold=-1
!      ibold=-1
!      jskip=0
!      !     do jorb=1,nlo+neh2
!      do jp=1,neh2
!         !            if(jorb.gt.nlo)then
!         !               jp=jorb-nlo
!         !               ilo=1
!         !            else if(jorb.le.nlo)then
!         !               jp=jorb
!         !               ilo=2
!         !     end if
!         !            j1=ir(jp,1)
!         !     j2=ir(jp,2)
!         ilo=1
!         j1=ir(jp,1)
!         j2=ir(jp,2)
!         j3=ir(jp,3)

!         !            if(j3.ne.0)then
!         !               thmax=-99d0
!         !               thmin=100d0
!         !               include "Steepest.F"
!         !            end if

!         if(l_table(j1)==lold .AND. ib_table(j1)==ibold)then
!            if(nlo==0)goto 2020
!            !               if(nlo.ne.0.and.jorb.ne.nlo+1)goto 2020
!            if(nlo /= 0 .AND. jp /= nlo+1)goto 2020
!         end if

!         do loop=0,rcut
!            if(dbg)etr=0.0d0
!            dth=0d0
!            !               thmin=(dble(loop)/dble(lcut)*2d0-1d0)*pi*0.5+thLB
!            !     !thmin=dble(loop)/dble(lcut)*pi*0.5-pi*0.5d0
!            !     thmin=(thUB(ilo)-thLB(ilo))*dble(loop)/dble(lcut)+thLB(ilo)
!            !               thmax=(thUB(ilo)-thLB(ilo))*dble(rcut-loop)/dble(rcut)+thLB(ilo)
!            thmax=(UB-LB)*dble(rcut-loop)/dble(rcut)+LB
!            !               stop

!            do rq=1,nqs      ! q-vector * spin
!               LD=llist(rq)
!               jsp=jspl(rq)

!               vM(:,:)=totvec(:,:)
!               vM(jp,ilo)=thmax
!               call mkdiff(jp,-ilo,neh2,nlo,ir,vM,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
!               dth=dth+dble(dt(1,1))
!            end do
!            dth=dth/dble(nqs)

!            r=dsin(thmax)
!            if(dbg)write(*,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
!            if(dbg)write(idum,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1

!            !               write(*,*) loop,thmax,sin(thmax),dth
!            !               write(*,*) loop,thmin,dth
!            !     if(dth.gt.0d0) exit
!            if(dbg)write(*,*) "thmax",loop,thmax,dth
!            if(dth <= 0d0) exit
!         end do
!         if(loop==rcut) stop "error! thmax was not determined appropriately! Denser rcut is needed!"
!         if(dbg)write(*,*) "================================"
!         if(dbg)write(idum,*) "================================"

!         !            stop

!         do loop=1,lcut
!            if(dbg)etr=0.0d0
!            dth=0.0d0

!            !             thmax=dble(loop)/dble(rcut)*(pi*0.5d0+thmin)-thmin
!            !     thmax=dble(loop)/dble(rcut)*pi*0.5d0+thmin
!            !               thmax=(thUB(ilo)-thmin)*dble(loop)/dble(rcut)+thmin
!            !     write(*,*) loop,thmax,1,thUB,thmin,rcut
!            !               thmin=thmax-(thmax-thLB(ilo))*dble(loop)/dble(lcut)
!            thmin=thmax-(thmax-LB)*dble(loop)/dble(lcut)
!            do rq=1,nqs      ! q-vector * spin
!               jsp=jspl(rq)

!               vM(:,:)=totvec(:,:)
!               vM(jp,ilo)=thmin
!               call mkdiff(jp,-ilo,neh2,nlo,ir,vM,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
!               !                  if(dbg)etr=etr+et/dble(nqs) !*weight(iq)
!               dth=dth+dble(dt(1,1))
!            end do
!            dth=dth/dble(nqs)

!            r=dsin(thmin)
!            if(dbg)write(*,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
!            if(dbg)write(idum,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1

!            !     if(dth.le.0d0) exit
!            if(dbg)write(*,*) "min",loop,thmin,dth
!            if(dth >= 0d0) exit
!         end do
!         if(loop==lcut) stop "error! thmin was not determined appropriately! Denser lcut is needed!"
!         if(dbg)write(*,*) "================================"
!         if(dbg)write(idum,*) "================================"

!         thR=thmax
!         thL=thmin
!         !            write(*,*) dsin(thL)
!         !            write(*,*) dsin(thR)
!         !            stop

!         do loop=1,nbsec
!            etr=0.0d0
!            etrL=0.0d0
!            etrR=0.0d0

!            dth=0.0d0
!            dthL=0.0d0
!            dthR=0.0d0

!            thM=(thR+thL)/2d0

!            do rq=1,nqs      ! q-vector * spin
!               jsp=jspl(rq)

!               vM(:,:)=totvec(:,:)
!               vM(jp,ilo)=thM
!               call mkdiff(jp,-ilo,neh2,nlo,ir,vM,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
!               dth=dth+dble(dt(1,1))

!               vL(:,:)=totvec(:,:)
!               vL(jp,ilo)=thL
!               call mkdiff(jp,-ilo,neh2,nlo,ir,vL,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
!               dthL=dthL+dble(dt(1,1))

!               vR(:,:)=totvec(:,:)
!               vR(jp,ilo)=thR
!               call mkdiff(jp,-ilo,neh2,nlo,ir,vR,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
!               dthR=dthR+dble(dt(1,1))
!            end do

!            dth=dth/dble(nqs)
!            dthL=dthL/dble(nqs)
!            dthR=dthR/dble(nqs)

!            !               write(*,*) dthL,dth,dthR

!            if(dth*dthL > 0)then
!               thL=thM
!            else if(dth*dthR > 0)then
!               thR=thM
!            else
!               stop "error"
!            end if

!            r=dsin(thM)
!            if(dbg)write(*,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
!            if(dbg)write(idum,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
!            if(thR-thL <= cbsec) exit
!         end do

! 2020    totvec(jp,ilo)=thM
!         !     do jsp=1,nsp
!         !            totvec(jp,ilo)=thM
!         ! nd do
!         !     write(*,*) jorb,lold

!         write(*,"(i6,a,es10.2,a,es10.2,a,es10.2,a,es10.2)") &
!              jp,"orb, [",thmin,",",thmax,"] ",thM," diff=",thR-thL

! 2021    lold=l_table(j1)
!         ibold=ib_table(j1)

!         if(dbg)write(*,*) "================================"
!         if(dbg)write(idum,*) "================================"
!      end do

!      ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!      !         if(flg_osym==1)call orbavg(totvec,jsp)

!      write(*,*) "----------------------------------------------"
!      !         do jsp=1,nsp
!      do j=1,nlo
!         ilo=2
!         r=dsin(totvec(j,ilo))
!         absd=abs(totvec(j,ilo)-oldvec(j,ilo))
!         write(*,"(i4,f12.5,2x,es10.1,3x,a20,i4,a)") j,r,absd,orbsym(j),jsp,"spin"
!      end do

!      do j=1,neh2
!         ilo=1
!         r=dsin(totvec(j,ilo))
!         absd=abs(totvec(j,ilo)-oldvec(j,ilo))
!         write(*,"(i4,f12.5,2x,es10.1,3x,a20)") j,r,absd,orbsym(j)
!      end do
!      !         end do
!      if(maxval(abs(totvec(:,:)-oldvec(:,:))) <= convc)then
!         write(*,*) "converged at ",itr,"th iteration"
!         exit
!      end if
!      oldvec=totvec
!      write(*,*) "----------------------------------------------"
!      write(*,*) "**********************",itr,"th iteration *******************"
!   end do

!   allocate(rmatout(1:LL,1:LL))
!   rmatout=0d0
!   !      do jsp=1,nsp
!   do i=1,LL
!      rmatout(i,i)=1d0
!   end do

!   do i=1+nlo,neh2
!      i1=ir(i,1)
!      i2=ir(i,2)
!      rmatout(i1,i1)= dcos(totvec(i,1))
!      rmatout(i2,i1)= dsin(totvec(i,1))
!      rmatout(i1,i2)=-dsin(totvec(i,1))
!      rmatout(i2,i2)= dcos(totvec(i,1))
!   end do

!   if(nlo /= 0)then
!      do i=1,nlo
!         i1=ir(i,1)
!         i2=ir(i,2)
!         i3=ir(i,3)
!         if(i3 <= 0) stop "error(nlo)"
!         rmatout(i1,i1)=dcos(totvec(i,1))*dcos(totvec(i,2))
!         rmatout(i2,i1)=dsin(totvec(i,1))*dcos(totvec(i,2))
!         rmatout(i3,i1)=                  dsin(totvec(i,2))
!      end do
!   end if
!   !      end do

!   !      do jsp=1,nsp
!   open(newunit=ifout,file="rotmat")
!   !         if(jsp==2)open(newunit=ifout,file="rotmat2")
!   do i=1,NMTO
!      do j=1,NMTO
!         write(ifout,"(2i4,2x,2(f12.5,x))") i,j,rmatout(i,j)
!      end do
!   end do
!   close(ifout)
!   !      end do

!   !      do jsp=1,nsp
!   open(newunit=ifout,file="totvec1")
!   !         if(jsp==2)open(newunit=ifout,file="totvec2")
!   open(newunit=ifo2,file="totvec1-sin")
!   !         if(jsp==2)open(newunit=ifo2,file="totvec2-sin")
!   do i=1,neh2
!      r=dsin(totvec(i,1))
!      r2=dsin(totvec(i,2))
!      write(ifout,*) totvec(i,1),totvec(i,2)
!      write(ifo2,*) r,r2
!   end do
!   close(ifout)
!   close(ifo2)
!   !      end do

!   stop
! end program mkhamloc

