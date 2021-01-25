      subroutine mkdiff(j,S,M,L,vec,ir,Omat,Hmat,dth)
      implicit none
      intent(in)  :: j,S,M,L,vec,ir,Omat,Hmat
      intent(out) :: dth
      logical     :: dbg,debug=.false.
      integer     :: j,j1,j2,M,S,nr,L,i,i1,i2,ir(M-S,2),jj
      real(8)     :: diff,vec(M-S),dth
      complex(8)  :: Omat(L,L),Hmat(L,L),Omats(M,M),Hmats(M,M)
     .     ,Ot(S,S),Ht(S,S),Oinv(S,S),rmat(M,S),OdH(S,M),OHOd(S,M),tmat(S,M),OL(S,M),HL(S,M)

      NR=M-S      
      Omats(1:M,1:M)=Omat(1:M,1:M)
      Hmats(1:M,1:M)=Hmat(1:M,1:M)
      
      rmat=0d0
      do i=1,S
         rmat(i,i)=1d0
      end do
      
      do i=1,NR
         i1=ir(i,1)
         i2=ir(i,2)
         rmat(i1,i1) =  dcos(vec(i))
         rmat(i2,i1) =  dsin(vec(i))         
      end do
      
      OL=matmul(transpose(rmat),Omats)
      HL=matmul(transpose(rmat),Hmats)                 
      
      Ot=matmul(OL,rmat)
      Ht=matmul(HL,rmat)
      
      call invs(S,Ot,Oinv) !Oinv is matrix inversion of Ot      
                  
      OdH=matmul(Oinv,HL)               
      OHOd=matmul(matmul(Oinv,Ht),matmul(Oinv,OL))
      
      tmat(:,:)=OdH(:,:)-OHOd(:,:)

      j1=ir(j,1)
      j2=ir(j,2)                  
      dth=tmat(j1,j1)*(-dsin(vec(j)))+tmat(j1,j2)*(dcos(vec(j)))
      
      return
      end subroutine mkdiff

      subroutine mktr(j,S,L,vec,ir,Omat,Hmat,etr)
      implicit none
      intent(in)  :: j,S,L,vec,ir,Omat,Hmat
      intent(out) :: etr
      logical     :: dbg,debug=.false.
      integer     :: j,j1,j2,S,L,i,i1,i2,ir(S,2),jj
      real(8)     :: etr,diff,vec(S)
      complex(8)  :: Omat(L,L),Hmat(L,L)
     .     ,Ot(S,S),Ht(S,S),Oinv(S,S),emat(S,S),Omats(S*2,S*2),Hmats(S*2,S*2)
     .     ,rmat(S*2,S),OdH(S,S*2),OHOd(S,S*2),tmat(S,S*2),OL(S,S*2),HL(S,S*2)

      Omats(1:S*2,1:S*2)=Omat(1:S*2,1:S*2)
      Hmats(1:S*2,1:S*2)=Hmat(1:S*2,1:S*2)
      
      rmat=0d0
      
      do i=1,S
         i1=ir(i,1)
         i2=ir(i,2)
         rmat(i1,i1) =  dcos(vec(i))
         rmat(i2,i1) =  dsin(vec(i))         
      end do
      
      do i=1,L-S
         rmat(L+1-i,L-S+1-i)=1.0d0
      end do
         
      do jj=1,S
         do i=1,S*2
            write(*,*) i,jj,rmat(i,jj)
         end do
      end do

      emat=matmul(Oinv,Ht)
      etr=0d0
      do i=1,S
         etr=etr+emat(i,i)
      end do

      return
      end subroutine mktr

      
      subroutine updatermat(totvec,LL,nr,ir,nsp,rmat)
      implicit none
      intent(in)  :: totvec,LL,nr,ir,nsp
      intent(out) :: rmat
      integer     :: i,LL,nr,nsp,ir(nr,2),jsp
      real(8)     :: totvec(1:nr,nsp)
      complex(8)  :: rmat(LL,LL,nsp)
      rmat=0d0

      do i=1,LL
         rmat(i,i,:)=1.0d0
      end do

      do i=1,nr
         do jsp=1,nsp
            rmat(ir(i,1),ir(i,1),jsp) =  dcos(totvec(i,jsp))
            rmat(ir(i,1),ir(i,2),jsp) = -dsin(totvec(i,jsp))
            rmat(ir(i,2),ir(i,1),jsp) =  dsin(totvec(i,jsp))
            rmat(ir(i,2),ir(i,2),jsp) =  dcos(totvec(i,jsp))
         end do
      end do              
      return
      end subroutine updatermat
      
      subroutine output(LL,NF,nsp,rmat)
      implicit none
      intent(in) :: LL,NF,nsp,rmat
      integer :: LL,NF,nsp,jsp,i,j,ifout,ifout2,ifoutb,ifout3
      complex(8) :: rmat(LL,LL,nsp)      
      
      do jsp=1,nsp
         if(jsp==1)open(newunit=ifout,file="rotmat1")
         if(jsp==2)open(newunit=ifout,file="rotmat2")     

         do i=1,NF
            do j=1,NF            
               write(ifout,"(2i4,2x,2(f12.5,x))") i,j,rmat(i,j,jsp)
            end do
         end do
         
         close(ifout)
      end do
      return
      end subroutine output

      subroutine invs(len,mat,inv)
      implicit none
      intent(in)  :: len,mat
      intent(out) :: inv
      integer     :: len,ipv(len),info
      complex(8)  :: mat(len,len),inv(len,len),work(len*3)      
      inv=mat
      call zgetrf(len, len, inv, len, ipv(1:len), info) 
      call zgetri(len, inv, len, ipv(1:len), work, len*3, info)     
      return
      end subroutine invs      
      
      program mkhamloc     
      use m_hamMTO,only:ReadHamMTOInfo,qplist,npair,ib_table,l_table,k_table,m_table,nsp,nkp,NMTO,lso,orbch,
     .     ibzweight,nqs,qprs,weight,itor,niqs,orbsymlo,eferm,orbavg,startprint
      implicit none
      logical:: lprint=.true.,savez=.false.,getz=.false.,dbg=.false.,LLmode=.false.
      logical,allocatable::eflg(:),loflg(:)

      integer(4):: ifmto,ifqpl,ifinput,ifout,i,j,jj,n,m,t,l,k,access,loop,nq,Lmax,i1,i2,NF,LD,imax,itr
     &     ,ikpd,ikp,ifih,it,iq,rq,lold,kold,ibold,ifig=-999,ii,ifile_handle,id,jd,fflg,LLdum,idum
     &     ,jsp,nn,LL,LE,lso1,flg_osym,flg_ssym,lenMTO,ncut,nevd,j1,j2,ioncut,jorb,lcut,rcut,nbsec,maxitr
     &     ,ifir,NR,NS,NM
      integer(4),allocatable:: nev(:),jspl(:),isptab(:),llist(:),ir(:,:)

      real(8)::qp(3),convc,h,ral,rimg,ryd,emax,emin,pi,epsdum,iden,fac,prdN,prdF,prdM,norm1,norm2
     &  ,dum,r,ecut,etr,etr0,et,dtmp,tmp,adum,fD,dth,dthL,dthR,thmin
     &  , etrL,etrR,thL,thR,thM,thmax,cbsec
      real(8),allocatable,dimension(:)  :: epsovl,convi,evl_D,vM,vL,vR
      real(8),allocatable,dimension(:,:):: totvec,oldvec
      
      complex(8),allocatable,dimension(:,:)   :: hamF,ovlF,hamD,ovlD,c_D,lmat1,lmat2,dmat,bdb
      complex(8),allocatable,dimension(:,:,:) :: ovlM,hamM,hamE,ovlE,damM,damE,rmat1

      character(20),allocatable:: orbsym(:)
      character*200 iname,indexname,ifc

      pi=4d0*atan(1d0)
      ryd=13.605693d0
      
      open(newunit=ifinput,file="input_mkhamloc")
      read(ifinput,*) iname 
      close(ifinput)            

      open(newunit=ifmto,file="mtoprm")      
      read(ifmto,*) convc,ecut,flg_osym,flg_ssym,maxitr,nbsec,cbsec,lcut,rcut
      close(ifmto)
      ecut=ecut/ryd+eferm
      
      call startprint(1)
      call ReadHamMTOInfo(iname)      
      allocate(convi(nsp),lmat1(NMTO,NMTO),lmat2(NMTO,NMTO))
      LL=NMTO
      
      open(newunit=ifir,file="ir_list")
      read(ifir,*) NM,NS
      nr=NM-NS
      write(*,*) "NR,NM,NS=",nr,nm,ns
      if(nr.le.0) stop "error! wrong size of nr in ir_list!"
      allocate(ir(nr,2),orbsym(nr))
      do i=1,nr
         read(ifir,*) ir(i,1),ir(i,2),orbsym(i)
      end do
      close(ifir)
      
c      call readlmat(NMTO,nsp,lmat1,lmat2)
      
      call ibzweight()
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(6,*) "====================open .K file! ===================="
      allocate(jspl(nqs),llist(nqs))
      open(newunit=ifih,file=trim(adjustl(iname))//".K",form='unformatted')      
      i=0 ; j=0; rq=0           !sanity check
      do          
         read(ifih,end=2019) qp,LD,lso1,epsdum,jsp         
         if(lso.ne.lso1) stop "error! inconsistent lso!"         
         if(jsp==1) i=i+1       !for check
         rq=rq+1
         llist(rq)=LD
         jspl(rq)=jsp
         if((qp(1)-qprs(1,rq))**2d0+(qp(2)-qprs(2,rq))**2d0+(qp(3)-qprs(3,rq))**2d0.gt.1d-10) stop "error! inconsistent qp!"         
         allocate(ovlF(1:LD,1:LD),hamF(1:LD,1:LD))         
         read(ifih) ovlF(1:LD,1:LD) 
         read(ifih) hamF(1:LD,1:LD) 
         deallocate(ovlF,hamF)
      end do
 2019 close(ifih)
      write(*,*) "**************sanity check is ok************************"
      if(i.ne.nkp.or.rq.ne.nqs)then
         write(*,*) i,nq,rq,nqs
         stop "inconsistent nq between *.Info and *.K"
      end if
      
      write(6,*)'Read: total # of q for Ham=',nkp,i,"(from .Info and .K)"
      write(6,*) "====================close .K file! ===================="      
      open(newunit=ifih,file=trim(adjustl(iname))//".K",form='unformatted')

      Lmax=maxval(llist)
      if(LLmode) llist(:)=LL
      allocate(ovlM(1:Lmax,1:Lmax,1:nqs),hamM(1:Lmax,1:Lmax,1:nqs),damM(1:Lmax,1:Lmax,1:nqs))
c      allocate(eflg(nqs)),epsovl(nqs))    
cccccccccccccccccccccc diagonalization process! ccccccccccccccccccccccccccc      
      allocate(nev(nqs),isptab(niqs)) 

      damM=0d0
      write(6,*) "====================diagonalization!===================="
      open(newunit=ioncut,file="ncut.kpoint")
      do rq=1,nqs
         read(ifih) qp,LD,lso1,epsdum,jsp
         
         allocate(ovlF(1:LD,1:LD),hamF(1:LD,1:LD))         
         read(ifih) ovlF(1:LD,1:LD) 
         read(ifih) hamF(1:LD,1:LD)
         if(LLmode) LD=LL
         
         ovlM(1:LD,1:LD,rq)=ovlF(1:LD,1:LD)
         hamM(1:LD,1:LD,rq)=hamF(1:LD,1:LD)
         
         allocate(hamD(LD,LD),ovlD(LD,LD),c_D(LD,LD),evl_D(LD))
         ovlD(1:LD,1:LD)=ovlF(1:LD,1:LD)
         hamD(1:LD,1:LD)=hamF(1:LD,1:LD)
         if(epsdum  < 1.000001d-14) call zhev_tk2(LD,hamD,ovlD,LD,nev(rq),evl_D,c_D,lprint,savez,getz,ifig)
         if(epsdum >= 1.000001d-14)
     .   call zhev_tk3(LD,hamD,ovlD,LD,nev(rq),evl_D,c_D,lprint,savez,getz,ifig,epsdum)

         ncut=0
         do i=1,nev(rq)
            if(evl_D(i).gt.ecut) exit
            ncut=ncut+1
         end do
         write(ioncut,"(a,2i5,2f15.5)") "iq,ncut,ecut,maxeigenval=",iq,ncut,ecut,evl_D(i)

         allocate(bdb(1:LD,1:LD))
         bdb=matmul(c_D(1:LD,1:ncut),transpose(dconjg(c_D(1:LD,1:ncut))))         
c         damM(1:LD,1:LD,rq)=matmul(ovlD,matmul(bdb,ovlD))
         damM(1:LD,1:LD,rq)=-hamM(1:LD,1:LD,rq)
         deallocate(ovlF,hamF,ovlD,hamD,c_D,evl_D,bdb)
      end do
      close(ioncut)
      write(6,*) "====================end diagonalization!================"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      allocate(vM(nr),vR(nr),vL(nr),totvec(1:nr,nsp),oldvec(1:nr,nsp))      
      write(6,*) "====================start main loop of mkhamloc!================"
      open(newunit=idum,file="mkhamloc.log")
c      totvec=-asin(1.0d0)
c      totvec=-0.8d0
      totvec=0d0
      oldvec=totvec

      do itr=1,maxitr
         write(*,"(a,i6,i6)") "itr/max=",itr,maxitr
         lold=-1
         ibold=-1
         do jorb=1,nr
            j1=ir(jorb,1)
            j2=ir(jorb,2)            
            if(l_table(j1)==lold.and.ib_table(j1)==ibold)then
               goto 2020
            end if

            write(*,*) "lcut"
            do loop=0,lcut
               if(dbg)etr=0.0d0
               dth=0.0d0
            
               thmin=dble(loop)/dble(lcut)*pi*0.5-pi*0.5d0            
               do rq=1,nqs      ! q-vector
                  LD=llist(rq)
                  jsp=jspl(rq)                  
                  
                  vM(:)=totvec(:,jsp)
                  vM(jorb)=thmin
                  call mkdiff(jorb,NS,NM,LD,vM,ir,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),dtmp)
c                  if(dbg)etr=etr+et/dble(nqs) 
                  dth=dth+dtmp
               end do
               dth=dth/dble(nqs)               

               r=dsin(thmin)
               if(dbg)write(*,"(i4,3(f11.6,2x),a12,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(dbg)write(idum,"(i4,3(f11.6,2x),a12,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1

               if(dth.gt.0d0) exit
            end do
            if(dbg)write(*,*) "================================"
            if(dbg)write(idum,*) "================================"

            write(*,*) "rcut"
            do loop=1,rcut
               if(dbg)etr=0.0d0
               dth=0.0d0
               
               thmax=dble(loop)/dble(rcut)*pi*0.5d0+thmin            
               do rq=1,nqs      ! q-vector
                  jsp=jspl(rq)                  

                  vM(:)=totvec(:,jsp)
                  vM(jorb)=thmax               
                  call mkdiff(jorb,NS,NM,LD,vM,ir,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),dtmp)
c                  if(dbg)etr=etr+et/dble(nqs) !*weight(iq)
                  dth=dth+dtmp
               end do
               dth=dth/dble(nqs)
               
               r=dsin(thmax)
               if(dbg)write(*,"(i4,3(f11.6,2x),a12,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(dbg)write(idum,"(i4,3(f11.6,2x),a12,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1

               if(dth.le.0d0) exit
            end do
            if(dbg)write(*,*) "================================"
            if(dbg)write(idum,*) "================================"
                        
            thR=thmax
            thL=thmin

            write(*,*) "bisec"
            do loop=1,nbsec
               etr=0.0d0
               etrL=0.0d0
               etrR=0.0d0
               
               dth=0.0d0
               dthL=0.0d0
               dthR=0.0d0
               
               thM=(thR+thL)/2d0
               
               do rq=1,nqs      ! q-vector
                  jsp=jspl(rq)

                  vM(:)=totvec(:,jsp)
                  vM(jorb)=thM               
                  call mkdiff(jorb,NS,NM,LD,vM,ir,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),dtmp)
c                  if(dbg)etr=etr+et/dble(nqs) 
                  dth=dth+dtmp

                  vL(:)=totvec(:,jsp)
                  vL(jorb)=thL
                  call mkdiff(jorb,NS,NM,LD,vL,ir,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),dtmp)
c                  if(dbg)etrL=etrL+et/dble(nqs) 
                  dthL=dthL+dtmp
                  
                  vR(:)=totvec(:,jsp)
                  vR(jorb)=thR
                  call mkdiff(jorb,NS,NM,LD,vR,ir,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),dtmp)
c                  if(dbg)etrR=etrR+et/dble(nqs) 
                  dthR=dthR+dtmp
               end do

               dth=dth/dble(nqs)
               dthL=dthL/dble(nqs)
               dthR=dthR/dble(nqs)

               if(dth*dthL.gt.0)then
                  thL=thM
               else if(dth*dthR.gt.0)then
                  thR=thM            
               else
                  stop "error"
               end if               
               
               r=dsin(thM)
               if(dbg)write(*,"(i4,3(f11.6,2x),a12,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(dbg)write(idum,"(i4,3(f11.6,2x),a12,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(thR-thL.le.cbsec) exit               
            end do
            
 2020       totvec(jorb,jsp)=thM
            lold=l_table(j1)
            ibold=ib_table(j1)
c            write(*,*) jorb,lold
            write(*,"(i6,a,es10.2,a,es10.2,a,es10.2,a,es10.2)")
     .           jorb,"orb, [",dsin(thmin),",",dsin(thmax),"] ",dsin(thM)," diff=",thR-thL
            
            if(dbg)write(*,*) "================================"
            if(dbg)write(idum,*) "================================"         
         end do
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         

         if(flg_osym==1)call orbavg(totvec,jsp)

         write(*,*) "----------------------------------------------"
         do jsp=1,nsp            
            do i=1,nr               
               r=dsin(totvec(i,jsp))
         write(*,"(i4,f12.5,2x,es10.1,3x,a12,i4,a)") i,r,abs(totvec(i,jsp)-oldvec(i,jsp)),orbsym(i),jsp,"spin"
            end do
         end do
         if(maxval(abs(totvec(:,:)-oldvec(:,:))).le.convc)then
            write(*,*) "converged at ",itr,"th iteration"
            exit
         end if
         oldvec=totvec
         write(*,*) "----------------------------------------------"
         write(*,*) "**********************",itr,"th iteration *******************"         
      end do

      allocate(rmat1(1:LL,1:LL,1:nsp))
      call updatermat(totvec,LL,nr,ir,nsp,rmat1)
      call output(LL,LL,nsp,rmat1)
      stop
      end program mkhamloc
