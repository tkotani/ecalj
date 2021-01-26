      subroutine mkdiff(jorb,neh2,nlo,ir,vec,M,L,Omat,Hmat,dth)
      implicit none
      intent(in)  :: jorb,neh2,nlo,ir,vec,M,L,Omat,Hmat
      intent(out) :: dth
      logical     :: dbg,debug=.false.
      integer     :: jorb,neh2,nlo,ir(neh2,3),M,L,j1,j2,j3,jj,i,i1,i2,i3
      real(8)     :: vec(neh2,2),dth
      complex(8)  :: Omat(L,L),Hmat(L,L),Ot(M,M),Ht(M,M),Oinv(M,M)
     . ,rmat(L,M),e2mat(L,L),lomat(L,L),dum(L,L),OdL(M,L),OL(M,L),HL(M,L)
      
      rmat=0d0
      do i=1,M
         rmat(i,i)=1d0
      end do
           
      do i=1+nlo,neh2
         i1=ir(i,1)
         i2=ir(i,2)
         rmat(i1,i1)=dcos(vec(i,1))
         rmat(i2,i1)=dsin(vec(i,1))
      end do

      if(nlo.ne.0)then
         do i=1,nlo
            i1=ir(i,1)
            i2=ir(i,2)
            i3=ir(i,3)
            if(i3.le.0) stop "error(nlo)"
            rmat(i1,i1)=dcos(vec(i,1))*dcos(vec(i,2))
            rmat(i2,i1)=dsin(vec(i,1))*dcos(vec(i,2))
            rmat(i3,i1)=               dsin(vec(i,2))
         end do
      end if      
      
      OL=matmul(transpose(rmat),Omat)
      HL=matmul(transpose(rmat),Hmat)                 
      
      Ot=matmul(OL,rmat)
      Ht=matmul(HL,rmat)
      
      call invs(M,Ot,Oinv) !Oinv is matrix inversion of Ot      
                  
      OdL=matmul(Oinv,HL)-matmul(matmul(Oinv,Ht),matmul(Oinv,OL))
      j1=ir(jorb,1)
      j2=ir(jorb,2)
      j3=ir(jorb,3)      
                     
      if(jorb.le.nlo)then
         dth=OdL(j1,j1)*dcos(vec(jorb,1))*(-dsin(vec(jorb,2)))+OdL(j1,j2)*dsin(vec(jorb,1))*(-dsin(vec(jorb,2)))
     .    +OdL(j1,j3)*dcos(vec(jorb,2))
      else
         dth=OdL(j1,j1)*(-dsin(vec(jorb,1))*dcos(vec(jorb,2)))+OdL(j1,j2)*(dcos(vec(jorb,1))*dcos(vec(jorb,2))) 
      end if
      
      return
      end subroutine mkdiff

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
     .     ibzweight,nqs,qprs,weight,itor,niqs,eferm,orbavg,startprint
      implicit none
      logical:: lprint=.true.,savez=.false.,getz=.false.,dbg=.false.,LLmode=.false.
      logical,allocatable::eflg(:),loflg(:)

      integer(4):: ifmto,ifqpl,ifinput,ifout,i,j,jj,n,m,t,l,k,access,loop,nq,Lmax,i1,i2,NF,LD,imax,itr
     &     ,ikpd,ikp,ifih,it,iq,rq,lold,kold,ibold,ifig=-999,ii,ifile_handle,id,jd,fflg,LLdum,idum
     &     ,jsp,nn,LL,LE,lso1,flg_osym,flg_ssym,lenMTO,nmax,nmin,nevd,j1,j2,ioncut,jorb,lcut,rcut,nbsec,maxitr
     &     ,ifir,NR,NS,NM,nlo,neh2,NL,jp,ilo
      integer(4),allocatable:: nev(:),jspl(:),isptab(:),llist(:),ir(:,:),ilorb(:)

      real(8)::qp(3),convc,h,ral,rimg,ryd,emax,emin,pi,epsdum,iden,fac,prdN,prdF,prdM,norm1,norm2
     &  ,dum,r,ecut,etr,etr0,et,dtmp,tmp,adum,fD,dth,dthL,dthR,thmin
     &  , etrL,etrR,thL,thR,thM,thmax,cbsec,absd
      real(8),allocatable,dimension(:)    :: epsovl,convi,evl_D
      real(8),allocatable,dimension(:,:)  :: vM,vL,vR
      real(8),allocatable,dimension(:,:,:):: totvec,oldvec
      
      complex(8),allocatable,dimension(:,:)   :: hamF,ovlF,hamD,ovlD,c_D,lmat1,lmat2,dmat,bdb
      complex(8),allocatable,dimension(:,:,:) :: ovlM,hamM,hamE,ovlE,damM,damE,rmatout
      
      character(40),allocatable:: orbsym(:)
      character*200 iname,indexname,ifc

      pi=4d0*atan(1d0)
      ryd=13.605693d0
      
      open(newunit=ifinput,file="input_mkhamloc")
      read(ifinput,*) iname 
      close(ifinput)            

      open(newunit=ifmto,file="mtoprm")      
      read(ifmto,*) convc,emin,emax,flg_osym,flg_ssym,maxitr,nbsec,cbsec,lcut,rcut
      close(ifmto)
      emin=emin/ryd+eferm
      emax=emax/ryd+eferm
      
      call startprint(1)
      call ReadHamMTOInfo(iname)      
      allocate(convi(nsp),lmat1(NMTO,NMTO),lmat2(NMTO,NMTO))
      LL=NMTO
      
      open(newunit=ifir,file="ir_list")
      read(ifir,*) NM,neh2,nlo
      NR=neh2+nlo
      NL=NR+NM
      write(*,*) "N(model),NL(model+EH2+lo)=",NM,NL
      write(*,*) "NR(EH2+lo),N(EH2),N(lo)=",NR,neh2,nlo
      
      if(nr.le.0) stop "error! wrong size of nr in ir_list!"
      allocate(ir(neh2,3),orbsym(nr))
      do i=1,neh2
         read(ifir,*) ir(i,1),ir(i,2),ir(i,3),orbsym(i)
      end do
      close(ifir)

c      do i=1,neh2
c         write(*,*) i,ir(i,1),ir(i,2),ir(i,3)
c      end do      
      
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

         nmax=1
         do i=1,nev(rq)
            if(evl_D(i).gt.emax) exit
            nmax=nmax+1
         end do
         nmin=1
         do i=1,nev(rq)
            if(evl_D(i).gt.emin) exit
            nmin=nmin+1
         end do
         
         write(ioncut,"(a,2i5,3f15.5)") "rq,nmin,emin,maxeigen=",rq,nmin,emin*ryd,(evl_D(nmin)-eferm)*ryd
         write(ioncut,"(a,2i5,3f15.5)") "rq,nmax,emax,maxeigen=",rq,nmax,emax*ryd,(evl_D(nmax)-eferm)*ryd

         allocate(bdb(1:LD,1:LD))         
c         bdb=matmul(c_D(1:LD,1:ncut),transpose(dconjg(c_D(1:LD,1:ncut))))
         bdb=matmul(c_D(1:LD,nmin:nmax),transpose(dconjg(c_D(1:LD,nmin:nmax))))
         damM(1:LD,1:LD,rq)=matmul(ovlD,matmul(bdb,ovlD))
c         damM(1:LD,1:LD,rq)=-hamM(1:LD,1:LD,rq)
         deallocate(ovlF,hamF,ovlD,hamD,c_D,evl_D,bdb)
      end do
      close(ioncut)
      write(6,*) "====================end diagonalization!================"
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      allocate(vM(neh2,2),vR(neh2,2),vL(neh2,2),totvec(neh2,2,nsp),oldvec(neh2,2,nsp))
      write(6,*) "====================start main loop of mkhamloc!================"
      open(newunit=idum,file="mkhamloc.log")
c      totvec=-asin(1.0d0)
c      totvec=-0.8d0
      totvec=0d0
      
      totvec(nlo+1:neh2,2,:)=0d0  !reset for local orbital
      oldvec=totvec
      
      do itr=1,maxitr
         write(*,"(a,i6,i6)") "itr/max=",itr,maxitr
         lold=-1
         ibold=-1
         do jorb=1,nlo+neh2
            if(jorb.gt.nlo) jp=jorb-nlo
            if(jorb.le.nlo) jp=jorb
            
            j1=ir(jp,1)
            j2=ir(jp,2)            
            if(l_table(j1)==lold.and.ib_table(j1)==ibold)then
               if(nlo==0)goto 2020
               if(nlo.ne.0.and.jorb.ne.nlo+1)goto 2020
            end if

            do loop=0,lcut
               if(dbg)etr=0.0d0
               dth=0.0d0
            
                thmin=(dble(loop)/dble(lcut)*2d0-1d0)*pi*0.5
c               thmin=dble(loop)/dble(lcut)*pi*0.5-pi*0.5d0
               do rq=1,nqs      ! q-vector
                  LD=llist(rq)
                  jsp=jspl(rq)                  
                  
                  vM(:,:)=totvec(:,:,jsp)
                  if(jorb.le.nlo) vM(jorb,2)=thmin
                  if(jorb.gt.nlo) vM(jorb-nlo,1)=thmin                 
                  call mkdiff(jorb,neh2,nlo,ir,vM,NM,NL,ovlM(1:NL,1:NL,rq),damM(1:NL,1:NL,rq),dtmp)
                  dth=dth+dtmp
               end do
               dth=dth/dble(nqs)               

               r=dsin(thmin)
               if(dbg)write(*,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(dbg)write(idum,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1

               if(dth.gt.0d0) exit
            end do
            if(dbg)write(*,*) "================================"
            if(dbg)write(idum,*) "================================"

            do loop=1,rcut
               if(dbg)etr=0.0d0
               dth=0.0d0
               
             thmax=dble(loop)/dble(rcut)*(pi*0.5d0+thmin)-thmin
c               thmax=dble(loop)/dble(rcut)*pi*0.5d0+thmin
               do rq=1,nqs      ! q-vector
                  jsp=jspl(rq)                  
                  
                  vM(:,:)=totvec(:,:,jsp)
                  if(jorb.le.nlo) vM(jorb,2)=thmax
                  if(jorb.gt.nlo) vM(jorb-nlo,1)=thmax
                  call mkdiff(jorb,neh2,nlo,ir,vM,NM,NL,ovlM(1:NL,1:NL,rq),damM(1:NL,1:NL,rq),dtmp)
                                    
c                  if(dbg)etr=etr+et/dble(nqs) !*weight(iq)
                  dth=dth+dtmp
               end do
               dth=dth/dble(nqs)
               
               r=dsin(thmax)
               if(dbg)write(*,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(dbg)write(idum,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1

               if(dth.le.0d0) exit
            end do
            if(dbg)write(*,*) "================================"
            if(dbg)write(idum,*) "================================"
                        
            thR=thmax
            thL=thmin
c            write(*,*) dsin(thL)
c            write(*,*) dsin(thR)
c            stop

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
                  
                  vM(:,:)=totvec(:,:,jsp)
                  if(jorb.le.nlo) vM(jorb,2)=thM
                  if(jorb.gt.nlo) vM(jorb-nlo,1)=thM
                  call mkdiff(jorb,neh2,nlo,ir,vM,NM,NL,ovlM(1:NL,1:NL,rq),damM(1:NL,1:NL,rq),dtmp)
                  dth=dth+dtmp

                  vL(:,:)=totvec(:,:,jsp)
                  if(jorb.le.nlo) vL(jorb,2)=thL
                  if(jorb.gt.nlo) vL(jorb-nlo,1)=thL
                  call mkdiff(jorb,neh2,nlo,ir,vL,NM,NL,ovlM(1:NL,1:NL,rq),damM(1:NL,1:NL,rq),dtmp)
                  dthL=dthL+dtmp

                  vR(:,:)=totvec(:,:,jsp)
                  if(jorb.le.nlo) vR(jorb,2)=thR
                  if(jorb.gt.nlo) vR(jorb-nlo,1)=thR
                  call mkdiff(jorb,neh2,nlo,ir,vR,NM,NL,ovlM(1:NL,1:NL,rq),damM(1:NL,1:NL,rq),dtmp)
                  dthR=dthR+dtmp
               end do

               dth=dth/dble(nqs)
               dthL=dthL/dble(nqs)
               dthR=dthR/dble(nqs)

c               write(*,*) dthL,dth,dthR

               if(dth*dthL.gt.0)then
                  thL=thM
               else if(dth*dthR.gt.0)then
                  thR=thM            
               else
                  stop "error"
               end if               
               
               r=dsin(thM)
               if(dbg)write(*,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(dbg)write(idum,"(i4,3(f11.6,2x),a20,i4,a,i4)") loop,etr,dth,r,orbsym(j1),jsp,"(spin)",j1
               if(thR-thL.le.cbsec) exit               
            end do
            
 2020       if(jorb.le.nlo) totvec(jorb,2,jsp)=thM
            if(jorb.ge.nlo) totvec(jorb-nlo,1,jsp)=thM
            lold=l_table(j1)
            ibold=ib_table(j1)
c     write(*,*) jorb,lold
            
c            do i=nlo+1,neh2
c               write(*,*) i,totvec(i,2,1),totvec(i,2,2)
c               if(totvec(i,2,1).ne.0d0) stop "error! LO constrain was brocken!(jsp==1)"
c               if(totvec(i,2,2).ne.0d0) stop "error! LO constrain was brocken!(jsp==2)"
c            end do
c            stop
            
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
               if(i.le.nlo)then
                  ilo=2
                  j=i
               else
                  ilo=1
                  j=i-nlo
               end if
               r=dsin(totvec(j,ilo,jsp))
               absd=abs(totvec(j,ilo,jsp)-oldvec(j,ilo,jsp))
         write(*,"(2i4,f12.5,2x,es10.1,3x,a20,i4,a)") i,j,r,absd,orbsym(j),jsp,"spin"
            end do
         end do
         if(maxval(abs(totvec(:,:,:)-oldvec(:,:,:))).le.convc)then
            write(*,*) "converged at ",itr,"th iteration"
            exit
         end if
         oldvec=totvec
         write(*,*) "----------------------------------------------"
         write(*,*) "**********************",itr,"th iteration *******************"         
      end do

      allocate(rmatout(1:LL,1:LL,1:nsp))
      
      do jsp=1,nsp
         if(jsp==1)open(newunit=ifout,file="rotmat1")
         if(jsp==2)open(newunit=ifout,file="rotmat2")     

         do i=1,NF
            do j=1,NF            
               write(ifout,"(2i4,2x,2(f12.5,x))") i,j,rmatout(i,j,jsp)
            end do
         end do
         
         close(ifout)
      end do
      stop
      end program mkhamloc
