      module m_procar
      use m_lmfinit,only: nlmax,nsp,nbas,stdo,sspec=>v_sspec,ssite=>v_ssite,nlmax,nspc
      use m_suham,only: ndhamx=>ham_ndhamx,ndham=>ham_ndham,nspx=>ham_nspx 
      use m_igv2x,only: igv2x,napw
      use m_mkpot,only: ppnl_rv
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize,procid
      use m_qplist,only: nkp,xdatt,qplist

      public m_procar_init, m_procar_closeprocar, m_procar_writepdos, dwgtall,nchanp
      
      private
      real(8),allocatable,protected:: dwgtall(:,:,:,:,:)
      integer::  nchanp=25      !total of s,p,d,f
      
      logical,private:: isp1init=.true.,isp2init=.true.,init=.true.
      integer,private:: iprocar1,iprocar2
      logical,private:: cmdopt0,fullmesh,debug,procaron
      
      contains
      subroutine m_procar_closeprocar()
      logical:: nexist
      close(iprocar1)
      inquire(iprocar2,exist=nexist)
      if(nexist) close(iprocar2)
      end
      
      subroutine m_procar_init(iq,isp,ef0,evl,ndimh,jsp,qp,nev,t_zv,ndimhx,nmx)
      implicit none
      complex(8):: t_zv(ndimhx,nmx)
      character*1000::ccc
      real(8):: ef0
      complex(8):: au,as,az
      real(8):: s11,s22,s33,s12,s13,s23,dwgt(100),dwgtt(100),xdat,qold(3),qp(3)
      complex(8),allocatable:: auspp(:,:,:,:,:)
      integer:: iq,isp,jspini,jspend,jspp,iprocar,iband,is,ilm,ndimh,nspc,jsp,ib,
     & nev,i,m,l,ndimhx,nmx
      real(8):: rydberg=13.6058d0,evl(ndhamx,nspx)
      logical:: cmdopt0
!! Decompose to alm components. PROCAR generation.
!     ! ??? no nsites is used? ! Compare following codes with those of mkpdos.F (PDOS).
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') !fullmesh stop just after do 2010 iq loop.
      debug = cmdopt0('--debugbndfp')
      PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
      
      if(isp1init.and.isp==1) then
         open(newunit=iprocar1,file='PROCAR.UP.'//trim(strprocid))
         isp1init=.false.
      endif   
      if(isp2init.and.isp==2) then
         open(newunit=iprocar2,file='PROCAR.DN.'//trim(strprocid))
         isp2init=.false.
      endif   
      if(procaron.and.fullmesh.and.init) then
         allocate(dwgtall(nchanp,nbas,ndhamx,nsp,nkp))
         dwgtall=0d0
         init=.false.
      endif
      
      allocate( auspp(nlmax,ndhamx,3,nsp,nbas) )
      auspp = 0d0
      call makusq(1, nbas,0, nev,jsp,1,qp,t_zv, auspp )
      if(nspx==1) then !spinoffdiagmode 
         jspini=1
         jspend=2
      else
         jspini=isp
         jspend=isp
      endif
      do jspp= jspini, jspend
         if(jspp==1) iprocar=iprocar1
         if(jspp==2) iprocar=iprocar2
         if(debug) write(stdo,*) 'm_procar',iprocar1,iprocar2,jspp,iprocar,ef0,nlmax,ndham,nspc,nsp,nbas
                  ccc="ion        s       py       pz       px      dxy      dyz      dz2      dxz   dx2-y2"//
     &                 "      f-3      f-2      f-1       f0       f1       f2       f3"//
     &                 "                                                                                 tot"
         write(iprocar,*)
         write(iprocar,*)
         write(iprocar,"('k-point ',i4,' :    ',3f11.8,'     weight = -------  : x =',f15.8)")iq,qp,xdatt(iq)
         write(iprocar,*)
         do iband = 1, nev
            write(iprocar,*)
            write(iprocar,"('band ',i3,' # energy ',f13.8,' # occ. -----' )")iband,(evl(iband,jsp)-ef0)*rydberg
            write(iprocar,*)
            dwgtt=0d0
            do ib = 1, nbas
               is  = ssite(ib)%spec
               ilm = 0
               dwgt=0d0
               do  l = 0, sspec(is)%lmxa
                  do  m = -l, l
                     ilm = ilm+1
                     au = auspp(ilm,iband,1,jspp,ib)
                     as = auspp(ilm,iband,2,jspp,ib)
                     az = auspp(ilm,iband,3,jspp,ib)
                     s11 = dconjg(au)*au*ppnl_rv(2,l+1,jspp,ib)
                     s22 = dconjg(as)*as*ppnl_rv(7,l+1,jspp,ib)
                     s33 = dconjg(az)*az*ppnl_rv(8,l+1,jspp,ib)
                     s12 = 0
                     s13 = 2*dconjg(au)*az*ppnl_rv(9,l+1,jspp,ib)
                     s23 = 2*dconjg(as)*az*ppnl_rv(10,l+1,jspp,ib)
                     dwgt(ilm)= s11+s22+s33 + s12+s13+s23
                  enddo
               enddo
               dwgtt = dwgtt + dwgt(1:ilm)
               if(ib==1)  write(iprocar,"(a)") trim(ccc)
               write(iprocar,"(i3,100(x,f8.5))")ib,(dwgt(i),i=1,nchanp),sum(dwgt)
               if(ib==nbas) write(iprocar,"('tot',100(x,f8.5))")(dwgtt(i),i=1,nchanp),sum(dwgtt)
               if(fullmesh) dwgtall(1:nchanp,ib,iband,jspp,iq) = dwgt(1:nchanp)
            enddo
         enddo
      enddo
c      if(debug) print *,'mmmmmmmm procid sum dwgt check=',procid,iq,jsp,sum(dwgtall),nev
      deallocate( auspp )
      end subroutine
      
!!--------------------------------------------------------
      subroutine m_procar_writepdos(evlall,nev_,ef0,kpproc)
      use m_mkqp,only: nkabc=> bz_nabc
      use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
      use m_ext,only:sname
      real(8) evlall(:,:,:)
      integer:: nev_(:)
c     include "mpif.h"
      logical:: cmdopt0
      integer:: kpproc(*)
      integer,allocatable:: ipqe(:,:,:),idtete(:,:)
      integer:: ik1,ik2,ik3,iq,nkk1,nkk2,nkk3,ntete,ifip,nevmin
      real(8):: qx(3),ef0!,emin,emax
      nkk1=nkabc(1)
      nkk2=nkabc(2)
      nkk3=nkabc(3)
      nevmin=minval(nev_(1:nkp))
      
!! pdos mode (--mkprocar and --fullmesh). ===
      if(debug) print *,'mmmm procid sum dwgt check=',procid,sum(dwgtall)
c      if(PROCARon.and.fullmesh) then
c         call mpi_barrier(MPI_comm_world,ierr)  
         call xmpbnd(kpproc,nbas*nchanp*ndhamx,nkp,nspx,dwgtall) 
         if(master_mpi) then  
c            if(debug) print *,'mmmmm tot procid sum dwgt check=',procid,sum(dwgtall)
!! tetrahedron without symmetry
            allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
            iq=0
            do ik3 = 1, nkk3 
               do ik2 = 1, nkk2
                  do ik1 = 1, nkk1
                     iq = iq+1
                     ipqe(ik1,ik2,ik3)=iq
                     qx = (ik1-1)*qlat(:,1)/nkk1 + (ik2-1)*qlat(:,2)/nkk2+ (ik3-1)*qlat(:,3)/nkk3
                     if(abs(sum(qx-qplist(:,iq)))>1d-6) then !sanity check
                        call rx("bndfp: qx/=qplist something strang")
                     endif  
                  enddo
               enddo
            enddo
            call tetirr(qlat, nkk1,nkk2,nkk3, ipqe, ntete,idtete) !tetrahedron
            write(*,*)" ntete 6*nkk1*nkk2*nkk3 nkp=",ntete,6*nkk1*nkk2*nkk3,nkp
!! doswt(1:nchanp=16,ndhamx,nsp,nkp)
!! dwgtall(1:16,nbas,ndhamx,nsp,nkp)
!!   If we will have a memory problem for dwgtall, we may store dwgtall(i,...) for a file.
!!   and repeat nchanp times.
c            eminp=-25.0/rydberg
c            emaxp=30.0/rydberg
c            ndos=5500
            write(stdo,"(' pdosmode: ndhamx nsp nspx =',4i7)") ndhamx, nsp, nspx,nevmin !nevmx
c            print *,' pdosmode: nchanp nbas ndos emin emax ef0=', nchanp,nbas,ndos, emin,emax,ef0
            write(stdo,"(' pdosmode: nchanp nbas ef0=',2i7,3f12.5)") nchanp,nbas,ef0 !mar2015
            open(newunit=ifip,file='pdosdata.'//trim(sname),form='unformatted')
            write(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp !ndos,nkp mar2015
            write(ifip) idtete,ipqe ! allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
            write(ifip) evlall  !  allocate(evlall(ndhamx,nspx,nkp))
            write(ifip) dwgtall !  allocate(dwgtall(nchanp,nbas,ndhamx,nspx,nkp))
            write(ifip) ef0 !eminp,emaxp,ef0
            close(ifip)
            if(nkp/=nkk1*nkk2*nkk3) call rx('pdosmode but nkp/=nkk1*nkk2*nkk3')
            if( cmdopt0('--tetraw')) then
               open(newunit=ifip,file='tetradata.dat',form='unformatted')
               write(ifip) ndhamx,nkp,ntete !ndos,nkp mar2015
               write(ifip) idtete
               close(ifip)
               open(newunit=ifip,file='qplistf.dat')
               do iq=1,nkp
                  write(ifip,"(3f15.8)") qplist(:,iq)
               enddo   
               close(ifip)
            endif               !! we need evlall(ndhamx,nkp) for qplist
            deallocate(idtete,ipqe) !remove this if you uncomment above lines.
         endif  
c      endif
      end subroutine
      end module
