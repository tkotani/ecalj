      module m_procar
      real(8),allocatable,protected:: dwgtall(:,:,:,:,:)
      integer::  nchanp=25      !total of s,p,d,f
      
      logical,private:: isp1init=.true.,isp2init=.true.,init=.true.
      integer,private:: iprocar1,iprocar2
      contains
      subroutine m_procar_init(fullmesh,iq,isp,nspx,nspc,ef0,evl,ndimh,jsp,qp,nev,t_zv,ndimhx,nmx)
      use m_lmfinit,only: nlmax,nsp,nbas,stdo,sspec=>v_sspec,ssite=>v_ssite,nlmax
      use m_suham,only: ndhamx=>ham_ndhamx,ndham=>ham_ndham
      use m_igv2x,only: igv2x,napw
      use m_mkpot,only: ppnl_rv
      use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize
      use m_qplist,only: nkp,xdatt
      implicit none
      complex(8):: t_zv(ndimhx,nmx)
      character*1000::ccc
      real(8):: ef0
      complex(8):: au,as,az
      real(8):: s11,s22,s33,s12,s13,s23,dwgt(100),dwgtt(100),xdat,qold(3),qp(3)
      complex(8),allocatable:: auspp(:,:,:,:,:)
      integer:: iq,isp,jspini,jspend,jspp,iprocar,iband,is,ilm,ndimh,nspc,jsp,nspx,ib,
     & iprocar1,iprocar2,nev,i,m,l,ndimhx,nmx
      real(8):: rydberg=13.6058d0,evl(ndhamx,nspx)
!! Decompose to alm components. PROCAR generation.
!     ! ??? no nsites is used? ! Compare following codes with those of mkpdos.F (PDOS).
      logical:: cmdopt0,fullmesh,debug,procaron
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
      call makusq( 1 , ssite , sspec ,   nbas , nbas , !sham ,slat ,
     .     0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .     nspc , jsp , 1 , qp , t_zv , ppnl_rv , auspp )
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
      end module
