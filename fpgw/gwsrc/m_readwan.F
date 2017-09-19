      module m_readwan

      use m_keyvalue,only: getkeyvalue
      use m_iqindx_wan,only: iqindx2_wan
!      use m_hamindex, only: iqimap !!tmp use; need to modify
      
      implicit none
      logical,private:: init=.true.
      real(8),private:: xmx2(3),rlapt(3,3)
      real(8),private::alat,plat(3,3),ef
      integer,private::nwf,natom,nrws,n1,n2,n3
      
      integer,private,allocatable:: irws(:),ibaswf(:)
      real(8),private,allocatable:: pos(:,:),rws(:,:),drws(:)
      complex(8),private,allocatable:: hrotr(:,:,:), evecc(:,:)

      !!!!! wannier eigenvalue and eigenvector
c$$$      complex(8),allocatable:: evecc_w(:,:,:,:)
c$$$      real(8),allocatable:: eval_w(:,:,:)
      complex(8),private,allocatable:: evecc_w(:,:,:,:)
      real(8),private,allocatable:: eval_w(:,:,:)
      
      integer::iwf
      
      contains

c---------------------------------------------
      subroutine readhrotr(is)
      implicit none
      integer(4):: is, ifh, n1, n2, n3, ifile_handle,nq
            
!!--- Required data set --------------------- (copied from htbplot.F)
!! q(:,1:nq) 
!! ef: fermi energy or VBM
!! alat: unit for primitive vector, atomic positions (in a.u.)
!! plat: primitive vector
!! rcut: real-space cutoff for tb parameters.
!! nwf:  # of Wannier function
!! nrws: # of R-R' pair
!! rws(3,i): R-R' vector, i=1,nrws
!! irws(i) : degeneracy, i=1,nrws
!! drws(i) : distance, i=1,nrws
!! hrotr(nwf,nwf,nrws) = Hmn(R) = <0m|H |Rn> (m=1,nwf; n=1,nwf, R=1,nrws)
!! 
!! natom: number of atoms in the primitive cell.
!! ibaswf(nwf): atomic position for the Wannier orbital.
!! pos(3,natom): atomic poisition in the cell.
      write (6,"('--- Reading HrotRS  :: isp=',I6)") is
      ifh=ifile_handle()
      if(is==1) open(ifh,file='HrotRS.up',form='unformatted')
      if(is==2) open(ifh,file='HrotRS.dn',form='unformatted')
      read(ifh)alat,plat,natom
      if (is==1) allocate(pos(3,natom))
      read(ifh)pos
      read(ifh)ef 
      read(ifh)nwf,nrws,n1,n2,n3

      write (6,"('nwf,nrws,n1,n2,n3=',I3,I6,3I4)") nwf,nrws,n1,n2,n3
      if (is==1) then
         allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),
     &        drws(n1*n2*n3*8),ibaswf(nwf), hrotr(nwf,nwf,nrws))
          !real space Hamiltonian in Wannier funciton basis
      endif
      read(ifh) irws,rws,drws,hrotr,ibaswf
      close(ifh)
      
      end subroutine readhrotr
c--------------------------------
      !!! generate eigenvalue list
      subroutine wan_readeigen(qbz,nqbz,is,nspin)
      implicit none

      integer,intent(in)::nqbz,nspin
      real(8),intent(in)::qbz(:,:)
      complex(8),allocatable:: hrotk(:,:,:),hrotkp(:,:),evec(:,:)
      real(8),allocatable:: eval(:)
      integer::iq,is

      !!! fat band plot test
      integer::iffb,ifile_handle,iband
      real(8):: rydberg
      logical::sw=.true.
      
      !! get HrotRS by file handling
      call readhrotr(is)

      if (init) then
         print *,"=== wan_readeigen::"
         allocate(hrotkp(nwf,nwf),evecc(nwf,nwf),eval(nwf))
         if (is==1) allocate(eval_w(nwf,nqbz,nspin),
     &        evecc_w(nwf,nwf,nqbz,nspin))
         do iq = 1,nqbz
!     write(6,*)' got get_hrotkp_ws iq =',iq
            call get_hrotkp_ws2(qbz(:,iq),hrotkp)
            call diag_hm2(hrotkp,nwf,eval,evecc)
            eval_w(1:nwf,iq,is)=eval
            evecc_w(1:nwf,1:nwf,iq,is)=evecc

ccccccccccc eigenvalue check
c$$$            if(iq==1) open(iwf,file="waneval_check.data")
c$$$!     if (kx==1 .and. iq > 100 ) open(iwf,file="wan_eval_check.data")
c$$$            if (qbz(2,iq)==1.0) then
c$$$               if (qbz(3,iq)==1.0) then
c$$$!     if (iq > 100) then
c$$$
c$$$                  write(iwf,"('q(1) ev_w1(:,kx)',6f9.4)") qbz(1,iq),eval_w(:,iq,1)
c$$$               endif
c$$$            endif
c$$$            if (iq==nqbz) close(iwf)
c$$$            if (iq==nqbz) call rx("check: end readwan...")

cccccccccccc   sumcheck for eigenvector  cccccccccccccccc
c$$$            if (iq==nqbz .and. is==2)  then
c$$$!     if (iq==1 .and. is==1) then
c$$$!     if (iq>50 .and. iq<70) then
c$$$               print *," iq:",iq
c$$$               do iwf=1,nwf 
c$$$                  write (6,"('sum evecc(up) evecc(dn) =',2f9.4)")
c$$$     &                 sum(abs(evecc_w(:,iwf,iq,1)**2)),
c$$$     &                 sum(abs(evecc_w(:,iwf,iq,2)**2))
c$$$               enddo
c$$$            endif
cccccccccccccccc  fat band plot (test: Gamma to H) ccccccccccccccccccc
c$$$            if (qbz(2,iq)/=0.0 .and. qbz(2,iq)/=1.0) cycle
c$$$            if (qbz(3,iq)/=0.0 .and. qbz(3,iq)/=1.0) cycle
c$$$            if (qbz(1,iq)>1.05) cycle
c$$$            write (6,"('SYML G-H  qbz',3f9.4)") qbz(:,iq)
c$$$
c$$$            !!! initial operation
c$$$            if (sw) then
c$$$               sw=.false.
c$$$               iffb=ifile_handle()
c$$$               if(is==1) open(iffb,file="fband.up.tmp")
c$$$               if(is==2) open(iffb,file="fband.dn.tmp")
c$$$               write(iffb,"('# ef nwf',f9.4,i4)") ef,nwf
c$$$               write(iffb,*)
c$$$            endif
c$$$            !!! write eigenvalue and eigenvector for wannier
c$$$            do iband=1,nwf
c$$$               write(iffb,"(i5,f13.5,' ',f13.6,i5,' ')",ADVANCE='NO')
c$$$     &          iq,qbz(1,iq),(eval_w(iband,iq,is)-ef)*rydberg(),iband
c$$$               do iwf=1,nwf
c$$$                  write(iffb,"(f13.6)",ADVANCE='NO')
c$$$     &                 (abs(evecc_w(iwf,iband,iq,is)))**2 
c$$$               enddo
c$$$               write(iffb,*)            
c$$$            enddo
c$$$
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
         enddo
         sw=.true.
         close(iffb)
      endif
      deallocate(hrotkp,evecc,eval)
      if (is==2) then
         init=.false.           !!!eigenvalue 
         print *,"=== wan_readeigen end"
         deallocate(pos,rws,drws,hrotr)
      endif
      if (allocated(hrotkp)) deallocate(hrotkp)
      end subroutine wan_readeigen
c--------------------------------
      subroutine wan_readeval(q,isp, ev_w, evc_w)
!      use m_hamindex, only: iqimap
      implicit none
      integer(4), intent(in) :: isp
      real(8), intent(in) :: q(3)
      real(8), intent(out) :: ev_w(nwf)
      complex(8), intent(out) :: evc_w(nwf,nwf)

      integer(4):: iq,iqindx,i
      real(8):: qu(3)
      
      if(init) call rx('wan_readeval: wan_readeigen should be called')
      call iqindx2_wan(q, iq, qu) !qu is used q. q-qu is a G vector.
c$$$      write (6,"('wan_readeval: iq=',I7,'  ---> iqimap, qu',I7,3f9.4)")
c$$$     &     iq,iqimap(iq),qu

c$$$      write (6,"('eval_w=',5f9.4)") eval_w(1:nwf,1,1)
c$$$      write (6,"('eval_w=',5f9.4)") eval_w(1:nwf,1,2)

      !!! return wannier eigenvalue : ev_w
      ev_w(1:nwf) = eval_w(1:nwf,iq,isp)
      !ev_w(1:nwf) = eval_w(1:nwf,iqimap(iq),isp)
      !!! return wannier eigenvalue : evc_w(orbital,band)
      evc_w(1:nwf,1:nwf) = evecc_w(1:nwf,1:nwf,iq,isp)
      end subroutine wan_readeval

c--------------------------------almpost same get_hrotkp_ws in maxloc3
      !!!! hrotr (q) => hrotkp
      subroutine get_hrotkp_ws2(q,hrotkp)

c see Ref.[2] eq.26
      implicit real*8(a-h,o-z)
      complex(8) :: hrotkp(nwf,nwf),
     &     ci,cikr,ceikr,ctmp
      real(8) :: q(3)
      integer::im,in,ir

      pi = 4d0* atan(1d0)
      ci = (0d0,1d0)

      hrotkp = (0d0,0d0)

      do ir = 1,nrws
         rk = sum(rws(:,ir)*q(:))
         cikr = ci * 2d0 * pi * rk
         ceikr = exp(cikr) / dble(irws(ir))
         do im = 1,nwf
            do in = 1,nwf
               hrotkp(im,in) = hrotkp(im,in) + 
     &              ceikr * hrotr(im,in,ir)
            enddo
         enddo
      enddo
      end subroutine get_hrotkp_ws2
c--------------------------------
      subroutine diag_hm2(zmat,ndim,eval,evecc_)
      implicit real*8(a-h,o-z)
      integer:: ndim,i,nev,nmx
      complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
      complex(8):: zmat(ndim,ndim),evecc_(ndim,ndim)
      real(8):: eval(ndim),wk(ndim,11)
      integer iwk(ndim)

      allocate(zmat2(ndim,ndim),ovlpc(ndim,ndim))
      nev  = ndim
      nmx  = ndim
      zmat2 = zmat
      ovlpc = (0d0,0d0)
      do i=1,ndim
         ovlpc(i,i) = (1d0,0d0)
      enddo
      evecc_ = (0d0,0d0)
      eval = 0d0
!      call diagno(ndim,zmat2,ovlpc,wk,iwk,evecc,eval)
      call diagcv(ovlpc,zmat2, evecc_, ndim, eval, nmx, 1d99, nev)
      deallocate(zmat2,ovlpc)
      end subroutine diag_hm2

c---------------------------------------
      subroutine write_qdata(ginv,nqtt_in,qtt_in)
      implicit none
      integer::ifwqb,ifile_handle
      integer::nqtt_in
      real(8)::qtt_in(3,nqtt_in)
      real(8)::ginv(3,3)
      
      open(ifwqb,file="wanqbz",form='unformatted')
      write(ifwqb) ginv
      write(ifwqb) nqtt_in
      write(ifwqb) qtt_in
      close(ifwqb)
      
      end subroutine write_qdata      
c---------------------------------------
      integer function ifile_handle()!! find open file handle
      implicit none
      integer:: i
      logical:: nexist
      do i=5001,9999
         inquire(unit=i,opened=nexist)
         if(.not.nexist) then
            ifile_handle=i
            return
         endif
      enddo
      stop 'ifile_handle: we did not find open file hundle'
      end
c---------------------------------------            
      real(8) function det33(am)
      implicit none
      real(8),intent(in) :: am(3,3)
      det33= am(1,1)*am(2,2)*am(3,3)
     &       -am(1,1)*am(3,2)*am(2,3)
     &       -am(2,1)*am(1,2)*am(3,3)
     &       +am(2,1)*am(3,2)*am(1,3)
     &       +am(3,1)*am(1,2)*am(2,3)
     &       -am(3,1)*am(2,2)*am(1,3)
      end function

c---------------------------------------      
      end module m_readwan


