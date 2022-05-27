!> lmf2gw() set variables to module variables by reading files from following input files.
!    gwb.head
!    NLAindx      
!    CLASS
!    lmfgw_kdivider
!    gwa.*
!    QGpsi,QGcou,Q0P,QIBZ
!      
c--   From QGpsi
Cr  QpGcut_psi  : Maxmum of |q+G| in a.u. for eigenfunction
Cr  ngp         : the number of the IPW (it is k-dependent)
Cr  ngvec
c--   From QGcou
Cr  QpGcut_cou  : Maxmum of |q+G| in a.u. for Coulomb matrix
Cr  ngc         : the number of the IPW (it is k-dependent)
c
!cr Radial mesh is specified from, nr, aa, and bb; mesh is r(i) = bb*(exp(aa*(i-1)) -1) ,i=1,nr
!!cr----NLAindx start---------------
!!cr     ldim2
!!cr     n l ibas for each line --- Repeat this line for the number of ldim2.
!!cr                            --- nphimx=maximum value of n.
!!     nsp,      ! =1 or 2, corresponding to para or ferro.
!!     nbas,     ! Number of atom in the primitive cell
!!     nclass,   ! Number of atomic class (or type) for in the primitive cell
!!     nrmx,     ! = maxval(nr(1:nclass))  Maximum number of nr
!!     ncoremx,  ! = maxval(ncore(1:nclass))
!!     lmxamx,   ! = maxval(lmxa(1:nclass))
!!     ngpmx,    ! Maximum number of G vector.
!!     nband,    ! Number of bands given by GWIN0
!!     nqbze,    ! = nqbz*(1+nq0i). Number of q points given by qg4gw
!!     iclass(nbas),   ! class is starting from 1.
!!     lmxa  (nclass), ! Maximum l number for each atom for augmentation.
!!     nr(nclass),     ! Size of radial mesh.
!!     konf(0:lmxamx,nclass),! Principle quantum numbers of valence electron.
!!c                              ! For example, 4s4p3d means konf(0)=4, konf(1)=4, konf(2)=3.
!!c                              ! Core orbitals are specified by
!!c                              !   1, 2,... konf(0)-1 for s
!!c                              !   2, 3,... konf(1)-1 for p
!!c                              !   3, 4,... konf(2)-1 for d, and so on.
!!c                              !
!!     ncore(nclass)   ! ncore = \sum_l  (konf(l)-1) - l
!!c                        ! Number of different core orbitals for each atom.
!!      integer:: nphi,          ! number of augmentation nRl channels, as distinct from:
!!     ldim2,         ! number of nRLm channels = 2*ldim for phi+phidot case.
!!     nphimx,        ! Maxmum number of phi for all l ---  2 for phi+phidot case.
!!     nindx(ldim2),  ! n    index
!!     lindx(ldim2),  ! l    index
!!     ibasindx(ldim2)! ibas index
!!c
!!      real(8) ::
!!     zz(nclass),     ! Atomic number.
!!     aa(nclass),bb(nclass),! Radial mesh are specified by these parameters with nr.
!!     bas(3,nbas),    ! Atomic posion in the Cartesian coordinate (alat unit),.
!!     alat,           ! Lattice constant in a.u.
!!     plat(3,3),   ! Primitive translation vectors plat(1:3,1) is 1st vector, ans so on.
!!     qbze(3,nqbze),  ! q points given by qg4gw
!!     efermi,         ! Fermi energy. It should be calculated for n1 n2 n3 given by GIWN0.
!!     ec(ncoremx, nclass, nsp),   ! Eigenvalues for core
!!     evl    (nband, nqirr, nsp), ! Eigenvalues
!!     vxclda (nband, nqirr, nsp), ! Lda XC potential <psi|Vxc(n_total)  |psi>
!!     gx (nrmx, 0:lmxamx, nphimx, nclass,nsp), ! Radial function.
!!         gx (nrmx, lmxamx, 1, nclass,nsp) is phi
!!         gx (nrmx, lmxamx, 2, nclass,nsp) is phidot
!!     gcore(nrmx, ncoremx, nclass,nsp)  ! Core radial function.
!!c    These radial functions are only the major part given by the scalar relativistic calculations.
!!c     gx and gcore = r \phi(r) = u(r), where \phi is the major part of the true radial functions.
!!c     gx is normalized as 1 = \int dr gx**2
!!c     gcore is the major part of the true radial functions.
!!      cphi(ldim2,  nband, nqirr,nsp), ! Coefficients of eigenfunction
!!      geig(ngpmx,  nband, nqirr,nsp) ! Coefficients of eigenfunction for IPW.
      
      module m_lmf2gw
      integer,allocatable,protected:: nindx(:),lindx(:),ibasindx(:),iantiferro(:),mnla(:,:)
      integer,protected :: nbandmx,nphimx,!n1q,n2q,n3q,
     &nsp,      ! =1 or 2, corresponding to para or ferro.
     &nbas,     ! Number of atom in the primitive cell
     &nclass,   ! Number of atomic class (or type) for in the primitive cell
     &nrmx,     ! = maxval(nr(1:nclass))  Maximum number of nr
     &ncoremx,  ! = maxval(ncore(1:nclass))
     &lmxamx,   ! = maxval(lmxa(1:nclass))
     &ngpmx,    ! Maximum number of G vector.
     &ldim2    ! = total number of augmentation functions nRlm
      character(8),allocatable,protected:: spid(:)
      integer,allocatable::
     &     iclass(:),lmxa_d  (:),nr(:),konf_d(:,:),ncore_d(:),ibasf(:)
      real(8),allocatable ::
     &     zz(:),aa(:),bb(:),ec_d (:,:,:),evl_d(:,:,:),
     &     gx_d(:,:,:,:,:),gcore_d(:,:,:,:), bas(:,:)
      real(8):: plat(3,3), alat, efermi
      complex(8),allocatable:: cphi_d(:,:,:,:)
      complex(8),allocatable:: geig_d(:,:,:,:)
c      real(8),allocatable :: qirr(:,:)
      logical,protected:: laf   !! - laf: antiferro switch
      integer,protected:: ngcmx,nqnum,nqnumc,nqtt,nq0i,iq0pin,nq0iadd,nqbz,nqibz,nqbzx
      real(8):: QpGcut_psi,QpGcut_cou
      real(8),allocatable :: wt(:),q0i(:,:)
      integer,allocatable,target:: ngvecptt(:,:,:),ngvecctt(:,:,:),ngptt(:),ngctt(:)
      real(8),allocatable:: qtt(:,:)

      integer,allocatable:: ngplist(:),ndimhall(:),iqindex(:)
      real(8),allocatable:: qplist(:,:)

c      character*256,allocatable:: extp(:)
c      integer,allocatable:: iprocq(:,:)
      integer:: nqirr     ! = Number of q points for irr=1 (see m_qplist, output of qg4gw).
      real(8),allocatable,protected :: qibz(:,:)

      contains
      subroutine lmf2gw()
!! Files, gwa.* gwb.* gw1.* gw2.* are converted to DATA4GW and CphiGeig.
!! Input files
!!   gwa.* : atomic data
!!   gwb.* : band data
!!   gw1.* : <psi|H|psi>
!!   gw2.* : <psi|H(without Vxc)|psi>
!!   CLASS, lmfgw_kdivider, NLAindx
c------------------------------------------------------------------
      use m_keyvalue,only: getkeyvalue
      use m_hamindex0,only: readhamindex0,nclass_in=>nclass,iclass_in=>iclasst,
     &   nindx_in=>nindx,lindx_in=>lindx,ibasindx_in=>ibasindx,nphimx_in=>nphimx
      implicit none
      integer:: iq0p
      integer::
     &ldim,     ! = sum ( (lmxa(1:nbas)+1)**2 )
     &nband    ! Number of bands given by GWIN0
      integer icor1,icorex,i,i1,i2,ibas,ibasx,ibx,ic,icore,
     .ifichkv,ifigw0,ifigwa,ifigwb,ifigwx1,ifigwx2,ifigwx3,isp,ispx,
     .ispxx,ix,kkk,kkkdummy,l,ldummy,lmxa,lxx,nclassx,m,n,
     .ncore,ndimh,ngp,nnc,nspdummy,IKP,NR_A !takao feb2012 ngc,ngcmx,
      real(8) ovv(20),ef0,z,a_a,b_a,rofi_anr
      character*120  ext0*256, ext*256
      real(8):: qqq(3),qxx(3), vvvv(18)
      integer:: ifi,ifefclass,icors(2)
      complex(8),allocatable:: zegf(:,:) ,geig(:,:)
      complex(8),allocatable:: cphi(:,:)
      real(8),allocatable:: evl(:), vvv1(:),vvv2(:),vvv3(:),rofi_A(:)
     &,gcore_A(:,:), ec_A(:)
      integer,allocatable:: konf(:,:),nncx(:,:),ngvecp(:,:),lmxaa(:)
      real(8),parameter ::  rydberg=13.6058d0
      ! nocore is obtained by inquire(file='NoCore',exist=nocore) in the upper routine.
      ! If nocore exist. you have to supply
      !  <psi|Vxc(n_valence)|psi>  to  vxclda (nband, nqirr).
      ! If not, <psi|Vxc(n_total)|psi>  to  vxclda.
c----------------------------------------------
      integer:: ificg
      integer:: procid,nrank,ifigwb_,ifigwx1_,ifigwx2_
      integer:: iq,iqq,iqqx,nxxx,ifibz
      character*256:: extn,aaa,fname
      integer,parameter :: nsize= 1 000 000
      integer:: ifiproc,nqixx,nspxx,numprocxx,ixxx,ifiqibz
      logical:: iopen
      integer::  id,nsizex,iqqxx,ifile_handle,ib,ii,ipqn,nn,nnn(3),ifiqg,ifiqgc,irr,irrq,iqibz
!! =================================================================
      call setprint0(30)
      open(newunit=ifigwb,file='gwb.head',form='unformatted')
      read (ifigwb) nbas,nsp,ldim2,nbandmx,lmxamx,ncoremx,nrmx,plat,alat,nqirr
      allocate(bas(3,nbas),lmxaa(nbas),qplist(3,nqirr),ngplist(nqirr),ndimhall(nqirr))
      read(ifigwb) bas,lmxaa,qplist,ngplist,ndimhall
      close(ifigwb)
     
c$$$!! NLAindx augmentation indexing !  call rdibasindx2(ldim2,nphi,nphimx,nindx,lindx,ibasindx,ifi)
c$$$      open(newunit=ifi,file='NLAindx')
c$$$      allocate(nindx(ldim2),lindx(ldim2),ibasindx(ldim2))
c$$$      read(ifi,*)
c$$$      read(ifi,*) i
c$$$      if (ldim2 /= i) call rx('lmf2gw: ldim2/=i')
c$$$      nphimx = 0
c$$$      nn = 0
c$$$      do
c$$$          read(ifi,*,err=101,end=101) ipqn,l,ib,ii
c$$$          do m=1,2*l+1
c$$$            nn = nn+1
c$$$            nindx   (nn) = ipqn ! principle
c$$$            lindx   (nn) = l    ! l 
c$$$            ibasindx(nn) = ib   ! ibas
c$$$          enddo
c$$$          nphimx = max(nphimx,ipqn)
c$$$       enddo
c$$$ 101   continue
c$$$       close(ifi)
c$$$!! readin CLASS
c$$$      open(newunit=ifefclass,file='CLASS')
c$$$      nclass=0
c$$$      do ibas = 1,nbas
c$$$        read(ifefclass,*)  ibasx,  nclassx
c$$$        if(nclassx>nclass) nclass= nclassx
c$$$  enddo
      call readhamindex0()
      nclass=nclass_in
      allocate(nindx(ldim2),lindx(ldim2),ibasindx(ldim2))
      nindx=nindx_in
      lindx=lindx_in
      ibasindx=ibasindx_in
      nphimx=nphimx_in
      allocate( iclass(nbas) )
      iclass=iclass_in
      allocate(lmxa_d(nclass), nr(nclass), ncore_d(nclass), konf_d(0:lmxamx,nclass), zz(nclass),
     &  aa(nclass), bb(nclass) )
      do ibas=1,nbas
        lmxa_d(iclass(ibas)) = lmxaa(ibas)
      enddo
!! ATOMIC PART ic = ibas scheme ==,  GET nrxx and ncoremx ----------------------
      open(newunit=ifigwa,file='gwa',form='unformatted')
      allocate(nncx(0:lmxamx,nbas),konf(lmxamx+1,nbas),spid(nbas))
      allocate(ec_d(ncoremx, nclass, nsp),gx_d (nrmx, 0:lmxamx, nphimx,  nclass,nsp),
     & gcore_d(nrmx, ncoremx, nclass,nsp)  )
      do 3001 ibas = 1, nbas
        read(ifigwa) z, nr_A, a_A, b_A, rofi_Anr,lmxa,nspdummy,ncore,spid(ibas)
        allocate(rofi_A(nr_A), gcore_A(nr_A,ncore),ec_A(ncore))
        read(ifigwa) konf(1:lmxa+1,ibas)
        read(ifigwa) rofi_A(1:nr_A)
        write(6,"(' site',i3,'  z=',f5.1,'  rmax=',f8.5,'  lmax=',i1,'  konf=',10i1)")
     &       ibas,z,rofi_A(nr_A),lmxa,konf(1:lmxa+1,ibas)
        ic = iclass(ibas)
        zz(ic)= z
        aa(ic)= a_A
        bb(ic)= b_A
        nr(ic)= nr_A
        ncore_d(ic) = ncore/nsp
        konf_d(0:lmxa,ic) = konf(1:lmxa+1,ibas)
        write(6,"('  l    g(rmax)    gp(rmax)',4x,'<g g>',9x,'<gp gp>',9x,'<g gp>')")
        do  l = 0, lmxa
        do  isp = 1, nsp
           read(ifigwa) lxx,ispxx  
           if(lxx /= l .or. isp /=ispxx) call rx('lmf2gw:lxx or isp wrong')
           read(ifigwa) gx_d(1:nr_A,l,1,ic,isp) !phi
           read(ifigwa) gx_d(1:nr_A,l,2,ic,isp) !phidot
           if (konf_d(l,ic) .ge. 10) read(ifigwa) gx_d(1:nr_A,l,3,ic,isp) !phiz
c$$$           m = 2
c$$$           if (konf_d(l,ic) .ge. 10) m = 3
c$$$           n = 0
c$$$           do  i1 = 1, m
c$$$           do  i2 = 1, i1
c$$$             n = n+1
c$$$             call gintxx(gx_d(1,l,i1,ic,isp),gx_d(1,l,i2,ic,isp),a_A,b_A,nr_A,ovv(n))
c$$$           enddo
c$$$           enddo
c$$$           write(6,'(i3,2f11.5,2f15.10,f12.6)')l,gx_d(nr_A,l,1,ic,isp),gx_d(nr_A,l,2,ic,isp),
c$$$     &          ovv(1),ovv(3),ovv(2)
c$$$           if (m .eq. 3) write(6,"('  gz(rmax); <gz gz> <gz g> <gz gp> = ',f8.5,f15.10,2f12.7)")
c$$$     &          gx_d(nr_A,l,3,ic,isp), ovv(6),ovv(4),ovv(5)
        enddo
        enddo
!! core part
        if (ncore .ne. 0) write(6,'(''  l  k isp       ecore      gc(rmax)     <gc gc>'')')
        icore = 0
        icors = 0
        do isp = 1, nsp
        do l = 0, lmxa
            nncx(l,ibas) = mod(konf(l+1,ibas),10)-1 -(l+1) +1
            nnc          = max(nnc,nncx(l,ibas))
            do kkk = l+1, mod(konf(l+1,ibas),10)-1
              icore = icore+1
              icors(isp) = icors(isp) +1
              icor1=icors(isp)
              read(ifigwa) icorex,ldummy,ispx,kkkdummy,ec_A(icore)
              if(icore/=icorex)  call rx('lmf2gw:icore/=icorex')
              read(ifigwa) gcore_A(1:nr_A,icore) ! gcore
              ec_d(icor1, ic, isp) = ec_A(icore)
              gcore_d(1:nr_A,icor1,ic,isp)  = gcore_A(1:nr_A,icore)
            enddo
        enddo
        enddo
        deallocate(rofi_A,gcore_A,ec_A)
 3001 continue
!!      
      allocate(iantiferro(nbas))
      read(ifigwa)iantiferro(1:nbas) !iantiferro may2015
      close(ifigwa)
      allocate(ibasf(nbas)) !AF pair
      ibasf=-999
      do ibas=1,nbas
        do ibasx=ibas+1,nbas !is this fine?
        if(abs(iantiferro(ibas))/=0 .and. iantiferro(ibas)+iantiferro(ibasx)==0) then
           ibasf(ibas)=ibasx
           exit
        endif
        enddo
        if(ibasf(ibas)/=-999) write(6,"(a,2i5)")' AF pair: ibas ibasf(ibas)=',ibas,ibasf(ibas)
      enddo
      laf=.false.
      if(sum(abs(iantiferro))/=0) laf=.true. 
      write(6,"(a,100i4)") ' antiferro index=',iantiferro(1:nbas)
      do ibas=1,nbas
        write(6,"(a,i4,a)") ' i spid=',ibas,' '//trim(spid(ibas))
      enddo
!------
      open(newunit=ifiqg ,file='QGpsi',form='unformatted')
      open(newunit=ifiqgc,file='QGcou',form='unformatted')
      read(ifiqg ) nqnum, ngpmx,QpGcut_psi,nqbz,nqirr
      read(ifiqgc) nqnumc,ngcmx,QpGcut_cou
      nqtt = nqnum
      allocate(qtt(3,nqtt),ngvecptt(3,ngpmx,nqtt),ngvecctt(3,ngpmx,nqtt),ngptt(nqtt),ngctt(nqtt),iqindex(nqtt))
      irrq=0
      do iq=1,nqtt
         read (ifiqg)  qtt(1:3,iq), ngptt(iq) , irr
         read (ifiqg)  ngvecptt(1:3,1:ngptt(iq),iq)
         read (ifiqgc) qxx, ngctt(iq)
         read (ifiqgc) ngvecctt(1:3,1:ngctt(iq),iq)
         if(irr==1) then
            irrq=irrq+1
            iqindex(irrq)=iq
         endif   
         if(sum(abs(qtt(:,iq)-qxx))>1d-8) call rx('QGpsi QGcou q/=q')
         write(*,"(' qtt =',i4,3f9.5)")iq, qtt(1:3,iq)
      enddo
      close(ifiqg)
      close(ifiqgc)
      write(6,*)'QpGcut_psi QpGcutCou =',QpGcut_psi,QpGcut_Cou
      
      set_mnla :block
      integer :: ix, ibas,lx,nx,mx,ic,nvmax(0:lmxamx,nclass)
      write(*,*) '--- set_mnla ---'
      nvmax=0
      do ix =1,ldim2
        nx = nindx(ix)  
        lx   = lindx(ix)
        ibas =ibasindx(ix)
        ic = iclass(ibas)
        if( nx> nvmax(lx,ic) ) nvmax(lx,ic) = nx
      enddo
      write(6,"(5a5)") 'm','n','l','atom','ix'
      allocate(mnla(4,ldim2))
      ix=0
      do nx = 1,3
        do ibas = 1,nbas
          ic = iclass(ibas)
          do lx = 0,lmxa_d(ic)
            if (nx.gt.nvmax(lx,ic)) then
              cycle
            endif
            do mx = -lx,lx
              ix=ix+1
              mnla(1,ix)=mx
              mnla(2,ix)=nx
              mnla(3,ix)=lx
              mnla(4,ix)=ibas
              write(6,"(5i5)") mnla(1:4,ix),ix
            enddo
          enddo
        enddo
      enddo
      if(ix.ne.ldim2) then
        write(6,*) 'Error in set_mnla: ix!=ldim2'
        write(6,*) 'ix2,ldim2=',ix,ldim2
        call rx('lmf2gw: Error in set_mnla: ix!=ldim2')
      endif
      endblock set_mnla

      end subroutine lmf2gw
      end module m_lmf2gw
