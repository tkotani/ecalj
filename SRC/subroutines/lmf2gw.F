!!cr Radial mesh is specified from, nr, aa, and bb; mesh is r(i) = bb*(exp(aa*(i-1)) -1) ,i=1,nr
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
!!     evl    (nband, nqbze, nsp), ! Eigenvalues
!!     vxclda (nband, nqbze, nsp), ! Lda XC potential <psi|Vxc(n_total)  |psi>
!!     gx (nrmx, 0:lmxamx, nphimx, nclass,nsp), ! Radial function.
!!         gx (nrmx, lmxamx, 1, nclass,nsp) is phi
!!         gx (nrmx, lmxamx, 2, nclass,nsp) is phidot
!!     gcore(nrmx, ncoremx, nclass,nsp)  ! Core radial function.
!!c    These radial functions are only the major part given by the scalar relativistic calculations.
!!c     gx and gcore = r \phi(r) = u(r), where \phi is the major part of the true radial functions.
!!c     gx is normalized as 1 = \int dr gx**2
!!c     gcore is the major part of the true radial functions.
!!      cphi(ldim2,  nband, nqbze,nsp), ! Coefficients of eigenfunction
!!      geig(ngpmx,  nband, nqbze,nsp) ! Coefficients of eigenfunction for IPW.
      
      module m_lmf2gw
      integer,allocatable,protected:: nindx(:),lindx(:),ibasindx(:),iantiferro(:)
      integer,protected :: nbandmx,nphimx,
     &nsp,      ! =1 or 2, corresponding to para or ferro.
     &nbas,     ! Number of atom in the primitive cell
     &nclass,   ! Number of atomic class (or type) for in the primitive cell
     &nrmx,     ! = maxval(nr(1:nclass))  Maximum number of nr
     &ncoremx,  ! = maxval(ncore(1:nclass))
     &lmxamx,   ! = maxval(lmxa(1:nclass))
     &ngpmx,    ! Maximum number of G vector.
     &ldim2,    ! = total number of augmentation functions nRlm
c     &nphi,     ! = total number of augmentation functions nRl
     &nqbze     ! = nqbz*(1+nq0i). Number of q points given by qg4gw
      character(8),allocatable,protected:: spid(:)
      integer,allocatable::
     &     iclass(:),lmxa_d  (:),nr(:),konf_d(:,:),ncore_d(:)
      real(8),allocatable ::
     &     zz(:),aa(:),bb(:),ec_d (:,:,:),evl_d(:,:,:),vxclda (:,:,:),
     &     gx_d(:,:,:,:,:),gcore_d(:,:,:,:)
      real(8):: plat(3,3), alat, efermi
      complex(8),allocatable:: cphi_d(:,:,:,:)
      complex(8),allocatable:: geig_d(:,:,:,:)
      real(8),allocatable :: qbze(:,:), bas(:,:)
c      logical :: nocore         ! Flag file NoCore exist or not.

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
      implicit none
      integer::
     &ldim,     ! = sum ( (lmxa(1:nbas)+1)**2 )
     &nband    ! Number of bands given by GWIN0
      integer icor1,icorex,i,i1,i2,ibas,ibasx,ibx,ic,icore,
     .ifichkv,ifigw0,ifigwa,ifigwb,ifigwx1,ifigwx2,ifigwx3,isp,ispx,
     .ispxx,ix,kkk,kkkdummy,l,ldummy,lmxa,lxx,n1q,n2q,n3q,nclassx,m,n,
     .ncore,ndimh,ngp,nnc,nq0i,nqbz,nspdummy,IKP,NR_A !takao feb2012 ngc,ngcmx,
      real(8) qlat(3,3),ovv(20),ef0,qpgcut_phi,qpgcut_cou,z,
     .a_a,b_a,rofi_anr
      character*120  ext0*256, ext*256
      real(8):: qqq(3)
      integer:: ifi,ifefclass,icors(2)
      complex(8),allocatable:: zegf(:,:) ,geig(:,:)
      complex(8),allocatable:: cphi(:,:)
      real(8),allocatable:: evl(:), vvv1(:),vvv2(:),vvv3(:),rofi_A(:)
     &,gcore_A(:,:), ec_A(:)
      integer,allocatable:: konf(:,:),lmxaa(:),nncx(:,:),ngvecp(:,:)
      real(8),parameter ::  rydberg=13.6058d0
      ! nocore is obtained by inquire(file='NoCore',exist=nocore) in the upper routine.
      ! If nocore exist. you have to supply
      !  <psi|Vxc(n_valence)|psi>  to  vxclda (nband, nqbze).
      ! If not, <psi|Vxc(n_total)|psi>  to  vxclda.
      real(8),allocatable :: wt(:),q0i(:,:)
c----------------------------------------------
      integer:: ificg
      integer:: procid,nrank,ifigwb_,ifigwx1_,ifigwx2_
      integer:: iq,iqq,iqqx,nqirr
      character*256:: extn,aaa,fname
      character*256,allocatable:: extp(:)
      integer,parameter :: nsize= 1 000 000
      integer:: ifiproc,nqixx,nspxx,numprocxx,ixxx
      integer,allocatable:: iprocq(:,:)
      logical:: iopen
      integer  id,nsizex,iqqxx,ifile_handle,ib,ii,ipqn,nn
!! =================================================================
      call setpr(30)
      open(newunit=ifigwb,file='gwb.head',form='unformatted')
      read (ifigwb) nbas,nsp,ldim2,ndimh,alat,qlat,ef0,nqbz,plat(1:3,1:3),nqbze, nqirr !feb2012takao
      nqbze=nqirr    !nqbze means number of calculated bands (irreducible only)         feb2012takao
      allocate(lmxaa(1:nbas), bas(1:3,1:nbas))
      read (ifigwb) lmxaa(1:nbas), bas(1:3,1:nbas)
      read (ifigwb) ncoremx,nrmx
      close(ifigwb)
      lmxamx = maxval(lmxaa(1:nbas))
      allocate( iclass(nbas) )
!! NLAindx augmentation indexing !  call rdibasindx2(ldim2,nphi,nphimx,nindx,lindx,ibasindx,ifi)
      open(newunit=ifi,file='NLAindx')
      allocate(nindx(ldim2),lindx(ldim2),ibasindx(ldim2))
      read(ifi,*)
      read(ifi,*) i
      if (ldim2 /= i) call rx('lmf2gw: ldim2/=i')
      nphimx = 0
      nn = 0
      do
          read(ifi,*,err=101,end=101) ipqn,l,ib,ii
          do m=1,2*l+1
            nn = nn+1
            nindx   (nn) = ipqn ! principle
            lindx   (nn) = l    ! l 
            ibasindx(nn) = ib   ! ibas
          enddo
          nphimx = max(nphimx,ipqn)
       enddo
 101   continue
       close(ifi)
!! readin CLASS
      open(newunit=ifefclass,file='CLASS')
      nclass=0
      do ibas = 1,nbas
        read(ifefclass,*)  ibasx,  nclassx
        if(nclassx>nclass) nclass= nclassx
      enddo
      write(6,'('' File CLASS:  nclass ='',i4)') nclass
      allocate( lmxa_d(nclass), nr(nclass),
     &  ncore_d(nclass), konf_d(0:lmxamx,nclass), zz(nclass),
     &  aa(nclass), bb(nclass) )
      rewind ifefclass
      do ibas=1,nbas
        read(ifefclass,*) ix,iclass(ibas)
        lmxa_d(iclass(ibas)) =lmxaa(ibas)
      enddo
      close(ifefclass)
!! Band Part
      open(newunit=ifiproc,file= 'lmfgw_kdivider')
      read(ifiproc,*) ext  !extension is read here
      read(ifiproc,*) nqixx, nspxx, nrank
      allocate(iprocq(nqixx,nspxx))
      do isp=1,nspxx
      do iqq=1,nqixx
        read(ifiproc,*) iqqxx, ispxx, ixxx
        if(iqqxx/=iqq) stop 'iqqxx/=iqq'
        if(ispxx/=isp) stop 'ispxx/=isp'
        iprocq(iqq,isp) = ixxx
        write(6,"('iqq isp irank=',i8,i2,i6)") iqq,isp, iprocq(iqq,isp)
      enddo
      enddo
      close(ifiproc)
      allocate(extp(0:nrank-1))
      extp(0) = trim(ext)
      write(6,"('  0 ext= ',a,a)") trim(extp(0)),' ----------'
      do procid=1,nrank-1
        write(extn,"(i10)") procid
        extp(procid)=trim(adjustl(ext))//'_'//trim(adjustl(extn))
        write(6,"(i3,' ext= ',a,a)") procid,trim(extp(procid)),' ----------'
      enddo
!!      
      open(newunit=ifichkv,file='VXCFP.chk')
      write(ifichkv,"(a)") '### LDA exchange correlation ###'
      write(ifichkv,"(a)")'#   qvec                 ikp iband    eigen VXC(ntotal)     VXC(nvalence)'
      open(newunit=ificg,file='CphiGeig',form='unformatted')
c      mrecb = 2*ldim2*nband *ndble/nwordr()
c      open(newunit=ifcphi,file='CPHI',form='unformatted',access='direct',recl=mrecb)
      nbandmx = ndimh
      write(6,'(" nqbze nqpz ndimh(max)=",3i5)') nqbze,nqbz,ndimh
      allocate( qbze(3,nqbze),evl_d(nbandmx, nqbze, nsp),vxclda(nbandmx, nqbze, nsp),
     &  vvv1(nbandmx),vvv2(nbandmx),vvv3(nbandmx),cphi(ldim2, nbandmx), evl(nbandmx) )
      evl_d =1d20
      vxclda=1d20
      iqq=0
      do ikp = 1,nqbze
        iqq=iqq+1
        do isp = 1, nsp
          fname='gwb'//trim(extp(iprocq(iqq,isp)))
          inquire(file=trim(fname),opened=iopen,number=ifigwb_)
          if(.not.iopen) open(newunit=ifigwb_,file=trim(fname),form='unformatted')
          read (ifigwb_) qbze(:,ikp),nband,iq
          write(6,"(a,3i6,3x,3f11.6)") 'ikp isp nband q=',ikp,isp,nband, qbze(:,ikp)
          evl =1d20
          cphi=1d20
          read (ifigwb_) evl(1:nband), cphi(:,1:nband)
          write(ificg) cphi(1:ldim2,1:nbandmx)
          evl_d (1:nband, ikp,isp) = evl (1:nband)
!! Plane wave part
          read (ifigwb_) ngp !,ngc
          allocate( geig(ngp, nbandmx), ngvecp(3,ngp) )
          geig=1d20
          read(ifigwb_) geig(:,1:nband) !ngvecp,
          write(ificg)  ngp
          if(ngp>0) write(ificg)  geig(1:ngp,1:nbandmx)
          deallocate(geig,ngvecp) !takao feb2012 ,ngvecc )
          fname='gw1'//trim(extp(iprocq(iqq,isp)))
          inquire(file=trim(fname),opened=iopen,number=ifigwx1_)
          if(.not.iopen) open(newunit=ifigwx1_,file=trim(fname),form='unformatted')
          fname='gw2'//trim(extp(iprocq(iqq,isp)))
          inquire(file=trim(fname),opened=iopen,number=ifigwx2_)
          if(.not.iopen) open(newunit=ifigwx2_,file=trim(fname),form='unformatted')
          read(ifigwx1_ ) qqq, (vvv1(i),i=1,nband)
          read(ifigwx2_ ) qqq, (vvv2(i),i=1,nband)
          vxclda (1:nband, ikp,isp)= vvv1(1:nband)-vvv2(1:nband)
          vvv3=vvv2
          do ibx =1,nband
             write(ifichkv,"(3f8.4,2i3,6f16.8)") qqq,ikp,ibx,
     &            vvv1(ibx),vvv1(ibx)-vvv2(ibx),vvv3(ibx)-vvv2(ibx),
     &            rydberg*vvv1(ibx),rydberg*(vvv1(ibx)-vvv2(ibx)),
     &            rydberg*(vvv3(ibx)-vvv2(ibx))
          enddo
 888      continue
c we can not close here because extp(iprocq(iqq,isp)) can be the same for different iqq,isp
c          close(ifigwx1_),  close(ifigwx2_),   close(ifigwb_)
        enddo
      enddo
      close(ifichkv)
      close(ificg)
      deallocate( cphi,evl, vvv1,vvv2,vvv3)
!! ATOMIC PART ic = ibas scheme ==,  GET nrxx and ncoremx ----------------------
c      nrmx    = 0
c      ncoremx = 0
      open(newunit=ifigwa,file='gwa'//trim(ext),form='unformatted')
      allocate(nncx(0:lmxamx,nbas),konf(lmxamx+1,nbas),spid(nbas))
      allocate(ec_d(ncoremx, nclass, nsp),gx_d (nrmx, 0:lmxamx, nphimx,  nclass,nsp),
     & gcore_d(nrmx, ncoremx, nclass,nsp)  )
      do ibas = 1, nbas
        read(ifigwa) z, nr_A, a_A, b_A, rofi_Anr
     &  ,lmxa,nspdummy,ncore,spid(ibas)
        allocate(rofi_A(nr_A),ec_A(ncore))
        read(ifigwa) konf(1:lmxa+1,ibas)
        read(ifigwa) rofi_A(1:nr_A)
        do  l = 0, lmxa
          do  isp = 1, nsp
            read(ifigwa) lxx,ispx  !!!ibas==ic,l,n=1
            read(ifigwa)
            read(ifigwa)
            if (konf(l+1,ibas) .ge. 10) read(ifigwa)
          enddo
        enddo
        do  l = 0, lmxa
          do  isp = 1, nsp
            nncx(l,ibas) = mod(konf(l+1,ibas),10)-1 -(l+1) +1
            nnc          = max(nnc,nncx(l,ibas))
            do  kkk = l+1, mod(konf(l+1,ibas),10)-1
              icore = icore+1
              read(ifigwa) icore,ldummy,ispx,kkkdummy,ec_A(icore)
              read(ifigwa)
            enddo
          enddo
        enddo
        deallocate(rofi_A,ec_A)
c        if(nr_A   >nrmx   ) nrmx    = nr_A
c        if(ncore  >ncoremx) ncoremx = ncore
      enddo
c      nphi = ldim2
c      ncoremx = ncoremx/nsp
      gx_d = 0
      rewind ifigwa
      nnc =0
      konf =0
      write(6,'(/'' File gwa: reading atom data'')')
      do 3001 ibas = 1, nbas
        read(ifigwa) z, nr_A, a_A, b_A, rofi_Anr,lmxa,nspdummy,ncore
        allocate(rofi_A(nr_A),
     &  gcore_A(nr_A,ncore),ec_A(ncore))
        read(ifigwa) konf(1:lmxa+1,ibas)
        read(ifigwa) rofi_A(1:nr_A)
        write(6,369) ibas,z,rofi_A(nr_A),lmxa,konf(1:lmxa+1,ibas)
  369   format(' site',i3,'  z=',f5.1,'  rmax=',f8.5,'  lmax=',i1,'  konf=',10i1)
!!       A bit rough treatment
        ic    = iclass(ibas)
        zz(ic)= z
        aa(ic)= a_A
        bb(ic)= b_A
        nr(ic)= nr_A
        ncore_d(ic) = ncore/nsp
        konf_d(0:lmxa,ic) = konf(1:lmxa+1,ibas)
        write(6,371)
  371   format('  l    g(rmax)    gp(rmax)',4x,'<g g>',9x,'<gp gp>',9x,'<g gp>')
!! --- PHIU
        do  l = 0, lmxa
          do  isp = 1, nsp
            read(ifigwa) lxx,ispxx  !!!ibas==ic,l,n=1
            if(lxx /= l .or. isp /=ispxx) stop 'hhbecfp:lxx or isp wrong'
            read(ifigwa) gx_d(1:nr_A,l,1,ic,isp)
            read(ifigwa) gx_d(1:nr_A,l,2,ic,isp)
            if (konf_d(l,ic) .ge. 10) read(ifigwa) gx_d(1:nr_A,l,3,ic,isp)
            m = 2
            if (konf_d(l,ic) .ge. 10) m = 3
            n = 0
            do  i1 = 1, m
              do  i2 = 1, i1
                n = n+1
                call gintxx(gx_d(1:nr_A,l,i1,ic,isp),
     .          gx_d(1:nr_A,l,i2,ic,isp),a_A,b_A,nr_A,ovv(n))
              enddo
            enddo
            write(6,'(i3,2f11.5,2f15.10,f12.6)')
     .      l,gx_d(nr_A,l,1,ic,isp),gx_d(nr_A,l,2,ic,isp),
     .      ovv(1),ovv(3),ovv(2)
            if (m .eq. 3) write(6,321) gx_d(nr_A,l,3,ic,isp),
     .      ovv(6),ovv(4),ovv(5)
  321       format('  gz(rmax); <gz gz> <gz g> <gz gp> = ',f8.5,
     .      f15.10,2f12.7)

          enddo
        enddo
!! --- ECORE PHICU
        if (ncore .ne. 0)
     .  write(6,'(''  l  k isp       ecore      gc(rmax)     <gc gc>'')')
        icore = 0
        icors = 0
        do isp = 1, nsp
          do l = 0, lmxa
            nncx(l,ibas) = mod(konf(l+1,ibas),10)-1 -(l+1) +1
            nnc          = max(nnc,nncx(l,ibas))
            do  kkk = l+1, mod(konf(l+1,ibas),10)-1
              icore = icore+1
              icors(isp) = icors(isp) +1
              icor1=icors(isp)
              read(ifigwa) icorex,ldummy,ispx,kkkdummy,ec_A(icore)
              if(icore/=icorex)  stop 'lmf2gw:icore/=icorex'
              read(ifigwa) gcore_A(1:nr_A,icore)
              ec_d(icor1, ic, isp) = ec_A(icore)
              gcore_d(1:nr_A,icor1,ic,isp)  = gcore_A(1:nr_A,icore)
              call gintxx(gcore_A(1,icore),gcore_A(1,icore), a_A,b_A,nr_A,ovv(1))
              write(6,'(3i3,f15.6,f12.6,f15.10,3i3)')
     &        l,kkk,isp,ec_A(icore),gcore_A(nr_A,icore),ovv(1),icore,icor1
            enddo
          enddo
        enddo
        deallocate(rofi_A,gcore_A,ec_A)
 3001 continue
      allocate(iantiferro(nbas))
      read(ifigwa)iantiferro(1:nbas) !iantiferro may2015
      close(ifigwa)
      write(6,"(a,100i4)") ' antiferro index=',iantiferro(1:nbas)
      do ibas=1,nbas
        write(6,"(a,i4,a)") ' i spid=',ibas,' '//trim(spid(ibas))
      enddo
      close(ifi)
      end subroutine lmf2gw
      end module m_lmf2gw
