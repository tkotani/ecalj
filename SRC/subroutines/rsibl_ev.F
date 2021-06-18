Ckino from rsibl write eigenvectors
      module m_rsibl_ev
      contains
      subroutine rsibl_ev(ssite,sspec,nbas,isp,q,iq,ndimh,nspc, !,slat
     .napw,igapw,iprmb,nevec,evec,k1,k2,k3, n_eiglist,eiglist)
      use m_struc_def 
      use m_w_psir
      use m_lmfinit,only:
     &     lat_alat,nspec
      use m_lattic,only: lat_qlat
      use m_lattic,only: lat_vol
      use m_supot,only: lat_nabc
      use m_supot,only: lat_ng
      use m_supot,only: lat_gmax
      use m_lattic,only: lat_plat,rv_a_opos
C- Add smooth part of output density into smrho and forces.
C ----------------------------------------------------------------------
Ci Inputs
Ci  x lfrce :if nonzero, accumulate contribution to force
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: rsibl1
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: ngcut
Ci     Stored:    *
Ci     Passed to: tbhsi rsibl1 uspecb
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat plat qlat gmax nabc ng ogv okv vol
Ci     Stored:    *
Ci     Passed to: *
Ci   lfrce :1 calculate contribution to forces
Ci   nbas  :size of basis
Ci   q     :Bloch vector
Ci   iq    :index to current k-point
Ci   ndimh :dimension of hamiltonian
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   napw  :number of augmented PWs in basis
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nevec :number of eigenvectors with nonzero weights
Ci   evec  :eigenvectors
Ci   ewgt  :eigenvector weights
Ci   k1..3 :dimensions smpot,smrho
Ci   smpot :smooth potential on uniform mesh, needed for forces
Co Outputs
Co   smrho :smooth density accumulated for this qp
Co   f     :force contribution accumulated for this qp
Cr Remarks
Cm MPI
Cm   Parallelise over the eigenvector loop. The vector block size is
Cm   chosen (in the range 6-16, by dstrbp.f) so as to distribute the
Cm   work optimally across processes. Two work arrays of the size of
Cm   smrho are allocated from the heap as buffers. Only one will be
Cm   needed under MPI-2. See comments in hsibl.
Cb Bugs
Cb    replace call to gvgvcomp and pass ipv as input
Cb    The non-F90 version should work, but it is no longer tested
Cu Updates
Cu   29 Dec 08 Unsuccessful attempt to make work with openmp
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu   10 Sep 06 Added MPI parallelization in the spin-coupled case
Cu   23 Dec 04 Extended to spin-coupled case
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu   12 Oct 00 Use q-dependent list of G vectors
Cu    6 Jul 00 attempt to vectorize by grouping eigenvectors in blocks
Cu   17 Jun 00 Spin polarized
Cu   23 May 00 Adapted from nfp rsif_q.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer procid, master, nproc, mpipid
      integer isp,k1,k2,k3,ndimh,nevec,iprmb(*),iq,nspc
      integer napw,igapw(3,napw)
      real(8):: q(3) 
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
c      type(s_lat)::slat
      double complex evec(ndimh,nspc,nevec)
      integer:: n_eiglist
      integer:: eiglist(n_eiglist)
C ... Local parameters
      integer n0,nkap0,nermx,npmx,nblk,nlmto
      parameter (n0=10,nkap0=3,nermx=100,npmx=128)
        integer:: ngabc(3) , n1 , n2 , n3 , nrt , net , ng , 
     .  nglob , ltop , nlmtop , ogq , og2 , ohe , ohr , oyl , oylw , 
     .  oiv , iprint 
       integer,allocatable :: iv_a_okv(:) 
       real(8),allocatable :: rv_a_ogv(:) 
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx)
      double precision alat,qlat(3,3),plat(3,3),q0(3),gmax,xx
      double precision vol
      double precision etab(nermx),rtab(nermx)
      integer ivec,nvec
      integer,allocatable:: ivp(:)
      complex(8),allocatable::psi(:,:,:),psir(:,:,:),vpsi(:,:,:), wk(:,:,:)
      real(8),allocatable:: cosi(:),sini(:),wk2(:)
      integer:: ivecini,ivecend,nbas
      integer,allocatable:: w_oiv(:)
      real(8),allocatable:: w_ogq(:),w_oyl(:),w_oylw(:),w_og2(:),w_ohe(:),w_ohr(:)
      complex(8),allocatable:: w_osmbuf(:)
      real(8),allocatable:: w_ofrbuf(:)
      real(8),allocatable:: ewgt(:)
      complex(8),allocatable:: smrho(:,:,:,:), smpot(:,:,:,:)
      logical:: cmdopt
      real(8):: w(1)!dummy
      nproc  = mpipid(0)
      procid = mpipid(1)
      if (nevec .le. 0) return
      nlmto = ndimh-napw
      alat = lat_alat
      plat = lat_plat
      qlat = lat_qlat
      gmax = lat_gmax
      ngabc= lat_nabc
      ng = lat_ng
      vol= lat_vol
C ... Setup for q-dependent gv ... also makes kv, gv+q and iv
C     NB: gv generated by gvlst2 has q already added to it!
      call pshpr(iprint()-30)
c      call gvlist(alat,plat,q,n1,n2,n3,gmax,500,0,ng,xx,xx,xx,xx)
      call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,0,500,0,ng,xx,xx,xx,xx)
      allocate(rv_a_ogv(abs(ng*3)))
      rv_a_ogv(:)=0.0d0
      allocate(iv_a_okv(abs(ng*3)))
      iv_a_okv(:)=0
      allocate(smrho(k1,k2,k3,isp))
      allocate(smpot(k1,k2,3,isp))
      allocate(w_oiv(ng*3))
c      call gvlist(alat, plat, q, n1, n2, n3, gmax,509, ng, ng, iv_a_okv, rv_a_ogv , w_oiv , w_oiv )
      call gvlst2(alat, plat, q, n1, n2, n3, 0d0,gmax,0,509, ng, ng, iv_a_okv, rv_a_ogv, w_oiv, w_oiv )
      call poppr
C     For PW basis ... for now.
      if (napw .gt. 0) then
         allocate(ivp(napw))
         call gvgvcomp(ng,w_oiv,napw,igapw,ivp)
      else
         allocate(ivp(1))
      endif
C --- Tables of energies, rsm, indices to them ---
      call tbhsi(sspec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)
C --- Allocate and occupy arrays for yl, energy factors, rsm factors ---
      nlmtop = (ltop+1)**2
      allocate(w_ogq(ng*3), w_oyl(ng*nlmtop), w_oylw(ng*nlmtop), w_og2(ng), w_ohe(ng*net), w_ohr(ng*nrt))

Ckino H_L(G)= \frac{-4 pi}{e-G^2} {cal Y}_L(-iG) exp(gamma(e-G^2))
Ckino hsibl1 calculaets he=1/(e-G^2) and hr=exp(-gamma G^2)
Ckino the other parts are calculated in rsibl5.
      call dpzero(q0,3)
      if (nlmto .gt. 0) then
         call hsibl1 ( net , etab , nrt , rtab , ltop , alat , q0 , ng 
     .   , rv_a_ogv , w_ogq , w_og2 , w_oyl , w_ohe ,  w_ohr ) 
      endif
      deallocate(w_og2)
      nblk = 16 !nblk is for futre parallelization
C  --- Loop over eigenstates ---
      allocate(psi(ng,nspc,nblk),vpsi(ng,nspc,nblk),wk(ng,nspc,nblk))
      allocate(psir(k1,k2,k3),cosi(ng),sini(ng),wk2(ng))
      ivecini= 1
      ivecend= nevec
      allocate(ewgt(nevec))
      ewgt=1.0d0
      write(*,'(a,9i5)') 'ivecinic,ivecend,nblk,nspc,nspec=',ivecini,ivecend,nblk,nspc,nspec
      ivecloop: do  ivec = ivecini,ivecend, nblk
      nvec = min(nblk, nevec-ivec+1) 
C   ... Add together total smooth wavefunction
ckino   rsibl1 calculates fourier transformed smoothed Hankel, H(G), to psi
        call rsibl1(0,ssite,sspec,q,nbas,iprmb,ng,w_ogq,w_oiv,n1,n2,
     .  n3,qlat,cosi,sini,w_oyl,w_oylw,w_ohe,w_ohr,wk,
     .  wk2,vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc,
     .  ewgt,ivec,nvec,evec,w,psi,w)
Ckino    rsiblp adds PW(G) to psi 
        call rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,dsqrt(vol),
     .  evec(1,1,ivec),psi)
Ckino now psi= H(G) + PW(G)
Ckino write psi=F0 part
         if (.true.) then
             call w_psir(ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2
     i   , k3 , iv_a_okv , isp, q ,iq 
     i   ,  n_eiglist,eiglist
     i   , plat,alat, nbas, rv_a_opos, sspec(ssite(1:nbas)%spec)%z
     w  , psir )
         endif
      enddo ivecloop
      deallocate(psi,vpsi,wk,psir,cosi,sini,wk2)
       if (allocated(ewgt)) deallocate(ewgt)
       if (allocated(rv_a_ogv)) deallocate(rv_a_ogv)
       if (allocated(iv_a_okv)) deallocate(iv_a_okv) 
      deallocate(ivp)
      end subroutine rsibl_ev
      end module m_rsibl_ev
