! ino from rsibl write eigenvectors
module m_rsibl_ev
  public rsibl_ev
  private
contains
  subroutine rsibl_ev(ssite,sspec,nbas,isp,q,iq,ndimh,nspc,&
    napw,igapw,iprmb,nevec,evec,k1,k2,k3, n_eiglist,eiglist)
    use m_struc_def
    use m_w_psir
    use m_lmfinit,only:  lat_alat,nspec
    use m_lattic,only: lat_qlat
    use m_lattic,only: lat_vol
    use m_supot,only: lat_nabc
    use m_supot,only: lat_ng
    use m_supot,only: lat_gmax
    use m_lattic,only: lat_plat,rv_a_opos
    !- Add smooth part of output density into smrho and forces.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i  x lfrce :if nonzero, accumulate contribution to force
    !i   nbas  :size of basis
    !i   ssite :struct for site-specific information; see routine usite
    !i     Elts read: spec pos
    !i     Stored:    *
    !i     Passed to: rsibl1
    !i   sspec :struct for species-specific information; see routine uspec
    !i     Elts read: ngcut
    !i     Stored:    *
    !i     Passed to: tbhsi rsibl1 uspecb
    !i   slat  :struct for lattice information; see routine ulat
    !i     Elts read: alat plat qlat gmax nabc ng ogv okv vol
    !i     Stored:    *
    !i     Passed to: *
    !i   lfrce :1 calculate contribution to forces
    !i   nbas  :size of basis
    !i   q     :Bloch vector
    !i   iq    :index to current k-point
    !i   ndimh :dimension of hamiltonian
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   napw  :number of augmented PWs in basis
    !i   igapw :vector of APWs, in units of reciprocal lattice vectors
    !i   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
    !i   nevec :number of eigenvectors with nonzero weights
    !i   evec  :eigenvectors
    !i   ewgt  :eigenvector weights
    !i   k1..3 :dimensions smpot,smrho
    !i   smpot :smooth potential on uniform mesh, needed for forces
    !o Outputs
    !o   smrho :smooth density accumulated for this qp
    !o   f     :force contribution accumulated for this qp
    !r Remarks
    !m MPI
    !m   Parallelise over the eigenvector loop. The vector block size is
    !m   chosen (in the range 6-16, by dstrbp.f) so as to distribute the
    !m   work optimally across processes. Two work arrays of the size of
    !m   smrho are allocated from the heap as buffers. Only one will be
    !m   needed under MPI-2. See comments in hsibl.
    !b Bugs
    !b    replace call to gvgvcomp and pass ipv as input
    !b    The non-F90 version should work, but it is no longer tested
    !u Updates
    !u   29 Dec 08 Unsuccessful attempt to make work with openmp
    !u   05 Jul 08 (T. Kotani) output density for new PW part
    !u   10 Sep 06 Added MPI parallelization in the spin-coupled case
    !u   23 Dec 04 Extended to spin-coupled case
    !u   25 Aug 04 Adapted to extended local orbitals
    !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
    !u   15 Feb 02 (ATP) Added MPI parallelization
    !u   27 Aug 01 Extended to local orbitals.
    !u   12 Oct 00 Use q-dependent list of G vectors
    !u    6 Jul 00 attempt to vectorize by grouping eigenvectors in blocks
    !u   17 Jun 00 Spin polarized
    !u   23 May 00 Adapted from nfp rsif_q.f
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: procid, master, nproc, mpipid
    integer :: isp,k1,k2,k3,ndimh,nevec,iprmb(*),iq,nspc
    integer :: napw,igapw(3,napw)
    real(8):: q(3)
    type(s_site)::ssite(*)
    type(s_spec)::sspec(*)
    !      type(s_lat)::slat
    double complex evec(ndimh,nspc,nevec)
    integer:: n_eiglist
    integer:: eiglist(n_eiglist)
    ! ... Local parameters
    integer :: n0,nkap0,nermx,npmx,nblk,nlmto
    parameter (n0=10,nkap0=3,nermx=100,npmx=128)
    integer:: ngabc(3) , n1 , n2 , n3 , nrt , net , ng , &
         nglob , ltop , nlmtop , ogq , og2 , ohe , ohr , oyl , oylw , &
         oiv , iprint
    integer,allocatable :: iv_a_okv(:)
    real(8),allocatable :: rv_a_ogv(:)
    equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
    integer :: iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx)
    double precision :: alat,qlat(3,3),plat(3,3),q0(3),gmax,xx
    double precision :: vol
    double precision :: etab(nermx),rtab(nermx)
    integer :: ivec,nvec
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
    if (nevec <= 0) return
    nlmto = ndimh-napw
    alat = lat_alat
    plat = lat_plat
    qlat = lat_qlat
    gmax = lat_gmax
    ngabc= lat_nabc
    ng = lat_ng
    vol= lat_vol
    ! ... Setup for q-dependent gv ... also makes kv, gv+q and iv
    !     NB: gv generated by gvlst2 has q already added to it!
    call pshpr(iprint()-30)
    !      call gvlist(alat,plat,q,n1,n2,n3,gmax,500,0,ng,xx,xx,xx,xx)
    call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,0,500,0,ng,xx,xx,xx,xx)
    allocate(rv_a_ogv(abs(ng*3)))
    rv_a_ogv(:)=0.0d0
    allocate(iv_a_okv(abs(ng*3)))
    iv_a_okv(:)=0
    allocate(smrho(k1,k2,k3,isp))
    allocate(smpot(k1,k2,3,isp))
    allocate(w_oiv(ng*3))
    !      call gvlist(alat, plat, q, n1, n2, n3, gmax,509, ng, ng, iv_a_okv, rv_a_ogv , w_oiv , w_oiv )
    call gvlst2(alat, plat, q, n1, n2, n3, 0d0,gmax,0,509, ng, ng, iv_a_okv, rv_a_ogv, w_oiv, w_oiv )
    call poppr
    !     For PW basis ... for now.
    if (napw > 0) then
       allocate(ivp(napw))
       call gvgvcomp(ng,w_oiv,napw,igapw,ivp)
    else
       allocate(ivp(1))
    endif
    ! --- Tables of energies, rsm, indices to them ---
    call tbhsi(sspec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)
    ! --- Allocate and occupy arrays for yl, energy factors, rsm factors ---
    nlmtop = (ltop+1)**2
    allocate(w_ogq(ng*3), w_oyl(ng*nlmtop), w_oylw(ng*nlmtop), w_og2(ng), w_ohe(ng*net), w_ohr(ng*nrt))

    ! ino H_L(G)= \frac{-4 pi}{e-G^2} {cal Y}_L(-iG) exp(gamma(e-G^2))
    ! ino hsibl1 calculaets he=1/(e-G^2) and hr=exp(-gamma G^2)
    ! ino the other parts are calculated in rsibl5.
    call dpzero(q0,3)
    if (nlmto > 0) then
       call hsibl1 ( net , etab , nrt , rtab , ltop , alat , q0 , ng &
            , rv_a_ogv , w_ogq , w_og2 , w_oyl , w_ohe ,  w_ohr )
    endif
    deallocate(w_og2)
    nblk = 16 !nblk is for futre parallelization
    !  --- Loop over eigenstates ---
    allocate(psi(ng,nspc,nblk),vpsi(ng,nspc,nblk),wk(ng,nspc,nblk))
    allocate(psir(k1,k2,k3),cosi(ng),sini(ng),wk2(ng))
    ivecini= 1
    ivecend= nevec
    allocate(ewgt(nevec))
    ewgt=1.0d0
    write(*,'(a,9i5)') 'ivecinic,ivecend,nblk,nspc,nspec=',ivecini,ivecend,nblk,nspc,nspec
    ivecloop: do  ivec = ivecini,ivecend, nblk
       nvec = min(nblk, nevec-ivec+1)
       !   ... Add together total smooth wavefunction
       ! ino   rsibl1 calculates fourier transformed smoothed Hankel, H(G), to psi
       call rsibl1(0,ssite,sspec,q,nbas,iprmb,ng,w_ogq,w_oiv,n1,n2, &
            n3,qlat,cosi,sini,w_oyl,w_oylw,w_ohe,w_ohr,wk, &
            wk2,vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc, &
            ewgt,ivec,nvec,evec,w,psi,w)
       ! ino    rsiblp adds PW(G) to psi
       call rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,dsqrt(vol), &
            evec(1,1,ivec),psi)
       ! ino now psi= H(G) + PW(G)
       ! ino write psi=F0 part
       if ( .TRUE. ) then
          call w_psir(ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2 &
               , k3 , iv_a_okv , isp, q ,iq &
               ,  n_eiglist,eiglist &
               , plat,alat, nbas, rv_a_opos, sspec(ssite(1:nbas)%spec)%z &
               , psir )
       endif
    enddo ivecloop
    deallocate(psi,vpsi,wk,psir,cosi,sini,wk2)
    if (allocated(ewgt)) deallocate(ewgt)
    if (allocated(rv_a_ogv)) deallocate(rv_a_ogv)
    if (allocated(iv_a_okv)) deallocate(iv_a_okv)
    deallocate(ivp)
  end subroutine rsibl_ev
end module m_rsibl_ev
