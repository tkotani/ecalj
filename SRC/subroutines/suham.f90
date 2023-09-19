module m_suham
  use m_MPItk,only:master_mpi
  use m_ftox
  integer,protected,public:: ham_ndham, ham_ndhamx,ham_nspx  !only these integers returned.
  public m_suham_init
  private
contains
  subroutine m_suham_init() !Get Hamiltonian dimension ham_ndham= Hamiltonian dimensition= nAPW +nlmto
    use m_struc_def
    use m_lmfinit,only:lso,nspc,nsp,nlmto, pwmode=>ham_pwmode,pwemax &!,pwemin &
         ,alat=>lat_alat, lat_tolft,nkaph, pot_nlma, pot_nlml ,stdo,nlmto
    use m_supot,only: lat_ng, rv_a_ogv
    use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
    !Main output is ndham = estimate for upper dimension of hamiltonian, including possible APW part
    !l   ndim  :total number of lmto orbitals = nl**2 * nbas
    !l   npwmin: lower limit to number of PWs, estimated by sweeping
    !l         :over a fine mesh of qp.
    !l   npwmax: upper limit to number of PWs, estimated by sweeping
    !l         : over a fine mesh of qp.
    !l   nqdiv: loop over mesh of q-points to estimate npwmax
    !l        : nqdiv is fineness of q-mesh.
    !l  npwpad: a 'safety' padding to npwmax in case npwmax
    !l        : underestimates actual upper limit !required???
    implicit none
    integer :: hord,i,i1,i2,iprint, lidim,lihdim,nclasp,ndim,neul, &
         nlspcp,nspx,nttab,partok,igets,nvi,nvl,ib,is,lmxa,lmxl,isw,nbf,lmxax, j1,j2,j3,m,npw,npwmin,npwmax
    integer:: ndham=0
    double precision :: q(3),Gmin,Gmax,xx
    integer,parameter:: nqdiv=12
    integer ::npwpad
    double precision :: ckbas,cksumf,kap2(10)
    integer:: oalph , oiax , ontab , os , ng
    real(8),allocatable :: rv_a_ogvx(:)
    integer,allocatable :: kv_iv(:), igv2_iv(:)
    integer:: obas , oo , opp ,  opti , obs , omagf,i_copy_size
    call tcn('m_suham_init')
    ! --- Hamiltonian offsets, orbital permutation table ---
    if (mod(pwmode,10)==2 .AND. master_mpi) write(stdo,"(a)")' suham: no MTO basis'
    !   ... PW setup : estimate upper bound to number of G vectors to set up upper bound to hamiltonian dimension
    ham_ndham = nlmto
    ndham     = nlmto
    if (pwemax>0 .AND. mod(pwmode,10)>0) then !       Gmin = dsqrt(pwemin)
       Gmax = dsqrt(pwemax)
       if (mod(pwmode/10,10) == 1) then
          if(master_mpi)write(stdo,*)'Estimate max size of PW basis from combinations of recip. lattice vectors ...'
          npwmax = -1
          npwmin = 99999
          do  j1 = 0, nqdiv
             do  j2 = 0, nqdiv
                do  j3 = 0, nqdiv
                   do   m = 1, 3
                      q(m) = (qlat(m,1)/nqdiv)*j1 +(qlat(m,2)/nqdiv)*j2 +(qlat(m,3)/nqdiv)*j3
                   enddo
                   call pshpr(iprint()-40)
                   call getgv2(alat,plat,qlat,q,Gmax,1,npw,xx) 
                   call poppr
                   npwmin = min(npwmin,npw)
                   npwmax = max(npwmax,npw)
                enddo
             enddo
          enddo
          npwpad = max(nint((npwmax-npwmin)*0.2d0),3) !we need this for safe?
       else
          q=0d0
          call pshpr(0)
          call getgv2(alat,plat,qlat,q,Gmax,1, npw,xx) 
          call poppr
          npwmin = npw
          npwmax = npw
          npwpad = 0
       endif
       ndham = npwmax + npwpad
       if (mod(pwmode,10) /= 2) ndham = nlmto + npwmax + npwpad
       ham_ndham = ndham
       !! Dimension maximum of Hamiltonian is ndhamx (ispx=1,nspx)
       !! Note spinoffdiag=T case: we use nspc,nsp,nspx,ndhamx (a little complicated, I think).
       !     nspx*nspc=nsp
       if(master_mpi)write(stdo,ftox)'suham: PW basis Emax=',ftof(pwemax,3),'npw,ndham=',npw,ndham
    endif
    if(nspc==2 .AND. nsp==1) call rx('suham: nspc==2 but nsp=1')
    ham_ndhamx= ndham*nspc !nspc=2 for spin-coupled case
    ham_nspx  = nsp/nspc
    call tcx('m_suham_init')
  end subroutine m_suham_init
end module m_suham

