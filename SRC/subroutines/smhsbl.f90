subroutine smhsbl(ssite,sspec,vavg,q,ndimh,iprmb, napw,igapw, h,s)
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg
  use m_lmfinit,only: lat_alat,nbas,nkaphh,lhh
  use m_lattic,only: lat_vol
  use m_uspecb,only:uspecb
  use m_struc_def
  use m_lattic,only:lat_plat
  use m_orbl,only: Orblib1,Orblib2,ktab1,ltab1,offl1,norb1,ktab2,ltab2,offl2,norb2
  use m_ropyln,only: ropyln
  use m_smhankel,only:hhibl
  !- Smoothed Bloch Hamiltonian (constant potential) and overlap matrix
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 compute both hamiltonian and overlap
  !i         :  otherwise, compute overlap only.
  !i         :  In this case, vavg is not used
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: *
  !i     Stored:    *
  !i     Passed to: uspecb
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: ocg ojcg oidxcg ocy
  !i     Stored:    *
  !i     Passed to: hhibl
  !i   vavg  :constant potential (MT zero) to be added to h
  !i   q     :Bloch wave vector
  !i   ndimh :dimension of hamiltonian
  !i   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
  !o Outputs
  !o   h     :smooth Bloch hamiltonian added to h (= s * vavg)
  !o   s     :smooth Bloch overlap added to s
  !b Bugs
  !b   Use of qpgv(1,1:napw) is poor design ... rewrite
  !r Remarks
  !r  *How orbital information is extracted and deployed.
  !r   Orbital specification requires the following information:
  !r     1. orbital type, specified by the triplet (l,rsmh,eh,pz)
  !r        (pz is relevant only for local orbitals)
  !r        Information is stored in spec->orbp, unpacked by uspecb
  !r        For local orbitals spec->pz is also required.
  !r        Note that each orbital type has 2*l+1 orbitals.
  !r     3. Location of orbitals (offsets) in hamiltonian
  !r        Information is stored in iprmb (passed array)
  !r
  !r   Routines needing this information unpack the data as follows:
  !r
  !r      call    Purpose, Data and format
  !r     (Input)
  !r     -------  ---------------------------------
  !r
  !r     uspecb    Extract orbital parameters for all orbital types.
  !r    (sspec)    rsmh(l,ik),eh(l,ik), l=0,lh(ik), ik=1..nkapi.
  !r               Entries for which rsmh(l,ik)>0 have envelopes.
  !r               Entries for which ik=nkapi and pz(l)>0 are local orbitals.
  !r               Note that these possiblities can simultaneously occur.
  !r
  !r     orbl      norb : (total number of orbital types)
  !r     (iprmb)   ltab(io),ktab(io) : specifes l and index ik to
  !r               rsmh(l,ik),eh(l,ik) for orbital type io (io=1..norb).
  !r               Thus, the parameters
  !r                 l=ltab(io),ik=ktab(io),rsmh(l,ik),eh(l,ik)
  !r               completely specify envelope information for type io
  !r               offl(io): hamiltonian offset for orbital type io
  !r
  !r     gtbsl1    Serves two purposes.
  !r     (norb,    1. Blocks orbitals with common (rsmh,eh) and
  !r     ltab,     consecutive l so that mesh tabulations are more
  !r     ktab,     efficiently generated.
  !r     rsmh,eh)  2. Marks each orbital type as being :
  !r                  a. valence function with envelope;
  !r                  b. local orbital with envelope;
  !r                  c. local orbital without envelope.
  !r               blks(io) = size of contiguous block of orbitals.  Note
  !r               that if io+1, io+2, ... are contiguous to io,
  !r               blks(io+1)=blks(io+2)=...=0
  !r
  !r     A typical loop that does NOT block orbitals therefore looks like
  !r
  !r     do  io= 1, norb
  !r       l  = ltab(io)
  !r       ik = ktab(io)
  !r   C   off = orbital index in iprmb order
  !r       off = offl(io)
  !r       do  ilm1 = l**2+1, (l+1)**2
  !r         off = off+1
  !r         ...
  !r       enddo
  !r     enddo
  !r
  !r     A typical loop that DOES block orbitals therefore looks like
  !r
  !r     do  io= 1, norb
  !r     if (blks(io) .ne. 0) then
  !r       l   = ltab(io)
  !r       ik  = ktab(io)
  !r       nlm1 = l**2+1
  !r       nlm2 = nlm1 + blks1(io1)-1
  !r   C   off = orbital index in iprmb order
  !r       off = offl(io)
  !r       do  ilm1 = nlm1, nlm2
  !r         off = off+1
  !r         ...
  !r       enddo
  !r     endif
  !r     enddo
  !r
  !r  *The Bloch sum of s(q) is computed as s = sum_T s_T exp(iq.T)
  !r   where T are the set of lattice vectors.
  !r
  !r  *NB: The matrix subblock in s for sites (i,j) is computed with
  !r       the connecting vector pi-pj.  This convention is different
  !r       from the tight-binding convention, where s is computed for
  !r       pj-pi.  Thus this (i,j) matrix subblock is the hermitian
  !r       conjugate of the same subblock in the tight-binding case.
  !r       tight-binding case.  Note, however, that the entire s
  !r       would not correspond to the hermitian conjugate of the tight
  !r       binding case.
  !r
  !r  *MPI
  !r   See remarks to hsibl. Buffers for h and s are F90 ALLOCATEd as
  !r   they need to be used locally as dimensioned. A buffer is taken
  !r   from the heap for ALLREDUCE.
  !u Updates
  !u   05 Jul 08 (T. Kotani) output density for new PW part
  !u   12 Aug 04 First implementation of extended local orbitals
  !u   14 Aug 02 Added overlap-only option
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   27 Aug 01 Extended to local orbitals.
  !u   02 Mar 01 Bug fix for multiple-kappa case
  !u   19 May 00 Adapted from nfp smhs_q.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: procid,master
  integer :: mode,ndimh,iprmb(1),napw,igapw(3,napw)
  real(8):: q(3) , vavg
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double complex h(ndimh,ndimh),s(ndimh,ndimh)
  ! ... Local parameters
  integer :: nlms,kdim,n0,nkap0
  parameter (nlms=25, kdim=1, n0=10, nkap0=3)
  integer :: nlmto
  integer :: i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,is1,is2,nlm1,nlm2,l1,l2,ig
  integer :: lh1(nkap0),lh2(nkap0),nkap1,nkap2
  !      integer ltab1(n0*nkap0),ktab1(n0*nkap0),offl1(n0*nkap0),norb1,ik1,
  integer:: ik1,blks1(n0*nkap0),ntab1(n0*nkap0)
  !      integer ltab2(n0*nkap0),ktab2(n0*nkap0),offl2(n0*nkap0),norb2,ik2,
  integer:: ik2,blks2(n0*nkap0),ntab2(n0*nkap0)
  double precision :: e1(n0,nkap0),rsm1(n0,nkap0),p1(3), &
       e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx
  double complex s0(nlms,nlms,0:kdim,nkap0,nkap0)
  ! ... for PW part
  integer :: lmxax,lmxa,nlmax
  double precision :: qpg2,alat,plat(3,3),qlat(3,3),vol,srvol,tpiba,pi, &
       denom,gam,fpi,ddot
  real(8),allocatable:: yl(:),ylv(:,:),qpgv(:,:),qpg2v(:)
  complex(8),allocatable:: srm1l(:)
  complex(8):: ovl,srm1,phase,fach
  parameter (srm1=(0d0,1d0))
  integer:: iloop,iloopmx
  call tcn('smhsbl')
  procid = 0
  master = 0
  ! --- Setup ---
  nlmto = ndimh-napw

  ! ... Setup below needed for PWs
  if (napw > 0) then
     alat=lat_alat
     plat=lat_plat
     vol=lat_vol
     pi = 4d0*datan(1d0)
     tpiba = 2d0*pi/alat
     srvol = dsqrt(vol)
     fpi = 4*pi
     call dinv33(plat,1,qlat,vol)
     vol = dabs(vol)*(alat**3)
     !     Find largest lmxa ... should be made elsewhere
     lmxax = -1
     do  ib1 = 1, nbas
        is1=ssite(ib1)%spec
        lmxa=sspec(is1)%lmxa
        lmxax = max(lmxax,lmxa)
     enddo
     nlmax=(lmxax+1)**2
     allocate(ylv(napw,nlmax),yl(nlmax),qpgv(3,napw),qpg2v(napw))
     !     Note: rearrange order
     do  ig = 1, napw
        qpgv(:,ig) = tpiba * (q + matmul(qlat, igapw(:,ig)))
     enddo
     call ropyln(napw,qpgv(1,1:napw),qpgv(2,1:napw),qpgv(3,1:napw), &
          lmxax,napw,ylv,qpg2v)
     allocate(srm1l(0:lmxax))
     srm1l(0) = 1d0
     do  l1 = 1, lmxax
        srm1l(l1) = (srm1)**l1
     enddo
  endif

  if (nlmto >0 ) then
     do ib1=1,nbas
        is1=ssite(ib1)%spec
        p1=ssite(ib1)%pos
        call uspecb(is1,rsm1,e1)
        !     Row info telling smhsbl where to poke s0 made by hhibl
        call orblib1(ib1)!,0,nlmto,iprmb,norb1,ltab1,ktab1,xx,offl1,xx)
        call gtbsl1(8+16,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)
        do  ib2 = ib1, nbas
           is2=ssite(ib2)%spec
           p2=ssite(ib2)%pos
           call uspecb(is2,rsm2,e2)
           !     Column info telling smhsbl where to poke s0 made by hhibl
           call orblib2(ib2) !,0,nlmto,iprmb,norb2,ltab2,ktab2,xx,offl2,xx)
           call gtbsl1(8+16,norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)
           !     ... M.E. <1> and <T> between all envelopes connecting ib1 and ib2
           do  i1 = 1, nkaphh(is1) !nkap1
              do  i2 = 1, nkaphh(is2) !nkap2
                 nlm1 = (lhh(i1,is1)+1)**2
                 nlm2 = (lhh(i2,is2)+1)**2
                 if (nlm1 > nlms .OR. nlm2 > nlms) &
                      call rx('smhsbl: increase nlms')
                 call hhibl ( 11 , p1 , p2 , q , rsm1 ( 1 , i1 ) , rsm2 ( 1 , &
                      i2 ) , e1 ( 1 , i1 ) , e2 ( 1 , i2 ) , nlm1 , nlm2 , 1 , nlms &
                      , nlms , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy &
                      , s0 ( 1 , 1 , 0 , i1 , i2 ) )
              enddo
           enddo
           !     ... Loop over orbital indices, poke block of integrals into s,h
           do  io2 = 1, norb2
              if (blks2(io2) /= 0) then
                 !     l2,ik2 = l and kaph indices, needed to locate block in s0
                 l2  = ltab2(io2)
                 ik2 = ktab2(io2)
                 !     i2 = orbital index in iprmb order
                 i2 = offl2(io2)
                 do  ilm2 = l2**2+1, (l2+1)**2
                    i2 = i2+1
                    !     if (mode .eq. 0) then
                    do  io1 = 1, norb1
                       if (blks1(io1) /= 0) then
                          !     l1,ik1 = l and kaph indices, needed to locate block in s0
                          l1  = ltab1(io1)
                          ik1 = ktab1(io1)
                          !     i1 = orbital index in iprmb order
                          i1 = offl1(io1)
                          do  ilm1 = l1**2+1, (l1+1)**2
                             i1 = i1+1
                             s(i1,i2) = s(i1,i2) + s0(ilm1,ilm2,0,ik1,ik2)
                             h(i1,i2) = h(i1,i2) - s0(ilm1,ilm2,1,ik1,ik2) &
                                  + vavg*s0(ilm1,ilm2,0,ik1,ik2)
                          enddo
                       endif
                    enddo
                 enddo
              endif
           enddo
           !     ... end loop over ib2
        enddo
        !     ... Hsm (i1) \times PW (i2)  Takao. Similar logic in fsmbl
        do  ig = 1, napw
           i2 = ig + nlmto
           qpg2 = qpg2v(ig)
           phase = exp(srm1*alat*ddot(3,qpgv(1,ig),1,p1,1))
           !     phase = exp( srm1 * sum(qpgv(:,ig)*p1)*alat  )
           do  io1 = 1, norb1
              if (blks1(io1) == 0) cycle
              l1  = ltab1(io1)
              ik1 = ktab1(io1)
              i1  = offl1(io1)
              denom = e1(l1+1,ik1) - qpg2
              gam   = 1d0/4d0*rsm1(l1+1,ik1)**2
              !     Note: fach depends on l.
              fach  = -fpi/denom * phase * srm1l(l1) * exp(gam*denom)
              do  ilm1 = l1**2+1, (l1+1)**2
                 i1 = i1+1
                 !     <Hsm^bloch | exp(i q+G r)/srvol > and
                 !     <Hsm^bloch | nabla + vavg |  exp(i q+G r)/srvol > in a cell.
                 ovl = fach * ylv(ig,ilm1)/srvol ! JMP Eq.(9.4)
                 !     gradient of ovl
                 !     ovl = srm1*qpgv(1,ig) * ovl
                 !     #if MPI
                 !     sbuf(i1,i2) = sbuf(i1,i2) + ovl
                 !     hbuf(i1,i2) = hbuf(i1,i2) + qpg2*ovl + vavg*ovl
                 s(i1,i2) = s(i1,i2) + ovl
                 h(i1,i2) = h(i1,i2) + qpg2*ovl + vavg*ovl
              enddo
           enddo
        enddo
        !     ... end loop over ib1
     enddo
  endif

  ! ... PW x PW part (diagonal matrix)
  do  ig = 1, napw
     i2 = ig + nlmto
     s(i2,i2) = s(i2,i2) + 1
     h(i2,i2) = h(i2,i2) + qpg2v(ig) + vavg
  enddo
  if (napw > 0) then
     deallocate(yl,ylv,qpgv,qpg2v,srm1l)
  endif
  ! ... Occupy second half of matrix
  do  i1 = 1, ndimh
     do  i2 = i1, ndimh
        h(i2,i1) = dconjg(h(i1,i2))
        s(i2,i1) = dconjg(s(i1,i2))
     enddo
  enddo
  call tcx('smhsbl')
end subroutine smhsbl


