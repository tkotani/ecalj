module m_addrbl !- Adds to the smooth and local output density and to eigval sum
  real(8),allocatable,protected,public:: swtk(:,:,:) !spin weight
  public:: addrbl, swtkzero, m_addrbl_allocate_swtk
  private
contains
  subroutine m_addrbl_allocate_swtk(ndham,nsp,nkp)
    integer::ndham,nsp,nkp
    if(allocated(swtk)) deallocate(swtk)
    allocate(swtk(ndham,nsp,nkp))
  end subroutine m_addrbl_allocate_swtk
  subroutine swtkzero()
    swtk=0d0
  end subroutine swtkzero
  subroutine addrbl(isp,q,iq,smpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,evec,evl,nevl,&
       smrho,sumqv,sumev,sv_p_oqkkl,sv_p_oeqkkl,f)
    use m_struc_def
    use m_suham,only: &
         ndham=>ham_ndham,ndhamx=>ham_ndhamx,nspx=>ham_nspx
    use m_lmfinit,only:alat=>lat_alat,nbas, ispec,sspec=>v_sspec,nsp,nspc,lmet=>bz_lmet,&
         lekkl, zbak ,lfrce=>ctrl_lfrce
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol
    use m_supot,only: lat_nabc,k1,k2,k3
    use m_igv2x,only: napw,ndimh,ndimhx,igapw=>igv2x
    use m_subzi, only: lwtkb,nevmx, lswtk,wtkb=>rv_a_owtkb
    use m_mkqp,only: wtkp=>rv_a_owtkp
    use m_mkpot,only: qval_=>qval
    use m_ropyln,only: ropyln
    use m_rsibl,only:rsibl
    use m_rlocbl,only: rlocbl
    !- Adds to the smooth and local output density and to eigval sum
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   sspec :struct for species-specific information; see routine uspec
    !i   isp   :current spin channel
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   q     :Bloch vector
    !i   wtkp  :q-point weights from symmetry operations
    !i   ndham :leading dimension of evl
    !i   ndimh :dimension of hamiltonian
    !i   napw  :number of PWs in APW part of basis
    !i   igapw :PWs in units of reciprocal lattice vectors
    !i   lmet  :See Remarks in subzi
    !i         :0 assume insulator
    !i         :1 save eigenvectors to disk
    !i         :2 read weights from file, if they are needed
    !i         :3 always make two band passes; weights never needed a priori
    !i         :4 BZ integration with 3-point scheme
    !i   lwtkb :0 set of weights is not given; use 3-point interpolation
    !i         :1 or 2 given a set of weights
    !i         :-1 needs weights, but not yet available
    !i   wtkb  :integration weights, needed if lwtkb is 1 or 2
    !i   lswtk :<1 do nothing
    !i         :1 given a set of weights, make 'spin weights' swtk
    !i   iq    :index to current k-point
    !i   lfrce :if nonzero, accumulate contribution to force
    !i   ldos  :if nonzero, accumulate density-of-states
    !i   k1,k2,k3 dimensions of smpot,smrho
    !i   smpot :smooth potential on uniform mesh (mkpot.f), for forces
    !i   vconst:additional constant potential
    !i   qval  :total valence charge
    !i   evec  :eigenvectors
    !i   evl   :eigenvalues
    !i   nev
    ! xxx   ef0   :estimate for fermi level
    ! xx   def   :When lmet=4, charge also accmulated for ef0+def and ef0-def
    ! xx   esmear:(sampling integration) gaussian broadening
    ! xx         :sign and integer part are extensions; see mkewgt.f
    !i   emin  :energy lower bound when adding to sampling dos
    !i   emax  :energy upper bound when adding to sampling dos
    !i   ndos  :number of energy mesh points
    !i   osig,otau,oppi  augmentation matrices
    ! xxx   lcplxp=1 only now !lcplxp:0 if ppi is real; 1 if ppi is complex
    !i   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
    !o Outputs
    !o   sumqv :integrated charge, resolved by spin
    !o   sumev :sum of eigenvalues
    !o   dos   :sampling density of states, if ldos=1
    !o   smrho :smooth density on uniform mesh
    !o   oqkkl :local part of density matrix
    !o   oeqkkl:local part of energy-weighted density matrix
    !o   f     :eigenvalue contribution to forces
    !o   swtk  :'spin weights' to determine global magnetic moment, nspc=2
    !o         : swtk = diagonal part of  (z)^-1 sigmz z
    !r Remarks
    !u Updates
    !u   05 Jul 08 (T. Kotani) output density for new PW part
    !u             Option to accumulate energy-weighted output density
    !u   09 Jun 07 Makes spin weights (noncollinear case)
    !u   02 Jan 06 sumqv resolved by spin
    !u   17 Jan 05 Extension of esmear to Fermi distribution
    !u   23 Dec 04 Extended to spin-coupled case
    !u   18 Nov 04 Sampling integration properly handles Fermi distribtion
    !u    1 Sep 04 Adapted to handle complex ppi
    !u   23 Jan 01 Added lrout switch
    !u   17 Jun 00 spin polarized
    !u   22 May 00 Adapted from nfp add_densw.f
    ! ----------------------------------------------------------------------
    implicit none
    integer :: isp,iq!,lfrce
    integer :: nevl
    double precision :: emax,emin,qval,vconst
    real(8):: q(3), evl(ndham,nsp), sumev(2,3), sumqv(3,2), f(3,*)
    type(s_cv1) :: sv_p_oppi(3,1)
    type(s_rv1) :: sv_p_otau(3,1)
    type(s_rv1) :: sv_p_osig(3,1)
    type(s_rv1) :: sv_p_oeqkkl(3,1)
    type(s_rv1) :: sv_p_oqkkl(3,1)
    double complex evec(ndimh,nspc,ndimh,nspc),smrho(k1,k2,k3,isp),smpot(k1,k2,k3,isp)
    integer:: i , k , nevec , lmxax , lmxa , nlmax , nlmto , ig
    double precision :: vavg,wgt,tpiba
    integer :: ipiv(ndimh*2),i_copy_size
    complex(8),allocatable:: evecc(:,:,:,:),work(:,:,:,:)
    real(8),allocatable:: qpgv(:,:),qpg2v(:),ylv(:,:)
    double precision :: ewgt(ndimh*nspc),epsnevec
    qval= qval_- zbak
    if (lwtkb < 0) return
    call tcn('addrbl')
    nlmto = ndimh-napw
    lmxax = -1
    do  i = 1, nbas
       k=ispec(i) 
       lmxa=sspec(k)%lmxa
       lmxax = max(lmxax,lmxa)
    enddo
    nlmax=(lmxax+1)**2
    if (napw>0) then
       allocate(ylv(napw,nlmax),qpgv(3,napw),qpg2v(napw))
       tpiba = 2d0*4d0*datan(1d0)/alat
       do  ig = 1, napw
          qpgv(:,ig) = tpiba * ( q + matmul(qlat,igapw(:,ig)) )
       enddo
       call ropyln(napw,qpgv(1,1:napw),qpgv(2,1:napw),qpgv(3,1:napw), lmxax,napw,ylv,qpg2v)
    else
       allocate(ylv(1,1),qpgv(1,1),qpg2v(1))
    endif
    ! --- Decide how many states to include and make their weights ---
    ! ... Case band weights not passed: make sampling weights
    if (lwtkb == 0) then
       call rxx(nspc.ne.1,'lwtkb=0 not implemented in noncoll case')
       wgt = abs(wtkp(iq))/nsp
       call mkewgt(lmet,wgt,qval,ndimh,evl(1,isp),nevec,ewgt,sumev,sumqv(1,isp))
       call dscal(nevec,wgt,ewgt,1)
    else! ... Case band weights are passed
       call dcopy(nevl,wtkb(1,isp,iq),1,ewgt,1)
       do  10  i = nevl, 1, -1
          nevec = i
          if (abs(wtkb(i,isp,iq)) > epsnevec()) goto 12
10     enddo
12     continue
    endif
    ! ... Force from smooth analytic hamiltonian and overlap
    if (lfrce > 0 ) then
       call rxx(nspc.ne.1,'forces not implemented in noncoll case')
       vavg = vconst
       if (nlmto>0) then
          call fsmbl(vavg , q , ndimh , nlmto, nevec, evl(1,isp) , evec ,  ewgt , f )
       endif
       if (napw>0) then
          call fsmbpw (vavg, ndimh , nlmto, nevec , evl(1,isp), evec , ewgt , napw, qpgv, &
               qpg2v , ylv , nlmax , lmxax , alat , dsqrt ( vol ) , f )
       endif
    endif
    ! ... Add to smooth density
    call rsibl( lfrce, isp , q , &
         iq , ndimh , nspc , napw , igapw , nevec &
         , evec , ewgt , k1 , k2 , k3 , smpot , smrho , f )
    ! ... Add to local density coefficients
    call rlocbl (lfrce , nbas , isp , q , &
         ndham , ndimh , nspc , napw , igapw , nevec &
         , evec , ewgt , evl , sv_p_osig , sv_p_otau , sv_p_oppi &
         , lekkl , sv_p_oqkkl , sv_p_oeqkkl , f )
    ! ... Weights for spin moments
    if (lswtk>0 .AND. nspc==2) then
       allocate(evecc(ndimh,2,ndimh,2),work(ndimh,2,ndimh,2))
       call zcopy(ndimhx**2,evec,1,evecc,1)
       call zgetrf(nevl,nevl,evecc,ndimhx,ipiv,i)
       if (i /= 0) call rx('addrbl: failed to generate overlap')
       call zgetri(nevl,evecc,ndimhx,ipiv,work,ndimhx**2,i) !evecc is evec^-1
       do  i = 1, ndimh
          do  k = 1, ndimh
             swtk(i,1,iq)= swtk(i,1,iq) + evecc(i,1,k,1)*evec(k,1,i,1) - evecc(i,1,k,2)*evec(k,2,i,1)
             swtk(i,2,iq)= swtk(i,2,iq) + evecc(i,2,k,1)*evec(k,1,i,2) - evecc(i,2,k,2)*evec(k,2,i,2)
          enddo
       enddo
       deallocate(evecc,work)
    endif
    deallocate(qpgv,qpg2v,ylv)
    call tcx('addrbl')
  end subroutine addrbl

  subroutine addsds(ndimh,evl,wgt,emin,emax,ndos,dos) !,esmear
    use m_lmfinit,only: bz_w,bz_n
    !- Add to sampling dos
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   ndimh :hamiltonian dimension
    !i   evl   :eigenvalues
    !i   wgt   :eigenvalue weights
    !i   emin  :lower bound for DOS
    !i   emax  :upper bound for DOS
    !i   esmear:Parameter that describes gaussian broadening.
    !i         :Integer part >0 for for generalized gaussian broadening
    !i         :and is the the Methfessel-Paxton integration order
    !i         :Fractional part is the broadening width.
    !i         :Integer part <0 => Fermi-Dirac broadening used
    !i         :Fractional part is the temperature
    !i         :(see delstp.f)
    !i         :integer part above 100's digit is stripped.
    !i   ndos  :dimensions dos
    !o Outputs
    !o   dos   :DOS accumulated for these eigenvalues
    !l Local variables
    !l         :
    !r Remarks
    !r
    !u Updates
    !u   17 Jan 05 Extension of esmear to Fermi distribution
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: ndimh,ndos
    double precision :: evl(1),dos(ndos,2),wgt,emin,emax,esmear
    ! ... Local parameters
    integer :: nord,ie,i1,i2,i
    double precision :: width,de,eigval,ei,sn,dn,fn,x
    !      width = dabs(esmear) - int(dabs(esmear))
    !     nord = dsign(1d0,esmear) * mod(int(dabs(esmear)),100)
    width= abs(bz_w)
    nord = bz_n
    de = (emax-emin)/(ndos-1)
    do  ie = 1, ndimh
       eigval = evl(ie)
       i1 = (eigval-emin-width*5d0)/de + 1
       i2 = (eigval-emin+width*5d0)/de + 2
       i1 = max0(i1,1)
       i1 = min0(i1,ndos)
       i2 = max0(i2,1)
       i2 = min0(i2,ndos)
       do  i = i1, i2
          ei = emin + (i-1)*de
          x = (eigval-ei)/width
          call delstp(nord,x,dn,fn,sn)
          dos(i,2) = dos(i,2) + wgt*fn
          dos(i,1) = dos(i,1) + (wgt/width)*dn
       enddo
       do i = i2+1,ndos
          dos(i,2) = dos(i,2) + wgt
       enddo
    enddo
  end subroutine addsds
  subroutine mkewgt(lmet,wgt,qval,ndimh,evl, nevec,ewgt,sumev,sumqv)
    use m_lmfinit,only: bz_w,bz_n,stdo
    !- State-dependent weights for sampling BZ integration.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   lmet  :0 assume insulator
    !i         :1 save eigenvectors to disk
    !i         :2 read weights from file, if they are needed
    !i         :3 always make two band passes; weights never needed a priori
    !i         :4 BZ integration with 3-point scheme
    !i   wgt   :symmetry weight for this qp
    !i   qval  :valence charge; number of occupied states
    !i   ndimh :dimensions evl
    !i   evl   :eigenvalues
    ! xxx   ef0   :trial Fermi level
    ! xxx   def   :uncertainty in Fermi level; quantities are accumulated for
    ! xxx         :ef0-def, ef0, ef0+def
    !i   esmear:Parameter that describes gaussian broadening.
    !i         :Integer part >0 for for generalized gaussian broadening
    !i         :and is the the Methfessel-Paxton integration order
    !i         :Fractional part is the broadening width.
    !i         :Integer part <0 => Fermi-Dirac broadening used
    !i         :Fractional part is the temperature
    !i         :(see delstp.f)
    !i         :integer part above 100's digit is stripped.
    !o Outputs
    !o   nevec :number of states with nonzero weights
    !o   ewgt  :weights
    !o   sumev :sumev(*) are sumev(1..2) at Fermi levels
    !o         :sumev(1) = sum of eigenvalues
    !o         :sumev(2) = entropy
    !o   sumqv :sum of charges
    !r Remarks
    !r   Weights are made for generalized Gaussian broadening, where a
    !r   delta function is expressed as a hermite polynomial*gaussian,
    !r   or for Fermi-Dirac broadening (see delstp.f)
    !u Updates
    !u   17 Jan 05 Extension of esmear to Fermi distribution
    !u   31 May 00 Adapted from nfp mkewgt
    ! ----------------------------------------------------------------------
    implicit none
    integer :: lmet,ndimh,nevec!,numq
    double precision :: qval,wgt !esmear,
    double precision :: evl(ndimh),ewgt(1),sumev(2),sumqv
    integer :: i,i1,i2,ie,iprint,iq,k,nord
    double precision :: dn,eff,fevec,fn,sn,width,wtop,wx,x
    if (qval/2 > ndimh) call fexit2(-1,111,'%N Exit -1 : MKEWGT: ' &
         //'Basis with ndimh=%i is insufficient to carry q/2=%d states', &
         ndimh,qval/2)
    ! --- Nonmetal case ---
    if (lmet == 0) then
       !        if (numq .ne. 1) call rxi('mkewgt: nonmetal but numq=',numq)
       fevec = qval/2d0
       nevec = fevec + 0.500001d0
       wtop = fevec-(nevec-1)
       if (dabs(wtop-1d0) > 1d-6) write(stdo,"(' mkewgt (warning):'// &
            ' uneven occupation for nonmet, top wgt=',d15.7)")wtop
       do  i = 1, nevec
          ewgt(i) = 1d0
          if (i == nevec) ewgt(i) = wtop
          sumev(1) = sumev(1) + wgt*ewgt(i)*evl(i)
          sumev(2) = 0d0
          sumqv    = sumqv   + wgt*ewgt(i)
       enddo
    endif
  end subroutine mkewgt
  subroutine fsmbl(vavg,q,ndimh,nlmto,nevec,evl,evec,ewgt, f) 
    !- Force from smoothed hamiltonian (constant potential) and overlap
    use m_lmfinit,only: lhh,nkaphh,ispec,nbas,sspec=>v_sspec
    use m_uspecb,only:uspecb
    use m_struc_def
    use m_orbl,only: Orblib1,Orblib2,ktab1,ltab1,offl1,norb1,ktab2,ltab2,offl2,norb2
    use m_smhankel,only: hhigbl
    use m_lattic,only: rv_a_opos
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   sspec :struct for species-specific information; see routine uspec
    !i     Elts read: *
    !i     Stored:    *
    !i     Passed to: uspecb
    !i   slat  :struct for lattice information; see routine ulat
    !i     Stored:    *
    !i     Passed to: hhigbl
    !i   vavg  :constant potential (MT zero) to be added to h
    !i   q     :Bloch wave vector
    !i   ndimh :dimension of hamiltonian
    !i   nevec :number of occupied eigenvectors
    !i   evl   :eigenvalues
    !i   evec  :eigenvectors
    !i   ewgt  :eigenvector weights
    !o Outputs
    !o   f
    !r Remarks
    !u Updates
    !u   05 Jul 08 Decouple ndimh from nlmto, for PW basis
    !u   10 Apr 02 Redimensioned eh,rsmh to accommodate larger lmax
    !u   15 Feb 02 (ATP) Added MPI parallelization
    !u   10 Sep 01 Extended to local orbitals.
    !u   23 May 00 Adapted from nfp fsm_q.f
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ndimh,nlmto,nevec
    real(8):: q(3)
    double precision :: evl(ndimh),f(3,nbas),ewgt(nevec),vavg
    double complex evec(ndimh,ndimh)
    integer :: nlms,k0,n0,nkap0
    parameter (nlms=25, k0=1, n0=10, nkap0=3)
    integer:: i1,i2,ib1,ib2,ilm1,ilm2,io1,io2,iq,is1,is2,l1,l2,ik1,ik2,ivec,m,nlm1,nlm2
    integer :: lh1(nkap0),lh2(nkap0),nkap1,nkap2,nlm21,nlm22,nlm11,nlm12
    integer:: blks1(n0*nkap0),ntab1(n0*nkap0)
    integer:: blks2(n0*nkap0),ntab2(n0*nkap0)
    real(8) :: e1(n0,nkap0),rsm1(n0,nkap0),p1(3),wt,e2(n0,nkap0),rsm2(n0,nkap0),p2(3),xx(n0)
    complex(8):: s(nlms,nlms,0:k0,nkap0,nkap0), ds(nlms,nlms,0:k0,3,nkap0,nkap0),ccc(3),ssum(3)
    integer:: iloopend,ibl1,ibl2,ol1,oi1,ol2,oi2
    if (nevec <= 0) return
    call tcn ('fsmbl')
    do ib1=1,nbas
       is1=ispec(ib1) !ssite(ib1)%spec
       p1 =rv_a_opos(:,ib1) !ssite(ib1)%pos
       call uspecb(is1,rsm1,e1)
       call orblib1(ib1) !norb1,ltab1,ktab1,offl1
       call gtbsl1(0,norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)
       do ib2=ib1+1,nbas
          is2=ispec(ib2) !ssite(ib2)%spec
          p2 =rv_a_opos(:,ib2) !ssite(ib2)%pos
          call uspecb(is2,rsm2,e2)
          call orblib2(ib2)
          call gtbsl1(0,norb2,ltab2,ktab2,rsm2,e2,ntab2,blks2)
          !     ... M.E. <1> and <KE> between all envelopes connecting ib1 and ib2
          do i1 = 1, nkaphh(is1) !nkap1
             do i2 = 1, nkaphh(is2) !nkap2
                nlm1 = (lhh(i1,is1)+1)**2
                nlm2 = (lhh(i2,is2)+1)**2
                if (nlm1 > nlms .OR. nlm2 > nlms) call rx('fsmbl: increase nlms')
                call hhigbl(11,p1,p2,q,rsm1(1,i1),rsm2(1,i2),e1(1,i1),e2(1,i2),nlm1,nlm2,1,nlms,nlms,k0,&
                     s(1,1,0,i1,i2), ds(1,1,0,1,i1,i2))
             enddo
          enddo
          do io2= 1, norb2
             ik2= ktab2(io2)
             ol2= ltab2(io2)**2
             oi2= offl2(io2)
             do io1= 1, norb1
                ik1= ktab1(io1)
                ol1= ltab1(io1)**2
                oi1= offl1(io1)
                do  ivec = 1, nevec
                   ssum = 0d0
                   do ibl2=1,blks2(io2) 
                      do ibl1=1,blks1(io1) 
                         ccc = [(vavg*ds(ol1+ibl1,ol2+ibl2,0,m,ik1,ik2) &
                              -       ds(ol1+ibl1,ol2+ibl2,1,m,ik1,ik2) &
                              - evl(ivec)*ds(ol1+ibl1,ol2+ibl2,0,m,ik1,ik2),m=1,3)]
                         ssum = ssum + dconjg(evec(oi1+ibl1,ivec))*ccc*evec(oi2+ibl2,ivec)
                      enddo
                   enddo
                   wt = ewgt(ivec)
                   f(:,ib1) = f(:,ib1) - 2*wt*ssum
                   f(:,ib2) = f(:,ib2) + 2*wt*ssum
                enddo
             enddo
          enddo
       enddo
    enddo
    call tcx ('fsmbl')
  end subroutine fsmbl
  subroutine fsmbpw(vavg,ndimh,nlmto,nevec,evl,evec,ewgt,napw,qpgv,qpg2v,ylv,nlmax,lmxax,alat,sqv, f)
    use m_struc_def
    use m_uspecb,only:uspecb
    use m_orbl,only: Orblib1,ktab1,ltab1,offl1,norb1
    use m_lattic,only: rv_a_opos
    use m_lmfinit,only: nbas,sspec=>v_sspec,ispec
    !- Force from smoothed hamiltonian (constant potential), PW contribution
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   ssite :struct for site-specific information; see routine usite
    !i     Elts read: spec pos
    !i     Stored:    *
    !i     Passed to: *
    !i   sspec :struct for species-specific information; see routine uspec
    !i     Elts read: *
    !i     Stored:    *
    !i     Passed to: uspecb
    !i   vavg  :constant potential (MT zero) to be added to h
    !i   ndimh :dimension of hamiltonian
    !i   nlmto :dimension of lmto part of hamiltonian
    !i   nevec :number of occupied eigenvectors
    !i   evl   :eigenvalues
    !i   evec  :eigenvectors
    !i   ewgt  :eigenvector weights
    !i   napw  :number of augmented PWs in basis
    !i   qpgv
    !i   qpg2v
    !i   ylv
    !i   nlmax
    !i   lmxax
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   sqv   :square root of volume
    !o Outputs
    !o   f     :PW contribution to force is added to f
    !r Remarks
    !u Updates
    !u   04 Jul 08 (T. Kotani) first created
    ! ----------------------------------------------------------------------
    implicit none
    integer :: ndimh,napw,nlmax,nlmto,nevec,lmxax
    integer,parameter:: n0=10, nkap0=3
    real(8),parameter:: pi = 4d0*datan(1d0), fpi = 4*pi
    real(8):: evl(ndimh),f(3,nbas),ewgt(nevec),vavg,qpgv(3,napw),qpg2v(napw),qpg2,alat,sqv
    real(8):: gam,denom, e1(n0,nkap0),rsm1(n0,nkap0),p1(3),xx(n0),wt,ylv(napw,nlmax),ssum(3)
    integer :: i1,i2,ib1,ilm1,io1,iq,is1,l1,ik1,ig,ivec,nglob, nlm11,nlm12,m
    integer:: blks1(n0*nkap0),ntab1(n0*nkap0),lh1(nkap0),nkap1
    complex(8):: phase,fach,ovl,ccc(3),sum, srm1l(0:n0),evec(ndimh,ndimh),img=(0d0,1d0)
    integer:: ibl1,oi1,ol1
    if (nevec <= 0) return
    call tcn ('fsmbpw')
    srm1l(0)=1d0
    do  l1 = 1, lmxax
       srm1l(l1) = img**l1
    enddo
    ib1loop: do 1000 ib1=1,nbas
       is1=ispec(ib1) !ssite(ib1)%spec
       p1=rv_a_opos(:,ib1) !ssite(ib1)%pos
       call uspecb(is1,rsm1,e1)
       call orblib1(ib1) !norb1,ltab1,ktab1,xx,offl1,xx)
       call gtbsl8(norb1,ltab1,ktab1,rsm1,e1,ntab1,blks1)
       !   ... Hsm (i1) \times i(q+G)[(q+G)**2+const] PW (i2) Takao. Taken from smhsbl.f
       !       i1--> Hsm, i2--> PW
       igloop: do 2000 ig = 1, napw
          i2 = ig + nlmto
          qpg2 = qpg2v(ig)
          phase = exp(img*alat*sum(qpgv(:,ig)*p1))
          iorbloop: do 3000 io1 = 1, norb1
             l1  = ltab1(io1)
             ik1 = ktab1(io1)
             ol1 = ltab1(io1)**2
             oi1 = offl1(io1)
             denom = e1(l1+1,ik1) - qpg2
             gam   = 1d0/4d0*rsm1(l1+1,ik1)**2
             fach  = -fpi/denom * phase * srm1l(l1) * exp(gam*denom)
             iorbblock: do 3010 ibl1 = 1,blks1(io1) 
                !             s(i1,i2) = ovl
                !             h(i1,i2) = (qpg2 + vavg) * ovl
                ovl = fach * ylv(ig,ol1+ibl1)/sqv ! Eq. 9.4 in JMP39 3393
                do ivec = 1, nevec  !        gradient PW * (H - E S)
                   ccc = [(ovl * img*qpgv(m,ig) * (qpg2 + vavg - evl(ivec)),m=1,3)]
                   ssum = ewgt(ivec)*[(dconjg(evec(oi1+ibl1,ivec))*ccc(m)*evec(i2,ivec),m=1,3)]
                   f(:,ib1) = f(:,ib1) - 2d0*ssum
                enddo
3010         enddo iorbblock
3000      enddo iorbloop
2000   enddo igloop
1000 enddo ib1loop
    call tcx ('fsmbpw')
  end subroutine fsmbpw
end module m_addrbl
