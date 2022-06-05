module m_addrbl
  real(8),allocatable,protected:: swtk(:,:,:)
  public:: swtk, swtkzero, addrbl,m_addrbl_allocate_swtk
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

  subroutine addrbl ( isp,q, &
       iq,lfrce, smpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi, &
       evec,evl,nevl,  smrho, sumqv,sumev, sv_p_oqkkl,sv_p_oeqkkl, f)
    use m_struc_def
    use m_suham,only: &
         ndham=>ham_ndham,ndhamx=>ham_ndhamx,nspx=>ham_nspx
    use m_lmfinit,only:alat=>lat_alat,nbas, ssite=>v_ssite,sspec=>v_sspec,nsp,nspc,lmet=>bz_lmet,&
         lekkl, zbak !,ldos=>ctrl_ldos!,bz_dosw
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol
    use m_supot,only: lat_nabc,k1,k2,k3
    use m_igv2x,only: napw,ndimh,ndimhx,igapw=>igv2x
    use m_lmfinit,only: iprmb
    use m_subzi, only: lwtkb,nevmx, lswtk,wtkb=>rv_a_owtkb
    use m_mkqp,only: wtkp=>rv_a_owtkp
    use m_mkpot,only: qval_=>qval
    use m_ropyln,only: ropyln
    use m_rsibl,only:rsibl
    !- Adds to the smooth and local output density and to eigval sum
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   ssite :struct for site-specific information; see routine usite
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
    integer :: isp,iq,lfrce
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
       k=ssite(i)%spec
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
       ! ... Case band weights are passed
    else
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
       !        print *, '!! avg=1=smpot=1'; vavg=1; smpot=1 !; nevec=1
       !        print *, '!! smpot=vavg'; smpot=vavg
       !       call zprm('evecs',2,evec,ndimh,ndimh,nevec)
       if (nlmto>0) then
          call fsmbl ( nbas , ssite , sspec ,  vavg , q , ndimh ,&
          nlmto , iprmb , nevec , evl ( 1 , isp ) , evec , &
               ewgt , f )
       endif
       if (napw>0) then
          call fsmbpw ( nbas , ssite , sspec , vavg , ndimh , nlmto , iprmb &
               , nevec , evl ( 1 , isp ) , evec , ewgt , napw , qpgv &
               , qpg2v , ylv , nlmax , lmxax , alat , dsqrt ( vol ) , f )
       endif
    endif
    ! ... Add to smooth density
    call rsibl ( ssite , sspec , lfrce , nbas , isp , q , &
         iq , ndimh , nspc , napw , igapw , iprmb , nevec &
         , evec , ewgt , k1 , k2 , k3 , smpot , smrho , f )
    ! ... Add to local density coefficients
    call rlocbl ( ssite , sspec , lfrce , nbas , isp , q , &
         ndham , ndimh , nspc , napw , igapw , iprmb ,  nevec &
         , evec , ewgt , evl , sv_p_osig , sv_p_otau , sv_p_oppi &
    , lekkl , sv_p_oqkkl , sv_p_oeqkkl , f )
    ! ... Weights for spin moments
    if (lswtk>0 .AND. nspc==2) then
       allocate(evecc(ndimh,2,ndimh,2),work(ndimh,2,ndimh,2))
       call zcopy(ndimhx**2,evec,1,evecc,1)
       call zgetrf(nevl,nevl,evecc,ndimhx,ipiv,i)
       if (i /= 0) call rx('addrbl: failed to generate overlap')
       call zgetri(nevl,evecc,ndimhx,ipiv,work,ndimhx**2,i)
       do  i = 1, ndimh
          do  k = 1, ndimh
             swtk(i,1,iq) = swtk(i,1,iq) + evecc(i,1,k,1)*evec(k,1,i,1) &
                  - evecc(i,1,k,2)*evec(k,2,i,1)
             swtk(i,2,iq) = swtk(i,2,iq) + evecc(i,2,k,1)*evec(k,1,i,2) &
                  - evecc(i,2,k,2)*evec(k,2,i,2)
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


end module m_addrbl
