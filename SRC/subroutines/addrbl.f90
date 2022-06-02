      module m_addrbl
      real(8),allocatable,protected:: swtk(:,:,:)
      public:: swtk, swtkzero, addrbl,m_addrbl_allocate_swtk
      private
      contains
      subroutine m_addrbl_allocate_swtk(ndham,nsp,nkp)
      integer::ndham,nsp,nkp
      if(allocated(swtk)) deallocate(swtk)
      allocate(swtk(ndham,nsp,nkp))
      end
      
      subroutine swtkzero()
      swtk=0d0
      end
      
      subroutine addrbl ( isp,q, 
     .     iq,lfrce, smpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,
     .     evec,evl,nevl,  !ef0,def, !emin,emax, !ndos,dos,
     o      smrho, sumqv,sumev, sv_p_oqkkl,sv_p_oeqkkl, f)
      use m_struc_def 
      use m_suham,only: 
     &     ndham=>ham_ndham,ndhamx=>ham_ndhamx,nspx=>ham_nspx
      use m_lmfinit,only:alat=>lat_alat,nbas, ssite=>v_ssite,sspec=>v_sspec,nsp,nspc,lmet=>bz_lmet,lekkl, !lcplxp,
     &     zbak !,ldos=>ctrl_ldos!,bz_dosw
      use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol
      use m_supot,only: lat_nabc,k1,k2,k3
      use m_igv2x,only: napw,ndimh,ndimhx,igapw=>igv2x
      use m_lmfinit,only: iprmb
      use m_subzi, only: lwtkb,nevmx, lswtk,wtkb=>rv_a_owtkb
      use m_mkqp,only: wtkp=>rv_a_owtkp
      use m_mkpot,only: qval_=>qval
C- Adds to the smooth and local output density and to eigval sum
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci   sspec :struct for species-specific information; see routine uspec
Ci   isp   :current spin channel
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   q     :Bloch vector
Ci   wtkp  :q-point weights from symmetry operations
Ci   ndham :leading dimension of evl
Ci   ndimh :dimension of hamiltonian
Ci   napw  :number of PWs in APW part of basis
Ci   igapw :PWs in units of reciprocal lattice vectors
Ci   lmet  :See Remarks in subzi
Ci         :0 assume insulator
Ci         :1 save eigenvectors to disk
Ci         :2 read weights from file, if they are needed
Ci         :3 always make two band passes; weights never needed a priori
Ci         :4 BZ integration with 3-point scheme
Ci   lwtkb :0 set of weights is not given; use 3-point interpolation
Ci         :1 or 2 given a set of weights
Ci         :-1 needs weights, but not yet available
Ci   wtkb  :integration weights, needed if lwtkb is 1 or 2
Ci   lswtk :<1 do nothing
Ci         :1 given a set of weights, make 'spin weights' swtk
Ci   iq    :index to current k-point
Ci   lfrce :if nonzero, accumulate contribution to force
Ci   ldos  :if nonzero, accumulate density-of-states
Ci   k1,k2,k3 dimensions of smpot,smrho
Ci   smpot :smooth potential on uniform mesh (mkpot.f), for forces
Ci   vconst:additional constant potential
Ci   qval  :total valence charge
Ci   evec  :eigenvectors
Ci   evl   :eigenvalues
Ci   nev
Cixxx   ef0   :estimate for fermi level
Cixx   def   :When lmet=4, charge also accmulated for ef0+def and ef0-def
Cixx   esmear:(sampling integration) gaussian broadening
Cixx         :sign and integer part are extensions; see mkewgt.f
Ci   emin  :energy lower bound when adding to sampling dos
Ci   emax  :energy upper bound when adding to sampling dos
Ci   ndos  :number of energy mesh points
Ci   osig,otau,oppi  augmentation matrices
cixxx   lcplxp=1 only now !lcplxp:0 if ppi is real; 1 if ppi is complex
Ci   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
Co Outputs
Co   sumqv :integrated charge, resolved by spin
Co   sumev :sum of eigenvalues
Co   dos   :sampling density of states, if ldos=1
Co   smrho :smooth density on uniform mesh
Co   oqkkl :local part of density matrix
Co   oeqkkl:local part of energy-weighted density matrix
Co   f     :eigenvalue contribution to forces
Co   swtk  :'spin weights' to determine global magnetic moment, nspc=2
Co         : swtk = diagonal part of  (z)^-1 sigmz z
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu             Option to accumulate energy-weighted output density
Cu   09 Jun 07 Makes spin weights (noncollinear case)
Cu   02 Jan 06 sumqv resolved by spin
Cu   17 Jan 05 Extension of esmear to Fermi distribution
Cu   23 Dec 04 Extended to spin-coupled case
Cu   18 Nov 04 Sampling integration properly handles Fermi distribtion
Cu    1 Sep 04 Adapted to handle complex ppi
Cu   23 Jan 01 Added lrout switch
Cu   17 Jun 00 spin polarized
Cu   22 May 00 Adapted from nfp add_densw.f
C ----------------------------------------------------------------------
      implicit none
      integer isp,iq,lfrce 
      integer nevl
      double precision emax,emin,qval,vconst 
      real(8):: q(3), evl(ndham,nsp), sumev(2,3), sumqv(3,2), f(3,*)
      type(s_cv1) :: sv_p_oppi(3,1)
      type(s_rv1) :: sv_p_otau(3,1)
      type(s_rv1) :: sv_p_osig(3,1)
      type(s_rv1) :: sv_p_oeqkkl(3,1)
      type(s_rv1) :: sv_p_oqkkl(3,1)
      double complex evec(ndimh,nspc,ndimh,nspc),smrho(k1,k2,k3,isp),smpot(k1,k2,k3,isp)
      integer:: i , k , nevec , lmxax , lmxa , nlmax , nlmto , ig 
      double precision vavg,wgt,tpiba
      integer ipiv(ndimh*2),i_copy_size
      complex(8),allocatable:: evecc(:,:,:,:),work(:,:,:,:)
      real(8),allocatable:: qpgv(:,:),qpg2v(:),ylv(:,:)
      double precision ewgt(ndimh*nspc),epsnevec

      qval= qval_- zbak 
      if (lwtkb .lt. 0) return
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
C --- Decide how many states to include and make their weights ---
C ... Case band weights not passed: make sampling weights
      if (lwtkb .eq. 0) then
        call rxx(nspc.ne.1,'lwtkb=0 not implemented in noncoll case')
        wgt = abs(wtkp(iq))/nsp
        call mkewgt(lmet,wgt,qval,ndimh,evl(1,isp),!ef0,def, !!esmear
     .  nevec,ewgt,sumev,sumqv(1,isp))
        call dscal(nevec,wgt,ewgt,1)
C ... Case band weights are passed
      else
        call dcopy(nevl,wtkb(1,isp,iq),1,ewgt,1)
        do  10  i = nevl, 1, -1
          nevec = i
          if (abs(wtkb(i,isp,iq)) .gt. epsnevec()) goto 12
   10   continue
   12   continue
      endif
C ... Force from smooth analytic hamiltonian and overlap
      if (lfrce .gt. 0 ) then
        call rxx(nspc.ne.1,'forces not implemented in noncoll case')
        vavg = vconst
c        print *, '!! avg=1=smpot=1'; vavg=1; smpot=1 !; nevec=1
C        print *, '!! smpot=vavg'; smpot=vavg
C       call zprm('evecs',2,evec,ndimh,ndimh,nevec)
        if (nlmto>0) then
          call fsmbl ( nbas , ssite , sspec ,  vavg , q , ndimh , !slat ,
     .     nlmto , iprmb , nevec , evl ( 1 , isp ) , evec , 
     .     ewgt , f )
        endif
        if (napw>0) then
          call fsmbpw ( nbas , ssite , sspec , vavg , ndimh , nlmto , iprmb 
     .     , nevec , evl ( 1 , isp ) , evec , ewgt , napw , qpgv 
     .     , qpg2v , ylv , nlmax , lmxax , alat , dsqrt ( vol ) , f )
        endif
      endif
C ... Add to smooth density
      call rsibl ( ssite , sspec , lfrce , nbas , isp , q , 
     .   iq , ndimh , nspc , napw , igapw , iprmb , nevec 
     .   , evec , ewgt , k1 , k2 , k3 , smpot , smrho , f )
C ... Add to local density coefficients
      call rlocbl ( ssite , sspec , lfrce , nbas , isp , q ,
     .        ndham , ndimh , nspc , napw , igapw , iprmb ,  nevec 
     .        , evec , ewgt , evl , sv_p_osig , sv_p_otau , sv_p_oppi  !,lcplxp 
     .        , lekkl , sv_p_oqkkl , sv_p_oeqkkl , f )
C ... Weights for spin moments
      if (lswtk>0 .and. nspc==2) then
          allocate(evecc(ndimh,2,ndimh,2),work(ndimh,2,ndimh,2))
          call zcopy(ndimhx**2,evec,1,evecc,1)
          call zgetrf(nevl,nevl,evecc,ndimhx,ipiv,i)
          if (i .ne. 0) call rx('addrbl: failed to generate overlap')
          call zgetri(nevl,evecc,ndimhx,ipiv,work,ndimhx**2,i)
          do  i = 1, ndimh
            do  k = 1, ndimh
              swtk(i,1,iq) = swtk(i,1,iq) + evecc(i,1,k,1)*evec(k,1,i,1)
     .        - evecc(i,1,k,2)*evec(k,2,i,1)
              swtk(i,2,iq) = swtk(i,2,iq) + evecc(i,2,k,1)*evec(k,1,i,2)
     .        - evecc(i,2,k,2)*evec(k,2,i,2)
            enddo
          enddo
          deallocate(evecc,work)
      endif
      deallocate(qpgv,qpg2v,ylv)
      call tcx('addrbl')
      end subroutine addrbl


      subroutine addsds(ndimh,evl,wgt,emin,emax,ndos,dos) !,esmear
      use m_lmfinit,only: bz_w,bz_n
C- Add to sampling dos
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndimh :hamiltonian dimension
Ci   evl   :eigenvalues
Ci   wgt   :eigenvalue weights
Ci   emin  :lower bound for DOS
Ci   emax  :upper bound for DOS
Ci   esmear:Parameter that describes gaussian broadening.
Ci         :Integer part >0 for for generalized gaussian broadening
Ci         :and is the the Methfessel-Paxton integration order
Ci         :Fractional part is the broadening width.
Ci         :Integer part <0 => Fermi-Dirac broadening used
Ci         :Fractional part is the temperature
Ci         :(see delstp.f)
Ci         :integer part above 100's digit is stripped.
Ci   ndos  :dimensions dos
Co Outputs
Co   dos   :DOS accumulated for these eigenvalues
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Jan 05 Extension of esmear to Fermi distribution
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ndimh,ndos
      double precision evl(1),dos(ndos,2),wgt,emin,emax,esmear
C ... Local parameters
      integer nord,ie,i1,i2,i
      double precision width,de,eigval,ei,sn,dn,fn,x
c      width = dabs(esmear) - int(dabs(esmear))
c     nord = dsign(1d0,esmear) * mod(int(dabs(esmear)),100)
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


      end module 
