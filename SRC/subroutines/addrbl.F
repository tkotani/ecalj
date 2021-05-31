CSFCPP#define F90 1
      subroutine addrbl ( isp,q, lrout,
     .     iq,lfrce, smpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,
     .     evec,evl,nevl,  ef0,def,emin,emax,ndos,
     o     dos, smrho, sumqv,sumev, sv_p_oqkkl,sv_p_oeqkkl, f)
      use m_struc_def 
      use m_suham,only: iv_a_oindxo
      use m_lmfinit,only:lat_alat,nbas, ssite=>v_ssite,sspec=>v_sspec,nsp,nspc,lmet=>bz_lmet,lekkl,lcplxp,
     &     ctrl_zbak!,bz_dosw
      use m_lattic,only: lat_qlat, lat_vol
      use m_supot,only: lat_nabc,k1,k2,k3
      use m_igv2x,only: napw,ndimh,ndimhx,igapw=>igv2x,m_Igv2x_set
      use m_suham,only: ndham=>ham_ndham,ham_ldham,ndhamx=>ham_ndhamx,nspx=>ham_nspx 
      use m_subzi, only: lwtkb,nevmx,numq,  lswtk,wtkb=>rv_a_owtkb,swtk=>rv_a_oswtk,ldos
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
Ci   lrout :0 no eigenvectors generated; exit after accumulating
Ci         :  eigenvalue weights and dos
Ci   wtkb  :integration weights, needed if lwtkb is 1 or 2
Ci   lswtk :<1 do nothing
Ci         :1 given a set of weights, make 'spin weights' swtk
Ci   iq    :index to current k-point
Ci   lfrce :if nonzero, accumulate contribution to force
Ci   ldos  :if nonzero, accumulate density-of-states
Ci   k1,k2,k3 dimensions of smpot,smrho
Ci   smpot :smooth potential on uniform mesh (mkpot.f), for forces
Ci   vconst:additional constant potential
Co   numq  :number of Fermi levels. Usu. 1, but numq=3 if lmet=4
Ci   qval  :total valence charge
Ci   evec  :eigenvectors
Ci   evl   :eigenvalues
Ci   nev
Ci   ef0   :estimate for fermi level
Ci   def   :When lmet=4, charge also accmulated for ef0+def and ef0-def
Cixx   esmear:(sampling integration) gaussian broadening
Cixx         :sign and integer part are extensions; see mkewgt.f
Ci   emin  :energy lower bound when adding to sampling dos
Ci   emax  :energy upper bound when adding to sampling dos
Ci   ndos  :number of energy mesh points
Ci   osig,otau,oppi  augmentation matrices
ci   lcplxp:0 if ppi is real; 1 if ppi is complex
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
C ... Passed parameters
      integer isp,iq,lrout,
     .lfrce,ndos !,numq,lekkl !,lcplxp,napw,igapw(3,napw)ndham,ndimh,
      integer nevl
      double precision def,ef0,emax,emin,qval,vconst !,wtkp(*)
      real(8):: q(3) ,  dos(ndos,2,isp) !wtkb(ndham,nsp,1) , swtk(ndham,nsp,1) ,
     ., evl(ndham,nsp) , sumev(2,3) , sumqv(3,2) , f(3,*)
c      type(s_site)::ssite(*)
c      type(s_spec)::sspec(*)
c      type(s_lat)::slat
c      type(s_ham)::sham
      type(s_rv1) :: sv_p_oppi(3,1)
      type(s_rv1) :: sv_p_otau(3,1)
      type(s_rv1) :: sv_p_osig(3,1)
      type(s_rv1) :: sv_p_oeqkkl(3,1)
      type(s_rv1) :: sv_p_oqkkl(3,1)
      double complex evec(ndimh,nspc,ndimh,nspc),smrho(k1,k2,k3,isp),
     .smpot(k1,k2,k3,isp)
C ... Local parameters
      integer:: i , k , nevec , nglob , ngabc(3) , n1 , n2 , 
     . n3 ,  lmxax , lmxa , nlmax , nlmto , ig !ndimhx ,
      double precision vavg,wgt,alat,qlat(3,3),tpiba,vol
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer ipiv(ndimh*2),i_copy_size
      complex(8),allocatable:: evecc(:,:,:,:),work(:,:,:,:)
      real(8),allocatable:: qpgv(:,:),qpg2v(:),ylv(:,:)
      double precision ewgt(ndimh*nspc*numq),epsnevec

      qval= qval_- ctrl_zbak(1) 
c      emin=bz_dosw(1)
c      emax=bz_dosw(2)
      
      if (lwtkb .lt. 0) return
      call tcn('addrbl')
      nlmto = ndimh-napw
c      ndimhx = ndimh*nspc
      i_copy_size=size(lat_nabc)
      call icopy(i_copy_size,lat_nabc,1,ngabc,1)
C     Setup for PW part of basis
      alat=lat_alat
      i_copy_size=size(lat_qlat)
      call dcopy(i_copy_size,lat_qlat,1,qlat,1)
C     Find largest lmxa ... should be made elsewhere
      lmxax = -1
      do  i = 1, nbas
        k=ssite(i)%spec
        lmxa=sspec(k)%lmxa
        lmxax = max(lmxax,lmxa)
      enddo
      nlmax=(lmxax+1)**2
CSFCPP#if F90
      if (napw .gt. 0) then
        allocate(ylv(napw,nlmax),qpgv(3,napw),qpg2v(napw))
        tpiba = 2d0*4d0*datan(1d0)/alat
        do  ig = 1, napw
          qpgv(:,ig) = tpiba * ( q + matmul(qlat,igapw(:,ig)) )
        enddo
        call ropyln(napw,qpgv(1,1:napw),qpgv(2,1:napw),qpgv(3,1:napw),
     .  lmxax,napw,ylv,qpg2v)
      else
        allocate(ylv(1,1),qpgv(1,1),qpg2v(1))
      endif
CSFCPP#else
CSFCPP      if (napw .gt. 0) then
CSFCPP        call rx('require F90 for PW code')
CSFCPP      endif
CSFCPP#endif

C --- Decide how many states to include and make their weights ---
C ... Case band weights not passed: make sampling weights
      if (lwtkb .eq. 0) then
        call rxx(nspc.ne.1,'lwtkb=0 not implemented in noncoll case')
        wgt = abs(wtkp(iq))/nsp
        call mkewgt(lmet,wgt,qval,ndimh,evl(1,isp),ef0,def,numq, !!esmear
     .  nevec,ewgt,sumev,sumqv(1,isp))
        call dscal(nevec*numq,wgt,ewgt,1)
C ... Case band weights are passed
      else
        if (numq .ne. 1) call rx('addbrl: incompatible numq')
        call dcopy(nevl,wtkb(1,isp,iq),1,ewgt,1)
        do  10  i = nevl, 1, -1
          nevec = i
          if (abs(wtkb(i,isp,iq)) .gt. epsnevec()) goto 12
   10   continue
   12   continue
      endif

C ... Add to sampling dos
      if (ldos .eq. 1) call addsds(nevl,evl(1,isp),abs(wtkp(iq))/nsp,
     .emin,emax,ndos,dos(1,1,isp)) !,esmear

C ... Force from smooth analytic hamiltonian and overlap
      if (lfrce .gt. 0 .and. lrout .gt. 0) then
        call rxx(nspc.ne.1,'forces not implemented in noncoll case')
        vavg = vconst
c        print *, '!! avg=1=smpot=1'; vavg=1; smpot=1 !; nevec=1
C        print *, '!! smpot=vavg'; smpot=vavg
C       call zprm('evecs',2,evec,ndimh,ndimh,nevec)

        if (nlmto .gt. 0) then
          call fsmbl ( nbas , ssite , sspec ,  vavg , q , ndimh , !slat ,
     .     nlmto , iv_a_oindxo , numq , nevec , evl ( 1 , isp ) , evec , 
     .     ewgt , f )
        endif
C        print *, 'after fsmbl'
C        do  i = 1, 3
C          print 543, f(:,i)
C  543     format(1p3e15.7)
C        enddo

        if (napw .gt. 0) then
          vol=lat_vol
          call fsmbpw ( nbas , ssite , sspec , vavg , ndimh , nlmto , iv_a_oindxo 
     .     , numq , nevec , evl ( 1 , isp ) , evec , ewgt , napw , qpgv 
     .     , qpg2v , ylv , nlmax , lmxax , alat , dsqrt ( vol ) , f )
C        print *, 'after fsmblpw'
C        do  i = 1, 3
C          print 543, f(:,i)
C        enddo
        endif
      endif

C ... Add to smooth density
      if (lrout .gt. 0) then
        call rsibl ( ssite , sspec , lfrce , nbas , isp , q , !slat , 
     .   iq , ndimh , nspc , napw , igapw , iv_a_oindxo , numq , nevec 
     .   , evec , ewgt , k1 , k2 , k3 , smpot , smrho , f )
C       call zprm3('smrho after rsibl',0,smrho,k1,k2,k3)
      endif
C      print *, 'after rsibl'
C      do  i = 1, 3
C        print 543, f(:,i)
C      enddo

C ... Add to local density coefficients
      if (lrout .gt. 0) then
        call rlocbl ( ssite , sspec , lfrce , nbas , isp , q , !slat , 
     .        ndham , ndimh , nspc , napw , igapw , iv_a_oindxo , numq , nevec 
     .        , evec , ewgt , evl , sv_p_osig , sv_p_otau , sv_p_oppi , lcplxp 
     .        , lekkl , sv_p_oqkkl , sv_p_oeqkkl , f )
C       call zprm3('smrho after rlocbl',0,smrho,k1,k2,k3)
      endif

C      print *, 'after rlocbl'
C      do  i = 1, 3
C        print 543, f(:,i)
C      enddo


C ... Weights for spin moments
      if (lswtk .gt. 0 .and. nspc .eq. 2) then
c this section is moved in bndfp.F
c        if (ndimhx .ne. nevl) then 
c          call info0(30,0,0,' addrbl: eigenvector matrix not square'
c     .    //' ... spin weights not evaluated')
c          lswtk = -2
c        else
CSFCPP#if F90
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
C        print 345,i,swtk(i,1,iq),swtk(i,2,iq), swtk(i,1,iq)+swtk(i,2,iq)
C  345   format(i4,3f14.8)
          enddo
          deallocate(evecc,work)
c        endif
CSFCPP#else
CSFCPP        lswtk = -2
CSFCPP#endif
      endif

CSFCPP#if F90
      deallocate(qpgv,qpg2v,ylv)
CSFCPP#endif

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


