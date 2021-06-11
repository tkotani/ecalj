CSFCPP#define F90 1

c#define USEOPTIONAL 1

      module m_prlcb1
      contains
      subroutine prlcb1(mode,ndimh,ik,norb,ltab,ktab,rsmh,offl,nlmbx,
     .nlma,kmax,b0,db0,b,db)

C- Poke strux and grads from b0,db0 to full arrays b,db
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 poke b only
Ci         :1 poke b and db
Ci   ndimh :dimension of hamiltonian
Ci   ib    :site of strux head
Ci   nlmbx :dimensions b,db
Ci   nlma  :augmentation L-cutoff
Ci   kmax  :Pkl polynomial cutoff
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   b0    :L-ordered strux for one ik block and pair of sites
Ci   db0   :gradient of b0
Co Outputs
Co   b     :subblock corresponding to b0 is poked into b
Co   db    :subblock corresponding to db0 is poked into db
Cr Remarks
Cr   b0,db0 have normal L ordering in both row and column dimensions.
Cr   b,db   have normal L ordering in rows but iprmb ordering in columns
Cr   This routine is identical in function to paugq1 (augmbl.f) except:
Cr     the gradient db of b can optionally be filled
Cr     array indices to b are ordered differently
Cu Updates
Cu   25 Aug 04 Adapted to extended local orbitals
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,kmax,ndimh,ik,nlma,nlmbx
      integer n0,nkap0
      parameter (n0=10,nkap0=3)
      integer norb,ltab(norb),ktab(norb),offl(norb)
      double precision rsmh(n0,nkap0)
c#ifdef USEOPTIONAL
      double complex b0(0:kmax,nlma,nlmbx),b(ndimh,nlma,0:kmax)
      complex(8),optional::db0(0:kmax,nlma,nlmbx,3),db(ndimh,nlma,0:kmax,3)
c#else
c      double complex b0(0:kmax,nlma,nlmbx),b(ndimh,nlma,0:kmax),
c     .db0(0:kmax,nlma,nlmbx,3),db(ndimh,nlma,0:kmax,3)
c#endif
C ... Local parameters
      integer i1,ik1,k,ilma,ilmb,iorb,l1,nlm1,nlm2
      integer blks(norb),ntab(norb)
      double precision xx

C     Block into groups of consecutive l
      call gtbsl1(4+16,norb,ltab,ktab,rsmh,xx,ntab,blks)

      do  iorb = 1, norb
        ik1 = ktab(iorb)
C       Loop only over orbitals belonging to this energy block
        if (ik1 .eq. ik .and. blks(iorb) .ne. 0) then
          l1  = ltab(iorb)
          nlm1 = l1**2+1
          nlm2 = nlm1 + blks(iorb)-1
C         i1 = index to hamiltonian offset
          i1 = offl(iorb)
C          print *, iorb,ik,l1,i1+1,i1+nlm2-nlm1+1
C  333     format('block',i4,' ik=',i4,' l=',i4,
C     .      ' start=',i4,' stop=',i4)
          do  ilmb = nlm1, nlm2
            i1 = i1+1
            if (mode .eq. 0) then
              do  ilma = 1, nlma
                do  k = 0, kmax
                  b(i1,ilma,k) = b0(k,ilma,ilmb)
                enddo
              enddo
            else
              do  ilma = 1, nlma
                do  k = 0, kmax
                  b(i1,ilma,k) =    b0(k,ilma,ilmb)
                  db(i1,ilma,k,1) = db0(k,ilma,ilmb,1)
                  db(i1,ilma,k,2) = db0(k,ilma,ilmb,2)
                  db(i1,ilma,k,3) = db0(k,ilma,ilmb,3)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      end subroutine prlcb1

      end module m_prlcb1

      subroutine bstrux(mode,ssite,sspec,cg,indxcg,jcg,cy,iprmb, !,slat
     .nbas,ia,pa,rsma,q,kmax,nlma,ndimh,napw,igapw,b,db)
      use m_struc_def  !Cgetarg
      use m_prlcb1
      use m_lmfinit,only:lat_alat
      use m_lattic,only: lat_qlat
      use m_lattic,only: lat_vol
      use m_uspecb,only:uspecb

C- Structure constants for P_kL expansion of Bloch lmto + PW around site ia
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :Whether b or both db are used, and determines ordering
Ci   mode  :0 Only b            b = b(0:kmax, nlma, ndimh)
Ci          1 Both b and db     b = b(ndimh, nlma, 0:kmax)
Ci          2 Only b            b = b(ndimh, nlma, 0:kmax)
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: alat qlat vol
Ci     Stored:
Ci     Passed to: hxpbl hxpgbl
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:
Ci     Stored:
Ci     Passed to: uspecb
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nbas  :size of basis
Ci   ia    :augmentation around site ia
Ci   pa    :position of site ia
Ci   rsma  :augmentation smoothing radius
Ci   q     :q-point for Bloch sum
Ci   kmax  :polynomial cutoff
Ci   nlma  :number of augmentation channels
Ci   ndimh :dimension of hamiltonian
Ci   napw  :number of PWs in basis
Ci   igapw :list of APW PWs
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   alat  :length scale of lattice and basis vectors, a.u.
Cl Local variables
Cl   nlmto :number of lmto's = ndimh - napw
Co Outputs
Co   b     : mode0  b(0:kmax,nlma, ndimh)
Co         : mode1  b(ndimh, nlma, 0:kmax)
Co         : mode2  b(ndimh, nlma, 0:kmax)
Co   db    : mode1  db(ndimh,nlma, 0:kmax, 3)
Co         :        Gradient is wrt head shift; use -db for grad wrt pa
Cr Remarks
Cr   Coefficients b are referred to as C_kL in the LMTO book.
Cr   bstrux requires an F90 compiler
Cu Updates
Cu   14 Jan 09 Bug fix, mode=2
Cu   05 Jul 08 (T. Kotani) adapted from augmbl; new PW part.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      real(8):: pa(3) , q(3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
c      type(s_lat)::slat

      double precision cg(1),cy(1),rsma
      integer kmax,nbas,ndimh,mode,ia,nlma,napw
      integer iprmb(1),igapw(3,napw),indxcg(1),jcg(1)
      double complex b(*)  !b(0:kmax,nlma,ndimh) mode 0;
      !b(ndimh,nlma,0:kmax) mode 1,2
      double complex db(*) !db(ndimh,nlma,0:kmax,3)
C ... Local parameters
      integer nlmto,nlmbx,ib,is,ik,n0,nkap0,nkapi,norb,nlmh
      parameter (nlmbx=25,n0=10,nkap0=3)
      integer lh(nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      double precision p(3),xx,alat,qlat(3,3),vol,srvol
CSFCPP#if F90
      complex(8),allocatable:: b0(:),db0(:),bos(:)

C --- Setup ---

      alat=lat_alat
      qlat=lat_qlat
      vol=lat_vol

      srvol = dsqrt(vol)
      nlmto = ndimh-napw

C     Zero out strux to eliminate contributions from local orbitals
      call dpzero(b,(kmax+1)*nlma*ndimh*2)
      if (mode .eq. 1) then
        call dpzero(db,(kmax+1)*nlma*ndimh*3*2)
      endif

C --- b for MTO  (Written as C_kl in LMTO book) ---
      if (nlmto .gt. 0) then
        allocate(b0((kmax+1)*nlma*nlmbx),bos((kmax+1)*nlmbx))
        call dpzero(b0,(kmax+1)*nlma*nlmbx*2)
        if (mode .eq. 1) then
          call dpzero(db,(kmax+1)*nlma*ndimh*3*2)
          allocate(db0((kmax+1)*nlma*nlmbx*3))
        endif
        do  ib = 1, nbas
          is=ssite(ib)%spec
          p=ssite(ib)%pos
C       This line augments no local orbitals
C       ik = 1
C       This line augments onsite extended local orbitals only
C       if (ia .eq. ib) ik = 2
C       This line augments extended local orbitals all sites
c          call uspecb(0,2,sspec,is,is,lh,rsmh,eh,nkapi)
          call uspecb(is,lh,rsmh,eh,nkapi)
          
          
C       Position in h; l,k indices for orbitals connected w/ ib
          call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C       Loop over blocks of envelope functions
          do  ik = 1, nkapi
            nlmh = (lh(ik)+1)**2
            if (nlmh .gt. nlmbx) call rxi('augmbl: need nlmbx',nlmh)
            if (nlmh .gt. nlma .and. ia.eq.ib)
     .      call rx('augmbl: nlmh.gt.nlma')
            if (mode .eq. 0.or.mode .eq. 2) then
              call hxpbl(p,pa,q,rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh,
     .        nlma,kmax,nlma,cg,indxcg,jcg,cy,b0)!,slat
            elseif (mode .eq. 1) then
              call hxpgbl(p,pa,q,rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh,
     .        nlma,kmax,nlmbx,nlma,cg,indxcg,jcg,cy,b0, !slat,
     .        db0)
            else
              call rxi('bfactor: bad mode',mode)
            endif

c#ifndef ALL3C
            if (ib .eq. ia) then
              call hxpos(rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh,kmax,bos)
              call paugq2(kmax,nlmh,nlma,bos,b0)
            endif
c#endif

C         Note: indices of b are ordered differently by mode (see Outputs)
            if (mode .eq. 0) then
              call paugq1(kmax,nlma,kmax,ik,norb,ltab,ktab,rsmh,offl,
     .        b0,b)
            elseif (mode .eq. 1) then
              call prlcb1(1,ndimh,ik,norb,ltab,ktab,rsmh,offl,nlmbx,
     .        nlma,kmax,b0,db0,b,db)
            elseif (mode .eq. 2) then
c#ifdef USEOPTIONAL
              call prlcb1(mode=0,ndimh=ndimh,ik=ik,norb=norb,ltab=ltab,
     .        ktab=ktab,rsmh=rsmh,offl=offl,nlmbx=nlmbx,
     .        nlma=nlma,kmax=kmax,b0=b0,b=b)
c#else
c              call prlcb1(0,ndimh,ik,norb,ltab,ktab,rsmh,offl,nlmbx,
c     .        nlma,kmax,b0,xx,b,xx)
c#endif
            endif
          enddo
        enddo
        deallocate(b0,bos)
        if (mode .eq. 1) deallocate(db0)
      endif

      call paugqp(mode,kmax,nlma,kmax,ndimh,napw,igapw,alat,qlat,srvol,
     .q,pa,rsma,b,b,db)

      if (mode .eq. 0) then
C       call zprm('b',2,b,(kmax+1)*nlma,(kmax+1)*nlma,ndimh)
      endif
      if (mode .eq. 1 .or. mode .eq. 2) then
C       call zprm('b mode 1 or 2',2,b,ndimh,ndimh,(kmax+1)*nlma)
      endif

CSFCPP#else
CSFCPP      call rx('bstrux requires F90')
CSFCPP#endif
      end subroutine bstrux


      subroutine paugq1(kmax,nlma,k0,ik,norb,ltab,ktab,rsmh,offl,b0,b)

C- Poke strux from b0 to full array b
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :Pkl polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   k0    :dimensions b0 and b
Ci   ik    :energy
Co   norb  :number of orbital types for ib; see Remarks
Co   ltab  :table of l-quantum numbers for each type
Co   ktab  :table of energy index for each type
Co   offl  :offl(norb) offset in h to this block of orbitals
Ci   b0    :L-ordered strux for one ik block and pair of sites
Co Outputs
Co   b     :subblock corresponding to b0 is poked into b
Cr Remarks
Cr   b0 has normal L ordering in both row and column dimensions.
Cr   b  has normal L ordering in row but iprmb ordering in col dim.
Cr   This routine is identical in function to prlcb1 (rlocbl.f) except:
Cr     no gradient db in this routine
Cr     array indices to b are ordered differently
Cu Updates
Cu   25 Aug 04 Adapted to extended local orbitals
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kmax,nlma,k0,norb,ltab(norb),ktab(norb),offl(norb),ik
      integer n0,nkap0
      parameter (n0=10,nkap0=3)
      double precision rsmh(n0,nkap0)
      double complex b0(0:k0,nlma,1),b(0:k0,nlma,1)
C ... Local parameters
      integer ilmb,ilma,k,iorb,l1,ik1,i1,nlm1,nlm2
      integer blks(norb),ntab(norb)
      double precision xx

C     Block into groups of consecutive l
      call gtbsl1(4+16,norb,ltab,ktab,rsmh,xx,ntab,blks)

      do  iorb = 1, norb
        ik1 = ktab(iorb)
C       Loop only over orbitals belonging to this energy block
        if (ik1 .eq. ik .and. blks(iorb) .ne. 0) then
          l1  = ltab(iorb)
          nlm1 = l1**2+1
          nlm2 = nlm1 + blks(iorb)-1
C         i1 = index to hamiltonian offset
          i1 = offl(iorb)
          do  ilmb = nlm1, nlm2
            i1 = i1+1
            do  ilma = 1, nlma
              do  k = 0, kmax
                b(k,ilma,i1) = b0(k,ilma,ilmb)
              enddo
              do  k = kmax+1, k0
                b(k,ilma,i1) = 0d0
              enddo
            enddo
          enddo
        endif
      enddo

      end subroutine paugq1


      subroutine paugqp(mode,kmax,nlma,k0,ndimh,napw,igapw,alat,qlat,
     .srvol,q,pa,rsma,b0,b1,db)
C- Make PW part of strux b
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 Make  b = b0(0:kmax,nlma,ndimh)
Ci          1 Make  b = b1(ndimh,nlma,0:kmax)
Ci            and  db = db(ndimh,nlma,0:kmax,3)
Ci          2 Make  b = b1(ndimh,nlma,0:kmax)
Ci   kmax  :Pkl polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   k0    :dimensionsb
Ci   ndimh :hamiltonian dimension
Ci   napw  :number of augmented PWs in basis
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   srvol :sqrt(vol)
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   q     :q-point for Bloch sum
Ci   pa    :position of site ia
Ci   rsma  :augmentation smoothing radius
Co Outputs
Co   b     :PW part of 1-center epansion is poked into b
Cr Remarks
Cu Updates
Cu   05 Jul 08 (T. Kotani) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,kmax,nlma,k0,napw,igapw(3,napw),ndimh
      double precision rsma,alat,qlat(3,3),q(3),srvol,pa(3)
      double complex b0(0:k0,nlma,ndimh),b1(ndimh,nlma,0:k0),
     .db(ndimh,nlma,0:k0,3)
C ... Local parameters
      integer k,lmxa,ll,l,ig,ilm,m,nlmto
      double precision gamma,qpg(3),pi,tpiba,qpg2,ddot,facexp,
     .rsmal,pgint,dfac(0:kmax),fac2l(0:nlma),yl(nlma),fpi,fac ,qk
C     complex(8),allocatable:: bpw(:,:)
      double complex srm1,srm1l,gfourier,phase,facilm,b
      parameter (srm1=(0d0,1d0))

      if (napw .eq. 0) return

cccccccccccccccccccccccccccccccccccccccccccc
c      b0=0d0
c      if(mode==0) return
c      b1=0d0
c      db=0d0
c      return
cccccccccccccccccccccccccccccccccccccccccccc

CSFCPP#if F90
      nlmto = ndimh - napw
      pi = 4d0*datan(1d0)
      fpi = 4*pi
      tpiba = 2d0*pi/alat
      gamma = rsma**2/4d0
      lmxa = ll(nlma)
C     allocate(bpw(0:kmax,nlma))
      fac2l(0) = 1d0
      do  l = 1, lmxa+1
        fac2l(l) = fac2l(l-1) * (2*l-1)
      enddo
c     fac2l(l)=(2l-1)!! data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/ See msrtx3.f
      dfac(0) = 1d0
      do  k = 1, kmax
        dfac(k) = dfac(k-1)*k
      enddo
      do  ig = 1, napw
        qpg = tpiba * ( q + matmul(qlat,igapw(1:3,ig)) )
        call ropyln(1,qpg(1),qpg(2),qpg(3),lmxa,1,yl,qpg2)
        phase = exp(srm1*alat*ddot(3,qpg,1,pa,1))
        facexp = exp(-gamma*qpg2)
        ilm = 0
        rsmal = 1d0 !takao 1 to 1d0 June2011 (may give little effects).
        srm1l = 1d0 !

        if (mode .eq. 0) then
          do  l = 0, lmxa
            do  m = 1, 2*l+1
              ilm = ilm + 1
              facilm = srm1l*yl(ilm)
              fac = fac2l(l+1)/rsmal/fpi
              qk = 1
              do  k = 0, kmax
                pgint =  dfac(k)*fac        ! Eq. 12.8 in JMP39 3393
                gfourier = qk*facilm*facexp ! Eq. 5.17
                b0(k,ilm,ig+nlmto) = gfourier/pgint/srvol*phase
                fac = fac * 4/rsma**2
                qk = -qpg2 * qk
              enddo
            enddo
            rsmal = rsmal*rsma
            srm1l = srm1l * srm1
          enddo

        elseif (mode .eq. 1 .or. mode .eq. 2) then
          do  l = 0, lmxa
            do  m = 1, 2*l+1
              ilm = ilm + 1
              facilm = srm1l*yl(ilm)
              fac = fac2l(l+1)/rsmal/fpi
              qk = 1
              if (mode .eq. 1) then
                do  k = 0, kmax
                  pgint =  dfac(k)*fac        ! Eq. 12.8 in JMP39 3393
                  gfourier = qk*facilm*facexp ! Eq.5.17
                  b = gfourier/pgint/srvol*phase
                  b1(ig+nlmto,ilm,k) = b
                  db(ig+nlmto,ilm,k,1) = -srm1*qpg(1) * b
                  db(ig+nlmto,ilm,k,2) = -srm1*qpg(2) * b
                  db(ig+nlmto,ilm,k,3) = -srm1*qpg(3) * b
                  fac = fac * 4/rsma**2
                  qk = -qpg2 * qk
                enddo
              else
                do  k = 0, kmax
                  pgint =  dfac(k)*fac        ! Eq. 12.8 in JMP39 3393
                  gfourier = qk*facilm*facexp ! Eq.5.17
                  b = gfourier/pgint/srvol*phase
                  b1(ig+nlmto,ilm,k) = b
                  fac = fac * 4/rsma**2
                  qk = -qpg2 * qk
                enddo
              endif
            enddo
            rsmal = rsmal*rsma
            srm1l = srm1l * srm1
          enddo

        else
          call rxi('paugqp: bad mode',mode)
        endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C        print *,' --- test code: Pkl fitting vs. true (bessel) ---'
C        call upack('spec rmt',sspec,isa,rmt,0,0,0)
C        ilm = 2                 !any lm for test
C        ndiv= 30
C        l = ll(ilm)
C        print *,' ilm l ig=',ilm,l,ig
C        allocate( fi(0:lmxa),gi(0:lmxa), pkl(0:kmax,0:lmxa) )
C        do ix=0,ndiv
C          if(ix==0) then
C            rr=1d-4
C          else
C            rr= ix/dble(ndiv) *rmt
C          endif
CC ... rf2: exact solution.  jl(|q+G| r)\times 4*pi* i**l *YL(q+G) [expansion of exp(i q+G r)]
C          absqpg = sqrt(qpg2)
C          call bessl((absqpg*rr)**2, lmxa, fi, gi)
C          bess = fi(l) *(absqpg*rr)**l !bessel function jl(qbsqpg*rr)
C          rf2 = bess * 4*pi*srm1**l * cy(ilm)*yl(ilm)/absqpg**l
C
CC ... rf1: sum of Pkl
C          call radpkl(rr,rsma,kmax,lmxa,kmax,pkl)
C          rf1 = 0d0
C          do k=0,kmax
Cc              print *,' k b pkl=',k, b(k,ilm,ig+nlmto), pkl(k,l)
C            rf1 = rf1 +  b(k,ilm,ig+nlmto) * pkl(k,l)*rr**l
C          enddo
C
C          zz= rr*absqpg
C          write(6,"(f6.3,3x,2d12.4,'    ratio=',2d12.4 ,3x,12d12.4)")
C     .      rr, rf2, rf1/rf2
C     .      ,(sin(zz)/zz)       !j0
C     .      ,(sin(zz)/zz**2 - cos(zz)/zz) !j1
C     .      ,(3*sin(zz)/zz**3 - 3*cos(zz)/zz**2-sin(zz)/zz) !j2
C        enddo
C        stop 'xxxxxxxxx test end xxxxxxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   ... Phase factor from atomic position
c       b(:,:, nlmto+ig) = bpw(:,:) *exp( srm1*alat*sum(qpg*pa))
C       Note:  construction wasteful of memory!  clean up
C        if (mode .eq. 0) then   ! b(ik,ilm,nlmto+ig)
C          call zcopy( (kmax+1)*nlma,
C     .              bpw* exp( srm1*alat*sum(qpg*pa)), 1,
C     .              b( (kmax+1)*nlma*(nlmto+ig-1)+1),1)
C        elseif (mode .eq. 1) then
C          do  ilm = 1,nlma
C            do  k = 0, kmax
C              idx = nlmto+ig + (ilm-1)*ndimh + k*ndimh*nlma ! b(nlmto+ig,ilm,ik)
C              b(idx) = exp(srm1*alat*sum(qpg*pa)) * bpw(k,ilm)
C            enddo
C          enddo
C        endif

      enddo
C     deallocate(fac2l,yl,bpw)

CSFCPP#else
CSFCPP      call rx('bstrux requires F90')
CSFCPP#endif

      end subroutine paugqp


      subroutine paugq2(kmax,nlmh,nlma,bos,b0)
C- Subtract on-site strux for ib=ia, leaving tail expansion
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :polynomial cutoff in PkL expansion
Ci   nlmh  :L-cutoff for this (site,energy) block
Ci   nlma  :dimensions b0
Ci   bos   :on-site strux
Co Outputs
Co   b0    :On-site part of strux subtracted
Cr Remarks
Cr   b0 has normal L ordering in both row and column dimensions.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kmax,nlmh,nlma
      double precision bos(0:kmax,nlmh)
      double complex b0(0:kmax,nlma,1)
C ... Local parameters
      integer ilm,k
      do  ilm = 1, nlmh
        do  k = 0, kmax
          b0(k,ilm,ilm) = b0(k,ilm,ilm)-bos(k,ilm)
        enddo
      enddo
      end subroutine paugq2



