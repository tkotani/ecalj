module m_bstrux
  !! Structure constants for P_kL expansion of Bloch lmto + PW around site ia
  use m_struc_def,only: s_cv3,s_cv4
  use m_ftox
  use m_lgunit,only:stdo
  use m_MPItk,only:procid
  public:: bstrux_set, bstr, dbstr, m_bstrux_init
  complex(8),pointer,protected::  bstr(:,:,:)
  complex(8),pointer,protected:: dbstr(:,:,:,:)
  type(s_cv3),allocatable,protected,private,target:: p_bstr(:,:)
  type(s_cv4),allocatable,protected,private,target:: p_dbstr(:,:)
  real(8),allocatable,private:: qall(:,:)
  private
contains

  subroutine bstrux_set(ia,qin)
    use m_qplist,only: qplist,iqini,iqend
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat
    use m_lmfinit,only: lfrce=>ctrl_lfrce
    use m_ftox
    use m_lgunit,only:stdo
    implicit none
    real(8):: qin(3),q(3),eps=1d-10
    integer:: iq,iqx,ia
    !sss    call shorbz(qin,q,qlat,plat)
    q=qin !sss
    iq=-999
    do iqx=iqini,iqend
       if( sum( (q-qall(:,iqx))**2 )<eps) then
          iq=iqx
          exit
       endif
    enddo
    if(iq==-999) then
       write(stdo,ftox)'qin=',ftof(qin)
       write(stdo,ftox)'  q=',ftof(q)
       do iqx=iqini,iqend
          write(stdo,ftox)'bstrux_set iq q=',procid,iqx,ftof(qall(:,iqx))
       enddo   
       call rx('err:bstrux_set')
    endif   
    bstr => p_bstr(ia,iq)%cv3
    if(lfrce/=0) dbstr=> p_dbstr(ia,iq)%cv4
  end subroutine bstrux_set
  
  subroutine m_bstrux_init() !q for qplist --> not yet for sugw.
    use m_qplist,only: qplist,iqini,iqend
    use m_lmfinit,only: lfrce=>ctrl_lfrce,nlmax,kmxt,nspec,nbas,ispec,sspec=>v_sspec,rsma
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat,rv_a_opos
    use m_igv2x,only: napw, igvapw=>igv2x, ndimh,m_Igv2x_setiq !igvapwin=>igv2x,
    integer:: kmaxx,ia,isa,lmxa,lmxb,kmax,nlmb,nlma,mode,inn(3),ig,iq,ndimhmax
    real(8):: pa(3),qin(3),q(3),qlatinv(3,3),qss(3)
    !    integer,allocatable:: igvapw(:,:)
    call tcn('m_bstrux_init')
    if(allocated(qall)) deallocate(qall,p_bstr,p_dbstr)
    allocate(qall(3,iqini:iqend),p_bstr(nbas,iqini:iqend),p_dbstr(nbas,iqini:iqend))
                    !iqini:iqend for each rank
    do 1200 iq = iqini, iqend !This is a big iq loop
       qin = qplist(:,iq)
       call m_Igv2x_setiq(iq) ! Get napw,ndimh, and so on for given iq
       !allocate(igvapw(3,napw))
       !! See hambl.F calling augmbl, calling bstrux_set
       !! input and output
       !!   qpg(ig) = tpiba * ( qin + matmul(qlat,igapwin(1:3,ig))) for h,o,hso
       !! internal
       !!   qpg(ig) = tpiba * ( q  + matmul(qlat,igapw(1:3,ig)))
       !!  NOTE: both qpg are the same for given ig.
       !! qlat*igapw = qlat*igqwin + (qin-q) ---> igvapw = igvapwin + matmul(qlatinv,qin-q)
       !sss call shorbz(qin,q,qlat,plat) !Get q. Is this fine?
       q=qin
       !qlatinv = transpose(plat)
       !inn = nint(matmul(qlatinv,qin-q))
       !       do ig=1,napw
       !          igvapw(:,ig) = inn + igvapwin(:,ig)
       !       enddo
       !igvapw(1:3,1:napw) = igvapwin(1:3,1:napw)
       do ia=1,nbas
          isa=ispec(ia) 
          pa=rv_a_opos(:,ia) 
          lmxa=sspec(isa)%lmxa !max l of augmentation
          kmax=sspec(isa)%kmxt !max of radial k
!          rsma=sspec(isa)%rsma
          nlma = (lmxa+1)**2
          if (lmxa == -1) cycle
          allocate(                p_bstr(ia,iq)%cv3(ndimh,nlma,0:kmax) )
          if(lfrce/=0) allocate(  p_dbstr(ia,iq)%cv4(ndimh,nlma,0:kmax,3) )
          !  --- Make strux to expand all orbitals at site ia ---
          mode = 2
          if(lfrce/=0) mode=1
          qss = q+ [1d-8,2d-8,3d-8] !for stabilizing deneracy ordering (this works well?)
          call bstrux (mode,ia,pa,rsma(isa),qss,kmax,nlma,ndimh,napw,igvapw, p_bstr(ia,iq)%cv3,p_dbstr(ia,iq)%cv4)
       enddo
       !       deallocate(igvapw)
       qall(:,iq)=q
       !write(stdo,ftox)'m_bstrux_init qin',procid,iq,ftof(qin)
       !write(stdo,ftox)'m_bstrux_init q  ',procid,iq,ftof(qall(:,iq))
1200 enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do iq=iqini,iqend
!       write(stdo,ftox)'m_bstrux_init',procid,iq,ftof(qall(:,iq))
!    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    call tcx('m_bstrux_init')
  end subroutine m_bstrux_init
  ! sssssssssssssssssssssssssssssssssssssss
  subroutine bstrux(mode,ia,pa,rsma,q,kmax,nlma,ndimh,napw,igapw,  b, db)
    use m_smhankel,only: hxpbl,hxpgbl
    use m_struc_def
    use m_lmfinit,only:alat=>lat_alat,lhh,nkaphh,nkapii,ispec,sspec=>v_sspec,nbas
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol,rv_a_opos
    use m_uspecb,only: uspecb
    use m_orbl,only: Orblib, norb,ltab,ktab,offl
    use m_smhankel,only: hxpos
    !- Structure constants for P_kL expansion of Bloch lmto + PW around site ia
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :Whether b or both db are used, and determines ordering
    !i   mode  :0 Only b            b = b(0:kmax, nlma, ndimh)
    !i          1 Both b and db     b = b(ndimh, nlma, 0:kmax)
    !i          2 Only b            b = b(ndimh, nlma, 0:kmax)
    !i   sspec :struct for species-specific information; see routine uspec
    !i     Elts read:
    !i     Stored:
    !i     Passed to: uspecb
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   cy    :Normalization constants for spherical harmonics
    !i   nbas  :size of basis
    !i   ia    :augmentation around site ia
    !i   pa    :position of site ia
    !i   rsma  :augmentation smoothing radius
    !i   q     :q-point for Bloch sum
    !i   kmax  :polynomial cutoff
    !i   nlma  :number of augmentation channels
    !i   ndimh :dimension of hamiltonian
    !i   napw  :number of PWs in basis
    !i   igapw :list of APW PWs
    !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !l Local variables
    !l   nlmto :number of lmto's = ndimh - napw
    !o Outputs
    !o   b     : mode0  b(0:kmax,nlma, ndimh)
    !o         : mode1  b(ndimh, nlma, 0:kmax)
    !o         : mode2  b(ndimh, nlma, 0:kmax)
    !o   db    : mode1  db(ndimh,nlma, 0:kmax, 3)
    !o         :        Gradient is wrt head shift; use -db for grad wrt pa
    !r Remarks
    !r   Coefficients b are referred to as C_kL in the LMTO book.
    !r   bstrux requires an F90 compiler
    !u Updates
    !u   14 Jan 09 Bug fix, mode=2
    !u   05 Jul 08 (T. Kotani) adapted from augmbl; new PW part.
    ! ----------------------------------------------------------------------
    implicit none
    real(8):: pa(3) , q(3)
    double precision :: rsma
    integer :: kmax,ndimh,mode,ia,nlma,napw
    integer :: igapw(3,napw)
    double complex b(*)  !b(0:kmax,nlma,ndimh) mode 0;
    !b(ndimh,nlma,0:kmax) mode 1,2
    double complex db(*) !db(ndimh,nlma,0:kmax,3)
    ! ... Local parameters
    integer :: nlmto,nlmbx,ib,is,ik,n0,nkap0,nkapi,nlmh
    parameter (nlmbx=25,n0=10,nkap0=3)
    integer :: lh(nkap0)
    double precision :: eh(n0,nkap0),rsmh(n0,nkap0)
    !      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
    double precision :: p(3),xx,srvol
    complex(8),allocatable:: b0(:),db0(:)
    real(8),allocatable:: bos(:)
    call tcn('bstrux')
    srvol = dsqrt(vol)
    nlmto = ndimh-napw
    !     Zero out strux to eliminate contributions from local orbitals
    call dpzero(b,(kmax+1)*nlma*ndimh*2)
    if (mode == 1) then
       call dpzero(db,(kmax+1)*nlma*ndimh*3*2)
    endif
    ! --- b for MTO  (Written as C_kl in LMTO book) ---
    if (nlmto > 0) then
       allocate(b0((kmax+1)*nlma*nlmbx),bos((kmax+1)*nlmbx))
       call dpzero(b0,(kmax+1)*nlma*nlmbx*2)
       if (mode == 1) then
          call dpzero(db,(kmax+1)*nlma*ndimh*3*2)
          allocate(db0((kmax+1)*nlma*nlmbx*3))
       endif
       do  ib = 1, nbas
          is= ispec(ib) 
          p = rv_a_opos(:,ib) 
          call uspecb(is,rsmh,eh)!  Position in h; l,k indices for orbitals connected w/ ib
          call orblib(ib) !return norb,ltab,ktab,offl
          !       Loop over blocks of envelope functions
          do  ik = 1, nkaphh(is)
             nlmh = (lhh(ik,is)+1)**2
             if (nlmh > nlmbx) call rxi('augmbl: need nlmbx',nlmh)
             if (nlmh > nlma .AND. ia == ib) call rx('augmbl: nlmh > nlma')
             if (mode == 0 .OR. mode == 2) then
                call hxpbl(p,pa,q,rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh, &
                     nlma,kmax,nlma,b0)  !,cg,indxcg,jcg,cy
             elseif (mode == 1) then
                call hxpgbl(p,pa,q,rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh, &
                     nlma,kmax,nlmbx,nlma,b0, db0) !cg,indxcg,jcg,cy,
             else
                call rxi('bfactor: bad mode',mode)
             endif
             if (ib == ia) then
                call hxpos(rsmh(1,ik),rsma,eh(1,ik),kmax,nlmh,kmax,bos)
                call paugq2(kmax,nlmh,nlma,bos,b0)
             endif
             !         Note: indices of b are ordered differently by mode (see Outputs)
             if (mode == 0) then
                call paugq1(kmax,nlma,kmax,ik,norb,ltab,ktab,rsmh,offl, b0,b)
             elseif (mode == 1) then
                call prlcb1(1,ndimh,ik,norb,ltab,ktab,rsmh,offl,nlmbx,  nlma,kmax,b0,db0,b,db)
             elseif (mode == 2) then
                call prlcb1(mode=0,ndimh=ndimh,ik=ik,norb=norb,ltab=ltab, &
                     ktab=ktab,rsmh=rsmh,offl=offl,nlmbx=nlmbx, &
                     nlma=nlma,kmax=kmax,b0=b0,b=b)
             endif
          enddo
       enddo
       deallocate(b0,bos)
       if (mode == 1) deallocate(db0)
    endif
    call paugqp(mode,kmax,nlma,kmax,ndimh,napw,igapw,alat,qlat,srvol,q,pa,rsma,b,b,db)
    call tcx('bstrux')
  end subroutine bstrux
  subroutine prlcb1(mode,ndimh,ik,norb,ltab,ktab,rsmh,offl,nlmbx,nlma,kmax,b0,db0,b,db)
    !- Poke strux and grads from b0,db0 to full arrays b,db
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 poke b only
    !i         :1 poke b and db
    !i   ndimh :dimension of hamiltonian
    !i   ib    :site of strux head
    !i   nlmbx :dimensions b,db
    !i   nlma  :augmentation L-cutoff
    !i   kmax  :Pkl polynomial cutoff
    !i   b0    :L-ordered strux for one ik block and pair of sites
    !i   db0   :gradient of b0
    !o Outputs
    !o   b     :subblock corresponding to b0 is poked into b
    !o   db    :subblock corresponding to db0 is poked into db
    !r Remarks
    !r   b0,db0 have normal L ordering in both row and column dimensions.
    !r   b,db   have normal L ordering in rows 
    !r   This routine is identical in function to paugq1 (augmbl.f) except:
    !r     the gradient db of b can optionally be filled
    !r     array indices to b are ordered differently
    !u Updates
    !u   25 Aug 04 Adapted to extended local orbitals
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: mode,kmax,ndimh,ik,nlma,nlmbx
    integer :: n0,nkap0
    parameter (n0=10,nkap0=3)
    integer :: norb,ltab(norb),ktab(norb),offl(norb)
    double precision :: rsmh(n0,nkap0)
    double complex b0(0:kmax,nlma,nlmbx),b(ndimh,nlma,0:kmax)
    complex(8),optional::db0(0:kmax,nlma,nlmbx,3),db(ndimh,nlma,0:kmax,3)
    ! ... Local parameters
    integer :: i1,ik1,k,ilma,ilmb,iorb,l1,nlm1,nlm2
    integer :: blks(norb),ntab(norb),ol,oi,iblk
    double precision :: xx
    !     Block into groups of consecutive l
    call gtbsl1(4+16,norb,ltab,ktab,rsmh,xx,ntab,blks)
    do  iorb = 1, norb
       ik1 = ktab(iorb)
       if(ik1 /= ik) cycle
       ol = ltab(iorb)**2
       oi = offl(iorb)
       do iblk = 1, blks(iorb)
          if (mode == 0) then
             b(oi+iblk,1:nlma,:) = transpose(b0(:,1:nlma, ol+iblk))
          else
             b(oi+iblk,:,0:kmax) = transpose(b0(0:kmax,:,ol+iblk))
             do  k = 0, kmax
                db(oi+iblk,:,k,:) = db0(k,:,ol+iblk,:)
             enddo
          endif
       enddo
    enddo
  end subroutine prlcb1
  subroutine paugq1(kmax,nlma,k0,ik,norb,ltab,ktab,rsmh,offl,b0, b)
    !- Poke strux from b0 to full array b
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   kmax  :Pkl polynomial cutoff
    !i   nlma  :augmentation L-cutoff
    !i   k0    :dimensions b0 and b
    !i   ik    :energy
    !o   norb  :number of orbital types for ib; see Remarks
    !o   ltab  :table of l-quantum numbers for each type
    !o   ktab  :table of energy index for each type
    !o   offl  :offl(norb) offset in h to this block of orbitals
    !i   b0    :L-ordered strux for one ik block and pair of sites
    !o Outputs
    !o   b     :subblock corresponding to b0 is poked into b
    !r Remarks
    !r   b0 has normal L ordering in both row and column dimensions.
    !r   b  has normal L ordering in row 
    !r   This routine is identical in function to prlcb1 (rlocbl.f) except:
    !r     no gradient db in this routine
    !r     array indices to b are ordered differently
    !u Updates
    !u   25 Aug 04 Adapted to extended local orbitals
    ! ----------------------------------------------------------------------
    implicit none
    integer :: kmax,nlma,k0,norb,ltab(norb),ktab(norb),offl(norb),ik
    integer :: n0,nkap0
    parameter (n0=10,nkap0=3)
    double precision :: rsmh(n0,nkap0)
    double complex b0(0:k0,nlma,1),b(0:k0,nlma,1)
    integer :: ilmb,ilma,k,iorb,l1,ik1,i1,nlm1,nlm2
    integer :: blks(norb),ntab(norb),oi,ol,nn
    double precision :: xx
    call gtbsl1(4+16,norb,ltab,ktab,rsmh,xx,ntab,blks) !! Block into groups of consecutive l
    do  iorb = 1, norb
       ik1 = ktab(iorb) ! Loop only over orbitals belonging to this energy block
       nn = blks(iorb)
       if (ik1==ik .AND. nn/= 0) then !blks(iorb): size of block of iorb
          l1 = ltab(iorb)
          oi = offl(iorb)
          ol = l1**2
          b(0:kmax,   1:nlma,oi+1:oi+nn) = b0(0:kmax,1:nlma,ol+1:ol+nn)
          b(kmax+1:k0,1:nlma,oi+1:oi+nn) = 0d0
       endif
    enddo
  end subroutine paugq1
  subroutine paugqp(mode,kmax,nlma,k0,ndimh,napw,igapw,alat,qlat,srvol,q,pa,rsma, b0,b1,db)
    use m_ropyln,only: ropyln
    !- Make PW part of strux b
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 Make  b = b0(0:kmax,nlma,ndimh)
    !i          1 Make  b = b1(ndimh,nlma,0:kmax)
    !i            and  db = db(ndimh,nlma,0:kmax,3)
    !i          2 Make  b = b1(ndimh,nlma,0:kmax)
    !i   kmax  :Pkl polynomial cutoff
    !i   nlma  :augmentation L-cutoff
    !i   k0    :dimensionsb
    !i   ndimh :hamiltonian dimension
    !i   napw  :number of augmented PWs in basis
    !i   igapw :vector of APWs, in units of reciprocal lattice vectors
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   srvol :sqrt(vol)
    !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
    !i   q     :q-point for Bloch sum
    !i   pa    :position of site ia
    !i   rsma  :augmentation smoothing radius
    !o Outputs
    !o   b     :PW part of 1-center epansion is poked into b
    !r Remarks
    !u Updates
    !u   05 Jul 08 (T. Kotani) first created
    ! ----------------------------------------------------------------------
    implicit none
    integer :: mode,kmax,nlma,k0,napw,igapw(3,napw),ndimh
    double precision :: rsma,alat,qlat(3,3),q(3),srvol,pa(3)
    complex(8):: b0(0:k0,nlma,ndimh),b1(ndimh,nlma,0:k0), db(ndimh,nlma,0:k0,3)
    integer :: k,lmxa,ll,l,ig,ilm,m,nlmto
    double precision :: gamma,qpg(3),pi,tpiba,qpg2(1),ddot,facexp, &
         rsmal,pgint,dfac(0:kmax),fac2l(0:nlma),yl(nlma),fpi,fac ,qk
    double complex srm1,srm1l,gfourier,phase,facilm,b
    parameter (srm1=(0d0,1d0))
    if (napw == 0) return
    nlmto = ndimh - napw
    pi = 4d0*datan(1d0)
    fpi = 4*pi
    tpiba = 2d0*pi/alat
    gamma = rsma**2/4d0
    lmxa = ll(nlma)
    !     fac2l(l)=(2l-1)!! data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/ See msrtx3.f
    fac2l(0) = 1d0
    do  l = 1, lmxa+1
       fac2l(l) = fac2l(l-1) * (2*l-1)
    enddo
    dfac(0) = 1d0
    do  k = 1, kmax
       dfac(k) = dfac(k-1)*k
    enddo
    do  ig = 1, napw
       qpg = tpiba * ( q + matmul(qlat,igapw(1:3,ig)) )
       call ropyln(1,qpg(1),qpg(2),qpg(3),lmxa,1,yl,qpg2)
       phase = exp(srm1*alat*ddot(3,qpg,1,pa,1))
       facexp = exp(-gamma*qpg2(1))
       ilm = 0
       rsmal = 1d0 !takao 1 to 1d0 June2011 (may give little effects).
       srm1l = 1d0 !

       if (mode == 0) then
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
                   qk = -qpg2(1) * qk
                enddo
             enddo
             rsmal = rsmal*rsma
             srm1l = srm1l * srm1
          enddo
       elseif (mode == 1 .OR. mode == 2) then
          do  l = 0, lmxa
             do  m = 1, 2*l+1
                ilm = ilm + 1
                facilm = srm1l*yl(ilm)
                fac = fac2l(l+1)/rsmal/fpi
                qk = 1
                if (mode == 1) then
                   do  k = 0, kmax
                      pgint =  dfac(k)*fac        ! Eq. 12.8 in JMP39 3393
                      gfourier = qk*facilm*facexp ! Eq.5.17
                      b = gfourier/pgint/srvol*phase
                      b1(ig+nlmto,ilm,k) = b
                      db(ig+nlmto,ilm,k,1) = -srm1*qpg(1) * b
                      db(ig+nlmto,ilm,k,2) = -srm1*qpg(2) * b
                      db(ig+nlmto,ilm,k,3) = -srm1*qpg(3) * b
                      fac = fac * 4/rsma**2
                      qk = -qpg2(1) * qk
                   enddo
                else
                   do  k = 0, kmax
                      pgint =  dfac(k)*fac        ! Eq. 12.8 in JMP39 3393
                      gfourier = qk*facilm*facexp ! Eq.5.17
                      b = gfourier/pgint/srvol*phase
                      b1(ig+nlmto,ilm,k) = b
                      fac = fac * 4/rsma**2
                      qk = -qpg2(1) * qk
                   enddo
                endif
             enddo
             rsmal = rsmal*rsma
             srm1l = srm1l * srm1
          enddo
       else
          call rxi('paugqp: bad mode',mode)
       endif
       ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !        print *,' --- test code: Pkl fitting vs. true (bessel) ---'
       !        call upack('spec rmt',sspec,isa,rmt,0,0,0)
       !        ilm = 2                 !any lm for test
       !        ndiv= 30
       !        l = ll(ilm)
       !        print *,' ilm l ig=',ilm,l,ig
       !        allocate( fi(0:lmxa),gi(0:lmxa), pkl(0:kmax,0:lmxa) )
       !        do ix=0,ndiv
       !          if(ix==0) then
       !            rr=1d-4
       !          else
       !            rr= ix/dble(ndiv) *rmt
       !          endif
       !C ... rf2: exact solution.  jl(|q+G| r)\times 4*pi* i**l *YL(q+G) [expansion of exp(i q+G r)]
       !          absqpg = sqrt(qpg2)
       !          call bessl((absqpg*rr)**2, lmxa, fi, gi)
       !          bess = fi(l) *(absqpg*rr)**l !bessel function jl(qbsqpg*rr)
       !          rf2 = bess * 4*pi*srm1**l * cy(ilm)*yl(ilm)/absqpg**l

       !C ... rf1: sum of Pkl
       !          call radpkl(rr,rsma,kmax,lmxa,kmax,pkl)
       !          rf1 = 0d0
       !          do k=0,kmax
       !c              print *,' k b pkl=',k, b(k,ilm,ig+nlmto), pkl(k,l)
       !            rf1 = rf1 +  b(k,ilm,ig+nlmto) * pkl(k,l)*rr**l
       !          enddo

       !          zz= rr*absqpg
       !          write(6,"(f6.3,3x,2d12.4,'    ratio=',2d12.4 ,3x,12d12.4)")
       !     .      rr, rf2, rf1/rf2
       !     .      ,(sin(zz)/zz)       !j0
       !     .      ,(sin(zz)/zz**2 - cos(zz)/zz) !j1
       !     .      ,(3*sin(zz)/zz**3 - 3*cos(zz)/zz**2-sin(zz)/zz) !j2
       !        enddo
       !        stop 'xxxxxxxxx test end xxxxxxxxxxxxxxxxxxxx'
       ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       !   ... Phase factor from atomic position
       !       b(:,:, nlmto+ig) = bpw(:,:) *exp( srm1*alat*sum(qpg*pa))
       !       Note:  construction wasteful of memory!  clean up
       !        if (mode .eq. 0) then   ! b(ik,ilm,nlmto+ig)
       !          call zcopy( (kmax+1)*nlma,
       !     .              bpw* exp( srm1*alat*sum(qpg*pa)), 1,
       !     .              b( (kmax+1)*nlma*(nlmto+ig-1)+1),1)
       !        elseif (mode .eq. 1) then
       !          do  ilm = 1,nlma
       !            do  k = 0, kmax
       !              idx = nlmto+ig + (ilm-1)*ndimh + k*ndimh*nlma ! b(nlmto+ig,ilm,ik)
       !              b(idx) = exp(srm1*alat*sum(qpg*pa)) * bpw(k,ilm)
       !            enddo
       !          enddo
       !        endif
    enddo
  end subroutine paugqp
  subroutine paugq2(kmax,nlmh,nlma,bos,b0)
    !- Subtract on-site strux for ib=ia, leaving tail expansion
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   kmax  :polynomial cutoff in PkL expansion
    !i   nlmh  :L-cutoff for this (site,energy) block
    !i   nlma  :dimensions b0
    !i   bos   :on-site strux
    !o Outputs
    !o   b0    :On-site part of strux subtracted
    !r Remarks
    !r   b0 has normal L ordering in both row and column dimensions.
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    integer :: kmax,nlmh,nlma
    real(8):: bos(0:kmax,nlmh)
    complex(8):: b0(0:kmax,nlma,1)
    integer :: ilm,k
    do  ilm = 1, nlmh
       b0(0:kmax,ilm,ilm) = b0(0:kmax,ilm,ilm)-bos(0:kmax,ilm)
    enddo
  end subroutine paugq2
end module m_bstrux

