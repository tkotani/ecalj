!>Structure constants for P_kL expansion of Bloch lmto + PW around site ia
module m_bstrux 
  ! bstr are stored in p_bstr(ia,iq)%cv3(ndimh,nlma,0:kmax) by m_bstrx_init
  ! "call bstrux_set(ia,iq)" rerurns  bstr(ndimh,nlma,0:kmax) and dbstr.
  use m_lmfinit,only: lmxa_i=>lmxa, kmxt_i=>kmxt,afsym
  use m_struc_def,only: s_cv3,s_cv4
  use m_lgunit,only:stdo
  use m_MPItk,only:procid
  use m_ll,only:ll
  use m_ftox
  use m_nvfortran,only:findloc
  public:: bstrux_set, bstr, dbstr, m_bstrux_init
  complex(8),pointer,protected::  bstr(:,:,:)
  complex(8),pointer,protected:: dbstr(:,:,:,:)
  type(s_cv3),allocatable,protected,private,target:: p_bstr(:,:)
  type(s_cv4),allocatable,protected,private,target:: p_dbstr(:,:)
  real(8),allocatable,private:: qall(:,:)
  private
  integer,private:: iqii,iqee
contains
  subroutine bstrux_set(ia,qin)!set bstr and dbstr for given ibas and q
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat
    implicit none
    real(8):: qin(3),q(3),eps=1d-10
    integer:: iq,iqx,ia !!!!! 2023-04-25 obatadebug    q=qin !    call shorbz(qin,q,qlat,plat) !Get q. Is this fine?
    logical:: lll(iqii:iqee)
    lll=[(sum( (qin-qall(:,iqx))**2 )<eps,iqx=iqii,iqee)]
    iq = findloc(lll,value=.true.,dim=1)+iqii-1
    Errorexit:if(iq<=0) then !error exit
       write(stdo,ftox)'qin=',ftof(qin)
       do iqx=iqii,iqee; write(stdo,ftox)'bstrux_set iq q=',procid,iqx,ftof(qall(:,iqx));  enddo
       call rx('err:bstrux_set can not find given qin')
    endif Errorexit
    bstr => p_bstr(ia,iq)%cv3
    dbstr=> p_dbstr(ia,iq)%cv4
  end subroutine bstrux_set
  subroutine m_bstrux_init() !q for qplist --> not yet for sugw.
    use m_qplist,only: qplist,iqini,iqend,nkp
    use m_lmfinit,only: nlmax,kmxt,nspec,nbas,ispec,rsma
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat,rv_a_opos
    use m_igv2x,only: napw, igvapw=>igv2x, ndimh,m_Igv2x_setiq !igvapwin=>igv2x,
    integer:: kmaxx,ia,isa,lmxa,lmxb,kmax,nlmb,nlma,inn(3),ig,iq,ndimhmax!,iqii,iqee
    real(8):: pa(3),qin(3),q(3),qlatinv(3,3),qss(3),qxx(3)
    logical:: cmdopt0
    call tcn('m_bstrux_init')
    if(allocated(qall)) deallocate(qall,p_bstr,p_dbstr)
    iqii=iqini
    iqee=iqend
    if(afsym) iqii=1   ! AF rotation map q points to another q point (rotwave); thus we need p_bstr(:,iq_mapped).
    if(afsym) iqee=nkp ! bugfix at 2023-11-18
    !                    For simplicity, we calculate p_bstr(:,iq) for all iq for all processes. This may be not efficients but works.
    allocate(qall(3,iqii:iqee),p_bstr(nbas,iqii:iqee),p_dbstr(nbas,iqii:iqee))
    iqloop: do 1200 iq = iqii, iqee ! iqii:iqee for each rank
       qin = qplist(:,iq)
       call m_Igv2x_setiq(iq) ! Get napw,ndimh,igvapw and so on for given iq !!! 2023-04-25 obatadebug call shorbz(qin,q,qlat,plat) !Get q. Is this fine?
       q=qin
       ibasloop: do ia=1,nbas !  --- Make strux to expand all orbitals at site ia ---
          isa=  ispec(ia) 
          pa =  rv_a_opos(:,ia) 
          lmxa= lmxa_i(isa) !max l of augmentation
          kmax= kmxt_i(isa) !max of radial k
          nlma = (lmxa+1)**2
          if (lmxa == -1) cycle
          allocate(   p_bstr(ia,iq)%cv3(ndimh,nlma,0:kmax)   )
          allocate(  p_dbstr(ia,iq)%cv4(ndimh,nlma,0:kmax,3) )
          qss = q !+ [1d-8,2d-8,3d-8] for stabilizing deneracy ordering (this works well?) --> this conflict with qshortn 
          call bstrux(ia,pa,rsma(isa),qss,kmax,nlma,ndimh,napw,igvapw, p_bstr(ia,iq)%cv3,p_dbstr(ia,iq)%cv4)
       enddo ibasloop
       qall(:,iq)=q  !write(stdo,ftox)'m_bstrux_init qin',procid,iq,ftof(qin),ftof(qall(:,iq))
1200 enddo iqloop
    call tcx('m_bstrux_init')
  end subroutine m_bstrux_init
  subroutine bstrux(ia,pa,rsma,q,kmax,nlma,ndimh,napw,igapw,  b, db) !Structure constants for P_kL expansion of Bloch lmto + PW around site ia
    use m_smhankel,only: hxpbl,hxpgbl
    use m_lmfinit,only:alat=>lat_alat,lhh,nkaphh,nkapii,ispec,nbas,n0,nkap0
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol,rv_a_opos
    use m_uspecb,only: uspecb
    use m_orbl,only: Orblib, norb,ltab,ktab,offl
    use m_smhankel,only: hxpos
    use m_ropyln,only: ropyln
    use m_ftox
    !i Inputs
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
    !o     b(ndimh,nlma, 0:kmax)
    !o    db(ndimh,nlma, 0:kmax, 3)     ! Gradient is wrt head shift; use -db for grad wrt pa
    !l Local variables
    !l   nlmto :number of lmto's = ndimh - napw
    !r Remarks  Coefficients b are referred to as C_kL in the LMTO book.
    implicit none
    integer :: kmax,ndimh,ia,nlma,napw, nlmto,ib,is,ik,nkapi,nlmh,k,lmxa,l,ig,ilm,m, igapw(3,napw),lh(nkap0),ikx
    complex(8):: b(ndimh,nlma,0:kmax), db(ndimh,nlma,0:kmax,3)
    real(8):: pa(3), q(3),rsma,eh(n0,nkap0),rsmh(n0,nkap0),p(3),xx,srvol
    real(8):: gamma,qpg(3),tpiba,qpg2(1),facexp,rsmal,pgint,dfac(0:kmax),fac2l(0:nlma),yl(nlma),fac
    complex(8):: srm1=(0d0,1d0),gfourier,phase,facilm
    real(8),parameter:: pi = 4d0*datan(1d0),fpi = 4*pi
    call tcn('bstrux')
    srvol = dsqrt(vol)
    nlmto = ndimh-napw
    b=0d0 
    db=0d0
    if(nlmto==0) goto 500 
    do ib= 1, nbas !MTO part 
       is= ispec(ib) 
       p = rv_a_opos(:,ib) 
       call uspecb(is,rsmh,eh)
       call orblib(ib) !Return norb,ltab,ktab,offl
       do ik = 1, nkaphh(is) ! Loop over blocks of envelope functions
          nlmh = (lhh(ik,is)+1)**2
          b0tob: block
            integer::ol,oi,ik1,iorb,iblk,k,ntab(norb),blks(norb)
            complex(8):: b0(0:kmax,nlma,nlmh),db0(0:kmax,nlma,nlmh,3)
            real(8):: bos(0:kmax,nlmh)
            call hxpgbl(p,pa,q,rsmh(:,ik),rsma,eh(:,ik),kmax,nlmh,nlma,kmax,nlmh,nlma,b0,db0)
            if (ib == ia) then
               call hxpos(rsmh(:,ik),rsma,eh(:,ik),kmax,nlmh,kmax,bos)    !Subtract on-site strux for ib=ia, leaving tail expansion
               forall(ilm = 1:nlmh) b0(0:kmax,ilm,ilm) = b0(0:kmax,ilm,ilm)-bos(0:kmax,ilm)
            endif
            call gtbsl1(0,norb,ltab,ktab,rsmh,xx,ntab,blks)
            do  iorb = 1, norb
               ik1 = ktab(iorb)
               if(ik1 /= ik) cycle
               ol = ltab(iorb)**2 !atomic  offset
               oi = offl(iorb)    !overall offset
               do iblk = 1, blks(iorb)
                  do  k = 0, kmax
                     b(oi+iblk,1:nlma,k)    = b0(k, 1:nlma,ol+iblk)  
                     db(oi+iblk,1:nlma,k,1:3) = db0(k,1:nlma,ol+iblk,1:3)
                  enddo
               enddo
            enddo
          endblock b0tob
       enddo
    enddo
500 continue
!    write(6,*)'xxxxxxxxx1111111'
    if(napw == 0) goto 1000
    tpiba = 2d0*pi/alat !APW part ! call paugqp(kmax,nlma,ndimh,napw,igapw,alat,qlat,srvol,q,pa,rsma,b,db)
    gamma = rsma**2/4d0
    lmxa = ll(nlma)
    fac2l(0) = 1d0
    do  l = 1, lmxa+1
       fac2l(l) = fac2l(l-1) * (2*l-1) !   fac2l(l)=(2l-1)!! data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/ See msrtx3.f
    enddo
    dfac(0) = 1d0
    do  k = 1, kmax
       dfac(k) = dfac(k-1)*k
    enddo
    do  ig = 1, napw
       qpg = tpiba * ( q + matmul(qlat,igapw(1:3,ig)) )
       call ropyln(1,qpg(1),qpg(2),qpg(3),lmxa,1,yl,qpg2)
       phase = exp(srm1*alat*sum(qpg*pa))
       facexp = exp(-gamma*qpg2(1))
       do ilm=1,(lmxa+1)**2
          l=ll(ilm)
          facilm = srm1**l*yl(ilm)
          fac = fac2l(l+1)/rsma**l/fpi
          do  k = 0, kmax
             pgint =  dfac(k)*fac*(4/rsma**2)**k    ! Eq. 12.8 in JMP39 3393
             gfourier = (-qpg2(1))**k*facilm*facexp ! Eq.5.17
             b(ig+nlmto,ilm,k)   = gfourier/pgint/srvol*phase
             db(ig+nlmto,ilm,k,:) = -srm1*qpg(:) * b(ig+nlmto,ilm,k)
          enddo
       enddo
    enddo 
1000 continue
!    write(6,*)'xxxxxxxxx111111122222222'
    call tcx('bstrux')
  end subroutine bstrux
end module m_bstrux
