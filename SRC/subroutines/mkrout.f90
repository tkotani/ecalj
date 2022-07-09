module m_mkrout
  use m_struc_def
  use m_uspecb,only:uspecb
  use m_lmfinit,only:nkapii,lhh
  use m_orbl,only: Orblib,ktab,ltab,offl,norb

  public:: m_Mkrout_init, orhoat_out, frcbandsym, qbyl_rv,hbyl_rv, sumec,sumtc,sumt0

  type(s_rv1),protected,allocatable :: orhoat_out(:,:)
  real(8),allocatable,protected::  frcbandsym(:,:), hbyl_rv(:,:,:), qbyl_rv(:,:,:)
  real(8),protected:: sumec,sumtc,sumt0

  private
contains

  subroutine m_mkrout_init()
    use m_lmfinit,only: sspec=>v_sspec,ispec,nsp,n0,ctrl_lfrce,lrout,nbas
    use m_bandcal,only: sv_p_oqkkl,sv_p_oeqkkl,frcband
    use m_mkpot,only: hab_rv,sab_rv
    use m_suham,only: ham_ndham
    integer:: ib,is,nr,lmxl,nlml
    !! contribution to force via band
    !     ! for force
    call tcn('m_mkrout_init')
    if(ctrl_lfrce>0 ) then
       if(allocated(frcbandsym)) deallocate(frcbandsym)
       allocate(frcbandsym(3,nbas)) !for forces
       frcbandsym=frcband         !part of atomic force from band part
    endif
    !!
    if( .NOT. allocated(orhoat_out)) allocate(orhoat_out(3,nbas))
    do  ib = 1, nbas
       is =  ispec(ib) !ssite(ib)%spec
       nr =  sspec(is)%nr
       lmxl= sspec(is)%lmxl
       nlml = (lmxl+1)**2
       if (lmxl > -1) then
          if(allocated(orhoat_out(1,ib)%v)) then
             deallocate(orhoat_out(1,ib)%v,orhoat_out(2,ib)%v,orhoat_out(3,ib)%v)
          endif
          allocate( orhoat_out(1,ib)%v(nr*nlml*nsp))
          allocate( orhoat_out(2,ib)%v(nr*nlml*nsp))
          allocate( orhoat_out(3,ib)%v(nr*nsp))
       endif
    enddo
    if( .NOT. allocated(qbyl_rv)) then
       allocate(qbyl_rv(n0,nsp,nbas))
       allocate(hbyl_rv(n0,nsp,nbas))
    endif
    call mkrout(sv_p_oqkkl,sv_p_oeqkkl, orhoat_out, hab_rv,sab_rv,qbyl_rv,hbyl_rv)
    if (lrout/=00) then
       !!  Symmetrize output atomic density and forces frcbandsym is symmetrized
       call symrhoat ( orhoat_out, qbyl_rv, hbyl_rv, frcbandsym)
    endif
    call tcx('m_mkrout_init')
  end subroutine m_mkrout_init
  !!-------------------------------
  subroutine mkrout( sv_p_oqkkl, sv_p_oeqkkl, orhoat_out, hab,sab, qbyl, hbyl)
    use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg,procid,master,nkaph &
         , sspec=>v_sspec,ispec,nbas,nsp,lekkl,lrout,n0,nab,nlmto,nmcore,rsma
    use m_lgunit,only: stdo
    use m_struc_def
    use m_elocp,only: rsmlss=>rsml, ehlss=>ehl
    use m_density,only: v0pot,v1pot,pnuall,pnzall !read
    !- Assembles local output densities out of the qkkl, and core states
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   ssite :struct for site-specific information; see routine usite
    !i     Elts read: spec pnu ov0 ov1 pz
    !i     Stored:    *
    !i     Passed to: *
    !i   sspec :struct for species-specific information; see routine uspec
    !i     Elts read: rsma lmxa lmxl kmxt a nr rmt lmxb stc orhoc
    !i     Stored:    *
    !i     Passed to: corprm uspecb gtpcor
    !i   slat  :struct for lattice information; see routine ulat
    !i     Elts read: ocg ojcg oidxcg ocy
    !i     Stored:    *
    !i     Passed to: *
    !i   sham  :struct for parameters defining hamiltonian; see routine uham
    !i     Elts read: oindxo eterms
    !i     Stored:    eterms (see Outputs)
    !i     Passed to: *
    !i   nbas  :size of basis
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlmto :dimension of lower block of LMTO part of hamiltonian
    !l   lekkl :1 Find band CG from eqkkl/qkkl
    !i   oqkkl :local density-matrix (rhocbl or comparable routine)
    !o   oeqkkl:local part of energy-weighted density matrix
    !i   hab   :hamiltonian matrix elements of radial wave functions at each
    !i         :site.  See potpus for their definition.
    !i   sab   :overlap matrix elements of radial wave functions at each
    !i         :site.  See potpus for their definition.
    !i         :hab and sab are used here to find band cg
    !i   lrout :0 calculate core states part only
    !o Outputs
    !o   orhoat:vector of offsets containing site density
    !o   sham->eterms various integrals for the total energy are stored:
    !o         :(8)  sumec = sum of foca=0 core eigenvalues
    !o         :(9)  sumtc = sum of all core kinetic energies
    !o         :(12) sumt0 = sum of frozen core kinetic energies
    !o   qbyl  :l-decomposed charge
    !o   hbyl  :l-decomposed eigenvalue sum
    !r Remarks
    !r   u and s are linear combinations of and phi,phidot defined as:
    !r   u has val=1, slo=1 at rmax, s has val=0, slo=1
    !u Updates
    !u   30 Jul 08 (T. Kotani) Use ekkl to set hbyl
    !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
    !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
    !u   28 Aug 01 Extended to local orbitals.
    !u   18 Jun 00 spin polarized
    !u   20 May 00 adapted from nfp mk_dens2
    ! ----------------------------------------------------------------------
    implicit none
    type(s_rv1) :: orhoat_out(3,nbas)
    type(s_rv1) :: sv_p_oeqkkl(3,nbas)
    type(s_rv1) :: sv_p_oqkkl(3,nbas)
    real(8):: qbyl(n0,nsp,nbas) , hbyl(n0,nsp,nbas) , sab(nab,n0,nsp,nbas) &
         , hab(nab,n0,nsp,nbas)
    integer:: ib , ipr , iprint , is , k , kcor , kmax , lcor , lfoc &
         , lgunit , lmxa , lmxh , lmxl , nlma , nlmh , nlml , nlml1 , &
         nr , ncore , igetss
    real(8) ,allocatable :: dmatl_rv(:)
    real(8) ,allocatable :: dh_rv(:)
    real(8) ,allocatable :: dp_rv(:)
    real(8) ,allocatable :: fh_rv(:)
    real(8) ,allocatable :: fp_rv(:)
    real(8) ,allocatable :: rofi_rv(:)
    real(8) ,allocatable :: rss_rv(:)
    real(8) ,allocatable :: rus_rv(:)
    real(8) ,allocatable :: ruu_rv(:)
    real(8) ,allocatable :: rwgt_rv(:)
    real(8) ,allocatable :: sl_rv(:)
    real(8) ,allocatable :: ul_rv(:)
    real(8) ,allocatable :: gz_rv(:)
    real(8) ,allocatable :: vh_rv(:)
    real(8) ,allocatable :: vp_rv(:)
    real(8) ,allocatable :: xh_rv(:)
    real(8) ,allocatable :: xp_rv(:)
    real(8) ,allocatable :: chh_rv(:)
    real(8) ,allocatable :: chp_rv(:)
    real(8) ,allocatable :: cpp_rv(:)
    double precision :: a,ceh,pi,qcor(2),rfoc,rmt,smec,smtc,stc0, &
         sum1,sum2,sums1,sums2,xx,y0,z,ddot,res,rsml(n0),ehl(n0)
    integer :: nkap0,nkapi,nkape,nglob
    parameter (nkap0=3)
    integer :: lh(nkap0)
    double precision :: eh(n0,nkap0),rsmh(n0,nkap0)
    !      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
    integer :: ntab(n0*nkap0),blks(n0*nkap0)
    double precision :: qcorg,qcorh,qsc,cofg,cofh !pnu(n0,2),pnz(n0,2),
    real(8),pointer:: pnu(:,:),pnz(:,:)
    integer ::iwdummy,ifx!,nlmto
    real(8):: dat(6,nbas)
    logical:: mmtargetx=.false.
    !      nlmto= ham_ldham(1)
    ipr  = iprint()
    pi   = 4d0*datan(1d0)
    y0   = 1d0/dsqrt(4d0*pi)
    sums1 = 0
    call tcn('mkrout')
    sumtc = 0d0
    sumec = 0d0
    sumt0 = 0d0
    if(procid==master) inquire(file='mmtarget.aftest',exist=mmtargetx)
    do  ib = 1, nbas
       is = ispec(ib) !int(ssite(ib)%spec)
       lmxa=sspec(is)%lmxa
       lmxl=sspec(is)%lmxl
       kmax=sspec(is)%kmxt
       if (lmxa == -1) goto 10
       a=sspec(is)%a
       nr=sspec(is)%nr
       rmt=sspec(is)%rmt
       lmxh=sspec(is)%lmxb
       stc0=sspec(is)%stc
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       call uspecb(is,rsmh,eh)
       nkapi=nkapii(is)
       call orblib(ib)!norb , ltab , ktab , xx , offl , xx )
       call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)
       is =ispec(ib) !ssite(ib)%spec
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       call gtpcor(is,kcor,lcor,qcor)
       nlml = (lmxl+1)**2
       nlma = (lmxa+1)**2
       nlmh = (lmxh+1)**2
       allocate(rofi_rv(nr))
       allocate(rwgt_rv(nr))
       call radmsh( rmt , a , nr , rofi_rv )
       call radwgt( rmt , a , nr , rwgt_rv )
       if (lrout /= 0) then
          orhoat_out(1, ib )%v=0d0
          orhoat_out(2, ib )%v=0d0
          !   --- Assemble rho1 = products of augmented functions ---
          !   ... Set up all radial head and tail envelope functions, and their bc's
          allocate(fh_rv(nr*(lmxh+1)*nkaph))
          allocate(xh_rv(nr*(lmxh+1)*nkaph))
          allocate(vh_rv((lmxh+1)*nkaph))
          allocate(dh_rv((lmxh+1)*nkaph))
          call fradhd(nkaph, eh, rsmh, lhh(:,is), lmxh , nr, rofi_rv, fh_rv , xh_rv , vh_rv , dh_rv )
          allocate(fp_rv(nr*(lmxa+1)*(kmax+1)))
          allocate(xp_rv(nr*(lmxa+1)*(kmax+1)))
          allocate(vp_rv((lmxa+1)*(kmax+1)))
          allocate(dp_rv((lmxa+1)*(kmax+1)))
          call fradpk(kmax, rsma(is), lmxa , nr , rofi_rv , fp_rv, xp_rv , vp_rv , dp_rv )
          !   ... Augmented wave functions
          allocate(ul_rv(nr*(lmxa+1)*nsp))
          allocate(sl_rv(nr*(lmxa+1)*nsp))
          allocate(gz_rv(nr*(lmxa+1)*nsp))
          allocate(ruu_rv(nr*(lmxa+1)*2*nsp))
          allocate(rus_rv(nr*(lmxa+1)*2*nsp))
          allocate(rss_rv(nr*(lmxa+1)*2*nsp))
          !     call uspecb(0,4,sspec,is,is,lh,rsml,ehl,k)
          rsml= rsmlss(:,is)
          ehl=  ehlss(:,is)
          call makusp(n0 , z , nsp , rmt , lmxa , v0pot(ib)%v, a , nr , &
               xx , xx , pnu , pnz , rsml , ehl , ul_rv , sl_rv , gz_rv , ruu_rv &
               , rus_rv , rss_rv )
          !   ... Contracted density matrix as coffs to products of (u,s,gz)
          k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml*nsp
          allocate(chh_rv(k))
          chh_rv=0d0
          allocate(chp_rv(k))
          chp_rv=0d0
          allocate(cpp_rv(k))
          cpp_rv=0d0
          allocate(dmatl_rv((lmxa+1)**2*nlml*nsp*9))
          dmatl_rv=0d0
          call mkrou1(nsp, nlmh , nlma , nlml , kmax , rv_a_ocg , iv_a_ojcg &
               , iv_a_oidxcg , nkaph , nkapi , norb , ltab , ktab , blks &
               , sv_p_oqkkl(3,ib)%v, sv_p_oqkkl(2,ib)%v, sv_p_oqkkl(1,ib)%v& !, OQHH, OQHP , OQPP ,
               , vh_rv , dh_rv , vp_rv , dp_rv , chh_rv &
               , chp_rv , cpp_rv , dmatl_rv )
          !   ... True local density for this sphere 1st component
          call mkrou2( nsp , lmxa , nlml , pnz , dmatl_rv , nr , ul_rv &
               , sl_rv , gz_rv , ruu_rv , rus_rv , rss_rv , orhoat_out( 1 , ib )%v )
          !   --- Assemble rho2 = unaugmented products Pkl*Pk'l' ---
          !       H H product.  Include local orbitals w/ envelopes (nkapi->nkaph)
          call mkrou5( nsp , nr , nlml , nkaph , nkaph , fh_rv , lmxh &
               , nkaph , nkaph , fh_rv , lmxh , chh_rv , orhoat_out( 2 , ib )%v  )
          !       H Pkl product
          call mkrou5( nsp , nr , nlml , nkaph , nkaph , fh_rv , lmxh &
               , kmax + 1 , kmax + 1 , fp_rv , lmxa , chp_rv , orhoat_out( 2 , ib )%v   )
          !       Pkl Pkl product
          call mkrou5( nsp , nr , nlml , kmax + 1 , kmax + 1 , fp_rv , &
               lmxa, kmax+1, kmax + 1 , fp_rv , lmxa , cpp_rv , orhoat_out( 2 , ib )%v   )
          !   ... Site charge and eigenvalue sum decomposed by l
          !! As I describe below, hbyl give here is wrong (bug). maybe because of wrong hab.
          call mkrou3(2, lmxa , nlml , nsp , pnz , dmatl_rv , hab &
               (1, 1,1,ib), sab(1, 1, 1, ib ) , qbyl (1 , 1 ,ib ) , hbyl ( 1 , 1 , ib ) )
          ! Remake hbyl from energy-weighted density matrix
          ! This is because above codes to give hbyl is wrong. So, we need to go though anothe pass with lekkl=1
          !! (OPTION_PFLOAT=1). In future, code above to generate hbyl will be removed.
          ! ino Jan.04.2012:           rv_p_oqhh => sv_p_oeqkkl(3,ib)%v
          if(lekkl /= 0) then
             nlml1 = 1
             call dpzero( dmatl_rv , ( lmxa + 1) ** 2 * nlml1 * nsp * 9 )
             k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml1*nsp
             call dpzero(chh_rv , k )
             call dpzero(chp_rv , k )
             call dpzero(cpp_rv , k )
             call mkrou1(nsp, nlmh, nlma, nlml1 , kmax , rv_a_ocg , iv_a_ojcg &
                  , iv_a_oidxcg , nkaph , nkapi , norb , ltab , ktab , blks &
                  , sv_p_oeqkkl(3,ib)%v, sv_p_oeqkkl(2,ib)%v,sv_p_oeqkkl(1,ib)%v& ! & , OQHH , OQHP , OQPP,
                  , vh_rv , dh_rv , vp_rv , dp_rv , chh_rv &
                  , chp_rv , cpp_rv , dmatl_rv ) !dmatl_rv is energy weighted. c.f. previous call to mkrou1.
             call mkrou3( mode=1 , lmxa=lmxa, nlml=nlml1 , nsp=nsp &
                  , pnz= pnz, dmatl=dmatl_rv, sab=sab( 1 , 1 , 1 , ib ), qsum=hbyl(1,1,ib )  )
          endif
          !   --- Print charges for information ---
          call radsum ( nr , nr , nlml , nsp , rwgt_rv , orhoat_out( 1 , ib )%v, sum1 )
          call radsum ( nr , nr , nlml , nsp , rwgt_rv , orhoat_out( 2 , ib )%v, sum2 )
          sum1 = sum1/y0
          sum2 = sum2/y0
          if (nsp == 2) then
             sums1 = sum1 - 2d0 * ddot ( nr , rwgt_rv , 1 , orhoat_out( 1 , ib )%v, 1 ) / y0
             sums2 = sum2 - 2d0 * ddot ( nr , rwgt_rv , 1 , orhoat_out( 2 , ib )%v, 1 ) / y0
          endif
          !          if (ipr.ge.30 .and. nsp .eq. 1) write(stdo,200) ib,sum1,sum2,sum1-sum2
          !          if (ipr.ge.30 .and. nsp .eq. 2)
          !     .      write(stdo,200) ib,sum1,sum2,sum1-sum2, -sums1,-sums2,-sums1+sums2
          ! 200      format(i4,3f12.6,2x,3f12.6)

          if(nsp==1) dat(1:3,ib) = [sum1,sum2,sum1-sum2]
          if(nsp==2) dat(1:6,ib) = [sum1,sum2,sum1-sum2, -sums1,-sums2,-sums1+sums2]
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !! experimental block to keep magnetic moment for AF.controlled by uhval.aftest file. See elsewhere.
          if( nsp==2 .AND. procid==master .AND. mmtargetx) then
             if(ib==1) then
                open(newunit=ifx,file='mmagfield.aftest',position='append')
                write(ifx,"('          ',f23.15)",advance='no') -sums1
             else
                write(ifx,"(2x,f23.15)",advance='no') -sums1
             endif
             if(ib==nbas) write(ifx,*)
             if(ib==nbas) close(ifx)
          endif
          ! ccccccccccccccccccccccccccccccccccccccccccccccccc
       endif
       !       Contribution to mag outsite rmt: extrapolate tail to infinity
       call mkrou6 ( rofi_rv , orhoat_out( 1 , ib )%v , nr , nlml , &
            nsp , xx , xx , res )
       if (ipr >= 30 .AND. res /= 0) then
          write(stdo,211) res,res-sums1
211       format(7x,'contr. to mm extrapolated for r>rmt:',f11.6, &
               ' est. true mm =',f9.6)
       endif
       !   --- Make new core density and core eigenvalue sum ---
       if (lfoc == 0) then
          call pshpr(ipr+11)
          !write(stdo,*)
          call getcor ( 0 , z , a , pnu , pnz , nr , lmxa , rofi_rv , v1pot(ib)%v & !ssite(ib)%rv_a_ov1 &
               , kcor , lcor , qcor , smec , smtc , orhoat_out( 3 , ib )%v &
               , ncore , 0d0 , 0d0,   nmcore(is))
          call poppr
          sumtc = sumtc + smtc
          sumec = sumec + smec
       else
          !          if (ipr .ge. 81) write(stdo,288) stc0
          !  288     format(' foca..  use smtc = ',f14.6)
          sumtc = sumtc + stc0
          sumt0 = sumt0 + stc0
          call dpcopy ( sspec(is)%rv_a_orhoc , orhoat_out( 3 , ib )%v , 1 , nr* nsp , 1d0 )
       endif
       !     !
       if (allocated(dmatl_rv)) deallocate(dmatl_rv)
       if (allocated(cpp_rv)) deallocate(cpp_rv)
       if (allocated(chp_rv)) deallocate(chp_rv)
       if (allocated(chh_rv)) deallocate(chh_rv)
       if (allocated(rss_rv)) deallocate(rss_rv)
       if (allocated(rus_rv)) deallocate(rus_rv)
       if (allocated(ruu_rv)) deallocate(ruu_rv)
       if (allocated(gz_rv)) deallocate(gz_rv)
       if (allocated(sl_rv)) deallocate(sl_rv)
       if (allocated(ul_rv)) deallocate(ul_rv)
       if (allocated(dp_rv)) deallocate(dp_rv)
       if (allocated(vp_rv)) deallocate(vp_rv)
       if (allocated(xp_rv)) deallocate(xp_rv)
       if (allocated(fp_rv)) deallocate(fp_rv)
       if (allocated(dh_rv)) deallocate(dh_rv)
       if (allocated(vh_rv)) deallocate(vh_rv)
       if (allocated(xh_rv)) deallocate(xh_rv)
       if (allocated(fh_rv)) deallocate(fh_rv)
       if (allocated(rwgt_rv)) deallocate(rwgt_rv)
       if (allocated(rofi_rv)) deallocate(rofi_rv)
10     continue
    enddo
    if (ipr >= 30 .AND. lrout > 0) then
       ! write(stdo,"(a)")' mkrout: site(class) decomposed charge and magnetic moment. class->lmchk'
       if (nsp == 1) write(stdo,201)
       if (nsp == 2) write(stdo,202)
201    format(/' mkrout:  Qtrue      sm,loc       local')
202    format(/' mkrout:  Qtrue      sm,loc       local',8x,'true mm   smooth mm    local mm')
       do ib=1,nbas
          write(stdo,200) ib,dat(1:3*nsp,ib)
200       format(i4,3f12.6,2x,3f12.6)
       enddo
    endif
    ! --- Put sumec,sumtc into etot struct ---
    !      ham_eterms=eterms
    call tcx('mkrout')
  end subroutine mkrout

  subroutine mkrou1(nsp,nlmh,nlma,nlml,kmax,cg,jcg,indxcg, &
       nkaph,nkapi,norb,ltab,ktab,blks,qhh,qhp,qpp,vh,dh,vp,dp, &
       chh,chp,cpp,dmatl)
    !- Contracted density matrix for one site
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlmh  :(lmxh+1)**2, where lmxh = L-cutoff in basis
    !i   nlma  :(lmxa+1)**2, where lmxa = L-cutoff in augmentation
    !i   nlml  :L-cutoff for charge density on radial mesh
    !i   kmax  :polynomial cutoff in augmentation
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   nkaph :number of types of one l-quantum number in the MTO basis
    !i   nkapi :number of the nkaph functions that are envelope functions
    !i   norb  :number of orbital types for this site; see orbl.f
    !i   ltab  :table of l-quantum numbers for each type
    !i   ktab  :table of energy index for each type
    !i   blks  :blks(iorb) = size of contiguous block of orbitals for
    !i         :orbital type iorb and possibly iorb+1,iorb+2...
    !i   qhh   :head-head component of site part of density matrix (addrbl.f)
    !i   qhp   :head-tail component of site part of density matrix (addrbl.f)
    !i   qpp   :tail-tail component of site part of density matrix (addrbl.f)
    !i   vh    :values of head functions on MT sphere
    !i   dh    :slopes of tail functions on MT sphere
    !i   vp    :values of Pkl  functions on MT sphere
    !i   dp    :slopes of Pkl  functions on MT sphere
    !o Outputs
    !o   chh   :head-head product function coefficients
    !o   chp   :head-tail product function coefficients
    !o   cpp   :tail-tail product function coefficients
    !o   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a Y_lm
    !o         :expansion of the function products f_i(l1) f_j(l2)
    !o         :where f_i,f_j form this table.
    !o         :    (uu  us  uz)
    !o         :    (su  ss  sz)
    !o         :    (zu  zs  zz)
    !r Remarks
    !r   Transforms density-matrix as generated by addrbl.f into
    !r   contracted density-matrix of wave function products.
    !r   Radial parts of chh,chp,cpp are
    !r     (psi1)_l (psi2)_l'
    !r   Radial parts of dmatl are:
    !r     (u, s, or gz)_l (u, s, or gz)_l'
    !u Updates
    !u   29 Jul 08 Adapted from mkrout.f
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: nsp,nlmh,nlma,nlml,kmax,norb,nkaph,nkapi
    integer :: jcg(1),indxcg(1)
    integer :: ltab(norb),ktab(norb),blks(norb)
    double precision :: qhh(nkaph,nkaph,nlmh,nlmh,nsp)
    double precision :: qhp(nkaph,kmax+1,nlmh,nlma,nsp)
    double precision :: qpp(kmax+1,kmax+1,nlmh,nlma,nsp)
    !     double precision dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp)
    !     double precision vh(0:lmxh,nkaph),dh(0:lmxh,nkaph),
    !    .                 vp(0:lmxa,kmxa+1),dp(0:lmxa,kmxa+1)
    !     chh = chh(nkaph,nkaph,0:lmxh,0:lmxh,nlml,nsp)
    !     chp = chp(nkaph,kmax+1,0:lmxh,0:lmxa,nlml,nsp)
    !     cpp = cpp(kmax+1,kmax+1,0:lmxa,0:lmxa,nlml,nsp)
    double precision :: dmatl(*),vh(*),dh(*),vp(*),dp(*), chh(*),chp(*),cpp(*)
    double precision :: cg(1)
    integer :: ll,k,lmxa,lmxh
    integer :: n0,nkap0
    parameter (n0=10,nkap0=3)
    integer :: ltba(n0*nkap0),ktba(n0*nkap0),blka(n0*nkap0)
    lmxa = ll(nlma)
    lmxh = ll(nlmh)
    do  k = 0, kmax
       ltba(k+1) = 0
       ktba(k+1) = k+1
       blka(k+1) = nlma
    enddo
    ! ... Contracted density-matrix chh,chp,cpp from qhh,qhp,qpp
    !     H H product
    call mkrou4(nsp,nlml,cg,jcg,indxcg, &
         nkaph,norb,ltab,ktab,blks,lmxh,nlmh, &
         nkaph,norb,ltab,ktab,blks,lmxh,nlmh, &
         qhh,chh)
    !     H Pkl product
    call mkrou4(nsp,nlml,cg,jcg,indxcg, &
         nkaph,norb,ltab,ktab,blks,lmxh,nlmh, &
         kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma, &
         qhp,chp)
    !     Pkl Pkl product
    call mkrou4(nsp,nlml,cg,jcg,indxcg, &
         kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma, &
         kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma, &
         qpp,cpp)

    ! ... Contracted density matrix as coffs to products of (u,s,gz)
    !     H H product
    call mkcfus(nsp,lmxa,nlml, &
         nkaph,nkapi,vh,dh,lmxh, &
         nkaph,nkapi,vh,dh,lmxh, &
         chh,dmatl)
    !     H Pkl product
    call mkcfus(nsp,lmxa,nlml, &
         nkaph,nkapi,vh,dh,lmxh, &
         kmax+1,kmax+1,vp,dp,lmxa, &
         chp,dmatl)
    !     Pkl Pkl product
    call mkcfus(nsp,lmxa,nlml, &
         kmax+1,kmax+1,vp,dp,lmxa, &
         kmax+1,kmax+1,vp,dp,lmxa, &
         cpp,dmatl)
  end subroutine mkrou1

  subroutine mkrou2(nsp,lmxa,nlml,pnz,dmatl,nr,ul,sl,gz,ruu,rus,rss, rho)
    !- Assemble true site density from product function coefficients
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmxa  :augmentation l-cutoff
    !i   nlml  :l-cutoff for charge density
    !i   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a y_lm
    !i         :expansion of the function products f_i(l1) f_j(l2)
    !i         :where f_i,f_j form this table.
    !i         :    (uu  us  uz)
    !i         :    (su  ss  sz)
    !i         :    (zu  zs  zz)
    !i   nr    :number of radial mesh points
    !i   ul    :r*radial wave functions; see remarks
    !i   sl    :r*radial wave functions; see remarks
    !i   gz    :r*semicore wave functions; see remarks
    !i   ruu   :l-diagonal (uu) product including small component (makusp)
    !i   rus   :l-diagonal (us) product including small component (makusp)
    !i   rss   :l-diagonal (ss) product including small component (makusp)
    !o Outputs
    !o   rho   :charge density assembled
    !r Remarks
    !r   This routine uses linear combinations (u,s) of phi,phidot
    !r   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
    !r   ul and sl are returned as r * u and r * s, respectively.
    !r
    !r   Let phi_z be the w.f. corresponding to pnu_z.
    !r   A local orbital of the first type is defined as follows.
    !r      gz = r * ( phi_z - phi_z(rmax) u - (phi_z)'(rmax) s )
    !r   By construction, gz/r has both value = 0 and slope = 0 at rmax.
    !r   A local orbital of the second type is defined as gz=r*phi_z;
    !r   for r>rmax a smooth Hankel tail (spec'd by ehl,rsml) is attached.
    !r
    !r   Spherical density assumbled using ruu,rus,rss which include
    !r   small comonent of relativistic wave function.
    !u Updates
    !u   28 aug 01 extended to local orbitals.  altered argument list.
    ! ----------------------------------------------------------------------
    implicit none
    integer :: nsp,lmxa,nlml,nr,n0
    parameter (n0=10)
    double precision :: rho(nr,nlml,nsp),pnz(n0,nsp), &
         dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp), &
         ul(nr,0:lmxa,nsp),ruu(nr,0:lmxa,2,nsp), &
         sl(nr,0:lmxa,nsp),rus(nr,0:lmxa,2,nsp), &
         gz(nr,0:lmxa,nsp),rss(nr,0:lmxa,2,nsp)
    ! ... local parameters
    logical :: lpz1,lpz2
    integer :: l,i,mlm,l1,l2,isp
    double precision :: xuu,xus,xsu,xss,xuz,xsz,xzu,xzs,xzz
    !     double precision top
    call tcn('mkrou2')
    ! ccccccccc
    !      print *,'ddd ', sum(abs(dmatl(:,:,:,1:2,1:2,:)))
    !      print *,'ddd ', sum(abs(ul))
    !      print *,'ddd ruu ', sum(abs(ruu))
    !      print *,'ddd ', sum(abs(rus))
    !      print *,'ddd ', sum(abs(rss))
    !      print *,'ddd ul ', sum(abs(ul))
    !      print *,'ddd ', sum(abs(sl))

    ! --- Full density as products (u,s) (u,s); no small component ---
    do  isp = 1, nsp
       do  mlm = 1, nlml
          do  l1 = 0, lmxa
             do  l2 = 0, lmxa
                lpz1 = pnz(l1+1,1) .ne. 0
                lpz2 = pnz(l2+1,1) .ne. 0

                xuu = dmatl(l1,l2,mlm,1,1,isp)
                xus = dmatl(l1,l2,mlm,1,2,isp)
                xsu = dmatl(l1,l2,mlm,2,1,isp)
                xss = dmatl(l1,l2,mlm,2,2,isp)

                !           top = xuu*xuu + xus*xus + xsu*xsu + xss*xss
                !           if (dsqrt(top).gt.1d-6)
                !    .        write (6,700) mlm,l1,l2,xuu,xus,xsu,xss
                ! 700       format(3i5,4f14.8)

                do  i = 1, nr
                   rho(i,mlm,isp) = rho(i,mlm,isp) &
                        + xuu * ul(i,l1,isp) * ul(i,l2,isp) &
                        + xsu * sl(i,l1,isp) * ul(i,l2,isp) &
                        + xus * ul(i,l1,isp) * sl(i,l2,isp) &
                        + xss * sl(i,l1,isp) * sl(i,l2,isp)
                enddo

                if (lpz1 .OR. lpz2) then

                   xuz = dmatl(l1,l2,mlm,1,3,isp)
                   xsz = dmatl(l1,l2,mlm,2,3,isp)
                   xzu = dmatl(l1,l2,mlm,3,1,isp)
                   xzs = dmatl(l1,l2,mlm,3,2,isp)
                   xzz = dmatl(l1,l2,mlm,3,3,isp)
                   if (xuz /= 0 .OR. xsz /= 0 .OR. xzu /= 0 .OR. &
                        xzs /= 0 .OR. xzz /= 0) then
                      do  i = 1, nr
                         rho(i,mlm,isp) = rho(i,mlm,isp) &
                              + xuz * ul(i,l1,isp) * gz(i,l2,isp) &
                              + xsz * sl(i,l1,isp) * gz(i,l2,isp) &
                              + xzu * gz(i,l1,isp) * ul(i,l2,isp) &
                              + xzs * gz(i,l1,isp) * sl(i,l2,isp) &
                              + xzz * gz(i,l1,isp) * gz(i,l2,isp)
                      enddo
                   endif
                endif

             enddo
          enddo
       enddo
    enddo

    ! --- Remake spherical density including small component ---
    !      print *, 'skip remaking spher. density'
    !      return
    call dpzero(rho,nr)
    call dpzero(rho(1,1,nsp),nr)
    do  isp = 1, nsp
       do  l = 0, lmxa
          xuu = dmatl(l,l,1,1,1,isp)
          xus = dmatl(l,l,1,1,2,isp) + dmatl(l,l,1,2,1,isp)
          xss = dmatl(l,l,1,2,2,isp)
          xuz = dmatl(l,l,1,1,3,isp) + dmatl(l,l,1,3,1,isp)
          xsz = dmatl(l,l,1,2,3,isp) + dmatl(l,l,1,3,2,isp)
          xzz = dmatl(l,l,1,3,3,isp)

          do  i = 1, nr
             rho(i,1,isp) = rho(i,1,isp) &
                  +xuu*ruu(i,l,1,isp)+xus*rus(i,l,1,isp)+xss*rss(i,l,1,isp)
          enddo

          if (pnz(l+1,1) /= 0) then
             do  i = 1, nr
                rho(i,1,isp) = rho(i,1,isp) &
                     +xuz*ruu(i,l,2,isp)+xsz*rus(i,l,2,isp)+xzz*rss(i,l,2,isp)
             enddo
          endif
       enddo
    enddo
    call tcx('mkrou2')
  end subroutine mkrou2

  subroutine mkrou4(nsp,nlml,cg,jcg,indxcg, &
       nk1,norb1,ltab1,ktab1,blks1,lmx1,ndim1, &
       nk2,norb2,ltab2,ktab2,blks2,lmx2,ndim2, &
       qkk12,ckk)

    !- Assemble contracted density-matrix (in Y_lm form)
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlml  :charge density L-cutoff
    !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
    !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
    !i   indxcg:index for Clebsch Gordon coefficients
    !i   nk1   :no. orbital types for a given L of first function
    !i         :nk1 merely dimensions qkk12,ckk
    !i   norb1 :total number of orbitals of the first type
    !i   ltab1 :table of l quantum numbers for the norb1 orbitals
    !i   ktab1 :table of k numbers (orbital type) for the norb1 orbitals
    !i   blks1 :block size for grouping orbitals into blocks (gtbls1)
    !i   lmx1  :dimensions ckk
    !i   ndim1 :dimensions qkk12; ndim1 should be (lmx1+1)**2
    !i   nk2   :no. orbital types for a given L of second function
    !i         :nk2 merely dimensions qkk12,ckk
    !i   norb2 :total number of orbitals of the first type
    !i   ltab2 :table of l quantum numbers for the norb1 orbitals
    !i   ktab2 :table of k numbers (orbital type) for the norb1 orbitals
    !i   blks2 :block size for grouping orbitals into blocks (gtbls1)
    !i   lmx2  :dimensions ckk
    !i   ndim2 :dimensions qkk12; ndim1 should be (lmx1+1)**2
    !i   qkk12 :density matrix between norb1 and norb2 orbitals
    !o Outputs
    !o   ckk   :product function coefficients
    !r Remarks
    !r   Wave function products have a radial part, written as products
    !r   of the radial functions (psi1)_l (psi2)_l'
    !r   The angular part Y_L Y_L' is contracted into a single index M
    !r   using Clebsch Gordan coefficients.  Thus ckk has indices
    !r      ckk(k1,k2,l1,l2,M),
    !r   where k1 ranges over the orbital types of the first function
    !r   and   k2 ranges over the orbital types of the second function.
    !u Updates
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: nsp,nlml,norb1,lmx1,ndim1,norb2,lmx2,ndim2,nk1,nk2
    integer :: jcg(1),indxcg(1)
    integer :: ltab1(norb1),ktab1(norb1),blks1(norb1), &
         ltab2(norb2),ktab2(norb2),blks2(norb2)
    double precision :: qkk12(nk1,nk2,ndim1,ndim2,nsp),cg(1), &
         ckk(nk1,nk2,0:lmx1,0:lmx2,nlml,nsp)
    ! ... Local parameters
    integer :: ilm1,io1,l1,nlm11,nlm12,k1, &
         ilm2,io2,l2,nlm21,nlm22,k2, &
         icg,ll,mlm,ix,isp
    double precision :: xx

    !     call tcn('mkrou4')

    do  isp = 1, nsp
       do  io2 = 1, norb2
          if (blks2(io2) /= 0) then
             !       k2,l2 = k and starting l index for this block
             l2 = ltab2(io2)
             k2 = ktab2(io2)
             nlm21 = l2**2+1
             nlm22 = nlm21 + blks2(io2)-1
             do  ilm2 = nlm21, nlm22
                l2 = ll(ilm2)
                do  io1 = 1, norb1
                   if (blks1(io1) /= 0) then
                      !           k1,l1 = k and starting l index for this block
                      l1 = ltab1(io1)
                      k1 = ktab1(io1)
                      nlm11 = l1**2+1
                      nlm12 = nlm11 + blks1(io1)-1
                      do  ilm1 = nlm11, nlm12
                         l1 = ll(ilm1)

                         ix = max0(ilm1,ilm2)
                         ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
                         do icg = indxcg(ix),indxcg(ix+1)-1
                            mlm = jcg(icg)
                            if (mlm <= nlml) then
                               xx = cg(icg)*qkk12(k1,k2,ilm1,ilm2,isp)
                               ckk(k1,k2,l1,l2,mlm,isp) = ckk(k1,k2,l1,l2,mlm,isp)+xx
                            endif
                         enddo
                      enddo
                   endif
                enddo
             enddo
          endif
       enddo
    enddo

    !     call tcx('mkrou4')
  end subroutine mkrou4


  subroutine mkcfus(nsp,lmxa,nlml,nf1,nf1s,val1,slo1,lmx1, &
       nf2,nf2s,val2,slo2,lmx2,ckk,dmatl)

    !- Assemble contracted density matrix as coffs to products of (u,s,gz)
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmxa  :leading dimension of dmatl
    !i   nlml  :charge density L-cutoff
    !i   nf1   :number of orbital types of first kind for each l
    !i   nf1s  :number of functions of first kind for each l, for which
    !i         :there is a smooth part to be subtracted (which also
    !i         :corresponds to the functions which connect to envelope
    !i         :functions)
    !i   val1  :function values at MT boundary, first function
    !i   slo1  :function slopes at MT boundary, first function
    !i   lmx1  :dimensions val1,slo1
    !i   nf2   :number of orbital types of second kind for each l
    !i   nf2s  :number of functions of second kind for each l, for which
    !i         :a smooth part is to be subtracted (which also
    !i         :corresponds to the functions which connect to envelope
    !i         :functions)
    !i   val2  :function values at MT boundary, second function
    !i   slo2  :function slopes at MT boundary, second function
    !i   lmx2  :dimensions val1,slo1
    !i   ckk   :density matrix between  and  orbitals
    !o Outputs
    !o   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a Y_lm
    !o         :expansion of the function products f_i(l1) f_j(l2)
    !o         :where f_i,f_j form this table.
    !o         :    (uu  us  uz)
    !o         :    (su  ss  sz)
    !o         :    (zu  zs  zz)
    !r Remarks
    !r   Transforms contracted density-function matrix into coefficients
    !r   of wave function products. Radial parts are
    !r     (u, s, or gz)_l (u, s, or gz)_l'
    !u Updates
    !u   28 Aug 01 Extended to local orbitals.  Altered argument list.
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: nsp,nlml,lmxa,lmx1,lmx2,nf1,nf2,nf1s,nf2s
    double precision :: ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nsp), &
         val1(0:lmx1,nf1s),slo1(0:lmx1,nf1s), &
         val2(0:lmx2,nf2s),slo2(0:lmx2,nf2s)
    double precision :: dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp)
    ! ... Local parameters
    integer :: l1,k1,l2,k2,mlm,isp
    double precision :: xx

    !     call tcn('mkcfus')

    do  isp = 1, nsp
       do  mlm = 1, nlml
          do  k2 = 1, nf2s
             do  k1 = 1, nf1s
                do  l2 = 0, lmx2
                   do  l1 = 0, lmx1
                      xx = ckk(k1,k2,l1,l2,mlm,isp)

                      dmatl(l1,l2,mlm,1,1,isp) = dmatl(l1,l2,mlm,1,1,isp) &
                           + xx * val1(l1,k1) * val2(l2,k2)
                      dmatl(l1,l2,mlm,1,2,isp) = dmatl(l1,l2,mlm,1,2,isp) &
                           + xx * val1(l1,k1) * slo2(l2,k2)
                      dmatl(l1,l2,mlm,2,1,isp) = dmatl(l1,l2,mlm,2,1,isp) &
                           + xx * slo1(l1,k1) * val2(l2,k2)
                      dmatl(l1,l2,mlm,2,2,isp) = dmatl(l1,l2,mlm,2,2,isp) &
                           + xx * slo1(l1,k1) * slo2(l2,k2)

                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    !     call tcx('mkcfus')

    ! --- Products involving local orbitals ---
    if (nf1s >= nf1 .AND. nf2s >= nf2) return
    !     call tcn('mkcfus')

    do  isp = 1, nsp
       do  mlm = 1, nlml
          do  k2 = 1, nf2
             do  k1 = 1, nf1
                if (k1 > nf1s .OR. k2 > nf2s) then
                   do  l2 = 0, lmx2
                      do  l1 = 0, lmx1
                         xx = ckk(k1,k2,l1,l2,mlm,isp)

                         !               sc-sc product
                         if (k1 > nf1s .AND. k2 > nf2s) then
                            dmatl(l1,l2,mlm,3,3,isp) = dmatl(l1,l2,mlm,3,3,isp)+xx

                            !               sc-valence product
                         elseif (k1 > nf1s) then
                            dmatl(l1,l2,mlm,3,1,isp) = dmatl(l1,l2,mlm,3,1,isp) &
                                 + xx * val2(l2,k2)
                            dmatl(l1,l2,mlm,3,2,isp) = dmatl(l1,l2,mlm,3,2,isp) &
                                 + xx * slo2(l2,k2)

                            !               valence-sc product
                         elseif (k2 > nf2s) then
                            dmatl(l1,l2,mlm,1,3,isp) = dmatl(l1,l2,mlm,1,3,isp) &
                                 + val1(l1,k1) * xx
                            dmatl(l1,l2,mlm,2,3,isp) = dmatl(l1,l2,mlm,2,3,isp) &
                                 + slo1(l1,k1) * xx

                         endif

                         !                if (l1.eq.2.and.l2.eq.2 .and. mlm.eq.1) then
                         !                  print 333, k1,k2,xx,
                         !     .              dmatl(l1,l2,mlm,1,3,isp),dmatl(l1,l2,mlm,3,1,isp)
                         !  333             format(2i4,3f12.6)
                         !                endif

                      enddo
                   enddo
                endif
             enddo
          enddo
       enddo
    enddo

    !     call tcx('mkcfus')

    !      isp = 1
    !      do  mlm = 1, nlml
    !        print *, ' '
    !        do  l2 = 0, lmx2
    !          do  l1 = 0, lmx1
    !            print 334, l1,l2,mlm,
    !     .        dmatl(l1,l2,mlm,1,2,isp),
    !     .        dmatl(l2,l1,mlm,2,1,isp),
    !     .        dmatl(l1,l2,mlm,1,2,isp)-dmatl(l2,l1,mlm,2,1,isp)
    !            print 334, l1,l2,mlm,
    !     .        dmatl(l1,l2,mlm,1,3,isp),
    !     .        dmatl(l2,l1,mlm,3,1,isp),
    !     .        dmatl(l1,l2,mlm,1,3,isp)-dmatl(l2,l1,mlm,3,1,isp)
    !            print 334, l1,l2,mlm,
    !     .        dmatl(l1,l2,mlm,2,3,isp),
    !     .        dmatl(l2,l1,mlm,3,2,isp),
    !     .        dmatl(l1,l2,mlm,2,3,isp)-dmatl(l2,l1,mlm,3,2,isp)

    !  334       format(3i4,2f12.6,f14.8)
    !          enddo
    !        enddo
    !      enddo
    !      pause

  end subroutine mkcfus


  subroutine mkrou5(nsp,nr,nlml,nf1,nf1s,f1,lmx1,nf2,nf2s,f2,lmx2, &
       ckk,rho)

    !- Assemble smooth site density from contracted density matrix
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlml  :charge density L-cutoff
    !i   nf1   :number of orbital types of first kind for each l
    !i         :nf1 merely dimensions ckk
    !i   nf1s  :number of functions of first kind for each l, for which
    !i         :there are functions of f1 type defined
    !i   f1    :first function on a radial mesh
    !i   lmx1  :l-cutoff for f1 functions
    !i   nf2   :number of orbital types of second kind for each l
    !i         :nf2 merely dimensions ckk
    !i   nf2s  :number of functions of second kind for each l, for which
    !i         :there are functions of f2 type defined
    !i   f2    :second function on a radial mesh
    !i   lmx2  :l-cutoff for f2 functions
    !i   ckk   :density matrix between f1 and f2 orbitals
    !o Outputs
    !o   rho   :density assembled on mesh
    !r Remarks
    !u Updates
    !u   28 Aug 01 Extended to local orbitals.  Altered argument list.
    !u             Envelopes f1,f2 must be zero for all channels that
    !u             have no smooth counterparts to subtract.
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: nsp,nr,nlml,lmx1,lmx2,nf1,nf2,nf1s,nf2s
    double precision :: ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nsp), &
         f1(nr,0:lmx1,nf1s),f2(nr,0:lmx2,nf2s),rho(nr,nlml,nsp)
    ! ... Local parameters
    integer :: l1,k1,l2,k2,mlm,isp,i
    double precision :: xx

    call tcn('mkrou5')
!    rho(:,:,:) =rho(:,:,:)+ ckk(k1,k2,l1,l2,:,:)*[(f1(ir,l1,k1)*f2(:,l2,k2),ir=1,nr)]
    do  isp = 1, nsp
       do  mlm = 1, nlml
          do  k2 = 1, nf2s
             do  k1 = 1, nf1s
                do  l2 = 0, lmx2
                   do  l1 = 0, lmx1
                      rho(:,mlm,isp)=rho(:,mlm,isp) +ckk(k1,k2,l1,l2,mlm,isp)*f1(:,l1,k1)*f2(:,l2,k2)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    call tcx('mkrou5')
  end subroutine mkrou5


  subroutine mkrou6(rofi,rho,nr,nlml,nsp,rho0,decay,res)

    !- Fit tail of spin density; integrate charge beyond MT sphere
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   rofi  :radial mesh points
    !i   rho   :spin-polarized charge density
    !i   nr    :number of radial mesh points
    !i   nlml  :L-cutoff for charge density on radial mesh
    !o Outputs
    !o   rho0  :fit density of form rho0*exp(-decay*r)
    !o   decay :fit density of form rho0*exp(-decay*r)
    !o   res   :integral of fit density from rofi(nr) to infinity
    !l Local variables
    !l         :
    !r Remarks
    !r
    !u Updates
    !u   27 Sep 03  First created
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: nr,nlml,nsp
    double precision :: rofi(nr),rho(nr,nlml,2),rho0,decay,res
    ! ... Local parameters
    integer :: ir
    double precision :: norm(2,2),tnorm(2,2),rhs(2),y,dy,fac,y0,r0,pi,a,b
    res = 0
    if (nr < 10 .OR. nsp == 1) return
    call dpzero(norm,4)
    call dpzero(rhs,2)
    fac = 1
    if (rho(nr,1,1) < rho(nr,1,2)) fac = -1
    do  ir = nr-5, nr
       y = fac*(rho(ir,1,1) - rho(ir,1,2))/rofi(ir)**2
       !       If the spin density changes sign, nonsensical to try and fit
       if (y <= 0) return
       dy = dlog(y)
       norm(1,1) = norm(1,1) + 1
       norm(1,2) = norm(1,2) + rofi(ir)
       norm(2,1) = norm(2,1) + rofi(ir)
       norm(2,2) = norm(2,2) + rofi(ir)**2
       rhs(1) = rhs(1) + dy
       rhs(2) = rhs(2) + rofi(ir)*dy
    enddo
    call dinv22(norm,tnorm)
    a = tnorm(1,1)*rhs(1) + tnorm(1,2)*rhs(2)
    b = tnorm(2,1)*rhs(1) + tnorm(2,2)*rhs(2)
    pi = 4d0*datan(1d0)
    y0 = 1d0/dsqrt(4d0*pi)
    a = fac*exp(a)/y0
    b = -b
    !     Nonsensical if density not decaying fast enough
    if (b < 1) return
    r0 = rofi(nr)
    !     Integral a*exp(-b*r)*r*r = (2+2*b*r0+b**2*r0*2)/b**3*exp(-b*r0)
    res = a*(2+2*b*r0+b**2*r0**2)/b**3*exp(-b*r0)
    decay = b
    rho0 = a
  end subroutine mkrou6

  subroutine mkrou3(mode,lmxa,nlml,nsp,pnz,dmatl,hab,sab,qsum,hsum)
    !- l-decomposed charges and eigenvalue sum
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 do nothing
    !i         :1 make qsum only
    !i         :2 make qsum and hsum both
    !i   lmxa  :augmentation l-cutoff
    !i   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a y_lm
    !i         :expansion of the function products f_i(l1) f_j(l2)
    !i         :where f_i,f_j form this table.
    !i         :    (uu  us  uz)
    !i         :    (su  ss  sz)
    !i         :    (zu  zs  zz)
    !i   hab   :<u,s | h | u,s> for each pair uu, us, su, ss; see remarks
    !i   sab   :<u,s | 1 | u,s>
    !o Outputs
    !o   qsum  :l-decomposed sphere charge
    !o   hsum  :l-decomposed one-electron energy
    !r Remarks
    !r   qsum and hsum are used to find the band centers of gravity
    !r   u and s ae linear combinations radial wave functions defined as:
    !r   u has val=1, slo=1 at rmax, s has val=0, slo=1
    !u Updates
    !u   28 Aug 01 Extended to local orbitals.
    ! ----------------------------------------------------------------------
    implicit none
    integer :: mode,lmxa,n0,nab,nlml,nsp
    parameter (n0=10,nab=9)
    real(8),optional:: hab(nab,n0,nsp), hsum(n0,nsp)
    real(8):: qsum(n0,nsp) ,sab(nab,n0,nsp), &
         dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp),pnz(n0,2)
    integer :: l,isp,m
    double precision :: pi,srfpi,qz,hz
    if (mode == 0) return
    !     call tcn('mkrou3')
    pi = 4d0*datan(1d0)
    srfpi = dsqrt(4d0*pi)
    do  isp = 1, nsp
       do  l = 0, lmxa
          if (mode >= 1) qsum(l+1,isp) = 0d0
          if (mode >= 2) hsum(l+1,isp) = 0d0
       enddo
    enddo
    do  isp = 1, nsp
       do  l = 0, lmxa
          m = l+1
          qsum(m,isp) = &
               + dmatl(l,l,1,1,1,isp)*sab(1,m,isp)*srfpi &
               + dmatl(l,l,1,1,2,isp)*sab(2,m,isp)*srfpi &
               + dmatl(l,l,1,2,1,isp)*sab(3,m,isp)*srfpi &
               + dmatl(l,l,1,2,2,isp)*sab(4,m,isp)*srfpi
          if (mode >= 2) then
             hsum(m,isp) = &
                  + dmatl(l,l,1,1,1,isp)*hab(1,m,isp)*srfpi &
                  + dmatl(l,l,1,1,2,isp)*hab(2,m,isp)*srfpi &
                  + dmatl(l,l,1,2,1,isp)*hab(3,m,isp)*srfpi &
                  + dmatl(l,l,1,2,2,isp)*hab(4,m,isp)*srfpi
          endif
          !         ... uz, sz, zu, zs, zz terms
          if (pnz(m,1) /= 0) then
             qz = &
                  + dmatl(l,l,1,1,3,isp)*sab(5,m,isp)*srfpi &
                  + dmatl(l,l,1,2,3,isp)*sab(6,m,isp)*srfpi &
                  + dmatl(l,l,1,3,1,isp)*sab(5,m,isp)*srfpi &
                  + dmatl(l,l,1,3,2,isp)*sab(6,m,isp)*srfpi &
                  + dmatl(l,l,1,3,3,isp)*sab(7,m,isp)*srfpi
             ! 13sep2012takao add qz and hz
             qsum(m,isp) = qsum(m,isp) + qz  !qsum is including local orbital 13sep2012takao
             if (mode >= 2) then
                hz = &
                     + dmatl(l,l,1,1,3,isp)*hab(5,m,isp)*srfpi &
                     + dmatl(l,l,1,2,3,isp)*hab(6,m,isp)*srfpi &
                     + dmatl(l,l,1,3,1,isp)*hab(5,m,isp)*srfpi &
                     + dmatl(l,l,1,3,2,isp)*hab(6,m,isp)*srfpi &
                     + dmatl(l,l,1,3,3,isp)*hab(7,m,isp)*srfpi
                ! 13sep2012takao add qz and hz
                hsum(m,isp)=hsum(m,isp)+hz  !hsum is including local orbital 13sep2012takao
             endif
          endif

       enddo
    enddo
    !     call tcx('mkrou3')
  end subroutine mkrou3

  subroutine dinv22(a,ainv) ! ainv  :Inverse of a,  A.Chantis
    double precision :: a(2,2), ainv(2,2)
    double precision :: det,aloc(2,2)
    det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    if (det == 0d0) call rx('INV22: vanishing determinant')
    ainv(1,1) = a(2,2)/det
    ainv(2,2) = a(1,1)/det
    ainv(1,2) = -a(1,2)/det
    ainv(2,1) = -a(2,1)/det
  end subroutine dinv22
end module m_mkrout
