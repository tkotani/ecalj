module m_mkrout
  use m_lmfinit,only: nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,lmxb_i=>lmxb,lmxl_i=>lmxl,spec_a,&
       kmxt_i=>kmxt
  use m_ll,only:ll
  use m_struc_def,only:s_rv1
  use m_uspecb,only:uspecb
  use m_lmfinit,only:nkapii,lhh,nkap0,n0
  use m_orbl,only: Orblib,ktab,ltab,offl,norb, ntab,blks
  use m_symrhoat,only:symrhoat
  public:: m_Mkrout_init, orhoat_out, frcbandsym, qbyl_rv,hbyl_rv, sumec,sumtc,sumt0
  type(s_rv1),protected,allocatable :: orhoat_out(:,:)
  real(8),allocatable,protected::  frcbandsym(:,:), hbyl_rv(:,:,:), qbyl_rv(:,:,:)
  real(8),protected:: sumec,sumtc,sumt0
  private
contains
  subroutine m_mkrout_init()
    use m_lmfinit,only: ispec,nsp,lfrce,lrout,nbas
    use m_bandcal,only: oqkkl,oeqkkl,frcband
    use m_mkpot,only: hab_rv,sab_rv
    use m_suham,only: ham_ndham
    use m_rhocor,only: getcor
    integer:: ib,is,nr,lmxl,nlml
    call tcn('m_mkrout_init')
    if(lfrce>0 ) then
       if(allocated(frcbandsym)) deallocate(frcbandsym)
       allocate(frcbandsym(3,nbas)) !for forces
       frcbandsym=frcband         !part of atomic force from band part
    endif
    if( .NOT. allocated(orhoat_out)) allocate(orhoat_out(3,nbas))
    do  ib = 1, nbas
       is =  ispec(ib) 
       nr =  nr_i(is)
       lmxl= lmxl_i(is)
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
    call mkrout(oqkkl,oeqkkl, orhoat_out, hab_rv,sab_rv,qbyl_rv,hbyl_rv)
    if (lrout/=00) then !  Symmetrize output atomic density and forces frcbandsym is symmetrized
       call symrhoat ( orhoat_out, qbyl_rv, hbyl_rv, frcbandsym)
    endif
    call tcx('m_mkrout_init')
  end subroutine m_mkrout_init
  subroutine mkrout(oqkkl,oeqkkl,orhoat_out,hab,sab, qbyl,hbyl)!Assembles local output densities out of the qkkl, and core states
    use m_lmfinit,only: nkaphh, ispec,nbas,nsp,lrout,n0,nlmto,nmcore,rsma !lekkl=1
    use m_lgunit,only: stdo
    use m_struc_def
    use m_elocp,only: rsmlss=>rsml, ehlss=>ehl
    use m_density,only: v0pot,v1pot,pnuall,pnzall !read
    use m_hansr,only:corprm
    use m_rhocor,only: getcor
    use m_makusp,only: makusp
    use m_MPItk,only: master_mpi
    use m_fatom,only: sspec
    
    !i   nbas  :size of basis, nsp: 2 for spin-polarized case, otherwise 1. nlmto:dimension of MTO basis.
    !l   lekkl=T for eqkkl/qkkl
    !i   oqkkl :local density-matrix (rhocbl or comparable routine)
    !o   oeqkkl:local part of energy-weighted density matrix
    !i   hab   :hamiltonian matrix elements of radial wave functions at each site.  See potpus for their definition.
    !i   sab   :overlap matrix elements of radial wave functions at each site.  See potpus for their definition.
    !i         :hab and sab are used here to find band cg
    !i   lrout :0 calculate core states part only
    !o Outputs
    !o   orhoat:vector of offsets containing site density
    !o   sumec = sum of foca=0 core eigenvalues
    !o   sumtc = sum of all core kinetic energies
    !o   sumt0 = sum of frozen core kinetic energies
    !o   qbyl  :l-decomposed charge
    !o   hbyl  :l-decomposed eigenvalue sum
    !r Remarks u and s are linear combinations of and phi,phidot defined as: u has val=1, slo=1 at rmax, s has val=0, slo=1
    implicit none
    type(s_rv1) :: orhoat_out(3,nbas)
    type(s_rv5) :: oeqkkl(3,nbas)
    type(s_rv5) :: oqkkl(3,nbas)
    real(8):: qbyl(n0,nsp,nbas) , hbyl(n0,nsp,nbas) , sab(3,3,n0,nsp,nbas), hab(3,3,n0,nsp,nbas)
    integer:: ib,ipr,iprint,is,k,kcor,kmax,lcor,lfoc,lgunit,lmxa,lmxh,lmxl,nlma,nlmh,nlml,nlml1,r,ncore,ifx
    real(8) :: a,ceh,pi,qcor(2),rfoc,rmt,smec,smtc,stc0, sum1,sum2,sums1,sums2,xx,y0,z,res,rsml(n0),ehl(n0)
    integer :: nkapi,nkape,nr, lh(nkap0),nkaph
    real(8) :: eh(n0,nkap0),rsmh(n0,nkap0),qcorg,qcorh,qsc,cofg,cofh
    real(8),pointer:: pnu(:,:),pnz(:,:)
    real(8):: dat(6,nbas)
    logical:: mmtargetx=.false.
    real(8),allocatable:: rofi_rv(:),rwgt_rv(:)
    ipr  = iprint()
    pi   = 4d0*datan(1d0)
    y0   = 1d0/dsqrt(4d0*pi)
    sums1 = 0
    call tcn('mkrout')
    sumtc = 0d0
    sumec = 0d0
    sumt0 = 0d0
    if(master_mpi) inquire(file='mmtarget.aftest',exist=mmtargetx)
    ibloop: do  ib = 1, nbas
       is = ispec(ib)
       lmxa=lmxa_i(is) 
       if (lmxa == -1) cycle 
       lmxl=lmxl_i(is)
       kmax=kmxt_i(is)
       nkaph=nkaphh(is)
       a=spec_a(is)
       nr=nr_i(is)
       rmt=rmt_i(is)
       lmxh=lmxb_i(is)
       stc0=sspec(is)%stc
       call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
       call uspecb(is,rsmh,eh)
       nkapi=nkapii(is)
       call orblib(ib)!norb,ltab,ktab,norb, ntab,blks
       is =ispec(ib)
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       call gtpcor(is,kcor,lcor,qcor)
       nlml = (lmxl+1)**2
       nlma = (lmxa+1)**2
       nlmh = (lmxh+1)**2
       allocate(rofi_rv(nr),rwgt_rv(nr))
       call radmsh(rmt,a,nr,rofi_rv )
       call radwgt(rmt,a,nr,rwgt_rv )
       if (lrout /= 0) then
          orhoat_out(1,ib)%v=0d0
          orhoat_out(2,ib)%v=0d0
          !   --- Assemble rho1 = products of augmented functions ---
          !  ... Set up all radial head and tail envelope functions, and their bc's
          k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml*nsp
          fhblock: block
            real(8):: chh_rv(k),chp_rv(k),cpp_rv(k),dmatl_rv((lmxa+1)**2*nlml*nsp*9), &
                 fh_rv(nr*(lmxh+1)*nkaph),xh_rv(nr*(lmxh+1)*nkaph),&  !h:head
                 vh_rv((lmxh+1)*nkaph),dh_rv((lmxh+1)*nkaph),&
                 fp_rv(nr*(lmxa+1)*(kmax+1)),xp_rv(nr*(lmxa+1)*(kmax+1)),& !p:tail
                 vp_rv((lmxa+1)*(kmax+1)),dp_rv((lmxa+1)*(kmax+1)),&
                 ul_rv(nr*(lmxa+1)*nsp),sl_rv(nr*(lmxa+1)*nsp),&     !u value function
                 gz_rv(nr*(lmxa+1)*nsp),ruu_rv(nr*(lmxa+1)*2*nsp),&  !s slope function
                 rus_rv(nr*(lmxa+1)*2*nsp),rss_rv(nr*(lmxa+1)*2*nsp) !gz local orbital
            call fradhd(nkaph,eh,rsmh,lhh(:,is),lmxh,nr,rofi_rv, fh_rv,xh_rv,vh_rv,dh_rv) !head
            call fradpk(kmax,rsma(is),lmxa,nr,rofi_rv, fp_rv,xp_rv,vp_rv,dp_rv) !tail
            rsml= rsmlss(:,is)
            ehl = ehlss(:,is)
            call makusp(n0,z,nsp,rmt,lmxa,v0pot(ib)%v,a,nr,&! Augmented wave functions
                 pnu,pnz,rsml,ehl,ul_rv,sl_rv,gz_rv,ruu_rv, rus_rv,rss_rv )
            ! ... Contracted density matrix as coffs to products of (u,s,gz)
            chh_rv=0d0
            chp_rv=0d0
            cpp_rv=0d0
            dmatl_rv=0d0
            call mkrou1(nsp, nlmh,nlma,nlml,kmax, nkaph,nkapi,norb,ltab,ktab,blks &
                 ,oqkkl(3,ib)%v,oqkkl(2,ib)%v,oqkkl(1,ib)%v,vh_rv,dh_rv,vp_rv,dp_rv,  chh_rv,chp_rv,cpp_rv,dmatl_rv )
            call mkrou2(nsp,lmxa,nlml,pnz, dmatl_rv,nr,ul_rv ,sl_rv,gz_rv,ruu_rv,rus_rv,rss_rv,     orhoat_out( 1,ib )%v)!True local density.1st component
            !   --- Assemble rho2 = unaugmented products Pkl*Pk'l' ---
            call mkrou5(nsp,nr,nlml,nkaph,nkaph,fh_rv,lmxh ,nkaph,nkaph,fh_rv,lmxh,           chh_rv,orhoat_out( 2,ib )%v) !H H product.  
            call mkrou5(nsp,nr,nlml,nkaph,nkaph,fh_rv,lmxh ,kmax + 1,kmax + 1,fp_rv,lmxa,     chp_rv,orhoat_out( 2,ib )%v) !H Pkl product
            call mkrou5(nsp,nr,nlml,kmax + 1,kmax + 1,fp_rv,lmxa, kmax+1, kmax + 1,fp_rv,lmxa,cpp_rv,orhoat_out( 2,ib )%v) !Pkl Pkl product
            call mkrou3(lmxa,nlml,nsp,pnz,dmatl_rv,sab(1,1,1,1,ib),qbyl(1,1,ib)) !  qbyl = site charge decomposed by l
            nlml1 = 1
            dmatl_rv=0d0 
            k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml1*nsp
            chh_rv=0d0
            chp_rv=0d0
            cpp_rv=0d0 ! lekkl /= 0 mode only now. !history: hab were/are wrong =>it had given wrong hbyl (in previous mkrou3) 
            call mkrou1(nsp, nlmh, nlma, nlml1,kmax,nkaph,nkapi,norb,ltab,ktab,blks &
                 ,oeqkkl(3,ib)%v, oeqkkl(2,ib)%v,oeqkkl(1,ib)%v,vh_rv,dh_rv,vp_rv,dp_rv,chh_rv ,chp_rv,cpp_rv, dmatl_rv) ! energy weighted dmatl_rv
            call mkrou3(lmxa,nlml=nlml1,nsp=nsp,pnz=pnz,dmatl=dmatl_rv,sab=sab(1,1,1,1,ib ), qsum=hbyl(1,1,ib ))! hbyl is energy-weighted density matrix
            call radsum(nr,nr,nlml,nsp,rwgt_rv,orhoat_out( 1,ib )%v, sum1 )
            call radsum(nr,nr,nlml,nsp,rwgt_rv,orhoat_out( 2,ib )%v, sum2 )
            sum1 = sum1/y0
            sum2 = sum2/y0
            if (nsp == 2) then
               sums1 = sum1 - 2d0 * sum( rwgt_rv(1:nr)*orhoat_out(1,ib)%v(1:nr))/y0
               sums2 = sum2 - 2d0 * sum( rwgt_rv(1:nr)*orhoat_out(2,ib)%v(1:nr))/y0
            endif
            if(nsp==1) dat(1:3,ib) = [sum1,sum2,sum1-sum2]
            if(nsp==2) dat(1:6,ib) = [sum1,sum2,sum1-sum2, -sums1,-sums2,-sums1+sums2]
            ! cccccccccccccccccccccccccccccccccccccccccccccccccccccc
            !!exper. block to keep mag.moment for AF.controlled by uhval.aftest file. See elsewhere.
            if( nsp==2 .AND. master_mpi .AND. mmtargetx) then
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
          endblock fhblock
       endif
       call mkrou6(rofi_rv,orhoat_out( 1,ib )%v,nr,nlml,nsp,xx,xx,res )!Contribution to mag outsite rmt: extrapolate tail to infinity
       if(ipr>=30.AND.res/=0) write(stdo,"(7x,'contr. to mm extrapolated for r>rmt:',f11.6,' est. true mm =',f9.6)")res,res-sums1
       if (lfoc == 0) then ! Make new core density and core eigenvalue sum ---
          call pshpr(ipr+11)
          call getcor(0,z,a,pnu,pnz,nr,lmxa,rofi_rv,v1pot(ib)%v & 
              ,kcor,lcor,qcor,smec,smtc,orhoat_out( 3,ib )%v,ncore,[0d0],[0d0],nmcore(is))
          call poppr
          sumtc = sumtc + smtc
          sumec = sumec + smec
       else
          sumtc = sumtc + stc0
          sumt0 = sumt0 + stc0
          orhoat_out(3,ib)%v(1:nsp*nr)= sspec(is)%rv_a_orhoc(1:nr*nsp) 
       endif
       deallocate(rofi_rv,rwgt_rv)
    enddo ibloop
    if (ipr >= 30 .AND. lrout > 0) then ! write(stdo,"(a)")' mkrout: site(class) decomposed charge and magnetic moment. class->lmchk'
       if (nsp == 1) write(stdo,"(/' mkrout:  Qtrue      sm,loc       local')")
       if (nsp == 2) write(stdo,"(/' mkrout:  Qtrue      sm,loc       local',8x,'true mm   smooth mm    local mm')")
       do ib=1,nbas
          write(stdo,"(i4,3f12.6,2x,3f12.6)") ib,dat(1:3*nsp,ib)
       enddo
    endif
    call tcx('mkrout')
  end subroutine mkrout
  subroutine mkrou1(nsp,nlmh,nlma,nlml,kmax,nkaph,nkapi,norb,ltab,ktab,blks,qhh,qhp,qpp,vh,dh,vp,dp, chh,chp,cpp,dmatl)!Contracted density mat. for a site
    use m_lmfinit,only: cg=>rv_a_ocg,indxcg=>iv_a_oidxcg,jcg=>iv_a_ojcg,cy=>rv_a_ocy
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
    integer :: nsp,nlmh,nlma,nlml,kmax,norb,nkaph,nkapi
    integer :: ltab(norb),ktab(norb),blks(norb)
    real(8) :: qhh(nkaph,nkaph,nlmh,nlmh,nsp)
    real(8) :: qhp(nkaph,kmax+1,nlmh,nlma,nsp)
    real(8) :: qpp(kmax+1,kmax+1,nlmh,nlma,nsp)
    real(8) :: dmatl(*),vh(*),dh(*),vp(*),dp(*), chh(*),chp(*),cpp(*)
    integer :: k,lmxa,lmxh
    integer :: ltba(n0*nkap0),ktba(n0*nkap0),blka(n0*nkap0)
    lmxa = ll(nlma)
    lmxh = ll(nlmh)
    ltba = 0
    do  k = 0, kmax
       ktba(k+1) = k+1
       blka(k+1) = nlma
    enddo
    ! ... Contracted density-matrix chh,chp,cpp from qhh,qhp,qpp
    call mkrou4(nsp,nlml,cg,jcg,indxcg, nkaph,norb,ltab,ktab,blks,lmxh,nlmh, nkaph,norb,ltab,ktab,blks,lmxh,nlmh, qhh,chh) !H H product
    call mkrou4(nsp,nlml,cg,jcg,indxcg, nkaph,norb,ltab,ktab,blks,lmxh,nlmh, kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma, qhp,chp) !H Pkl product
    call mkrou4(nsp,nlml,cg,jcg,indxcg, kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma, kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma, qpp,cpp)!Pkl Pkl product
    ! ... Contracted density matrix as coffs to products of (u,s,gz)
    call mkcfus(nsp,lmxa,nlml, nkaph,nkapi,vh,dh,lmxh, nkaph,nkapi,vh,dh,lmxh, chh,dmatl) !H H product
    call mkcfus(nsp,lmxa,nlml, nkaph,nkapi,vh,dh,lmxh, kmax+1,kmax+1,vp,dp,lmxa, chp,dmatl) !H Pkl product
    call mkcfus(nsp,lmxa,nlml, kmax+1,kmax+1,vp,dp,lmxa, kmax+1,kmax+1,vp,dp,lmxa, cpp,dmatl) !Pkl Pkl product
  end subroutine mkrou1
  subroutine mkrou2(nsp,lmxa,nlml,pnz,dmatl,nr,ul,sl,gz,ruu,rus,rss, rho) !- Assemble true site density from product function coefficients
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
    integer :: nsp,lmxa,nlml,nr
    real(8) :: rho(nr,nlml,nsp),pnz(n0,nsp), &
         dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp), &
         ul(nr,0:lmxa,nsp),ruu(nr,0:lmxa,2,nsp), &
         sl(nr,0:lmxa,nsp),rus(nr,0:lmxa,2,nsp), &
         gz(nr,0:lmxa,nsp),rss(nr,0:lmxa,2,nsp)
    logical :: lpz1,lpz2
    integer :: l,i,mlm,l1,l2,isp
    real(8) :: xuu,xus,xsu,xss,xuz,xsz,xzu,xzs,xzz
    call tcn('mkrou2')
    ! --- Full density as products (u,s) (u,s); no small component ---
    do  isp = 1, nsp
       do  mlm = 1, nlml
          do  l1 = 0, lmxa
             do  l2 = 0, lmxa
                lpz1 = pnz(l1+1,1) .ne. 0
                lpz2 = pnz(l2+1,1) .ne. 0
                rho(:,mlm,isp) = rho(:,mlm,isp) &
                     + dmatl(l1,l2,mlm,1,1,isp) * ul(:,l1,isp) * ul(:,l2,isp) &
                     + dmatl(l1,l2,mlm,2,1,isp) * sl(:,l1,isp) * ul(:,l2,isp) &
                     + dmatl(l1,l2,mlm,1,2,isp) * ul(:,l1,isp) * sl(:,l2,isp) &
                     + dmatl(l1,l2,mlm,2,2,isp) * sl(:,l1,isp) * sl(:,l2,isp)
                if (lpz1.OR.lpz2) then
!                   if (xuz /= 0 .OR. xsz /= 0 .OR. xzu /= 0 .OR. xzs /= 0 .OR. xzz /= 0) then !commented 2023-5-24
                      rho(:,mlm,isp) = rho(:,mlm,isp) &
                           + dmatl(l1,l2,mlm,1,3,isp) * ul(:,l1,isp) * gz(:,l2,isp) &
                           + dmatl(l1,l2,mlm,2,3,isp) * sl(:,l1,isp) * gz(:,l2,isp) &
                           + dmatl(l1,l2,mlm,3,1,isp) * gz(:,l1,isp) * ul(:,l2,isp) &
                           + dmatl(l1,l2,mlm,3,2,isp) * gz(:,l1,isp) * sl(:,l2,isp) &
                           + dmatl(l1,l2,mlm,3,3,isp) * gz(:,l1,isp) * gz(:,l2,isp)
!                   endif
                endif
             enddo
          enddo
       enddo
    enddo
    rho(1:nr,1,1:nsp)=0d0  ! --- Remake spherical density including small component ---
    do  isp = 1, nsp
       do  l = 0, lmxa
          xuu = dmatl(l,l,1,1,1,isp)
          xus = dmatl(l,l,1,1,2,isp) + dmatl(l,l,1,2,1,isp)
          xss = dmatl(l,l,1,2,2,isp)
          xuz = dmatl(l,l,1,1,3,isp) + dmatl(l,l,1,3,1,isp)
          xsz = dmatl(l,l,1,2,3,isp) + dmatl(l,l,1,3,2,isp)
          xzz = dmatl(l,l,1,3,3,isp)
          rho(:,1,isp) = rho(:,1,isp) +xuu*ruu(:,l,1,isp)+xus*rus(:,l,1,isp)+xss*rss(:,l,1,isp)
          if(pnz(l+1,1)/=0) rho(:,1,isp)=rho(:,1,isp) +xuz*ruu(:,l,2,isp)+xsz*rus(:,l,2,isp)+xzz*rss(:,l,2,isp)
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
    implicit none
    integer :: nsp,nlml,norb1,lmx1,ndim1,norb2,lmx2,ndim2,nk1,nk2
    integer :: jcg(1),indxcg(1)
    integer :: ltab1(norb1),ktab1(norb1),blks1(norb1), &
         ltab2(norb2),ktab2(norb2),blks2(norb2)
    real(8) :: qkk12(nk1,nk2,ndim1,ndim2,nsp),cg(1), &
         ckk(nk1,nk2,0:lmx1,0:lmx2,nlml,nsp)
    integer :: ilm1,io1,l1,nlm11,nlm12,k1, &
         ilm2,io2,l2,nlm21,nlm22,k2, icg,mlm,ix,isp
    !     call tcn('mkrou4')
    do  isp = 1, nsp
       do  io2 = 1, norb2
          if (blks2(io2) ==0) cycle 
          l2 = ltab2(io2)  ! k2,l2 = k and starting l index for this block
          k2 = ktab2(io2)
          nlm21 = l2**2+1
          nlm22 = nlm21 + blks2(io2)-1
          do  ilm2 = nlm21, nlm22
             l2 = ll(ilm2)
             do  io1 = 1, norb1
                if(blks1(io1) ==0) cycle
                l1 = ltab1(io1) !  k1,l1 = k and starting l index for this block
                k1 = ktab1(io1)
                nlm11 = l1**2+1
                nlm12 = nlm11 + blks1(io1)-1
                do  ilm1 = nlm11, nlm12
                   l1 = ll(ilm1)
                   ix = max0(ilm1,ilm2)
                   ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
                   do icg = indxcg(ix),indxcg(ix+1)-1
                      mlm = jcg(icg)
                      if(mlm<=nlml) then
                         ckk(k1,k2,l1,l2,mlm,isp) = ckk(k1,k2,l1,l2,mlm,isp)+cg(icg)*qkk12(k1,k2,ilm1,ilm2,isp)
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    !     call tcx('mkrou4')
  end subroutine mkrou4
  subroutine mkcfus(nsp,lmxa,nlml, nf1,nf1s,val1,slo1,lmx1, nf2,nf2s,val2,slo2,lmx2,ckk,dmatl)!Assemble contracted density mat. as coffs to products of (u,s,gz)
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
    implicit none
    integer :: nsp,nlml,lmxa,lmx1,lmx2,nf1,nf2,nf1s,nf2s
    real(8) :: ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nsp), &
         val1(0:lmx1,nf1s),slo1(0:lmx1,nf1s), &
         val2(0:lmx2,nf2s),slo2(0:lmx2,nf2s)
    real(8) :: dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp)
    integer :: l1,k1,l2,k2,mlm,isp
    real(8) :: xx
    !     call tcn('mkcfus')
    do  isp = 1, nsp
       do  mlm = 1, nlml
          do  k2 = 1, nf2s
             do  k1 = 1, nf1s
                do  l2 = 0, lmx2
                   do  l1 = 0, lmx1
                      xx = ckk(k1,k2,l1,l2,mlm,isp)
                      dmatl(l1,l2,mlm,1,1,isp)= dmatl(l1,l2,mlm,1,1,isp)+ xx*val1(l1,k1) * val2(l2,k2)
                      dmatl(l1,l2,mlm,1,2,isp)= dmatl(l1,l2,mlm,1,2,isp)+ xx*val1(l1,k1) * slo2(l2,k2)
                      dmatl(l1,l2,mlm,2,1,isp)= dmatl(l1,l2,mlm,2,1,isp)+ xx*slo1(l1,k1) * val2(l2,k2)
                      dmatl(l1,l2,mlm,2,2,isp)= dmatl(l1,l2,mlm,2,2,isp)+ xx*slo1(l1,k1) * slo2(l2,k2)
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
                         if (k1 > nf1s .AND. k2 > nf2s) then!               sc-sc product
                            dmatl(l1,l2,mlm,3,3,isp) = dmatl(l1,l2,mlm,3,3,isp)+xx 
                         elseif (k1 > nf1s) then ! sc-valence product
                            dmatl(l1,l2,mlm,3,1,isp) = dmatl(l1,l2,mlm,3,1,isp) + xx*val2(l2,k2)
                            dmatl(l1,l2,mlm,3,2,isp) = dmatl(l1,l2,mlm,3,2,isp) + xx*slo2(l2,k2)
                         elseif (k2 > nf2s) then!  valence-sc product
                            dmatl(l1,l2,mlm,1,3,isp) = dmatl(l1,l2,mlm,1,3,isp) + xx*val1(l1,k1) 
                            dmatl(l1,l2,mlm,2,3,isp) = dmatl(l1,l2,mlm,2,3,isp) + xx*slo1(l1,k1) 
                         endif
                      enddo
                   enddo
                endif
             enddo
          enddo
       enddo
    enddo
    !     call tcx('mkcfus')
  end subroutine mkcfus
  subroutine mkrou5(nsp,nr,nlml,nf1,nf1s,f1,lmx1,nf2,nf2s,f2,lmx2, ckk,rho)
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
    implicit none
    integer :: nsp,nr,nlml,lmx1,lmx2,nf1,nf2,nf1s,nf2s,l1,k1,l2,k2,mlm,isp,i,ir
    real(8) :: ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nsp),f1(nr,0:lmx1,nf1s),f2(nr,0:lmx2,nf2s),rho(nr,nlml,nsp)
!    rho(:,:,:) =rho(:,:,:)+ ckk(k1,k2,l1,l2,:,:)*[(f1(ir,l1,k1)*f2(:,l2,k2),ir=1,nr)]
    do  k2 = 1, nf2s
       do  k1 = 1, nf1s
          do  l2 = 0, lmx2
             do  l1 = 0, lmx1
                do ir=1,nr
                   rho(ir,:,:)=rho(ir,:,:) + ckk(k1,k2,l1,l2,:,:)*f1(ir,l1,k1)*f2(ir,l2,k2)
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine mkrou5
  subroutine mkrou6(rofi,rho,nr,nlml,nsp,rho0,decay,res) !- Fit tail of spin density; integrate charge beyond MT sphere
    !i Inputs
    !i   rofi  :radial mesh points
    !i   rho   :spin-polarized charge density
    !i   nr    :number of radial mesh points
    !i   nlml  :L-cutoff for charge density on radial mesh
    !o Outputs
    !o   rho0  :fit density of form rho0*exp(-decay*r)
    !o   decay :fit density of form rho0*exp(-decay*r)
    !o   res   :integral of fit density from rofi(nr) to infinity
    implicit none
    integer :: nr,nlml,nsp
    real(8) :: rofi(nr),rho(nr,nlml,2),rho0,decay,res
    integer :: ir
    real(8) :: norm(2,2),tnorm(2,2),rhs(2),y,dy,fac,y0,r0,pi,a,b
    res = 0
    if (nr < 10 .OR. nsp == 1) return
    call dpzero(norm,4)
    call dpzero(rhs,2)
    fac = 1
    if (rho(nr,1,1) < rho(nr,1,2)) fac = -1
    do  ir = nr-5, nr
       y = fac*(rho(ir,1,1) - rho(ir,1,2))/rofi(ir)**2
       if (y <= 0) return ! If the spin density changes sign, nonsensical to try and fit
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
    if (b < 1) return !   Nonsensical if density not decaying fast enough
    r0 = rofi(nr)
    res = a*(2+2*b*r0+b**2*r0**2)/b**3*exp(-b*r0) ! Integral a*exp(-b*r)*r*r = (2+2*b*r0+b**2*r0*2)/b**3*exp(-b*r0)
    decay = b
    rho0 = a
  end subroutine mkrou6
  subroutine mkrou3(lmxa,nlml,nsp,pnz,dmatl,sab,qsum) ! l-decomposed charges and eigenvalue sum
    !i   lmxa  :augmentation l-cutoff
    !i   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a y_lm
    !i         :expansion of the function products f_i(l1) f_j(l2)
    !i         :where f_i,f_j form this table.
    !i         :    (uu  us  uz)
    !i         :    (su  ss  sz)
    !i         :    (zu  zs  zz)
    !i   sab   :<u,s | 1 | u,s>
    !o Outputs qsum  :l-decomposed sphere charge
    !r u and s are linear combinations radial wave functions defined as: u has val=1, slo=1 at rmax, s has val=0, slo=1
    implicit none
    integer ::lmxa,nlml,nsp,l,isp,m
    real(8):: qsum(n0,nsp) ,sab(3,3,n0,nsp), dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp),pnz(n0,2),qz,hz
    real(8),parameter:: pi=4d0*datan(1d0),srfpi = dsqrt(4d0*pi)
    qsum=0d0
    do  isp = 1, nsp
       do  l = 0, lmxa
          m = l+1
          qsum(m,isp) = &
               + dmatl(l,l,1,1,1,isp)*sab(1,1,m,isp)*srfpi &
               + dmatl(l,l,1,1,2,isp)*sab(1,2,m,isp)*srfpi &
               + dmatl(l,l,1,2,1,isp)*sab(2,1,m,isp)*srfpi &
               + dmatl(l,l,1,2,2,isp)*sab(2,2,m,isp)*srfpi
          if (pnz(m,1) /= 0) then!         ... uz, sz, zu, zs, zz terms
             qz =   dmatl(l,l,1,1,3,isp)*sab(1,3,m,isp)*srfpi &
                  + dmatl(l,l,1,2,3,isp)*sab(2,3,m,isp)*srfpi &
                  + dmatl(l,l,1,3,1,isp)*sab(1,3,m,isp)*srfpi &
                  + dmatl(l,l,1,3,2,isp)*sab(2,3,m,isp)*srfpi &
                  + dmatl(l,l,1,3,3,isp)*sab(3,3,m,isp)*srfpi
             qsum(m,isp) = qsum(m,isp) + qz  !qsum is including local orbital
          endif
       enddo
    enddo
  end subroutine mkrou3
  subroutine dinv22(a,ainv) ! ainv  :Inverse of a,  A.Chantis
    real(8) :: a(2,2), ainv(2,2),det,aloc(2,2)
    det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    if (det == 0d0) call rx('INV22: vanishing determinant')
    ainv(1,1) = a(2,2)/det
    ainv(2,2) = a(1,1)/det
    ainv(1,2) = -a(1,2)/det
    ainv(2,1) = -a(2,1)/det
  end subroutine dinv22
end module m_mkrout
