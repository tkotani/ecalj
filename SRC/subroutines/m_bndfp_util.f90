!>utils called in bndfp
module m_bndfp_util
  use m_MPItk,only:master_mpi
  use m_ll,only:ll
  use m_lgunit,only:stdo
  public mkekin,makdos,phispinsym_ssite_set,iorbtm
contains
  subroutine mkekin(osig,otau,oppi,oqkkl,vconst,smpot,smrho,sumev, ekinval) !- Evaluate the valence kinetic energy
    use m_struc_def
    use m_lmfinit,only:nsp,nspc,nbas,ispec,nlmto
    use m_lmfinit,only: lmxa_i=>lmxa,lmxb_i=>lmxb,kmxt_i=>kmxt
    use m_lattic,only: lat_vol
    use m_supot,only: n1,n2,n3
    use m_orbl,only: Orblib,ktab,ltab,offl,norb,ntab,blks
    !NOTE: When SOC included, ek contains HSO, because ek= Eband- V*n where Eband contains SO contr.
    !i Inputs
    !i   osig  :augmentation overlap integrals
    !i   otau  :augmentation kinetic energy integrals
    !i   oppi  :augmentation kinetic + potential integrals
    !i   oqkkl :local density-matrix; see rlocbl.f
    !i   vconst:constant potential added to hamiltonian
    !i   smpot :smooth input potential on uniform mesh (mkpot.f)
    !i         :smpot = Ves~ + vxc = Ves(n0 + compensating gaussians) + vxc
    !i   smrho :smooth output density n0 (no gaussians n0~-n0)
    !i   sumev :sum of eigenvalues
    !o Outputs
    !o   ekinval :kinetic energy
    !l Local variables
    !l   sraugm:sum_ib q * (tau+ppi-tau) : corresponds to valftr in mkpot.f
    !l   smresh:sm rho * sm V ; corresponds to valfsm in mkpot.f
    !r Remarks
    !r   The valence kinetic energy is evaluated in the usual way as
    !r        ekinval = sumev - srhov
    !r   where sumev is the sum-of-eigenvalues and srhov is the integral
    !r   of the output density and input potential.
    !r   Integrals of the density with the xc potential are folded into the
    !r   electrostatic parts, that is:
    !r     V0 = V0es + V0xc  V1 = V1es + V1xc  V2 = V2es + V2xc
    !r   and are not discussed here.
    !r
    !r   mkekin make the electrostatic integral
    !r     int (n0~ Ves~ + n1 Ves1 - n2 Ves2~)                         (40)
    !r   as described in
    !r      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
    !r      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
    !r      ed. (Springer-Verlag, Berlin) 2000.
    !r   for an output density defined through (smrho,qkkl) and an input
    !r   potential defined through smpot and the matrix elements (ppi-tau).

    !r   Consider only one atom/cell for simplicity. The output density is
    !r     n = sum_ij Dij F~i F~j with Dij = {sum_n w_n z*_in z_jn}   (32)
    !r   i.e. the contraction of the density matrix with partial densities
    !r     F~i F~j = Fi Fj +
    !r               sum_kLk'L' C^i_kL (P~kLP~k'L' - PkLPk'L') C^j_k'L' (33)
    !r             = n0_ij + n1_ij - n2_ij
    !r   Note that local parts of the partial densities have two `levels'
    !r   of decomposition, namely at the `ij' level as in Eq. 33, or
    !r   at a still finer level in which the (kLk'L') indices are not
    !r   summed over.  Thus
    !r     n{1,2} = sum_ij D_ij n{1,2}_ij
    !r     n{1,2}_ij = sum_kLk'L' C^i_kL n{1,2}_kL,k'L' C^j_k'L'
    !r     n{1,2}_kL,k'L' = PkL Pk'L'
    !r   Note also that the 'k index' may be a sum over polynomials, or when
    !r   function `heads' are dealt with, the function itself, as described
    !r   in augmat.f.  As in making the matrix elements, we have to deal
    !r   with three cases, HH; HP; PP, but this is a inessential detail
    !r   needed only because representing H with sums of polynomials tends
    !r   to be ill-conditioned, and in the description below we ignore it.
    !r
    !r   Densities n0 and n2 have corresponding n0~ and n2~ which include
    !r   the additional multipole terms that guarantee n1 and n2~ have
    !r   the same multipole moments.  Thus:
    !r     n0~ = n0 + sum_M q_M G_M
    !r     n2~ = n2 + sum_M q_M G_M
    !r   where q_M are the difference in multipole moments between n1 and n2
    !r     q_M = int dr Y_M r^m (n1 - n2)
    !r   We can define partial densities for multipole contributions as well
    !r     n2~-n2 = sum_ij D_ij (n2~-n2)_ij
    !r     (n2~-n2)_ij = sum_M Q_ijM G_M
    !r                 = sum_kLk'L'M C^i_kL Q_kkLL'M G_M C^j_k'L'
    !r   with the two forms decomposing q_M into two levels:
    !r     q_M = sum_ij D_ij Q_ijM
    !r     Q_ijM = sum_kLk'L' C^i_kL Q_kkLL'M C^j_k'L'
    !r     Q_kkLL'M = int dr Y_M r^m (P~kL P~k'L' - PkL Pk'L')         (27)
    !r
    !r   Using the identity
    !r     n2~ - n2 = n0~ - n0 = sum_M q_M G_M
    !r   Eq. 40 is evaluated as
    !r     int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + sum_M q_M G_M (Ves0~-Ves2~)
    !r   The first term is evaluated on the mesh and stored in srmesh
    !r   The remaining terms amount to products of the density-matrix
    !r   and the ppi matrix elements.  Thus:
    !r     int n1 Ves1 = sum_ij D_ij int n1_ij Ves1
    !r     int n1_ij Ves1 = sum_kLk'L' C^i_kL int n1_kL,k'L' Ves1 C^j_k'L'
    !r                    = sum_kLk'L' C^i_kL pi1_kk'LL' C^j_k'L'
    !r   where pi1 is the first term of the pi matrix element, Eq. 29:
    !r     pi1_kk'LL' = P~kL V1 P~k'L'
    !r   Similarly for the second term, substituting n2 for n1 and
    !r   Ves2~ for Ves1.
    !r     int n2 Ves2~ = sum_ij D_ij int n2_ij Ves2~
    !r     int n2_ij Ves2~ = sum_kLk'L' C^i_kL int n2_kL,k'L' Ves2~ C^j_k'L'
    !r                    = sum_kLk'L' C^i_kL pi2_kk'LL' C^j_k'L'
    !r     pi2_kk'LL' = P~kL V1 P~k'L'
    !r   The last term just amounts to products of the density-matrix and
    !r   the remaining parts of the ppi matrix element:
    !r     pi_kk'LL'  = pi1_kk'LL' - pi2_kk'LL' + pi3_kk'LL'
    !r     pi3_kk'LL' = sum_M Q_kkLL'M int G_M (Ves0~ - Ves2~)
    !r   Evaluating the last term in the electrostatic integral we have
    !r     rhoV_MP = int sum_M q_M G_M (Ves0~ - Ves2~)
    !r             = int sum_ij D_ij sum_M Q_ijM G_M (Ves0~ - Ves2~)
    !r             = sum_ij D_ij sum_kLk'L'M C^i_kL pi3_kk'LL' C^j_k'L'
    !r   which follows using the relationship between Q_kkLL'M and Q_ijM
    !r   Using the definition of the local density-matrix (see rlocbl.f)
    !r      qpp_kLk'L' = sum_ij D_ij C^i_kL C^j_k'L'
    !r   the electrostatic integral then becomes
    !r     int rhoVes = int n0 Ves0~ + n1 Ves1 - n2 Ves2~ + rhoV_MP
    !r                = int n0 Ves0~
    !r                + sum_ij D_ij sum_kLk'L' C^i_kL pi_kk'LL' C^j_k'L'
    !r                = int n0 Ves0~ + sum_kLk'L' qpp'LL' pi_kk'LL'
    !r
    implicit none
    type(s_cv5),target :: oppi(3,nbas)
    type(s_rv4),target :: otau(3,1)
    type(s_rv4),target :: osig(3,1)
    type(s_rv5),target :: oqkkl(3,nbas)
    real(8),pointer:: QPP(:,:,:,:,:), QHP(:,:,:,:,:),QHH(:,:,:,:,:), &
         SIGHH(:,:,:,:),SIGHP(:,:,:,:),SIGPP(:,:,:,:), &
         TAUHH(:,:,:,:),TAUHP(:,:,:,:),TAUPP(:,:,:,:)
    complex(8):: smpot(n1,n2,n3,nsp),smrho(n1,n2,n3,nsp)
    complex(8),pointer:: PPIHH(:,:,:,:,:),PPIHP(:,:,:,:,:),PPIPP(:,:,:,:,:)
    real(8):: sumev, ekinval , vconst
    integer :: ib,igetss,ipr,is,kmax,lgunit,lmxa,lmxh,n0, ngabc(3),nglob,nkap0,nlma,nlmh,ilm1
    logical :: lgors
    parameter (n0=10,nkap0=3)
    double precision :: qum1,qum2,sraugm,srhov,srmesh,sum1,sum2,sumh,sumq,sumt,vol,xx
    integer:: procid=0,ier=0,io,iorb
    integer,parameter::master=0
    logical:: iprx
!    call mpi_comm_rank(mpi_comm_world,procid,ier)
    iprx= master_mpi
    call tcn('mkekin')
    call getpr(ipr)
!    ngabc=lat_nabc
    vol=lat_vol
    srmesh = (dreal(sum(smpot*smrho)) + vconst*dreal(sum(smrho))) *vol/(n1*n2*n3) ! Integral n0(out) (Ves0~ + Vxc0), contribution from mesh
    ! Note that it does not include the term (n0~-n0) Ves0~
    sraugm = 0d0
    do  ib = 1, nbas !Integral rhout*veff, part from augmentation ---
       is = ispec(ib) 
       lmxa=lmxa_i(is)
       kmax=kmxt_i(is)
       lmxh=lmxb_i(is)
       if (lmxa == -1) cycle
       call orblib(ib) ! norb , ltab , ktab , offl ,ntab,blks
       nlma = (lmxa+1)**2
       nlmh = (lmxh+1)**2
       QHH   => oqkkl(3,ib)%v !head x head index=3
       QHP   => oqkkl(2,ib)%v !head x tail index=2
       QPP   => oqkkl(1,ib)%v !tail x tail index=1
       TAUHH => otau(3,ib)%v
       TAUHP => otau(2,ib)%v
       TAUPP => otau(1,ib)%v
       SIGHH => osig(3,ib)%v
       SIGHP => osig(2,ib)%v
       SIGPP => osig(1,ib)%v
       PPIHH => oppi(3,ib)%cv
       PPIHP => oppi(2,ib)%cv
       PPIPP => oppi(1,ib)%cv
       ! real(8):: &
       !      tauhh(nkaph,nkaph,0:lmxh,nsp),  sighh(nkaph,nkaph,0:lmxh,nsp), &
       !      tauhp(nkaph,0:kmax,0:lmxh,nsp), sighp(nkaph,0:kmax,0:lmxh,nsp), &
       !      taupp(0:kmax,0:kmax,0:lmxa,nsp),sigpp(0:kmax,0:kmax,0:lmxa,nsp), &
       !      qhh(nkaph,nkaph,nlmh,nlmh,nsp), &
       !      qhp(nkaph,0:kmax,nlmh,nlma,nsp), &
       !      qpp(0:kmax,0:kmax,nlma,nlma,nsp)
       ! complex(8)::&
       !      ppihh(nkaph,nkaph,nlmh,nlmh,nsp), &
       !      ppihp(nkaph,0:kmax,nlmh,nlma,nsp), &
       !      ppipp(0:kmax,0:kmax,nlma,nlma,nsp)
       !i   kmax  :cutoff in PkL expansion
       !i   lmxa  :dimensions sigpp, taupp
       !i   nlma  :L cutoff in PkL expansion
       !i   nkaph :dimensions augmentation matrices
       !i   norb  :number of orbitals for this site
       !i   ltab  :table of l quantum numbers for the orbitals
       !i   ktab  :table of k numbers (orbital type) for the orbitals
       !i   blks  :block size for grouping orbitals into blocks (gtbls1)
       !i   lmxh  :dimensions sighh, sighp, tauhh, tauhp
       !i   nlmh  :dimensions heads ppi and qhh and qhp
       !i   tauhh :head-head kinetic energy integrals (augmat.f)
       !i   sighh :head-head overlap integrals (augmat.f)
       !i   ppihh :head-head kinetic + potential integrals (augmat.f)
       !i   tauhp :head-tail kinetic energy integrals (augmat.f)
       !i   sighp :head-tail overlap integrals (augmat.f)
       !i   ppihp :head-tail kinetic + potential integrals (augmat.f)
       !i   taupp :tail-tail kinetic energy integrals (augmat.f)
       !i   sigpp :tail-tail overlap integrals (augmat.f)
       !i   ppipp :tail-tail potential integrals (augmat.f)
       !i   qhh   :head-head density matrix for this site
       !i   qhp   :head-tail density matrix for this site
       !i   qpp   :tail-tail density matrix for this site
       !i   nsp   :number of spin channels
       !i   nspc  :2 for coupled spins; otherwise 1
       pvgtknblock: block !Local contribution to kinetic energy for one site
         integer :: ilm1,ilm2,k1,k2,nlm11,nlm12,nlm21,nlm22,i, io1,io2,l1,l2
         sumt = sum([(sum(qpp(:,:,ilm1,ilm1,:)*taupp(:,:,ll(ilm1),:)),ilm1=1,nlma)]) !Pkl*Pkl
         sumq = sum([(sum(qpp(:,:,ilm1,ilm1,:)*sigpp(:,:,ll(ilm1),:)),ilm1=1,nlma)])
         sumh = sum(qpp(:,:,:,:,:)*ppipp(:,:,:,:,:))
         do  io2 = 1, norb;    if(blks(io2)==0) cycle 
            k2 = ktab(io2)
            nlm21 = ltab(io2)**2+1
            nlm22 = nlm21 + blks(io2)-1
            do  io1 = 1, norb; if(blks(io1)==0) cycle 
               associate(k1=>ktab(io1),nlm11 => ltab(io1)**2+1)
                 sumh=sumh+sum([( sum(qhh(k1,k2,ilm1,nlm21:nlm22,:)*ppihh(k1,k2,ilm1,nlm21:nlm22,:)),ilm1=nlm11,nlm11+blks(io1)-1)])
                 do  ilm1 = nlm11, nlm11+ blks(io1)-1
                    if( nlm21<= ilm1 .and. ilm1<=nlm21+blks(io2)-1) then !Hsm*Hsm
                       sumt = sumt + sum(qhh(k1,k2,ilm1,ilm1,:)*tauhh(k1,k2,ll(ilm1),:))
                       sumq = sumq + sum(qhh(k1,k2,ilm1,ilm1,:)*sighh(k1,k2,ll(ilm1),:))
                    endif
                 enddo
               endassociate
            enddo
         enddo
         do  io1 = 1, norb; if (blks(io1)==0) cycle !Hsm*Pkl
            associate(k1=>ktab(io1), nlm11=>ltab(io1)**2+1, nlm11e=>min(ltab(io1)**2+1 +blks(io1)-1,nlma) )
              sumh=sumh+sum([(sum(qhp(k1,:,ilm1,:,:)*ppihp(k1,:,ilm1,:,:)),ilm1= nlm11,nlm11+blks(io1)-1)])!site contribution to kinetic energy +potential
              sumt=sumt+sum([(sum(qhp(k1,:,ilm1,ilm1,:)*tauhp(k1,:,ll(ilm1),:)),ilm1= nlm11,nlm11e)]) !site contribution to kinetic energy
              sumq=sumq+sum([(sum(qhp(k1,:,ilm1,ilm1,:)*sighp(k1,:,ll(ilm1),:)),ilm1= nlm11,nlm11e)]) !site contribution to overlap (charge ?)
            endassociate
         enddo
       endblock pvgtknblock
       sraugm = sraugm + sumh - sumt
    enddo
    srhov   = srmesh + sraugm != n_out*Vin
    ekinval = sumev - srhov   != Eband - nout*Vin (V do not include SO term)
    if(ipr>= 30) write(stdo,"(/a)")' mkekin:'
    if(ipr>= 30) write(stdo,340) srmesh,sraugm,srhov,sumev,ekinval
340 format('   nout*Vin = smpart,onsite,total=:',3f14.6,/'    E_B(band energy sum)=',f12.6,'  E_B-nout*Vin=',f12.6)
    call tcx('mkekin')
  end subroutine mkekin
  subroutine makdos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax, ndos,dos)! Make density of states from bands
    !i  Input
    !i   nqp   :number of q-points
    !i   nband :number of bands
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   wgts  :band weights
    !i   evl   :band eigenvalues
    !i   n,w   :Methfessel-Paxton order and broadening parameters
    !i   tol   :(tol>0) allowed error in DOS due to truncating the gaussian,
    !i         :        to a finite energy range (number of bins)
    !i         :(tol<0) dimensionless energy window specifying truncation
    !i         :        of gaussian.  Energy window for which gaussian is
    !i         :        taken to be nonzero is set to -tol*w
    !i   emin, emax, ndos: energy range and number of energy mesh points
    !i   nbmx  :leading dimension of evl
    !o  Ouput
    !o    dos: density of states
    implicit none
    integer :: nqp,nband,nbmx,nsp,n,ndos
    double precision :: wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp),w,emin,emax,tol,wt,emesh
    integer :: i,isp,iband,iq,meshpt,mesh1,mesh2,mrange=999999,iprint
    double precision :: e,x,range,test,step,d,s,xx
    dos=0d0
    step = (emax - emin) / (ndos - 1)
    if( tol > 0d0 ) then
       do   i = 0, ndos-1
          x = i * step / w
          call delstp(0,x,test,s,xx)
          if ( test < tol ) then
             mrange = i + 1
             goto 3
          endif
       enddo
       if(iprint()> 30) print *,'makdos (warning) : tol too small'
3      continue
       range = 2 * mrange * step
       test = tol
    else
       range = -tol * w
       mrange = range / ( 2 * step )
       call delstp(0,-tol/2,test,s,xx)
    endif
    if (iprint() > 30) write (*,100) range/w,2*mrange,test
    iqloop: do  7  iq = 1, nqp
       wt = abs(wgts(iq)) / nsp
       ibandloop: do  61  iband = 1, nband
          isploop: do   isp = 1, nsp
             e = evl(iband,isp,iq)
             meshpt = (e - emin) / step
             mesh1 = meshpt - mrange
             mesh2 = meshpt + mrange
             if (mesh2 >= ndos) mesh2 = ndos-1
             if (mesh1 < 0) mesh1 = 0
             do   meshpt = mesh1, mesh2
                emesh = emin + meshpt * step
                x = (emesh - e) / w
                call delstp(n,x,d,s,xx)
                dos(meshpt,isp) = dos(meshpt,isp) + wt * d / w
             enddo
          enddo isploop
61     enddo ibandloop
7   enddo iqloop
100 format(/1x,'MAKDOS :  range of gaussians is ',f5.2,'W (',i4,' bins).'/11x,'Error estimate in DOS : ',1pe9.2,' per state.')
  end subroutine makdos
  subroutine phispinsym_ssite_set()
    use m_lmfinit,only: nbas,nsp,n0,ispec,lmxa_i=>lmxa
    use m_struc_def
    use m_lgunit,only:stdo
    use m_density,only: pnuall,pnzall
    implicit none
    integer:: ib,is,lmxa,l
    real(8):: pmean
    real(8),pointer:: pnu(:,:),pnz(:,:)
    if(master_mpi)write(6,*)' --phispinsym use spin-averaged potential for phi and phidot'
    do ib = 1,nbas
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       is   = ispec(ib) 
       lmxa = lmxa_i(is)
       do l=0,lmxa
          pmean = sum(pnu(l+1,1:nsp))/nsp
          if(master_mpi .AND. nsp==2) write(stdo,"('  ibas l=',2i3,' pnu=',2f10.5,' -->',f10.5)") &
               ib,l,pnu(l+1,1:nsp),pmean
          pnu(l+1,1:nsp) = pmean
          if(pnz(l+1,1)/=0 ) then
             pmean = sum(pnz(l+1,1:nsp))/nsp
             if(master_mpi .AND. nsp==2)write(stdo,"('  ibas l=',2i3,' pnz=',2f10.5,' -->',f10.5)") &
                  ib,l,pnz(l+1,1:nsp),pmean
             pnz(l+1,1:nsp) = pmean
          endif
       enddo
       pnuall(:,1:nsp,ib)=pnu(:,1:nsp)
       pnzall(:,1:nsp,ib)=pnz(:,1:nsp)
    enddo
  end subroutine phispinsym_ssite_set
  subroutine iorbtm()! Printout of orbital moments
    use m_lmfinit,only: ispec,nsp,nbas,slabl,lmxax
    use m_bandcal,only: orbtm=>orbtm_rv
    use m_lgunit,only:stdo
    !i Inputs
    !i   ics   :species table: class ic belongs to species ics(ic)
    !i   nl    :(global maximum l) + 1
    !i   nclass:number of inequivalent classes
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   orbtm :orbital moments
    implicit none
    integer :: isp,l,im,lm,m,l1,ipr,is,ibas,nl
    double precision :: amom,orbl(10)
    nl=lmxax+1
    call getpr(ipr)
    if (ipr < 20) return
    write(stdo,332)
332 format(/'IORBTM:  orbital moments :'/' ibas  Spec        spin   Moment decomposed by l ...')
    do  ibas = 1,nbas
       amom = 0
       do  isp = 1, nsp
          orbl(1:lmxax+1) = orbtm(1:lmxax+1,isp,ibas)
          amom = amom + sum(orbtm(1:lmxax+1,isp,ibas))
          write(stdo,"(i5,4x,a8,i6,8f12.6)") ibas,slabl(ibas),isp,(orbl(l1),l1=1,nl)
       enddo
       write(stdo,"(' total orbital moment',i4,':',f12.6)") ibas, amom
    enddo
  end subroutine iorbtm
endmodule m_bndfp_util
