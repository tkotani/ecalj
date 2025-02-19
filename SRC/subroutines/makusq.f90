module m_makusq !Accumulate coefficients (u,s,z) in all augmentation spheres for evec(:,iq,isp)
  use m_lmfinit,only: nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,lmxl_i=>lmxl,spec_a,kmxt_i=>kmxt,lmxb_i=>lmxb
  use m_ll,only:ll
  public makusq
  private
contains
  subroutine makusq(nsites,isite,nev,ispin,iq,q,evec, auszall)!Accumulate coefficients (u,s,z) in all augmentation spheres for evec(:,iq,isp)
    !note. For SO=1, ispin is neglected (all spin components are calculated simultaneously).
    use m_lmfinit,only: ispec,nbas,nlmax,nsp,nspc,nkapii,lhh,rsma,lso
    use m_igv2x,only: ndimh,nbandmx 
    use m_uspecb,only:uspecb
    use m_orbl,only: Orblib,ktab,ltab,offl,norb,blks
    use m_bstrux,only: bstrux_set,bstr
    implicit none
    intent(in)::    nsites,isite,nev,ispin,iq,q,evec
    intent(out)::                                   auszall
    integer,parameter:: n0=10,nkap0=3
    integer:: ispin,iq,nev,nsites,isite(nsites),ib,nkapi,is,nr,kmax,lmxa,lmxl,lmxh,i,nlma
    real(8):: q(3),eh(n0,nkap0),rsmh(n0,nkap0),a,rmt
    complex(8):: evec(ndimh,nspc,nev) !ndimhx = ndimh*nspc (Hamiltonian dimension). 
    complex(8),target:: auszall(nlmax,nbandmx,3,nsp,nsites,iq)
    auszall=0d0
    call tcn ('makusq')
    do  i = 1, nsites
       if (nsites == nbas) ib = i
       if (nsites /= nbas) ib = isite(i)
       is = ispec(ib)
       lmxa=lmxa_i(is)
       lmxl=lmxl_i(is)
       kmax=kmxt_i(is)
       nr=  nr_i(is)
       lmxh=lmxb_i(is)
       rmt =rmt_i(is)
       a =  spec_a(is)
       if (lmxa == -1) cycle
       call uspecb(is,rsmh,eh)
       nkapi= nkapii(is)
       nlma = (lmxa+1)**2
       SetupAllRadialheadANDtailFunctionsANDtheirBCs: block
         real(8):: rofi_rv(nr)
         real(8):: fh_rv(nr*(lmxh+1)*nkapi),   xh_rv(nr*(lmxh+1)*nkapi)
         real(8):: vh(0:lmxh,nkapi), dh(0:lmxh,nkapi)
         real(8):: fp_rv(nr*(lmxa+1)*(kmax+1)),xp_rv(nr*(lmxa+1)*(kmax+1)),vp(0:lmxa,0:kmax),dp(0:lmxa,0:kmax)
         call radmsh(rmt,a,nr,rofi_rv )
         call fradhd(nkapi,eh,rsmh,lhh(:,is),lmxh,nr,rofi_rv,fh_rv,xh_rv, vh,dh) !head part
         call fradpk(kmax,rsma(is),lmxa,nr,rofi_rv,fp_rv,xp_rv, vp,dp)           !tail part
         puqs11:block
           integer:: ivec,io1,l1,ik1,nlm11,nlm12,ilm1,i1,ilma,k,l,ispi,ispe,ispc,isp
           complex(8) ::cPkl(0:kmax,nlma)
           complex(8),pointer:: ausz(:,:,:,:)
           logical:: s12
           ausz => auszall(:,:,:,:,i,iq) !the coefficient for the projection onto (u,s,gz) for this site
           call bstrux_set(ib,q) !Get structure constant bstr
           ispi = merge(1,ispin,lso==1) !when lso==1 ispin is dummy, isp runs 1 and 2
           ispe = merge(2,ispin,lso==1)
           ispcloop: do isp = ispi,ispe ! For so=1, we do loop isp=1,2, 
              ispc = merge(isp,1,lso==1)
              Eigenfunctionloop: do ivec = 1, nev 
                 forall(ilma=1:nlma) cPkL(:,ilma)=matmul(evec(1:ndimh,ispc,ivec),bstr(1:ndimh,ilma,:))
                 call orblib(ib)     !Return norb,ltab,ktab,offl
                 ContributionFromHeadPart: do  io1 = 1, norb
                    l1  = ltab(io1)
                    ik1 = ktab(io1)
                    nlm11 = l1**2+1
                    nlm12 = nlm11 + blks(io1)-1
                    i1 = offl(io1)-nlm11+1 !  i1 = hamiltonian offset for first orbital in block
                    s12= merge(.true.,.false.,ik1 <= nkapi)
                    do  ilm1 = nlm11, nlm12
                       l = ll(ilm1)
                       if(s12)      ausz(ilm1,ivec,1:2,isp)= ausz(ilm1,ivec,1:2,isp)+ [vh(l,ik1),dh(l,ik1)]*evec(ilm1+i1,ispc,ivec)
                       if(.not.s12) ausz(ilm1,ivec,3,isp)  = ausz(ilm1,ivec,3,isp)  + evec(ilm1+i1,ispc,ivec)
                    enddo
                 enddo ContributionFromHeadPart
                 ContributionFromTailPart: do ilma = 1, nlma 
                    l = ll(ilma)
                    ausz(ilma,ivec,1:2,isp)=ausz(ilma,ivec,1:2,isp)+[sum(vp(l,:)*cPkL(:,ilma)),sum(dp(l,:)*cPkL(:,ilma))]
                 enddo ContributionFromTailPart
              enddo Eigenfunctionloop
           enddo ispcloop
         endblock puqs11
       endblock SetupAllRadialheadANDtailFunctionsANDtheirBCs
    enddo
    call tcx('makusq')
  end subroutine makusq
end module m_makusq
! Roughly examined at 2023jan ---------------------------------------------
!i Inputs
!i   nbas  :number of basis atoms
!i   isite :sites at which to calculate coefficients; see nsites
!i   nlmax :1st dimension of ausz (maximum nlma over all sites)
!i   ndham :dimensions ausz
!i   ndimh :dimensions evec
!i   napw  :number of G vectors in PW basis (gvlst2.f)
!i   nev   :number of eigenvectors for which to accumulate ausz
!i   nsp   :2 for spin-polarized case, otherwise 1
!i   nspc  :2 for so=1; otherwise 1
!i   isp   :spin channel
!i   iq    :qp index
!i   q     :Bloch vector
!i   evec  :eigenvectors for this q
!o Outputs
!o   ausz   :val,slo of w.f. at MT sphere surface added to ausz; see Remarks
!
!l Local variables
!l   ispc  :the current spin index in the coupled spins case.
!l         :Some quantities have no separate address space for each
!l         :spin in the indepedent-spins case (evec,evl,ewgt) but do
!l         :in the coupled-spins case.  A separate loop ispc=1..nspc
!l         :must be added for the latter case
!l         :ispc is the appropriate index for objects which distinguish
!l         :spins in the spin-coupled case only
!l   isp   :isp  is the appropriate index for objects which distinguish spins in the spin-uncoupled case only
!l   isp   :the current spin index in both independent and coupled cases.
!l         :isp is appropriate spin index for quantities that have
!l         :separate address space for each spin in every case
!l         :(potential- and density-like objects).
!r Remarks
!r   Makes coefficients for projection of wave function onto
!r   augmented functions (u,s,gz) which is valid inside the MT spheres.
!r   u and s are linear combinations of and phi,phidot defined as:
!r   u has val=1, slo=0 at rmax, s has val=0, slo=1
!r
!r   For example, for EELS matrix elements <nk|r|core> we will need
!r    |nk> = \sum_L(au_nkL*u_l*Y_L + as_nkL*s_l*Y_L + az_nkl*gz_l*Y_L)
!r
!r   These are generated from the potential later (see vcdmel)
!r   makusq returns the au_nkL and as_nkL az_nkl at one spin and k-pt for
!r   each of the sites in the unit cell.
!r   If nsites/=nbas, coeffs are made just for the nsites sites listed in isite.
!   ! ----------------------------------------------------------------------
!   !     In noncollinear case, isp=1 always => need internal ispc=1..2
!   !     isp is the current spin index in both cases:
!   !     isp = isp  in the collinear case
!   !         = ispc in the noncollinear case
!   !     whereas ispc=1 for independent spins, and spin index when nspc=2
