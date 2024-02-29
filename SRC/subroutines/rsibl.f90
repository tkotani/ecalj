module  m_rsibl
  use m_nvfortran,only:findloc
  public rsibl
  private
contains
  subroutine rsibl(lfrce,isp,q,iq,ndimh,nspc,napw,igapw,nevec,evec,ewgt,k1,k2,k3,smpot,smrho,f) !Add smooth part of output density into smrho and forces.
    use m_ll,only:ll
    use m_lmfinit,only: alat=>lat_alat,nspec,nbas,ispec,n0,nkap0
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol,plat=>lat_plat, pos=>rv_a_opos
    use m_supot,only: n1,n2,n3, lat_ng, gmax=>lat_gmax
    use m_uspecb,only: uspecb
    use m_shortn3,only:gvlst2
    use m_orbl,only: Orblib,ktab,ltab,offl,norb
    use m_sugcut,only:ngcut
    use m_hsibl,only: hsibl,hsibl1
    !i   lfrce :if nonzero, accumulate contribution to force
    !i   nbas  :size of basis
    !i   lfrce :1 calculate contribution to forces
    !i   nbas  :size of basis
    !i   q     :Bloch vector
    !i   iq    :index to current k-point
    !i   ndimh :dimension of hamiltonian
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   napw  :number of augmented PWs in basis
    !i   igapw :vector of APWs, in units of reciprocal lattice vectors
    !i   nevec :number of eigenvectors with nonzero weights
    !i   evec  :eigenvectors
    !i   ewgt  :eigenvector weights
    !i   k1..3 :dimensions smpot,smrho
    !i   smpot :smooth potential on uniform mesh, needed for forces
    
    !i   q     :Bloch wave number
    !i   ng    :number of G-vectors
    !i   gq    :2*pi/alat * (q+G) for all G-vectors
    !i   iv    :g-vectors as integer multiples of qlat (suphs0)
    !i   n1..3 :size uniform mesh for smooth density and potential
    !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
    !i   yl    :spherical harmonics for ng vectors
    !i   he    :table of energy factors
    !i   hr    :table of smoothing radius factors
    !i   vol   :cell volume
    !o   iprt  :index to which entry in rt a given orbital belongs
    !i   ipet  :index to which entry in etab a given orbital belongs
    !i   etab  :table of all inequivalent energies
    !i   rtab  :table of all inequivalent smoothing radii
    !i   ndimh :dimensions evec
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   ewgt  :weights for each of the trial fermi levels
    !i   nevec  :number of eigenstates to generate
    !i   evec  :eigenevectors
    !i   vspi  :potential * wave function, needed only for mode=1
    !o Outputs
    !o   smrho :smooth density accumulated for this qp
    !o   f     :force contribution accumulated for this qp
    !  psi :wave function 
    !   f  : forces
    implicit none
    integer,parameter:: nermx=100
    real(8),parameter:: pi = 4d0*datan(1d0),tpi = 2d0*pi
    integer :: lfrce,isp,k1,k2,k3,ndimh,nevec,iq,nspc, napw,igapw(3,napw),ng
    real(8):: q(3), ewgt(nevec) , f(3,nbas), wgt1,p(3),f0(3),xx(1),q0(3),etab(nermx),rtab(nermx)
    complex(8):: evec(ndimh,nspc,nevec),smrho(k1,k2,k3,*), smpot(k1,k2,k3,*),img=(0d0,1d0)
    integer :: nlmto,nrt , net , ltop , nlmtop , ogq , og2 , ohe , ohr , oiv , iprint
    integer :: iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx),ispc,is,ib,kb, ig,ixx(1),jg,i
    integer,allocatable :: iv_a_okv(:), ivp(:), igv(:,:)
    real(8),allocatable :: ogv(:,:),w_ogq(:,:),yl(:,:),w_og2(:),he(:,:),hr(:,:)
    complex(8),allocatable::psi(:,:,:),psir(:,:,:),vpsi(:,:,:), psi0(:,:,:,:),phase(:)
    if (nevec <= 0) return
    call tcn('rsibl')
    nlmto = ndimh-napw
    ng=lat_ng
    call pshpr(iprint()-30)
    call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,[0],000, 0,ng, ixx,xx,ixx)
    allocate(ogv(ng,3), igv(ng,3), iv_a_okv(ng*3))
    call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,[0],509,ng,ng, iv_a_okv,ogv,igv) ! q-dependent gv, kv, gv+q and iv   ! Note ogv has q already added to it!
    call poppr
    if(napw > 0) then
       allocate(ivp(napw) )
       do ig=1,napw
          ivp(ig)= findloc( [( sum(abs(igv(jg,:)-igapw(:,ig)))==0,jg=1,ng)], value=.true.,dim=1)
       enddo   
       if(any(ivp==0)) call rx('bug in rsibl.  Cannot find ivp')
    endif
    call tbhsi(nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop) ! --- Tables of energies, rsm, indices to them ---
    nlmtop = (ltop+1)**2
    allocate(w_ogq(ng,3),yl(ng,nlmtop), w_og2(ng), he(ng,net), hr(ng,nrt))
    q0=0d0
    if(nlmto>0) call hsibl1(net,etab,nrt,rtab,ltop,alat,q0,ng,ogv,w_ogq,w_og2, yl,he,hr) !H_L(G)=\frac{-4 pi}{e-G^2} {cal Y}_L(-iG) exp(gamma(e-G^2))  ! hsibl1 calculates he=1/(e-G^2) and hr=exp(-gamma G^2).
    deallocate(w_og2)
    allocate(psi(ng,nspc,nevec),vpsi(ng,nspc,nevec),psi0(ng,nspc,nevec,nbas),psir(k1,k2,k3),phase(ng))
    psi0=0d0
    rsibl1: block !    ! Make wave function for a block of evecs, or add contr. to forces
      integer:: blks(n0*nkap0),ntab(n0*nkap0),ncut(n0,nkap0),nkapi,jo,l2,lt,kp,ie,ir,ioff,nlm1,nlm2,ilm,io,l,ncutt
      real(8) :: e,rsm,eh(n0,nkap0),rsmh(n0,nkap0),fac
      do ib = 1, nbas
         is=ispec(ib) 
         ncut=ngcut(:,:,is)
         phase = [(exp(-img*tpi*sum(pos(:,ib)*(q+matmul(qlat,igv(ig,:))))),ig=1,ng)]
         call orblib(ib) !Return norb,ltab,ktab,offl
         call uspecb(is,rsmh,eh)
         call gtbsl1(1,norb,ltab,ktab,rsmh,eh,ntab,blks)
         do io = 1, norb
            if (blks(io) == 0) cycle
            jo = ntab(io)
            l2 = ltab(io)
            lt = ltab(jo)
            kp = ktab(io)
            ie = ipet(l2+1,kp,is)
            ir = iprt(l2+1,kp,is)
            ioff = offl(io)
            nlm1 = l2**2+1
            nlm2 = nlm1 + blks(io)-1
            rsm = rtab(ir)
            e   = etab(ie)
            ncutt=ncut(lt+1,kp) 
            fac = 4d0*pi*dexp(e*rsm*rsm*0.25d0)/vol
            do  ilm = nlm1, nlm2 ! ... Make vector evec*phase ! ... Combine G-dependent energy, rsm and YL factors
               l=ll(ilm)
               do i = 1, min(ng,ncutt)
                  psi0(i,1:nspc,1:nevec,ib) = psi0(i,1:nspc,1:nevec,ib)&
                       + he(i,ie)*hr(i,ir)* yl(i,ilm) *(0d0,-1d0)**(l+2)* fac*phase(i) *evec(ilm-nlm1+ioff+1,1:nspc,1:nevec) 
               enddo
            enddo
         enddo
      enddo
    endblock rsibl1
    psi(:,:,:) = sum(psi0(:,:,:,1:nbas),dim=4)
    if(napw>0) psi(ivp(:),:,1:nevec) = psi(ivp(:),:,1:nevec) + evec(nlmto+1:nlmto+napw,:,1:nevec)/vol**.5 !add PW(G) to psi
    ! psi= H(G) + PW(G) ! Add to real-space mesh, optionally make smpot*psi for forces
    do  ispc = 1, nspc !nspc=2 for SO=1, nspc=1 otherwise. nspx=nspc/nps
       do  i = 1, nevec
          call gvputf(ng,1,iv_a_okv,k1,k2,k3,psi(1,ispc,i),psir)
          call fftz3(psir,n1,n2,n3,k1,k2,k3,1,0,1)
          wgt1 = ewgt(i)
          smrho(:,:,:,isp+ispc-1)=smrho(:,:,:,isp+ispc-1)+ wgt1*dconjg(psir(:,:,:))*psir(:,:,:) !realspace density
          if (lfrce /= 0) then
             psir(:,:,:) = psir(:,:,:)*smpot(:,:,:,isp+ispc-1)
             call fftz3(psir,n1,n2,n3,k1,k2,k3,1,0,-1)
             call gvgetf(ng,1,iv_a_okv,k1,k2,k3,psir, vpsi(1,ispc,i))
          endif 
       enddo
    enddo
    AddForce: if (lfrce /= 0) then
       do ib = 1, nbas
          is=ispec(ib) 
          phase = exp(-img*tpi*sum(p*q)) * exp(-img*tpi*matmul(pos(:,ib), matmul(qlat, transpose(igv))))
          do i = 1, nevec
             f0(:)=2d0*vol*matmul( -sum(dimag(dconjg(vpsi(:,:,i))*psi0(:,:,i,ib)),dim=2),  w_ogq(:,:))
             f(:,ib) = f(:,ib)   + ewgt(i)*f0(:) !Add 2*Re( (v psi+) grad(psi) ) to f
             do kb = 1, nbas !!             This shouldn't be necessary?
                f(:,kb) = f(:,kb) - ewgt(i)*f0(:)/nbas
             enddo
          enddo
       enddo
    endif AddForce
    call tcx('rsibl')
  end subroutine rsibl
end module m_rsibl
