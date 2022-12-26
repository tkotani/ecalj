module m_potpus
  public potpus
  private
contains
  subroutine potpus(z,rmax,lmxa,v,vdif,a,nr,nsp,lso,rofi,pnu,pnz, &
       ehl,rsml,rs3,vmtz,nab,n0,ppnl,hab,vab,sab,sodb)
    use m_lmfinit,only: stdo,lrel,cc
    use m_ftox
    use m_vxtrap,only: vxtrap,rwftai
    !- Potential parameters for potential and boundary condition
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   z     :nuclear charge
    !i   rmax  :augmentation radius, in a.u.
    !i   lmxa  :augmentation L-cutoff
    !i   lso   :if nonzero calculate radial spin-orbit integrals
    !i   v     :spherical potential used to construct wave functions
    !i   vdif  :perturbation potential added to v in matrix elements
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   rofi  :radial mesh points
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   pnz   :boundary conditions for second p.q.n (local orbital).
    !i         :10s digit controls how local orbital included in hamiltonian;
    !i         :See Remarks
    !i         :The following four parameters are used to extrapolate
    !i         :quantities outside the MT radius.  Extrapolation is used to
    !i         :correct matrix elements for local orbitals.
    !i   ehl   :energy of smoothed Hankel tail for local orbital
    !i   rsml  :smoothing radius for smoothed Hankel tail of local orbital
    !i   rs3   :minimum smoothing radius for extrapolation of MT potential
    !i   vmtz  :muffin-tin zero: subtracted from V in the fitting procedure.
    !i         :The asymptotic form of V-vmtz is taken to be zero.
    !i   nab,n0:first and second dimensions of hab, vab, sab
    !o Outputs
    !o   ppnl  :NMTO pot pars (but no backward extrapolation; <phi phi>=1)
    !ox         :(1)  = inverse potential function (not created)
    !o         :(2)  = normalization of phi (= 1 here for now)
    !   pvalue(:), pslope(:), ppovl(:,:)
    !o         :(3)  = rmax * log derivative of phi at rmax      
    !o         :(4)  = rmax * log derivative of phidot at rmax
    !o         :(5)  = phi at rmax
    !o         :(6)  = phidot at rmax
    !o         :(7)  = normalization of phidot
    !o         :(8)  = normalization <gz|gz> for local orbital
    !o         :(9)  = overlap <phi|gz>
    !o         :(10) = overlap <phidot|gz>
    !o         :(11) = phz at rmax for local orbitals defined, before (phi,phidot) subtracted
    !o         :(12) = dphz (slope of phz) at rmax for local orbitals defined
    !o   hab   :matrix elements of the hamiltonian (spher. v) for this site
    !o         :specified by the potentials v, vdif and the boundary
    !o         :condition pnu (or pnz).  See Remarks.
    !o   vab   :matrix elts of the spher. pot, w/ the same structure as hab
    !o   sab   :matrix elts of the overlap w/ the same structure as hab
    !o   sodb  :radial integrals <psi1_s|SO|psi2_s>, where psi1,psi2
    !o         :are one of the ul,sl,gz and s is the spin index.
    !o         :SO = 2/(c^2) dV/dr*(1/r), V(r)=-2*z/r+0.5*(v(r)_s1 + v(r)_s2)
    !o         :SO = 1/(2*m^2*c^2)*(dV/dr*1/r), m=.5, c=274 (at. Rydberg units)
    !o         :Notice that for now only spin-diagonal matrix
    !o         :elements are evaluated, since we plan to calculate
    !o         :only H_so = so*sz*lz   ===>SO calcualtion aughsoc
    !l Local variables
    !l   lpzi  :flags how local orbitals is to be treated in current channel
    !l         :0 no local orbital gz
    !l         :1 value and slope of gz constructed to be zero at rmax
    !l         :  by admixture of phi,phidot
    !l         :2 a smooth Hankel tail is attached (extended local orbital)
    !l         :  and it is included explicitly in the basis.
    !l         :3 a smooth Hankel tail is attached (extended local orbital)
    !l         :  and is coupled to the valence states in an extended atom
    !l         :  approximation.
    !r Remarks
    !r   This routine makes linear combinations of phi,phidot and the
    !r   matrix elements of them in the supplied spherical potential.
    !r   Linear combinations (u,s) of phi,phidot are defined as :
    !r     u has val=1, slo=1 at rmax;   s has val=0, slo=1
    !r   NB: potpus actually works with ul=r*u and sl=r*s respectively.
    !r
    !r   There can additionally be local orbitals, of the following types,
    !r   as specified by pnz.
    !r      pnz = 0      : no local orbital
    !r      10 > pnz > 0 : local orbital of the first type
    !r      pnz > 10     : local orbital of the second type
    !r   From the boundary conditions (1s digit+fractional part of pnz),
    !r   wave function phi_z can be generated for r<rmax.
    !r   A local orbital of the first type is defined as follows.
    !r      gz = r * ( phi_z - phi_z(rmax) u - (phi_z)'(rmax) s )
    !r   By construction, gz/r has both value = 0 and slope = 0 at rmax.
    !r   A local orbital of the second type is defined as gz=r*phi_z;
    !r   for r>rmax a smooth Hankel tail (spec'd by ehl,rsml) is attached.
    !r
    !r  *Documentation of hab,vab,sab generated by potpus.
    !r   hab,vab,sab are matrix elements of the true wave function (including
    !r   the small component). in the spherical part of the potential V.
    !r   Matrix elements of vdif are included perturbatively.
    !r   Array sodb holds matrix elements for L.S coupling.
    !r
    !r   The leading dimension of the (hab,vab,sab,sodb) arrays correspond
    !r   to the following matrix elements:
    !r      (1)=<ul|(h,v,1)|ul>   (5)=<ul|(h,v,1)|gz>   (8)=<gz|(h,v,1)|ul>
    !r      (2)=<ul|(h,v,1)|sl>   (6)=<sl|(h,v,1)|gz>   (9)=<gz|(h,v,1)|sl>
    !r      (3)=<sl|(h,v,1)|ul>   (7)=<gz|(h,v,1)|gz>
    !r      (4)=<sl|(h,v,1)|sl>
    !r   hab(2) and hab(3) are related by a Wronskian, as are hab(5) and
    !r   hab(7), and hab(6) and hab(8); but both are kept for convenience,
    !r   because of the nonhermiticity in h.
    implicit none
    integer :: lmxa,lso,nr,nsp,n0,nab,nppn
    parameter (nppn=12)
    double precision :: z,rmax,a,rofi(1),v(nr,nsp),ehl(n0),rsml(n0), &
         pnu(n0,nsp),pnz(n0,nsp),ppnl(nppn,n0,2), & !,pp(n0,2,5)
         hab(3,3,n0,nsp),sab(3,3,n0,nsp),vab(3,3,n0,nsp),vdif(nr,nsp), &
         sodb(3,3,n0,nsp,2),rs3,vmtz
    integer :: ipr,ir,i,j,k,l,lpz,lpzi(0:n0),nrbig
    double precision :: m21,m11,m22,m12,d00,d10,d11,det,vi, &
         h00,h01,h10,h11, &
         phmins,phplus,q,r,s00,s01,s10,s11,tmc, &
         umegam,umegap,v00,v01,v10,v11,vl,xx,xxx,yyy,zzz, &
         b,ghg,ghgp,gphgp !enu,c,srdel,qpar,ppar,
    
!    real(8):: ,gf11(nr),gf22(nr),gf12(nr)
!    real(8):: gf11,gf22,gf12
    double precision :: d0z,d1z,dzz,v0z,v1z,s0z,s1z,vz0,vz1,sz0, &
         sz1,h0z,hz0,h1z,hz1,suz,ssz,szu,szs,szz, &
         vuz,vsz,vzu,vzs,vzz,huz,hsz,hzu,hzs,hzz
    double precision :: g(nr,2),gp(nr,2*4),gz(nr,2)
    double precision :: ev,phi,dphi,phip,dphip,p, dlphi,dlphip
    double precision :: ez,phz,dphz,phzp,dphzp,pz,dlphz,dlphzp,phz2,dphz2
    !     double precision ul,sl,D
    integer :: nrx
    parameter (nrx=1501)
    double precision :: rwgtx(nrx) !,gzbig(nrx,2)
    real(8),allocatable:: gzbig(:,:)
    !     double precision gbig(nrx*2),gpbig(nrx*8),vbig(nrx,2)
    !     Spin-Orbit related parameters
    integer :: moda(0:lmxa)
    double precision :: vavg(nr),dv(nr),sop(0:lmxa,nsp,nsp,3), &
         sopz(0:lmxa,nsp,nsp,3), &
         psi(nr,0:lmxa,nsp),dpsi(nr,0:lmxa,nsp), &
         pzi(nr,0:lmxa,nsp),ezum(0:8,nsp), &
         enumx(0:8,nsp),wrk(nr,4),m(2,2,0:lmxa,nsp), &
         x21,x11,x22,x12,vx00,vx01,vx10,vx11, &
         vx0z,vxz0,vx1z,vxz1,vxzz,vxzu,vxuz,vxzs,vxsz, xxxx(1,1)
    real(8):: rwgt(nr)
    call getpr(ipr)
    call vxtrap(1,z,rmax,lmxa,v,a,nr,nsp,pnz,rs3,vmtz,nrx,lpz,nrbig,rofi,rwgtx,xxxx(1,1))
    rwgt= rwgtx(1:nr)
    b   = rmax/(dexp(a*nr - a) - 1d0)
    ! ... Gradient of average v (for spin-orbit)
    if (lso /= 0) then
       if (lrel == 0) call rx('spin-orbit requires REL=1')
       if (nsp == 1)  call rx('spin-orbit requires NSPIN=2')
       vavg = .5d0*(v(:,1)+v(:,2))
       call radgra(a,b,nr,rofi,vavg,dv)
    endif
    hab=0d0
    vab=0d0
    sab=0d0
    if(lso/=0) sodb=0d0
    isploop: do 80  i = 1, nsp
       if(ipr>=40)             write(stdo,ftox)' potpus spin=',i,'pnu=',ftof(pnu(1:lmxa+1,i),3)
       if(ipr>=40.and.lpz/= 0) write(stdo,ftox)' pnz=',ftof(pnz(1:lmxa+1,i),3)
       lloop: do  10  l = 0, lmxa
          k = l+1
          lpzi(l) = 0
          if (pnz(k,i) >  0)  lpzi(l) = 1
          if (pnz(k,i) >= 10) lpzi(l) = 2
          moda(l) = 5
          if (lpzi(l)/=0) then ! ... lo wf gz and its sphere boundary parameters
             moda(l) = 6 
             call makrwf(10,z,rmax,l,v(1,i),a,nr,rofi,pnz(1,i),4,gz,gp,ez,phz,dphz,phzp,dphzp,pz)
             if (lpzi(l)==2) then ! Extend local orbital to large mesh; match gz to envelope
                allocate(gzbig(nrbig,2))
                gzbig(1:nr,:) = gz(1:nr,:)
                call rwftai(5,rmax,a,nr,nrbig,rofi,phz,dphz,xx,l,ehl(k), rsml(k),gzbig)
                if (gzbig(nr,1) /= gz(nr,1)) then !   If rwftai scales gzbig, rescale phz,gz
                   xx = gzbig(nr,1)/gz(nr,1)
                   phz  = phz*xx
                   dphz = dphz*xx
                   phzp = phzp*xx
                   dphzp= dphzp*xx
                   gz   = gz*xx !call dscal(nr,xx,gz(1,1),1) call dscal(nr,xx,gz(1,2),1)
                endif
             endif
             if(lso/=0) ezum(l,i) = ez     !for SO
             if(lso/=0) pzi(:,l,i)= gz(:,1)!for SO
          else
             phz = 0
             dphz = 0
          endif
          ! ... Valence wf g,gp, and their sphere boundary parameters
          call makrwf(10,z,rmax,l,v(1,i),a,nr,rofi,pnu(1,i),2,g,gp,ev,phi,dphi,phip,dphip,p)
          ghg    = ev   ! <g H g> = e <g g> = e
          ghgp   = 1d0  ! <g H gp> = <g (H-e) gp> + e <g gp> = <g g> = 1
          gphgp  = ev*p ! <gp H gp> = <gp (H-e) gp> + e <gp gp> = <gp g> + e p = ep
          dlphi  = rmax*dphi/phi
          dlphip = rmax*dphip/phip
          ! xxx Integrate g and gp outward on extended mesh --->removed here at 2022-dec-25
          if(lso /= 0) then ! ... Keep local copies of phi and phidot for SO coupling
             psi(:,l,i) = g(:,1)
             dpsi(:,l,i)= gp(:,1)
             enumx(l,i) = ev
          endif
          ppnl(:,k,i) = 0d0 ! potential parameters, no backwards integration, <phi phi>=1
          ppnl(1,k,i) = 0
          ppnl(2,k,i) = 1d0
          ppnl(3,k,i) = rmax * dlphi
          ppnl(4,k,i) = rmax * dlphip
          ppnl(5,k,i) = phi
          ppnl(6,k,i) = phip
          ppnl(7,k,i) = p
          ppnl(11,k,i) = phz
          ppnl(12,k,i) = dphz
          intg: block ! --- Integrals of w.f. products with spherical potential ---
          ! ... This branch computes integrals with products of (g,gp,gz)
          !     Convention: 00 (phi,phi) 10 (dot,phi) 11 (dot,dot)
          !     Convention: 0z (phi,sc) 1z (dot,sc) zz (sc,sc)
            integer:: fllp1
            real(8):: vii(nr),tmcc(nr),gf11(nr),gf22(nr),gf12(nr),xxxw(nr),yyyw(nr),zzzw(nr)
            fllp1 = l*(l+1)
            vii(1)=0d0
            vii(2:nr) = v(2:nr,i) - 2d0*z/rofi(2:nr)
            tmcc = cc - (vii-ev)/cc
            gf11(1)=0d0
            gf11(2:nr) = 1d0 + fllp1/(tmcc(2:nr)*rofi(2:nr))**2
            d00= sum(rwgt*vdif(:,i)* (gf11* g(:,1)* g(:,1)+ g(:,2)*g(:,2)) ) 
            d10= sum(rwgt*vdif(:,i)* (gf11*gp(:,1)* g(:,1)+gp(:,2)*g(:,2)) ) 
            d11= sum(rwgt*vdif(:,i)* (gf11*gp(:,1)*gp(:,1)+gp(:,2)*gp(:,2)))
            v00= sum(rwgt*vii*       (gf11* g(:,1)* g(:,1)+ g(:,2)*g(:,2)) ) 
            v10= sum(rwgt*vii*       (gf11*gp(:,1)* g(:,1)+gp(:,2)*g(:,2)) )
            v11= sum(rwgt*vii*       (gf11*gp(:,1)*gp(:,1)+gp(:,2)*gp(:,2)))
            if(lpzi(l)/=0) then !computes integrals with products of (g,gp) x gz
               gf12(1)=0d0
               gf22(1)=0d0
               tmcc= cc - (vii-ez)/cc
               gf22(2:nr) = 1d0 + fllp1/(tmcc(2:nr)*rofi(2:nr))**2
               gf12 = (gf11 + gf22)/2
               d0z = sum(rwgt*vdif(:,i)*(gf12*g(:,1) *gz(:,1)+ g(:,2) *gz(:,2)))
               d1z = sum(rwgt*vdif(:,i)*(gf12*gp(:,1)*gz(:,1)+ gp(:,2)*gz(:,2)))
               dzz = sum(rwgt*vdif(:,i)*(gf22*gz(:,1)*gz(:,1)+ gz(:,2)*gz(:,2)))
               v0z = sum(rwgt*vii*(gf12*g(:,1) *gz(:,1)+ g(:,2) *gz(:,2)))
               v1z = sum(rwgt*vii*(gf12*gp(:,1)*gz(:,1)+ gp(:,2)*gz(:,2)))
               vzz = sum(rwgt*vii*(gf22*gz(:,1)*gz(:,1)+ gz(:,2)*gz(:,2)))
               s0z = sum(rwgt*(gf12*g(:,1) *gz(:,1)+ g(:,2) *gz(:,2)))
               s1z = sum(rwgt*(gf12*gp(:,1)*gz(:,1)+ gp(:,2)*gz(:,2)))
               szz = sum(rwgt*(gf22*gz(:,1)*gz(:,1)+ gz(:,2)*gz(:,2)))
            endif
          endblock intg
          v00 = v00 + d00
          v10 = v10 + d10
          v11 = v11 + d11
          v01 = v10
          s00 = 1d0
          s10 = 0d0
          s01 = 0d0
          s11 = p           ! h+d
          h00 = ghg   + d00 ! h00 = <g H g> = e <g g> = e                       
          h01 = ghgp  + d10 ! h01 = <g H gp> = <g (H-e) gp> + e <g gp> = <g g>  
          h10 = 0d0   + d10 ! h10 = <gp H g> = 0d0
          h11 = gphgp + d11 ! h11 = <gp H gp> = <gp (H-e) gp> + e <gp gp> = <gp g> + e p = ep     
          !     Should not be needed since Wronskian explicit in makrwf
          call pvpus1(rmax,phi,dphi,phip,dphip,h01,h10)
          if (lpzi(l) /= 0) then
             v0z = v0z + d0z
             v1z = v1z + d1z
             vzz = vzz + dzz
             vz0 = v0z
             vz1 = v1z
             sz0 = s0z
             sz1 = s1z
             h0z = ez*s0z + d0z!? h0z = ez*s0z + d0z - ev*(phz*m11 + dphz*m21) ! if gz has u,s subt.
             hz0 = ev*s0z + d0z
             h1z = ez*s1z + d1z!? h1z = ez*s1z + d1z - ev*(phz*m12 + dphz*m22) ! if gz has u,s subt.
             hz1 = ev*sz1 + sz0 + d1z
             hzz = ez*szz + dzz
             call pvpus1(rmax,phi,dphi,phz,dphz,h0z,hz0) !       Put in Wronskians explicitly
             call pvpus1(rmax,phip,dphip,phz,dphz,h1z,hz1)
          endif

          ! --- Integrals of u-s products from phi,phidot products ---
          !     Linear transformation between (u,s) and (phi,phidot)

          !       ( u(r) )        ( phi(r)    )
          !       (      ) =    M (           )
          !       ( s(r) )        ( phidot(r) )

          !     Conditions u=1,s=0 and u'=0,s'=1 at rmt = >

          !       (1  0)      (phi     dphi )                    ( dphip  -dphi )
          !       (    )  = M (             )  = >  M = (det)^-1 (              )
          !       (0  1)      (phip    dphip)                    (-phip    phi  )

          !     where det = phi*dphip - phip*dphi

          !     To compute matrix elements of (u,s) from elements (phi,phidot),
          !     let u_i = u or s and phi_i = phi or phidot; i=1 or 2
          !     <u_i|u_j> = <(M phi)_i|(M phi)_j> = sum_lm M_il <phi_l|phi_m> M_mj
          det = phi*dphip - dphi*phip
          m11 = dphip/det
          m12 = -dphi/det
          m21 = -phip/det
          m22 = phi/det
          ! ... Keep local copies of mij for SO coupling
          if (lso /= 0) then
             m(1,1,l,i) = m11
             m(1,2,l,i) = m12
             m(2,1,l,i) = m21
             m(2,2,l,i) = m22
          endif
          !        block
          !       real(8)::mm(2,2)
          !       mm=reshape([m11,m21,m12,m22],[2,2])
          !       endblock

          ! ... hab(1)=huu   hab(2)=hus   hab(3)=hsu   hab(4)=hss
          hab(1,1,k,i) = m11*h00*m11 + m11*h01*m12 + m12*h10*m11 + m12*h11*m12
          hab(2,1,k,i) = m11*h00*m21 + m11*h01*m22 + m12*h10*m21 + m12*h11*m22
          hab(1,2,k,i) = m21*h00*m11 + m21*h01*m12 + m22*h10*m11 + m22*h11*m12
          hab(2,2,k,i) = m21*h00*m21 + m21*h01*m22 + m22*h10*m21 + m22*h11*m22
          !     Next 2 lines put in Wronskian explicitly
          !     Should no longer be needed: Wronskian explicit in makrwf.
          hab(2,1,k,i) = 0.5d0 * (hab(2,1,k,i)+hab(1,2,k,i)-rmax**2)
          hab(1,2,k,i) = hab(2,1,k,i)+rmax**2

          sab(1,1,k,i) = m11*s00*m11 + m11*s01*m12 + m12*s10*m11 + m12*s11*m12
          sab(2,1,k,i) = m11*s00*m21 + m11*s01*m22 + m12*s10*m21 + m12*s11*m22
          sab(1,2,k,i) = m21*s00*m11 + m21*s01*m12 + m22*s10*m11 + m22*s11*m12
          sab(2,2,k,i) = m21*s00*m21 + m21*s01*m22 + m22*s10*m21 + m22*s11*m22

          vab(1,1,k,i) = m11*v00*m11 + m11*v01*m12 + m12*v10*m11 + m12*v11*m12
          vab(2,1,k,i) = m11*v00*m21 + m11*v01*m22 + m12*v10*m21 + m12*v11*m22
          vab(1,2,k,i) = m21*v00*m11 + m21*v01*m12 + m22*v10*m11 + m22*v11*m12
          vab(2,2,k,i) = m21*v00*m21 + m21*v01*m22 + m22*v10*m21 + m22*v11*m22

          ! --- Integrals of transformed gz with (new gz, u, s) ---
          !     New gz = (gz0 - gz0(rmax) u - r*(gz0/r)'(rmax) s)
          !            = (gz0 - phz u - dphz s)
          !     To compute <u or s | gz> = M < phi or phidot | gz> :
          !     let u_i = u or s and phi_i = phi or phidot; i=1 or 2
          !     <u_i|gz> = <(M phi)_i | gz> = sum_km M_ik <phi_k|gz>
          !     <u_i|gz> = <(M phi)_i | gz> = <(M phi)_i | (gz0 - phz u - dphz s)>
          !              = sum_km M_ik <phi_k|gz> - phz <u_i|u_1> - dphz <u_i|u_2>


          ! ... Setup for transformation on gz -> local orbital
          !     At this point, phz, dphz are amount of gz at rmt
          !     Reuse phz,dphz to project amount of (u,s) onto gz.  Projection
          !     only applies when local orbital is a true local orbital.
          if (lpzi(l) == 2) then
             phz = 0
             dphz = 0
          endif
          if (lpzi(l) /= 0) then
             suz = m11*s0z + m12*s1z
             ssz = m21*s0z + m22*s1z
             szu = sz0*m11 + sz1*m12
             szs = sz0*m21 + sz1*m22
             szz = szz - phz*(suz+szu) - dphz*(ssz+szs) + phz**2*sab(1,1,k,i) + &
                  phz*dphz*(sab(2,1,k,i)+sab(1,2,k,i)) + dphz**2*sab(2,2,k,i)
             suz = suz - phz*sab(1,1,k,i) - dphz*sab(2,1,k,i)
             ssz = ssz - phz*sab(1,2,k,i) - dphz*sab(2,2,k,i)
             szu = szu - phz*sab(1,1,k,i) - dphz*sab(1,2,k,i)
             szs = szs - phz*sab(2,1,k,i) - dphz*sab(2,2,k,i)

             vuz = m11*v0z + m12*v1z
             vsz = m21*v0z + m22*v1z
             vzu = vz0*m11 + vz1*m12
             vzs = vz0*m21 + vz1*m22
             vzz = vzz - phz*(vuz+vzu) - dphz*(vsz+vzs) + phz**2*vab(1,1,k,i) + &
                  phz*dphz*(vab(2,1,k,i)+vab(1,2,k,i)) + dphz**2*vab(2,2,k,i)
             vuz = vuz - phz*vab(1,1,k,i) - dphz*vab(2,1,k,i)
             vsz = vsz - phz*vab(1,2,k,i) - dphz*vab(2,2,k,i)
             vzu = vzu - phz*vab(1,1,k,i) - dphz*vab(1,2,k,i)
             vzs = vzs - phz*vab(2,1,k,i) - dphz*vab(2,2,k,i)

             huz = m11*h0z + m12*h1z
             hsz = m21*h0z + m22*h1z
             hzu = hz0*m11 + hz1*m12
             hzs = hz0*m21 + hz1*m22
             hzz = hzz - phz*(huz+hzu) - dphz*(hsz+hzs) + phz**2*hab(1,1,k,i) + &
                  phz*dphz*(hab(2,1,k,i)+hab(1,2,k,i)) + dphz**2*hab(2,2,k,i)
             huz = huz - phz*hab(1,1,k,i) - dphz*hab(2,1,k,i)
             hsz = hsz - phz*hab(1,2,k,i) - dphz*hab(2,2,k,i)
             hzu = hzu - phz*hab(1,1,k,i) - dphz*hab(1,2,k,i)
             hzs = hzs - phz*hab(2,1,k,i) - dphz*hab(2,2,k,i)
             ! ... New gz val,slo=0 => hamiltonian is hermitian
             if (lpzi(l) == 1) then
                hzu = (hzu + huz)/2
                huz = hzu
                hzs = (hzs + hsz)/2
                hsz = hzs
             endif
             ! ... hab(5)=uz    hab(2)=sz    hab(7)=zz
             !     print *, 'zero out potpus local orbitals'
             sab(1,3,k,i) = suz !5
             sab(2,3,k,i) = ssz !6
             sab(3,3,k,i) = szz !7
             sab(3,1,k,i) = szu !8
             sab(3,2,k,i) = szs !9

             vab(1,3,k,i) = vuz
             vab(2,3,k,i) = vsz
             vab(3,3,k,i) = vzz
             vab(3,1,k,i) = vzu
             vab(3,2,k,i) = vzs

             hab(1,3,k,i) = huz
             hab(2,3,k,i) = hsz
             hab(3,3,k,i) = hzz
             hab(3,1,k,i) = hzu
             hab(3,2,k,i) = hzs
             ! ... NMTO potential parameters for local orbitals
             !     Note that (s0z,s1z) = <(phi,phidot)|gz0>. We need <(phi,phidot)|gz>.
             !     Let i=1 or 2 and define s_iz = (<phi|gz>,<phidot|gz>) for i=1,2

             !     (s_1z)     (phi_1) |
             !     (    ) = < (     ) | gz0 - phz u - dphz s>
             !     (s_2z)     (phi_2) |

             !                (phi_1) |             (u_1) |
             !            = < (     ) | gz0 - M^-1< (   ) | phz u_1 - dphz u_2>
             !                (phi_2) |             (u_2) |

             !              The first term are matrix elements (s0z,s1z)
             ppnl(8,k,i)  = szz
             xxx = (phz*sab(1,1,k,i) + dphz*sab(2,1,k,i))
             yyy = (phz*sab(1,2,k,i) + dphz*sab(2,2,k,i))
             ppnl(9,k,i)  = s0z - phi *xxx - dphi *yyy
             ppnl(10,k,i) = s1z - phip*xxx - dphip*yyy
          endif
10     enddo lloop
80  enddo isploop
    ! ... Calculate spin-orbit parameters
    if (lso /= 0) then
       call soprm(5,lpzi,psi,dpsi,pzi,nr,nsp,lmxa,lmxa,v,dv,enumx, ezum,z,rofi,rwgt,wrk,sop,sopz)
       !   ... Make the spin diagonal radial integrals
       do i = 1, nsp
          do l = 0, lmxa
             phz = ppnl(11,l+1,i)
             dphz = ppnl(12,l+1,i)
             k = l + 1
             m11 = m(1,1,l,i)
             m12 = m(1,2,l,i)
             m21 = m(2,1,l,i)
             m22 = m(2,2,l,i)
             v00 = sop(l,i,i,1)
             v01 = sop(l,i,i,2)
             v10 = v01
             v11 = sop(l,i,i,3)
             sodb(1,1,k,i,1) =m11*v00*m11+m11*v01*m12+m12*v10*m11+m12*v11*m12
             sodb(2,1,k,i,1) =m11*v00*m21+m11*v01*m22+m12*v10*m21+m12*v11*m22
             sodb(1,2,k,i,1) =m21*v00*m11+m21*v01*m12+m22*v10*m11+m22*v11*m12
             sodb(2,2,k,i,1) =m21*v00*m21+m21*v01*m22+m22*v10*m21+m22*v11*m22
             !  ...  Make the local orbitals spin diagonal integrals
             if (moda(l) == 6) then
                vzz = sopz(l,i,i,1)
                v0z = sopz(l,i,i,2)
                vz0 = v0z
                v1z = sopz(l,i,i,3)
                vz1 = v1z
                vzu = vz0*m11 + vz1*m12
                vuz = m11*v0z + m12*v1z
                vzs = vz0*m21 + vz1*m22
                vsz = m21*v0z + m22*v1z
                vzz = vzz - phz*(vuz+vzu) - dphz*(vsz+vzs) + &
                     phz**2*sodb(1,1,k,i,1) + &
                     phz*dphz*(sodb(2,1,k,i,1)+sodb(1,2,k,i,1)) + &
                     dphz**2*sodb(2,2,k,i,1)
                vuz = vuz - phz*sodb(1,1,k,i,1) - dphz*sodb(2,1,k,i,1)
                vsz = vsz - phz*sodb(1,2,k,i,1) - dphz*sodb(2,2,k,i,1)
                vzu = vzu - phz*sodb(1,1,k,i,1) - dphz*sodb(1,2,k,i,1)
                vzs = vzs - phz*sodb(2,1,k,i,1) - dphz*sodb(2,2,k,i,1)
                sodb(1,3,k,i,1) = vuz
                sodb(2,3,k,i,1) = vsz
                sodb(3,3,k,i,1) = vzz
                sodb(3,1,k,i,1) = vzu
                sodb(3,2,k,i,1) = vzs
             endif
          enddo
       enddo
       !   ... Make the spin off-diagonal radial integrals
       do l = 0, lmxa
          k = l + 1
          phz   = ppnl(11,k,1)
          dphz  = ppnl(12,k,1)
          phz2  = ppnl(11,k,2)
          dphz2 = ppnl(12,k,2)
          m11 = m(1,1,l,1)
          m12 = m(1,2,l,1)
          m21 = m(2,1,l,1)
          m22 = m(2,2,l,1)
          x11 = m(1,1,l,2)
          x12 = m(1,2,l,2)
          x21 = m(2,1,l,2)
          x22 = m(2,2,l,2)
          v00 = sop(l,1,2,1)
          v01 = sop(l,1,2,2)
          v10 = v01
          v11 = sop(l,1,2,3)
          vx00 = sop(l,2,1,1)
          vx01 = sop(l,2,1,2)
          vx10 = vx01
          vx11 = sop(l,2,1,3)
          !     ... up-down block
          sodb(1,1,k,1,2)=m11*v00*x11+m11*v01*x12+m12*vx10*x11+m12*v11*x12
          sodb(2,1,k,1,2)=m11*v00*x21+m11*v01*x22+m12*vx10*x21+m12*v11*x22
          sodb(1,2,k,1,2)=m21*v00*x11+m21*v01*x12+m22*vx10*x11+m22*v11*x12
          sodb(2,2,k,1,2)=m21*v00*x21+m21*v01*x22+m22*vx10*x21+m22*v11*x22
          !     ... down-up block
          sodb(1,1,k,2,2) = x11*vx00*m11+x11*v01*m12+x12*vx10*m11 +x12*vx11*m12
          sodb(2,1,k,2,2) = x11*vx00*m21+x11*v01*m22+x12*vx10*m21   +x12*vx11*m22
          sodb(1,2,k,2,2) = x21*vx00*m11+x21*v01*m12+x22*vx10*m11   +x22*vx11*m12
          sodb(2,2,k,2,2) = x21*vx00*m21+x21*v01*m22+x22*vx10*m21   +x22*vx11*m22
          !  ...  Make the local orbitals spin off-diagonal radial integrals
          if (moda(l) == 6) then
             vzz = sopz(l,1,2,1)
             v0z = sopz(l,1,2,2)
             vz0 = v0z
             v1z = sopz(l,1,2,3)
             vz1 = v1z
             vxzz = sopz(l,2,1,1)
             vx0z = sopz(l,2,1,2)
             vxz0 = vx0z
             vx1z = sopz(l,2,1,3)
             vxz1 = vx1z
             vzu = vz0*x11 + vz1*x12
             vuz = m11*vx0z + m12*vx1z
             vzs = vz0*x21 + vz1*x22
             vsz = m21*vx0z + m22*vx1z
             vxzu = vxz0*m11 + vxz1*m12
             vxuz = x11*v0z + x12*v1z
             vxzs = vxz0*m21 + vxz1*m22
             vxsz = x21*v0z + x22*v1z
             vzz = vzz - vzu*phz2 - vzs*dphz2 - vuz*phz - vsz*dphz &
                  +  sodb(1,1,k,1,2)*phz*phz2 &
                  +  sodb(2,1,k,1,2)*phz*dphz2 &
                  +  sodb(1,2,k,1,2)*dphz*phz2 &
                  +  sodb(2,2,k,1,2)*dphz*dphz2
             vxzz = vxzz - vxzu*phz - vxzs*dphz - vxuz*phz2 - vxsz*dphz2 &
                  +  sodb(1,1,k,2,2)*phz2*phz &
                  +  sodb(2,1,k,2,2)*phz2*dphz &
                  +  sodb(1,2,k,2,2)*dphz2*phz &
                  +  sodb(2,2,k,2,2)*dphz2*dphz
             vuz  = vuz  - sodb(1,1,k,1,2)*phz2 - sodb(2,1,k,1,2)*dphz2
             vxuz = vxuz - sodb(1,1,k,2,2)*phz  - sodb(2,1,k,2,2)*dphz
             vzu  = vzu  - sodb(1,1,k,1,2)*phz  - sodb(1,2,k,1,2)*dphz
             vxzu = vxzu - sodb(1,1,k,2,2)*phz2 - sodb(1,2,k,2,2)*dphz2
             vsz  = vsz  - sodb(1,2,k,1,2)*phz2  -sodb(2,2,k,1,2)*dphz2
             vxsz = vxsz - sodb(1,2,k,2,2)*phz   -sodb(2,2,k,2,2)*dphz
             vzs  = vzs  - sodb(2,1,k,1,2)*phz  - sodb(2,2,k,1,2)*dphz
             vxzs = vxzs - sodb(2,1,k,2,2)*phz2 - sodb(2,2,k,2,2)*dphz2
             sodb(1,3,k,1,2) = vuz
             sodb(2,3,k,1,2) = vsz
             sodb(3,3,k,1,2) = vzz
             sodb(3,1,k,1,2) = vzu
             sodb(3,2,k,1,2) = vzs
             sodb(1,3,k,2,2) = vxuz
             sodb(2,3,k,2,2) = vxsz
             sodb(3,3,k,2,2) = vxzz
             sodb(3,1,k,2,2) = vxzu
             sodb(3,2,k,2,2) = vxzs
          endif
       enddo
    endif
  end subroutine potpus

  subroutine pvpus1(r,f,df,g,dg,Tfg,Tgf)
    !- Forces K.E. or hamiltonian matrix elements to satisfy Wronskian
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   r     :radius
    !i   f     :value of first function at r
    !i   df    :df/dr
    !i   g     :value of second function at r
    !i   dg    :dg/dr
    ! o Inputs/Outputs
    ! o  Tfg   :<f | T | g> or <f | H | g>
    ! o        :Input value is overwritten with one that satisfies Wronskian:
    ! o  Tgf   :<g | T | f> or <g | H | f>
    !r Remarks
    !r   The matrix elements of the Laplacian operator int_0^r f -nabla g
    !r   for any two analytic radial functions (f,g), which satisfy
    !r   have val=0 or slope=0 at r=0 must also satisfy the Wronskian
    !r     Tfg-Tgf = -W(f,g)
    !u Updates
    !u   21 Jul 04  First created
    ! ----------------------------------------------------------------------
    implicit none
    double precision :: Tfg,Tgf,r,df,g,f,dg
    double precision :: diff,avg
    diff = r*r*(df*g - f*dg)!     diff = Tfg-Tgf
    avg = Tfg+Tgf
    Tfg = (avg+diff)/2
    Tgf = (avg-diff)/2
  end subroutine pvpus1
  subroutine soprm(mode,lpzi,phi,phid,phiz,nr,nsp,lmxs,lmx,v,dv,enu, ez,z,ri,wi,wk,sop,sopz)
    use m_lgunit,only:stdo
    !- Radial matrix elements between orbitals of different spin
    ! ---------------------------------------------------------------------
    !i Inputs
    !i   mode  :1 make spin-orbit parameters, i.e. matrix elements
    !i         :  <phi|so|phi> <phi|so|phidot> <phidot|so|phidot>
    !i         :  and for local orbitals that are present
    !i         :  <phiz|so|phiz> <phiz|so|phi> <phiz|so|phidot>
    !i         :2 make matrix elements for input magnetic field B
    !i         :  <phi|B|phi> <phi|B|phidot> <phidot|B|phidot>
    !i         :4 orthonormalize phi,phidot in separate spin channels
    !i         :5 1+4 above
    !i         :6 2+4 above
    !i   lpzi  :flags which channels have local orbitals
    !i   phi   :radial wave function * r
    !i   phid  :energy derivative of phi
    !i   phiz  :local orbital wave function
    !i   nr    :number of radial mesh points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   lmxs  :phi,phid,phiz,sop,sopz are dimensioned 0:lmxs
    !i   lmx   :lmx(j) = maximum l for atom j
    !i   v     :electron part of potential: complete potential is v(r)-2*z/r.
    !i   dv    :(mode=1) radial derivative of potential, dV/dr with V=(V+ + V-)/2
    !i         :(mode=2) magnetic field
    !i         :(mode=4) not used
    !i   enu   :enu's for making charge density
    !i   ez    :enu's for local orbitals
    !i   z     :nuclear charge
    !i   ri    :radial mesh
    !i   wi    :weights for integration on mesh
    !i   wk    :work array of length nr*4
    !o  Outputs
    !o   sop   :sop(l,is1,is2,i=1..3) : matrix elements between orbitals
    !o         :of spin is1 and is2 for quantum number l.
    !o         :Three types of integrals are calculated for i=1..3:
    !o         :<phi pert phi>  <phi pert phidot>  <phidot pert phidot>
    !o   sopz  :sopz(l,is1,is2,i=1..3) : matrix elements between local orbitals
    !o         :and orbitals of spin is1 and is2 for quantum number l.
    !o         :Three types of integrals are calculated for i=1..3:
    !o         :<phiz SO phiz>  <phiz SO phi>  <phiz SO phidot>.
    !r  Remarks
    !r   so = 2/(c^2) dV/dr*(1/r), V(r)=-2*z/r+v(r)
    !r   Note: so=1/(2*m^2*c^2)*(dV/dr*1/r), m=.5, c=274 (at. Rydberg units)
    !r   H_so = so*s^ dot l^, s^=0.5d0*sigma (Pauli matrix).
    !u Updates
    !u   17 Jan 07 Set ME for l=0 to zero (weren't necesssarily calc before)
    !u   11 Jul 05 Merged with sofp to make one routine
    !u   25 Apr 05 A. Chantis added local orbitals
    !u   05 Jan 04 leading dimensions of sop distinct from lmx
    !u   07 Feb 03 Added ability to compute matrix elements of external
    !u             field B.  New arg list and definition of mode.
    ! ---------------------------------------------------------------------
    implicit none
    integer :: lmx,mode,lpzi(0:lmx),lmxs,nr,nsp
    double precision :: z, &
         phi(nr,0:lmxs,nsp),phid(nr,0:lmxs,nsp),phiz(nr,0:lmxs,nsp), &
         ri(nr),sop(0:lmxs,nsp,nsp,3),v(nr,nsp),wk(nr,4),dv(nr), &
         wi(nr),sopz(0:lmxs,nsp,nsp,3),enu(0:8,nsp),ez(0:8,nsp)
    integer :: l,ir,is,is1,is2,ipr,mode0,lmin
    double precision :: c,pa,r,r1,r2,dot3,vavg,eavg,eavgz,dva,xx,xxz,xxavg,wkz(nr,4)
    common /cc/ c ! c = 274.071979d0 or 1d10 !  Speed of light, or infinity in nonrelativistic case
    data pa /1d0/
    call getpr(ipr)
    mode0 = mod(mode,4)  
    lmin = 1
    if (mode0 == 2) lmin=0
    ! --- Orthonormalize phi, phidot, neglecting small component ---
    if (mode/4 /= 0) then
       if (ipr > 50) print '(/'' soprm: overlaps  phi*phi     phi*phidot'')'
       do   is = 1, nsp
          do   l = 0, lmx
             r1 = dot3(nr,phi(1,l,is),phi(1,l,is),wi)
             phi(:,l,is) = phi(:,l,is)/dsqrt(r1) ! call dscal(nr,1,,1)
             r2 = dot3(nr,phi(1,l,is),phid(1,l,is),wi)
             phid(:,l,is)=phid(:,l,is) -r2*phi(:,l,is)! daxpy(nr,-r2,phi(1,l,is),1,phid(1,l,is),1)
             if (ipr > 50) write(stdo,"('  spin',i2,'  l=',i1,2f13.6)") is,l,r1,r2
          enddo
       enddo
    endif
    if (mode0==0) return
    ! --- Matrix elements for each l ---
    wk(1,:) = 0d0
    wkz(1,:) = 0d0
    if(mode0==2) wk(2:nr,1) = wi(2:nr) * dv(2:nr)
    ! .. Initialize matrix elements for s orbitals, in case not calculated
    l = 0
    do  is2 = 1, nsp
       do  is1 = 1, nsp
          do is = 1, 3
             sop(l,is1,is2,is) = 0d0
             if (lpzi(l) /= 0) sopz(l,is1,is2,is) = 0d0
          enddo
       enddo
    enddo
    do  l = lmin, lmx
       eavg = (enu(l,1)+enu(l,nsp))/2
       if (lpzi(l) /= 0) eavgz = (ez(l,1)+ez(l,nsp))/2
       do  is2 = 1, nsp
          do  is1 = 1, nsp
             if (mode0 == 1) then
                do  ir = 2, nr
                   r = ri(ir)
                   vavg = (v(ir,1)+v(ir,nsp))/2 - 2*z/r
                   dva  = dv(ir) + 2*z/r**2
                   xx = 1/r/(1d0+pa*(eavg-vavg)/c**2)**2
                   if (lpzi(l) /= 0) then
                      xxz = 1/r/(1d0+pa*(eavgz-vavg)/c**2)**2
                      xxavg = 0.5d0*(xx+xxz)
                      wkz(ir,1) = phiz(ir,l,is1)*dva
                      wkz(ir,2) = phiz(ir,l,is1)*xxz
                      wkz(ir,3) = phi(ir,l,is2)*xxavg
                      wkz(ir,4) = phid(ir,l,is2)*xxavg
                   endif
                   wk(ir,1) = phi(ir,l,is1)*dva
                   wk(ir,3) = phid(ir,l,is1)*dva
                   wk(ir,2) = phi(ir,l,is2)*xx
                   wk(ir,4) = phid(ir,l,is2)*xx
                enddo
                sop(l,is1,is2,1) = dot3(nr,wk,wk(1,2),wi)*2d0/c**2
                sop(l,is1,is2,2) = dot3(nr,wk,wk(1,4),wi)*2d0/c**2
                sop(l,is1,is2,3) = dot3(nr,wk(1,3),wk(1,4),wi)*2d0/c**2
                if (lpzi(l) /= 0) then
                   sopz(l,is1,is2,1) = dot3(nr,wkz,wkz(1,2),wi)*2d0/c**2
                   sopz(l,is1,is2,2) = dot3(nr,wkz,wkz(1,3),wi)*2d0/c**2
                   sopz(l,is1,is2,3) = dot3(nr,wkz,wkz(1,4),wi)*2d0/c**2
                endif
             elseif (mode0 == 2) then
                sop(l,is1,is2,1) = dot3(nr,phi(1,l,is1),phi(1,l,is2),wk)
                sop(l,is1,is2,2) = dot3(nr,phi(1,l,is1),phid(1,l,is2),wk)
                sop(l,is1,is2,3) = dot3(nr,phid(1,l,is1),phid(1,l,is2),wk)
             else
                call rxi('soprm: bad mode:',mode)
             endif
          enddo
       enddo
    enddo
    ! --- Printout ---
    if (ipr <= 50) return
    if (mode0 == 1) write(stdo,332) 'spin-orbit coupling'
    if (mode0 == 2) write(stdo,332) 'external field'
332 format(' soprm:  matrix elements for perturbation from ',a/ &
         13x,'l',4x,'<phi || phi>',2x,'<dot || phi>',2x,'<dot || dot>')
    if (nsp == 1) then
       do  l = lmin, lmx
          write(stdo,333) '          ',  l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
          if(lpzi(l)/=0)write(stdo,333) '          ', l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
       enddo
    else
       do  l = lmin, lmx
          write(stdo,333) 'up   up   ', l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
          write(stdo,333) 'down down ', l,sop(l,2,2,1),sop(l,2,2,2),sop(l,2,2,3)
          write(stdo,333) 'up   down ', l,sop(l,1,2,1),sop(l,1,2,2),sop(l,1,2,3)
          write(stdo,333) 'down up   ', l,sop(l,2,1,1),sop(l,2,1,2),sop(l,2,1,3)
          write(stdo,333)
          if (lpzi(l) /= 0) then
             write(stdo,335) 'up   up   ', l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
             write(stdo,335) 'down down ', l,sopz(l,2,2,1),sopz(l,2,2,2),sopz(l,2,2,3)
             write(stdo,335) 'up   down ', l,sopz(l,1,2,1),sopz(l,1,2,2),sopz(l,1,2,3)
             write(stdo,335) 'down up   ', l,sopz(l,2,1,1),sopz(l,2,1,2),sopz(l,2,1,3)
             write(stdo,335)
          endif
       enddo
    endif
333 format(1x,a,i3,1x, 3f14.8)
335 format(1x,a,i3,'l',3f14.8)
  end subroutine soprm
end module m_potpus
