module m_potpus
  public potpus
!  real(8),public, allocatable:: dlphi, dlphip, phi,phip, phz,dphz !3,4,5,6,11,12
!  real(8),public, allocatable:: s00,s11,szz,s0z,s1z !2,7,8,9,10
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
    real(8):: z,rmax,a,rofi(1),v(nr,nsp),ehl(n0),rsml(n0), &
         pnu(n0,nsp),pnz(n0,nsp),ppnl(nppn,n0,2), & !,pp(n0,2,5)
         hab(3,3,n0,nsp),sab(3,3,n0,nsp),vab(3,3,n0,nsp),vdif(nr,nsp), &
         sodb(3,3,n0,nsp,2),rs3,vmtz
    integer:: ipr,ir,i,j,k,l,lpz,lpzi(0:n0),nrbig
    real(8):: dmat(3,3),det,vi, &
         hmat(3,3),phmins,phplus,q,r,smat(3,3),tmc, &
         umegam,umegap,vmat(3,3),vl,xx,xxx,yyy,zzz, &
         b,ghg,ghgp,gphgp, g(nr,2),gp(nr,2*4),gz(nr,2),&
         ev,phi,dphi,phip,dphip,p, dlphi,dlphip,&
         ez,phz,dphz,phzp,dphzp,pz,dlphz,dlphzp,phz2,dphz2
    integer :: nrx
    parameter (nrx=1501)
    double precision :: rwgtx(nrx) !,gzbig(nrx,2)
    real(8),allocatable:: gzbig(:,:)    !     double precision gbig(nrx*2),gpbig(nrx*8),vbig(nrx,2)
    real(8):: vavg(nr),dv(nr),sop(0:lmxa,nsp,nsp,3), &     !     Spin-Orbit related parameters
         sopz(0:lmxa,nsp,nsp,3), &
         psi(nr,0:lmxa,nsp),dpsi(nr,0:lmxa,nsp), &
         pzi(nr,0:lmxa,nsp),ezum(0:8,nsp), &
         enumx(0:8,nsp),wrk(nr,4), x21,x11,x22,x12,vx00,vx01,vx10,vx11, &
         vx0z,vxz0,vx1z,vxz1,vxzz,vxzu,vxuz,vxzs,vxsz, xxxx(1,1)
    real(8),target:: m(2,2,0:lmxa,nsp)
    real(8):: rwgt(nr),szz_,vzz_,hzz_,vx13,vx23
    call getpr(ipr)
    call vxtrap(1,z,rmax,lmxa,v,a,nr,nsp,pnz,rs3,vmtz,nrx,lpz,nrbig,rofi,rwgtx,xxxx(1,1))
    rwgt= rwgtx(1:nr)
    b   = rmax/(dexp(a*nr - a) - 1d0)
    if (lso /= 0) then !Gradient of average v (for spin-orbit)
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
          if (lpzi(l)/=0) then ! ... lo wf gz and its sphere boundary parameters
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
          endif
          ! ... Valence wf g,gp, and their sphere boundary parameters
          call makrwf(10,z,rmax,l,v(1,i),a,nr,rofi,pnu(1,i),2,g,gp,ev,phi,dphi,phip,dphip,p)
          ghg    = ev   ! <g H g> = e <g g> = e
          ghgp   = 1d0  ! <g H gp> = <g (H-e) gp> + e <g gp> = <g g> = 1
          gphgp  = ev*p ! <gp H gp> = <gp (H-e) gp> + e <gp gp> = <gp g> + e p = ep
          dlphi  = rmax*dphi/phi
          dlphip = rmax*dphip/phip
          if(lso /= 0) then ! ... Keep local copies of phi and phidot for SO coupling
             psi(:,l,i) = g(:,1)
             dpsi(:,l,i)= gp(:,1)
             enumx(l,i) = ev
          endif
          RadialIntegrals: block ! --- Integrals of w.f. products with spherical potential ---
          ! ... This branch computes integrals with products of (g,gp,gz)
          !     Convention: 11 (phi,phi) 21 (dot,phi) 22 (dot,dot), 13 (phi,lo) 23 (dot,lo) 33 (lo,lo)
            integer:: fllp1
            real(8):: vii(nr),tmcc(nr),gf11(nr),gf22(nr),gf12(nr),xxxw(nr),yyyw(nr),zzzw(nr)
            fllp1 = l*(l+1)
            vii(1)=0d0
            vii(2:nr) = v(2:nr,i) - 2d0*z/rofi(2:nr)
            tmcc = cc - (vii-ev)/cc
            gf11(1)=0d0
            gf11(2:nr) = 1d0 + fllp1/(tmcc(2:nr)*rofi(2:nr))**2
            dmat(1,1)= sum(rwgt*vdif(:,i)* (gf11* g(:,1)* g(:,1)+ g(:,2)*g(:,2)) ) 
            dmat(2,1)= sum(rwgt*vdif(:,i)* (gf11*gp(:,1)* g(:,1)+gp(:,2)*g(:,2)) ) 
            dmat(2,2)= sum(rwgt*vdif(:,i)* (gf11*gp(:,1)*gp(:,1)+gp(:,2)*gp(:,2)))
            vmat(1,1)= sum(rwgt*vii*       (gf11* g(:,1)* g(:,1)+ g(:,2)*g(:,2)) ) 
            vmat(2,1)= sum(rwgt*vii*       (gf11*gp(:,1)* g(:,1)+gp(:,2)*g(:,2)) )
            vmat(2,2)= sum(rwgt*vii*       (gf11*gp(:,1)*gp(:,1)+gp(:,2)*gp(:,2)))
            dmat(1,2)= dmat(2,1)
            vmat(1,2)= vmat(2,1)
            vmat(1:2,1:2) = vmat(1:2,1:2) + dmat(1:2,1:2)
            smat(1:2,1) = [1d0,0d0]
            smat(1:2,2) = [0d0, p]
            hmat(1:2,1) = [ghg,0d0]   ! hmat(1,1)=<g H g> = e <g g> = e                       
            hmat(1:2,2) = [ghgp,gphgp]! hmat(1,2)=<g H gp> = <g (H-e) gp> + e <g gp> = <g g>  
            hmat(1:2,1:2)=hmat(1:2,1:2)+dmat(1:2,1:2)
            ! hmat(2,1)=<gp H g>=0d0, hmat(2,2)=<gp H gp>=<gp (H-e) gp>+e <gp gp>=<gp g> + e p =ep
            call pvpus1(rmax,phi,dphi,phip,dphip,hmat(1,2),hmat(2,1))
            ! pvpus1 not needed? since Wronskian explicit in makrwf
            if(lpzi(l)/=0) then !computes integrals with products of (g,gp) x gz
               gf12(1)=0d0
               gf22(1)=0d0
               tmcc= cc - (vii-ez)/cc
               gf22(2:nr) = 1d0 + fllp1/(tmcc(2:nr)*rofi(2:nr))**2
               gf12 = (gf11 + gf22)/2
               dmat(1,3) = sum(rwgt*vdif(:,i)*(gf12*g(:,1) *gz(:,1)+ g(:,2) *gz(:,2)))
               dmat(2,3) = sum(rwgt*vdif(:,i)*(gf12*gp(:,1)*gz(:,1)+ gp(:,2)*gz(:,2)))
               dmat(3,3) = sum(rwgt*vdif(:,i)*(gf22*gz(:,1)*gz(:,1)+ gz(:,2)*gz(:,2)))
               vmat(1,3) = sum(rwgt*vii*(gf12*g(:,1) *gz(:,1)+ g(:,2) *gz(:,2)))
               vmat(2,3) = sum(rwgt*vii*(gf12*gp(:,1)*gz(:,1)+ gp(:,2)*gz(:,2)))
               vmat(3,3) = sum(rwgt*vii*(gf22*gz(:,1)*gz(:,1)+ gz(:,2)*gz(:,2)))
               smat(1,3) = sum(rwgt*(gf12*g(:,1) *gz(:,1)+ g(:,2) *gz(:,2)))
               smat(2,3) = sum(rwgt*(gf12*gp(:,1)*gz(:,1)+ gp(:,2)*gz(:,2)))
               smat(3,3) = sum(rwgt*(gf22*gz(:,1)*gz(:,1)+ gz(:,2)*gz(:,2)))
               dmat(3,1:2) = dmat(1:2,3)
               vmat(3,1:2) = vmat(1:2,3)
               smat(3,1:2) = smat(1:2,3)
               vmat(1:2,3) = vmat(1:2,3) + dmat(1:2,3)
               vmat(3,1:2) = vmat(3,1:2) + dmat(3,1:2)
               vmat(3,3)   = vmat(3,3) + dmat(3,3)
               hmat(1,3) = ez*smat(1,3) + dmat(1,3)
               hmat(2,3) = ez*smat(2,3) + dmat(2,3)
               hmat(3,3) = ez*smat(3,3) + dmat(3,3)
               hmat(3,1) = ev*smat(1,3) + dmat(1,3)
               hmat(3,2) = ev*smat(3,2) + smat(3,1) + dmat(2,3)
               call pvpus1(rmax,phi,dphi,phz,dphz,  hmat(1,3),hmat(3,1))!Put in Wronskians explicitly
               call pvpus1(rmax,phip,dphip,phz,dphz,hmat(2,3),hmat(3,2))
            endif
          endblock RadialIntegrals
          !  Linear transformation between (u,s) and (phi,phidot)
          !     u_i(r) = \sum_j M_ij phi_j(r), where u_i=(u,s), phi_i=(phi,phidot)
          !     u: u(rmt)=1,u'(rmt)=0, and s:s(rmt)=0 s'(rmt)=1
          !       (1  0)      (phi     dphi )                    ( dphip  -dphi )
          !       (    )  = M (             )  = >  M = (det)^-1 (              )
          !       (0  1)      (phip    dphip)                    (-phip    phi  )
          !     where det = phi*dphip - phip*dphi
          !     To compute matrix elements of (u,s) from elements (phi,phidot),
          !     let u_i = u or s and phi_i = phi or phidot; i=1 or 2
          !     <u_i|u_j> = <(M phi)_i|(M phi)_j> = sum_lm M_il <phi_l|phi_m> M_mj
          det = phi*dphip - dphi*phip
          m(1,1,l,i) = dphip/det  !M matrix for u_i = \sum_j M_ij phi_j 
          m(1,2,l,i) = -dphi/det
          m(2,1,l,i) = -phip/det
          m(2,2,l,i) =  phi/det  
          matmmm :block ! Fix confusing indexing (transposed hab...). 2022-12-26
            real(8),pointer:: mmm(:,:)
            real(8):: mmmt(2,2) ! ... (u,s) x (u,s)
            mmm => m(:,:,l,i)
            mmmt= transpose(mmm)
            hab(1:2,1:2,k,i) = matmul(mmm,matmul(hmat(1:2,1:2),mmmt)) !<u,s|h|u,s>
            sab(1:2,1:2,k,i) = matmul(mmm,matmul(smat(1:2,1:2),mmmt)) 
            vab(1:2,1:2,k,i) = matmul(mmm,matmul(vmat(1:2,1:2),mmmt))
          endblock matmmm
          if (lpzi(l) == 2) then
             phz = 0
             dphz = 0
          endif
          if (lpzi(l) /= 0) then  ! <(u,s,gz) |h| (u,s,gz)>
             lpzint: block   !  gz = (gz0 - gz0(rmax) u - r*(gz0/r)'(rmax) s) = (gz0 - phz u - dphz s)
               real(8):: mm0(3,3),mmz(3,3),mmm(3,3),mmmt(3,3)
               MM0(:,1)=[m(:,1,l,i), 0d0] !1st col (u s gz0)^t= MM0 t(phi phidot gz0)^t
               MM0(:,2)=[m(:,2,l,i), 0d0] 
               MM0(:,3)=[0d0,   0d0, 1d0]
               MMz(:,1)=[1d0,  0d0, -phz] !1st col (u s gz)^t= MMz t(u s gz0)^t
               MMz(:,2)=[0d0,  1d0,-dphz]
               MMz(:,3)=[0d0,  0d0,  1d0]
               mmm =matmul(MMz,MM0)
               mmmt=transpose(mmm)
               sab(1:3,1:3,k,i) = matmul(mmm,matmul(smat(1:3,1:3),mmmt)) !<(u,s,gz)|(u,s,gz)>
               vab(1:3,1:3,k,i) = matmul(mmm,matmul(vmat(1:3,1:3),mmmt)) !<(u,s,gz)|v|(u,s,gz)>
               hab(1:3,1:3,k,i) = matmul(mmm,matmul(hmat(1:3,1:3),mmmt)) !<(u,s,gz)|h|(u,s,gz)>
             endblock lpzint
          endif
          ppnl(1:7,k,i) = [0d0, 1d0, rmax*dlphi, rmax*dlphip, phi, phip, p]
          if (lpzi(l) /= 0) then  
             ppnl(8,k,i)  = sab(3,3,k,i)
             ppnl(9,k,i)  = smat(1,3) - sum([phi,  dphi]*matmul(sab(1:2,1:2,k,i),[phz,dphz]))
             ppnl(10,k,i) = smat(2,3) - sum([phip,dphip]*matmul(sab(1:2,1:2,k,i),[phz,dphz]))
             ppnl(11:12,k,i) = [phz, dphz]
          else
             ppnl(8:12,k,i) = 0d0 
          endif
10     enddo lloop
80  enddo isploop
    if(lso==0) return
    ! Get spin-orbit parameters
    call soprm(5,lpzi,psi,dpsi,pzi,nr,nsp,lmxa,lmxa,v,dv,enumx, ezum,z,rofi,rwgt,wrk,sop,sopz)
    do i = 1, nsp ! Make the spin diagonal radial integrals
       do l = 0, lmxa
          phz = ppnl(11,l+1,i)
          dphz = ppnl(12,l+1,i)
          k = l + 1
          if (lpzi(l)==0) then 
             sodbmat: block  !note original sodb is transposed to avoid confusing defition. 2022-12-26
               real(8),pointer:: mmm(:,:)
               real(8):: mmmt(2,2),somat(2,2) ! ... (u,s) x (u,s)
               mmm => m(:,:,l,i) !M matrix for u_i = \sum_j M_ij phi_j 
               mmmt= transpose(mmm)
               somat(:,1) = [sop(l,i,i,1),sop(l,i,i,2)]            
               somat(:,2) = [sop(l,i,i,2),sop(l,i,i,3)]
               sodb(1:2,1:2,k,i,1) = matmul(mmm,matmul(somat,mmmt))
             endblock sodbmat
          else !  ...  Make the local orbitals spin diagonal integrals
             sodbint: block !gz= (gz0 - gz0(rmax) u - r*(gz0/r)'(rmax) s) = (gz0 - phz u - dphz s)
               real(8):: mm0(3,3),mmz(3,3),mmm(3,3),mmmt(3,3),somat(3,3)
               mm0(:,1)=[m(:,1,l,i), 0d0] !1st col (u s gz0)^t= MM0 t(phi phidot gz0)^t
               mm0(:,2)=[m(:,2,l,i), 0d0] 
               mm0(:,3)=[0d0,  0d0,  1d0]
               mmz(:,1)=[1d0,  0d0, -phz] !1st col (u s gz)^t= MMz t(u s gz0)^t
               mmz(:,2)=[0d0,  1d0,-dphz]
               mmz(:,3)=[0d0,  0d0,  1d0]
               mmm =matmul(mmz,mm0)
               mmmt=transpose(mmm)
               somat(:,1) = [sop(l,i,i,1), sop(l,i,i,2), sopz(l,i,i,2)]            
               somat(:,2) = [sop(l,i,i,2), sop(l,i,i,3), sopz(l,i,i,3)]
               somat(:,3) = [sopz(l,i,i,2),sopz(l,i,i,3),sopz(l,i,i,1)]
               sodb(1:3,1:3,k,i,1) = matmul(mmm,matmul(somat,mmmt)) !<(u,s,gz)|SO_diag|(u,s,gz)>
             endblock sodbint
          endif
       enddo
    enddo
    do l = 0, lmxa !   ... Make the spin off-diagonal radial integrals
       k = l + 1
       phz   = ppnl(11,k,1)
       dphz  = ppnl(12,k,1)
       phz2  = ppnl(11,k,2)
       dphz2 = ppnl(12,k,2)
       if (lpzi(l)==0) then 
          sooff2: block
            real(8),pointer::mmm1(:,:),mmm2(:,:)
            real(8):: soud(2,2),sodu(2,2)
            mmm1=>m(:,:,l,1)
            mmm2=>m(:,:,l,2)
            soud(:,1)=[sop(l,1,2,1),sop(l,2,1,2)] 
            soud(:,2)=[sop(l,1,2,2),sop(l,1,2,3)]
            sodu(:,1)=[sop(l,2,1,1),sop(l,2,1,2)]
            sodu(:,2)=[sop(l,1,2,2),sop(l,2,1,3)]
            sodb(1:2,1:2,k,1,2)= matmul(mmm1,matmul(soud,transpose(mmm2)))! up-dn block 
            sodb(1:2,1:2,k,2,2)= matmul(mmm2,matmul(sodu,transpose(mmm1)))! dn-up block
          endblock sooff2
       else
          sooff3: block ! gz= (gz0 - gz0(rmax) u - r*(gz0/r)'(rmax) s) = (gz0 - phz u - dphz s)
            real(8):: mm1(3,3),mm2(3,3),mm1z(3,3),mm2z(3,3),mmm1(3,3),mmm2(3,3)
            real(8):: soud(3,3),sodu(3,3),ppp(3,3)
            MM1(:,1)=[m(:,1,l,1), 0d0] !L 1st col (u s gz0)^t= MM0 t(phi phidot gz0)^t
            MM1(:,2)=[m(:,2,l,1), 0d0] 
            MM1(:,3)=[0d0,   0d0, 1d0]
            MM2(:,1)=[m(:,1,l,2), 0d0] !R 1st col (u s gz0)^t= MM0 t(phi phidot gz0)^t
            MM2(:,2)=[m(:,2,l,2), 0d0] 
            MM2(:,3)=[0d0,   0d0, 1d0]
            MM1z(:,1)=[1d0,  0d0, -phz] !L 1st col (u s gz)^t= MMz t(u s gz0)^t
            MM1z(:,2)=[0d0,  1d0,-dphz]
            MM1z(:,3)=[0d0,  0d0,  1d0]
            MM2z(:,1)=[1d0,  0d0, -phz2] !R 1st col (u s gz)^t= MMz t(u s gz0)^t
            MM2z(:,2)=[0d0,  1d0,-dphz2]
            MM2z(:,3)=[0d0,  0d0,  1d0]
            mmm1 =matmul(MM1z,MM1)
            mmm2 =matmul(MM2z,MM2)
            soud(:,1)=[sop(l,1,2,1), sop(l,2,1,2),  sopz(l,1,2,2)] !note matrix indexing
            soud(:,2)=[sop(l,1,2,2), sop(l,1,2,3),  sopz(l,1,2,3)]
            soud(:,3)=[sopz(l,2,1,2),sopz(l,2,1,3), sopz(l,1,2,1)]
            sodb(:,:,k,1,2)= matmul(mmm1,matmul(soud,transpose(mmm2)))! up-dn block <1|SOoff|2>
            sodu(:,1)=[sop(l,2,1,1), sop(l,2,1,2), sopz(l,2,1,2)]
            sodu(:,2)=[sop(l,1,2,2), sop(l,2,1,3), sopz(l,2,1,3)]
            sodu(:,3)=[sopz(l,1,2,2),sopz(l,1,2,3),sopz(l,2,1,1)]
            sodb(:,:,k,2,2)= matmul(mmm2,matmul(sodu,transpose(mmm1)))! dn-up block 
          endblock sooff3
       endif
    enddo
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
