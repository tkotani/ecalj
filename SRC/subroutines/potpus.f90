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
    !ox   pp    :potential parameters pp(l+1,isp,j):
    !ox         :(1): enu   linearization energies corresponding to input pnu
    !ox         :(2): c     band center of gravity
    !ox         :(3): srdel band width parameter
    !ox         :(4): qpar  band distortion parameter
    !ox         :(5): ppar  small parameter = integral phidot**2
    !o   ppnl  :NMTO pot pars (but no backward extrapolation; <phi phi>=1)
    !ox         :(1)  = inverse potential function (not created)
    !ox         :(2)  = normalization of phi (= 1 here for now)
    !   pvalue(:), pslope(:), ppovl(:,:)
    !o         :(3)  = rmax * log derivative of phi at rmax      
    !o         :(4)  = rmax * log derivative of phidot at rmax
    !o         :(5)  = phi at rmax
    !o         :(6)  = phidot at rmax
    !o         :(7)  = normalization of phidot
    !o         :(8)  = normalization <gz|gz> for local orbital
    !o         :(9)  = overlap <phi|gz>
    !o         :(10) = overlap <phidot|gz>
    !o         :(11) = phz at rmax for local orbitals defined, before
    !o                 (phi,phidot) subtracted
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
    !l   lorthd:if 0, phi,phidot constructed by calling phidx (normal mode)
    !l         :if 1, orthonormalize phidot to phi
    !l   lorthz:if 0, phi,phidot constructed by calling phidx (normal mode)
    !l         :if 1, orthonormalize phi and phidot to phiz in large sphere
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
    !u Updates
    !u   17 Jan 07 Fixed dimension of ezum (bug fix)
    !u   26 Sep 05 (A. Chantis) Bug fix: local orbitals in conjunction w/ SO
    !u   09 Jul 05 (MvS) bug fix for local orbitals and SO coupling
    !u   08 Jun 05 (MvS) extended ppnl to return phz and dphz
    !u   03 Feb 05 (A. Chantis) Revised parameters for full L.S
    !u   24 Dec 04 (A. Chantis) radial matrix elements for full L.S
    !u   12 Aug 04 First implementation of extended local orbitals
    !u   29 Jun 04 (A. Chantis) spin-diagonal radial spin-orbit integrals
    !u   06 Mar 02 potpus can work with scaled local orb gz (see SCALEGZ)
    !u   13 Feb 02 New matrix elements ppnl(8..10) for local orbitals
    !u   24 Aug 01 Extended to local orbitals.  Altered argument list.
    !u   20 Feb 01 (mvs) returns potential parameters pp,ppn; new arg list
    !u   16.05.00  (mvs) adapted from mc,nfp potpsr.f
    !u   02.01.95  (mvs) spin polarized
    !u   02.09.94  (msm) Potential vdif is included by perturbation.
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: lmxa,lso,nr,nsp,n0,nab,nppn
    parameter (nppn=12)
    double precision :: z,rmax,a,rofi(1),v(nr,nsp),ehl(n0),rsml(n0), &
         pnu(n0,nsp),pnz(n0,nsp),pp(n0,2,5),ppnl(nppn,n0,2), &
         hab(nab,n0,nsp),sab(nab,n0,nsp),vab(nab,n0,nsp),vdif(nr,nsp), &
         sodb(nab,n0,nsp,2),rs3,vmtz
    ! ... Local parameters
    integer :: ipr,ir,i,j,k,l,lpz,lpzi(0:n0),nrbig
    !       lorthz,lrel
    !      character*80 outs
    double precision :: m21,m11,m22,m12,d00,d10,d11,det,vi, &
         fllp1,gfac,gf11,gf22,gf12,h00,h01,h10,h11, &
         phmins,phplus,q,r,s00,s01,s10,s11,tmc, &
         umegam,umegap,v00,v01,v10,v11,vl,xx,xxx,yyy,zzz, &
         enu,c,srdel,qpar,ppar,b,ghg,ghgp,gphgp
    double precision :: d0z,d1z,dzz,v0z,v1z,s0z,s1z,vz0,vz1,sz0, &
         sz1,h0z,hz0,h1z,hz1,suz,ssz,szu,szs,szz, &
         vuz,vsz,vzu,vzs,vzz,huz,hsz,hzu,hzs,hzz
    double precision :: g(nr,2),gp(nr,2*4),gz(nr,2)
    double precision :: ev,phi,dphi,phip,dphip,p, dlphi,dlphip
    double precision :: ez,phz,dphz,phzp,dphzp,pz,dlphz,dlphzp,phz2,dphz2
    !     double precision ul,sl,D
    integer :: nrx
    parameter (nrx=1501)
    double precision :: rwgtx(nrx),gzbig(nrx*2)
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
    isploop: do  80  i = 1, nsp
       if (ipr >= 40) then
          write(stdo,ftox)' potpus spin=',i,'pnu=',ftof(pnu(1:lmxa+1,i),3)
          if (lpz /= 0) write(stdo,ftox), 'pnz=',ftof(pnz(1:lmxa+1,i),3)
       endif
       if(ipr>=50) then
          write(stdo,"(' l',8x,'enu',9x,'v',11x,'c',10x,'srdel',8x,'qpar',8x,'ppar')")
       endif
       ! --- Loop over l ---
       lloop: do  10  l = 0, lmxa
          k = l+1
          lpzi(l) = 0
          if (pnz(k,i) >  0)  lpzi(l) = 1
          if (pnz(k,i) >= 10) lpzi(l) = 2
          !        if (pnz(k,i) >= 20) lpzi(l) = 3
          moda(l) = 5
          if (lpzi(l) /= 0) moda(l) = 6
          ! ... Semicore wf gz and its sphere boundary parameters
          if (lpzi(l) /= 0) then
             call makrwf(10,z,rmax,l,v(1,i),a,nr,rofi,pnz(1,i),4, &
                  gz,gp,ez,phz,dphz,phzp,dphzp,pz)
             !   ... Keep local copies of gz for SO coupling
             if (lso /= 0) then
                call dcopy(nr,gz,1,pzi(1,l,i),1)
                ezum(l,i) = ez
             endif
             !       Extend local orbital to large mesh; match gz to envelope
             if (lpzi(l) > 1) then
                call dcopy(nr,gz,1,gzbig,1)
                call dcopy(nr,gz(1,2),1,gzbig(1+nrbig),1)
                call rwftai(5,rmax,a,nr,nrbig,rofi,phz,dphz,xx,l,ehl(k), rsml(k),gzbig)
                !         If rwftai scales gzbig, rescale phz,gz
                if (gzbig(nr) /= gz(nr,1)) then
                   xx = gzbig(nr)/gz(nr,1)
                   phz   = phz*xx
                   dphz  = dphz*xx
                   phzp  = phzp*xx
                   dphzp = dphzp*xx
                   call dscal(nr,xx,gz(1,1),1)
                   call dscal(nr,xx,gz(1,2),1)
                endif
                !     ... Recopy gz for SO coupling
                if (lso /= 0) call dcopy(nr,gz,1,pzi(1,l,i),1)
             endif
          else
             phz = 0
             dphz = 0
          endif
          ! ... Valence w.f. g,gp, and their sphere boundary parameters
          call makrwf(10,z,rmax,l,v(1,i),a,nr,rofi,pnu(1,i),2,g,gp,ev,phi,dphi,phip,dphip,p)
          !     <g H g> = e <g g> = e
          !     <g H gp> = <g (H-e) gp> + e <g gp> = <g g> = 1
          !     <gp H gp> = <gp (H-e) gp> + e <gp gp> = <gp g> + e p = ep
          ghg = ev
          ghgp = 1d0
          gphgp = ev*p
          dlphi  = rmax*dphi/phi
          dlphip = rmax*dphip/phip
          ! ... Integrate g and gp outward on extended mesh
          !     large component gp done nonrelativistically ... ok for large r
          !     small component gp computed large component gp
          !      For this routine to work, uncomment call to vxtrap
          !      if (lpzi(l) .gt. 1 .and. lorthz .ne. 0) then
          !        call rx('lorthz not implemented')

          !        call dmcpy(g,nr,1,gbig,nrbig,1,nr,2)
          !        call dpzero(gpbig,2*nrx)
          !        call dmcpy(gp,nr,1,gpbig,nrbig,1,nr,2)
          !        call rsq1(nr+1,ev,l,z,vbig(1,i),nrbig,gbig,xxx,yyy,j,a,b,rofi,
          !     .    nrbig)
          !        call rsqnri(ev,l,vbig(1,i),z,a,rofi,nr+1,nrbig,gbig,j,xxx,yyy,
          !     .    gpbig,q)
          !        call rsqnrs(ev,l,vbig(1,i),z,rofi,nr+1,nrbig,gpbig)

          !        call prrmsh('g-big before orth.',rofi,gbig,nrbig,nrbig,2)
          !        call prrmsh('gp-big before orth',rofi,gpbig,nrbig,nrbig,2)
          !        call ortrwf(4,z,l,vbig(1,i),nrbig,nr,nrbig,rofi,rwgt,ev,ev,ez,
          !     .    gbig,gpbig,gzbig,D)
          !        call prrmsh('g-big after orth.',rofi,gbig,nrbig,nrbig,2)
          !        call prrmsh('gp-big after orth',rofi,gpbig,nrbig,nrbig,2)

          !        call rx('done')
          !      endif

          ! ... Keep local copies of phi and phidot for SO coupling
          if(lso /= 0) then
             do ir = 1, nr
                psi(ir,l,i) = g(ir,1)
                dpsi(ir,l,i) = gp(ir,1)
             enddo
             enumx(l,i) = ev
          endif
          !$$$C   ... Scale gz so that <|gz-P(g,gp)|^2> = 1
          !$$$#if SCALEGZ
          !$$$          if (lpzi(l) .ne. 0) then
          !$$$
          !$$$CC       debugging ... check that mode 0 works properly when g and gp
          !$$$CC       are not orthogonal.  The next line
          !$$$C        call ortrwf(10,z,l,v(1,i),nr,nr,nr,rofi,rwgt,ev,ev,ez,g,gp,gz,D)
          !$$$CC       Should produce the same D as if g,gp were orthogonalized
          !$$$CC       and then computed, as done in the next two lines
          !$$$C        call ortrwf(12,z,l,v(1,i),nr,nr,nr,rofi,rwgt,ev,ev,ez,g,gp,gz,D)
          !$$$C        call ortrwf(11,z,l,v(1,i),nr,nr,nr,rofi,rwgt,ev,ev,ez,g,gp,gz,D)
          !$$$
          !$$$            call ortrwf(0,z,l,v(1,i),nr,nr,nr,rofi,rwgt,ev,ev,ez,g,gp,gz,D)
          !$$$            call dscal(nr*2,1/D,gz,1)
          !$$$            phz = phz/D
          !$$$            dphz = dphz/D
          !$$$            phzp = phzp/D
          !$$$            dphzp = dphzp/D
          !$$$C       call prrmsh('gz',rofi,gz,nr,nr,2)
          !$$$
          !$$$          endif
          !$$$#endif

          ! ... 2nd generation potential parameters (Methfessel scaling)
          umegam = -(phi*ghgp/phip)*(-l-1-dlphi)/(-l-1-dlphip)
          umegap = -(phi*ghgp/phip)*(l-dlphi)/(l-dlphip)
          phplus = phi + umegap*phip/ghgp
          phmins = phi + umegam*phip/ghgp
          enu = ev
          c   = ev + umegam
          vl  = ev + umegap
          srdel = phmins*dsqrt(0.5d0*rmax)
          q   = phmins/(2*(2*l+1)*phplus)
          qpar = 1d0/q
          ppar = 1d0/dsqrt(p)

          if (ipr >= 50) write(stdo,699) l,'    ',ev,vl,c,srdel,qpar,ppar
699       format(i2,a,f10.6,3f12.6,f12.4,f12.6)
          pp(k,i,1) = enu
          pp(k,i,2) = c
          pp(k,i,3) = srdel
          pp(k,i,4) = qpar
          pp(k,i,5) = ppar

          ! ... NMTO potential parameters, no backwards integration, <phi phi>=1
          call dpzero(ppnl(1,k,i),nppn)
          ppnl(1,k,i) = 0
          ppnl(2,k,i) = 1d0
          ppnl(3,k,i) = rmax * dlphi
          ppnl(4,k,i) = rmax * dlphip
          ppnl(5,k,i) = phi
          ppnl(6,k,i) = phip
          ppnl(7,k,i) = p
          ppnl(11,k,i) = phz
          ppnl(12,k,i) = dphz
          ! ... 2nd generation potential sc parameters (Methfessel scaling)
          !     NB: overwrites c,vl,srdel,q,qpar
          !     ppars not saved; only printed out
          if (lpzi(l) /= 0 .AND. ipr >= 50) then
             dlphz  = rmax*dphz/phz
             dlphzp = rmax*dphzp/phzp
             umegam = -(phz/phzp)*(-l-1-dlphz)/(-l-1-dlphzp)
             umegap = -(phz/phzp)*(l-dlphz)/(l-dlphzp)
             phplus = phz + umegap*phzp
             phmins = phz + umegam*phzp
             enu = ez
             c   = ez + umegam
             vl  = ez + umegap
             srdel = phmins*dsqrt(0.5d0*rmax)
             q   = phmins/(2*(2*l+1)*phplus)
             qpar = 1d0/q
             ppar = 1d0/dsqrt(pz)
             write(stdo,699) l,'(sc)',ez,vl,c,srdel,qpar,ppar
          endif

          !$$$C --- Rescale gp parameters and normalization of gp ---
          !$$$#if ORTHPHD
          !$$$          if (lorthd .ne. 0) then
          !$$$            call ortrwf(2,z,l,v,nrbig,nr,nrbig,rofi,rwgt,ev,ev,ez,g,gp,xx,xx)
          !$$$C      <g H gp> = <g (H-e) gp> = 1/sqrt(p)
          !$$$C      <gp H gp> = <gp (H-e) gp> + e <gp gp> = <gp g> + e = e
          !$$$            ghgp = 1/sqrt(p)
          !$$$            gphgp = ev
          !$$$C      Scaling of gp -> phip, <gp gp> also scaled
          !$$$            phip = phip/sqrt(p)
          !$$$            dphip = dphip/sqrt(p)
          !$$$            p = 1
          !$$$          endif
          !$$$#endif

          ! --- Integrals of w.f. products with spherical potential ---
          fllp1 = l*(l+1)
          v00 = 0d0
          v10 = 0d0
          v11 = 0d0
          d00 = 0d0
          d10 = 0d0
          d11 = 0d0

          s0z = 0d0
          s1z = 0d0
          szz = 0d0
          v0z = 0d0
          v1z = 0d0
          vzz = 0d0
          d0z = 0d0
          d1z = 0d0
          dzz = 0d0

          ! ... This branch computes integrals with products of (g,gp)
          !     Convention: 00 (phi,phi) 10 (dot,phi) 11 (dot,dot)
          if (lpzi(l) == 0) then
             block
               real(8):: vii(nr),tmcc(nr),gfac(nr)
               vii(1)=0d0
               vii(2:nr) = v(2:nr,i) - 2d0*z/rofi(2:nr)
               tmcc = cc - (vii-ev)/cc
               gfac(1)=0d0
               gfac(2:nr) = 1d0 + fllp1/(tmcc(2:nr)*rofi(2:nr))**2
               d00= sum(rwgt*vdif(:,i)* (gfac* g(:,1)*g(:,1)+ g(:,2)*g(:,2)) ) 
               d10= sum(rwgt*vdif(:,i)* (gfac* gp(:,1)*g(:,1)+ gp(:,2)*g(:,2)) ) 
               d11= sum(rwgt*vdif(:,i)* (gfac*gp(:,1)*gp(:,1)+gp(:,2)*gp(:,2)) )
               v00= sum(rwgt*vii*       (gfac* g(:,1)* g(:,1)+ g(:,2)* g(:,2)) ) 
               v10= sum(rwgt*vii*       (gfac*gp(:,1)* g(:,1)+gp(:,2)* g(:,2)) )
               v11= sum(rwgt*vii*       (gfac*gp(:,1)*gp(:,1)+gp(:,2)*gp(:,2)) )
             endblock
             ! ... This branch computes integrals with products of (g,gp,gz)
             !     Convention: 0z (phi,sc) 1z (dot,sc) zz (sc,sc)
          else
             !       This section overwrite gz with (gz - gz(rmax) u - gz'(rmax) s)
             !       Instead, we leave gz untouched and correct integrals later.
             !       It helps to avoid ambiguities about the meaning of inner
             !       products of scalar-relativistic wave functions.
             !        det = phi*dphip - dphi*phip
             !        m11 = dphip/det
             !        m12 = -dphi/det
             !        m21 = -phip/det
             !        m22 = phi/det
             !        do  j  = 1, 2
             !        do  ir = 2, nr
             !          ul = m11*g(ir,j) + m12*gp(ir,j)
             !          sl = m21*g(ir,j) + m22*gp(ir,j)
             !          gz(ir,j) = gz(ir,j) - phz*ul - dphz*sl
             !        enddo
             !        enddo
             !        call prrmsh('phi',rofi,g,nr,nr,2)
             !        call prrmsh('phidot',rofi,gp,nr,nr,2)
             !        call prrmsh('gz',rofi,gz,nr,nr,2)
             do  ir = 2, nr
                r = rofi(ir)
                vi = v(ir,i) - 2d0*z/r
                tmc = cc - (vi-ev)/cc
                gf11 = 1d0 + fllp1/(tmc*r)**2
                tmc = cc - (vi-ez)/cc
                gf22 = 1d0 + fllp1/(tmc*r)**2
                xxx = rwgt(ir)*vi
                yyy = rwgt(ir)*vdif(ir,i)
                gf12 = (gf11 + gf22)/2
                zzz = rwgt(ir)

                d0z = d0z + yyy*(gf12*g(ir,1)*gz(ir,1) + g(ir,2)*gz(ir,2))
                d1z = d1z + yyy*(gf12*gp(ir,1)*gz(ir,1) + gp(ir,2)*gz(ir,2))
                dzz = dzz + yyy*(gf22*gz(ir,1)*gz(ir,1) + gz(ir,2)*gz(ir,2))
                v0z = v0z + xxx*(gf12*g(ir,1)*gz(ir,1) + g(ir,2)*gz(ir,2))
                v1z = v1z + xxx*(gf12*gp(ir,1)*gz(ir,1) + gp(ir,2)*gz(ir,2))
                vzz = vzz + xxx*(gf22*gz(ir,1)*gz(ir,1) + gz(ir,2)*gz(ir,2))

                d00 = d00 + yyy*(gf11*g(ir,1)*g(ir,1) + g(ir,2)*g(ir,2))
                d10 = d10 + yyy*(gf11*gp(ir,1)*g(ir,1) + gp(ir,2)*g(ir,2))
                d11 = d11 + yyy*(gf11*gp(ir,1)*gp(ir,1) + gp(ir,2)*gp(ir,2))
                v00 = v00 + xxx*(gf11*g(ir,1)*g(ir,1) + g(ir,2)*g(ir,2))
                v10 = v10 + xxx*(gf11*gp(ir,1)*g(ir,1) + gp(ir,2)*g(ir,2))
                v11 = v11 + xxx*(gf11*gp(ir,1)*gp(ir,1) + gp(ir,2)*gp(ir,2))

                s0z = s0z + zzz*(gf12*g(ir,1)*gz(ir,1)+g(ir,2)*gz(ir,2))
                s1z = s1z + zzz*(gf12*gp(ir,1)*gz(ir,1)+gp(ir,2)*gz(ir,2))
                szz = szz + zzz*(gf22*gz(ir,1)*gz(ir,1)+gz(ir,2)*gz(ir,2))
             enddo
          endif

          v00 = v00 + d00
          v10 = v10 + d10
          v11 = v11 + d11
          v01 = v10
          s00 = 1d0
          s10 = 0d0
          s01 = 0d0
          s11 = p
          !     h00 = <g H g> = e <g g> = e
          !     h01 = <g H gp> = <g (H-e) gp> + e <g gp> = <g g>
          !     h11 = <gp H gp> = <gp (H-e) gp> + e <gp gp>
          !         = <gp g> + e p = ep
          h00 = ghg   + d00
          h01 = ghgp  + d10
          h10 = 0d0   + d10
          h11 = gphgp + d11
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
             h0z = ez*s0z + d0z
             !??     h0z = ez*s0z + d0z - ev*(phz*m11 + dphz*m21) ! if gz has u,s subt.
             hz0 = ev*s0z + d0z
             h1z = ez*s1z + d1z
             !??     h1z = ez*s1z + d1z - ev*(phz*m12 + dphz*m22) ! if gz has u,s subt.
             hz1 = ev*sz1 + sz0 + d1z
             hzz = ez*szz + dzz
             !       Put in Wronskians explicitly
             call pvpus1(rmax,phi,dphi,phz,dphz,h0z,hz0)
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
          hab(1,k,i) = m11*h00*m11 + m11*h01*m12 + m12*h10*m11 + m12*h11*m12
          hab(2,k,i) = m11*h00*m21 + m11*h01*m22 + m12*h10*m21 + m12*h11*m22
          hab(3,k,i) = m21*h00*m11 + m21*h01*m12 + m22*h10*m11 + m22*h11*m12
          hab(4,k,i) = m21*h00*m21 + m21*h01*m22 + m22*h10*m21 + m22*h11*m22
          !     Next 2 lines put in Wronskian explicitly
          !     Should no longer be needed: Wronskian explicit in makrwf.
          hab(2,k,i) = 0.5d0 * (hab(2,k,i)+hab(3,k,i)-rmax**2)
          hab(3,k,i) = hab(2,k,i)+rmax**2

          sab(1,k,i) = m11*s00*m11 + m11*s01*m12 + m12*s10*m11 + m12*s11*m12
          sab(2,k,i) = m11*s00*m21 + m11*s01*m22 + m12*s10*m21 + m12*s11*m22
          sab(3,k,i) = m21*s00*m11 + m21*s01*m12 + m22*s10*m11 + m22*s11*m12
          sab(4,k,i) = m21*s00*m21 + m21*s01*m22 + m22*s10*m21 + m22*s11*m22

          vab(1,k,i) = m11*v00*m11 + m11*v01*m12 + m12*v10*m11 + m12*v11*m12
          vab(2,k,i) = m11*v00*m21 + m11*v01*m22 + m12*v10*m21 + m12*v11*m22
          vab(3,k,i) = m21*v00*m11 + m21*v01*m12 + m22*v10*m11 + m22*v11*m12
          vab(4,k,i) = m21*v00*m21 + m21*v01*m22 + m22*v10*m21 + m22*v11*m22

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
             szz = szz - phz*(suz+szu) - dphz*(ssz+szs) + phz**2*sab(1,k,i) + &
                  phz*dphz*(sab(2,k,i)+sab(3,k,i)) + dphz**2*sab(4,k,i)
             suz = suz - phz*sab(1,k,i) - dphz*sab(2,k,i)
             ssz = ssz - phz*sab(3,k,i) - dphz*sab(4,k,i)
             szu = szu - phz*sab(1,k,i) - dphz*sab(3,k,i)
             szs = szs - phz*sab(2,k,i) - dphz*sab(4,k,i)

             vuz = m11*v0z + m12*v1z
             vsz = m21*v0z + m22*v1z
             vzu = vz0*m11 + vz1*m12
             vzs = vz0*m21 + vz1*m22
             vzz = vzz - phz*(vuz+vzu) - dphz*(vsz+vzs) + phz**2*vab(1,k,i) + &
                  phz*dphz*(vab(2,k,i)+vab(3,k,i)) + dphz**2*vab(4,k,i)
             vuz = vuz - phz*vab(1,k,i) - dphz*vab(2,k,i)
             vsz = vsz - phz*vab(3,k,i) - dphz*vab(4,k,i)
             vzu = vzu - phz*vab(1,k,i) - dphz*vab(3,k,i)
             vzs = vzs - phz*vab(2,k,i) - dphz*vab(4,k,i)

             huz = m11*h0z + m12*h1z
             hsz = m21*h0z + m22*h1z
             hzu = hz0*m11 + hz1*m12
             hzs = hz0*m21 + hz1*m22
             hzz = hzz - phz*(huz+hzu) - dphz*(hsz+hzs) + phz**2*hab(1,k,i) + &
                  phz*dphz*(hab(2,k,i)+hab(3,k,i)) + dphz**2*hab(4,k,i)
             huz = huz - phz*hab(1,k,i) - dphz*hab(2,k,i)
             hsz = hsz - phz*hab(3,k,i) - dphz*hab(4,k,i)
             hzu = hzu - phz*hab(1,k,i) - dphz*hab(3,k,i)
             hzs = hzs - phz*hab(2,k,i) - dphz*hab(4,k,i)

             ! ... New gz val,slo=0 => hamiltonian is hermitian
             if (lpzi(l) == 1) then
                hzu = (hzu + huz)/2
                huz = hzu
                hzs = (hzs + hsz)/2
                hsz = hzs
             endif

             ! ... hab(5)=uz    hab(2)=sz    hab(7)=zz
             !     print *, 'zero out potpus local orbitals'
             sab(5,k,i) = suz
             sab(6,k,i) = ssz
             sab(7,k,i) = szz
             sab(8,k,i) = szu
             sab(9,k,i) = szs

             vab(5,k,i) = vuz
             vab(6,k,i) = vsz
             vab(7,k,i) = vzz
             vab(8,k,i) = vzu
             vab(9,k,i) = vzs

             hab(5,k,i) = huz
             hab(6,k,i) = hsz
             hab(7,k,i) = hzz
             hab(8,k,i) = hzu
             hab(9,k,i) = hzs

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
             xxx = (phz*sab(1,k,i) + dphz*sab(2,k,i))
             yyy = (phz*sab(3,k,i) + dphz*sab(4,k,i))
             ppnl(9,k,i)  = s0z - phi *xxx - dphi *yyy
             ppnl(10,k,i) = s1z - phip*xxx - dphip*yyy
          endif
10     enddo lloop
80  enddo isploop
    !         ! --- Printout ---
    !         if (ipr >= 40) then
    !            write(stdo,301)
    ! 301        format(/' l,E',6x,'  a      <phi phi>       a*D', &
    !                 '        a*Ddot      phi(a)      phip(a)')
    !            do  l = 0, lmxa
    !               write(stdo,311) l,rmax,(ppnl(k,l+1,i),k=2,6)
    ! 311           format(i2,2x,2f12.6,f13.6,10f12.6)
    !            enddo
    !         endif

    !         if (ipr >= 60) then
    !            write(stdo,810)
    ! 810        format(/'  l',9x,'val-val',5x,'val-slo',5x,'slo-val', &
    !                 5x,'slo-slo')
    !            do   l = 0, lmxa
    !               k  = l+1
    !               write(stdo,811) l,'s',(sab(ir,k,i),ir=1,4)
    !               write(stdo,812)   'h',(hab(ir,k,i),ir=1,4)
    !               write(stdo,812)   'v',(vab(ir,k,i),ir=1,4)
    !            enddo
    ! 811        format(i3,3x,a1,6f12.6)
    ! 812        format(3x,3x,a1,6f12.6)

    !            if (lpz /= 0) then
    !               write(stdo,815)
    ! 815           format(/'  l',9x,' val-sc',5x,' slo-sc',5x,' sc-sc ', &
    !                    5x,'sc-val',6x,'sc-slo')
    !               do   l = 0, lmxa
    !                  k  = l+1
    !                  if (pnz(k,1) > 0) then
    !                     write(stdo,811) l,'s',(sab(ir,k,i),ir=5,9)
    !                     write(stdo,812)   'h',(hab(ir,k,i),ir=5,9)
    !                     write(stdo,812)   'v',(vab(ir,k,i),ir=5,9)
    !                  endif
    !               enddo

    !            endif
    !         endif
    !80   enddo

    ! ... Calculate spin-orbit parameters
    if (lso /= 0) then
       call soprm(5,lpzi,psi,dpsi,pzi,nr,nsp,lmxa,lmxa,v,dv,enumx, &
            ezum,z,rofi,rwgt,wrk,sop,sopz)

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
             sodb(1,k,i,1) =m11*v00*m11+m11*v01*m12+m12*v10*m11+m12*v11*m12
             sodb(2,k,i,1) =m11*v00*m21+m11*v01*m22+m12*v10*m21+m12*v11*m22
             sodb(3,k,i,1) =m21*v00*m11+m21*v01*m12+m22*v10*m11+m22*v11*m12
             sodb(4,k,i,1) =m21*v00*m21+m21*v01*m22+m22*v10*m21+m22*v11*m22

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
                     phz**2*sodb(1,k,i,1) + &
                     phz*dphz*(sodb(2,k,i,1)+sodb(3,k,i,1)) + &
                     dphz**2*sodb(4,k,i,1)
                vuz = vuz - phz*sodb(1,k,i,1) - dphz*sodb(2,k,i,1)
                vsz = vsz - phz*sodb(3,k,i,1) - dphz*sodb(4,k,i,1)
                vzu = vzu - phz*sodb(1,k,i,1) - dphz*sodb(3,k,i,1)
                vzs = vzs - phz*sodb(2,k,i,1) - dphz*sodb(4,k,i,1)
                sodb(5,k,i,1) = vuz
                sodb(6,k,i,1) = vsz
                sodb(7,k,i,1) = vzz
                sodb(8,k,i,1) = vzu
                sodb(9,k,i,1) = vzs
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
          sodb(1,k,1,2)=m11*v00*x11+m11*v01*x12+m12*vx10*x11+m12*v11*x12
          sodb(2,k,1,2)=m11*v00*x21+m11*v01*x22+m12*vx10*x21+m12*v11*x22
          sodb(3,k,1,2)=m21*v00*x11+m21*v01*x12+m22*vx10*x11+m22*v11*x12
          sodb(4,k,1,2)=m21*v00*x21+m21*v01*x22+m22*vx10*x21+m22*v11*x22
          !     ... down-up block
          sodb(1,k,2,2) = x11*vx00*m11+x11*v01*m12+x12*vx10*m11 &
               +x12*vx11*m12
          sodb(2,k,2,2) = x11*vx00*m21+x11*v01*m22+x12*vx10*m21 &
               +x12*vx11*m22
          sodb(3,k,2,2) = x21*vx00*m11+x21*v01*m12+x22*vx10*m11 &
               +x22*vx11*m12
          sodb(4,k,2,2) = x21*vx00*m21+x21*v01*m22+x22*vx10*m21 &
               +x22*vx11*m22

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
                  +  sodb(1,k,1,2)*phz*phz2 &
                  +  sodb(2,k,1,2)*phz*dphz2 &
                  +  sodb(3,k,1,2)*dphz*phz2 &
                  +  sodb(4,k,1,2)*dphz*dphz2
             vxzz = vxzz - vxzu*phz - vxzs*dphz - vxuz*phz2 - vxsz*dphz2 &
                  +  sodb(1,k,2,2)*phz2*phz &
                  +  sodb(2,k,2,2)*phz2*dphz &
                  +  sodb(3,k,2,2)*dphz2*phz &
                  +  sodb(4,k,2,2)*dphz2*dphz
             vuz  = vuz  - sodb(1,k,1,2)*phz2 - sodb(2,k,1,2)*dphz2
             vxuz = vxuz - sodb(1,k,2,2)*phz  - sodb(2,k,2,2)*dphz
             vzu  = vzu  - sodb(1,k,1,2)*phz  - sodb(3,k,1,2)*dphz
             vxzu = vxzu - sodb(1,k,2,2)*phz2 - sodb(3,k,2,2)*dphz2
             vsz  = vsz  - sodb(3,k,1,2)*phz2 - sodb(4,k,1,2)*dphz2
             vxsz = vxsz - sodb(3,k,2,2)*phz  - sodb(4,k,2,2)*dphz
             vzs  = vzs  - sodb(2,k,1,2)*phz  - sodb(4,k,1,2)*dphz
             vxzs = vxzs - sodb(2,k,2,2)*phz2 - sodb(4,k,2,2)*dphz2
             sodb(5,k,1,2) = vuz
             sodb(6,k,1,2) = vsz
             sodb(7,k,1,2) = vzz
             sodb(8,k,1,2) = vzu
             sodb(9,k,1,2) = vzs
             sodb(5,k,2,2) = vxuz
             sodb(6,k,2,2) = vxsz
             sodb(7,k,2,2) = vxzz
             sodb(8,k,2,2) = vxzu
             sodb(9,k,2,2) = vxzs
          endif
       enddo
       !  ... Print the so radial integrals of phi,phidot
       !        do   l = 0, lmxa
       !          do i = 1, nsp
       !            do j = 1, nsp
       !              print*, l, i, j ,(sop(l,j,i,ir),ir=1,3)
       !              write(stdo,811) l,'o',(sop(ir,j,i,2),ir=1,3)
       !              write(stdo,811) l,'o',(sop(ir,j,i,3),ir=1,4)
       !            enddo
       !          enddo
       !        enddo
       !        stop

       !   ... Write the so radial integrals of (u,s)
       !        do i = 1, nsp
       !          do   l = 0, lmxa
       !            k  = l+1
       !            write(stdo,811) l,'d',(sodb(ir,k,i,1),ir=1,4)
       !            write(stdo,811) l,'o',(sodb(ir,k,i,2),ir=1,4)
       !          enddo
       !        enddo
       !        stop
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
    !      print 333, (avg+diff)/2,(avg+diff)/2-Tfg,
    !     .           (avg-diff)/2,(avg-diff)/2-Tgf
    !  333 format(' Tfg now',f12.6,'  change=',f12.6,2x,
    !     .       ' Tgf now',f12.6,'  change=',f12.6)
  end subroutine pvpus1

  subroutine soprm(mode,lpzi,phi,phid,phiz,nr,nsp,lmxs,lmx,v,dv,enu, &
       ez,z,ri,wi,wk,sop,sopz)
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
    !     implicit none
    ! ... Passed parameters
    integer :: lmx,mode,lpzi(0:lmx),lmxs,nr,nsp
    double precision :: z, &
         phi(nr,0:lmxs,nsp),phid(nr,0:lmxs,nsp),phiz(nr,0:lmxs,nsp), &
         ri(nr),sop(0:lmxs,nsp,nsp,3),v(nr,nsp),wk(nr,4),dv(nr), &
         wi(nr),sopz(0:lmxs,nsp,nsp,3),enu(0:8,nsp),ez(0:8,nsp)
    ! ... Local parameters
    integer :: l,ir,is,is1,is2,ipr,mode0,lmin
    double precision :: c,pa,r,r1,r2,dot3,vavg,eavg,eavgz,dva,xx, &
         xxz,xxavg,wkz(nr,4)
    ! ... External calls
    external daxpy,dscal,getpr
    !     Speed of light, or infinity in nonrelativistic case
    common /cc/ c
    data pa /1d0/

    ! --- Setup ---
    call getpr(ipr)
    !      stdo = lgunit(1)
    mode0 = mod(mode,4)
    !     c = 274.071979d0
    lmin = 1
    if (mode0 == 2) lmin=0

    ! --- Orthonormalize phi, phidot, neglecting small component ---
    if (mode/4 /= 0) then
       if (ipr > 50) &
            print '(/'' soprm: overlaps  phi*phi     phi*phidot'')'
       do   is = 1, nsp
          do   l = 0, lmx
             r1 = dot3(nr,phi(1,l,is),phi(1,l,is),wi)
             call dscal(nr,1/dsqrt(r1),phi(1,l,is),1)
             r2 = dot3(nr,phi(1,l,is),phid(1,l,is),wi)
             call daxpy(nr,-r2,phi(1,l,is),1,phid(1,l,is),1)
             if (ipr > 50) write(stdo,334) is,l,r1,r2
334          format('  spin',i2,'  l=',i1,2f13.6)
          enddo
       enddo
    endif
    if (mode0 == 0) return

    ! --- Matrix elements for each l ---
    do  is = 1, 4
       wk(1,is) = 0d0
       wkz(1,is) = 0d0
    enddo
    if (mode0 == 2) then
       do  ir = 2, nr
          wk(ir,1) = wi(ir) * dv(ir)
       enddo
    endif

    ! .. Initialize matrix elements for s orbitals, in case not calculated
    l = 0
    do  is2 = 1, nsp
       do  is1 = 1, nsp
          do is = 1, 3
             sop(l,is1,is2,is) = 0d0
             if (lpzi(l) /= 0) then
                sopz(l,is1,is2,is) = 0d0
             endif
          enddo
       enddo
    enddo

    do  l = lmin, lmx
       eavg = (enu(l,1)+enu(l,nsp))/2
       if (lpzi(l) /= 0) then
          eavgz = (ez(l,1)+ez(l,nsp))/2
       endif
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
          write(stdo,333) '          ', &
               l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
          if (lpzi(l) /= 0) then
             write(stdo,333) '          ', &
                  l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
          endif
       enddo
    else
       do  l = lmin, lmx
          write(stdo,333) 'up   up   ', &
               l,sop(l,1,1,1),sop(l,1,1,2),sop(l,1,1,3)
          write(stdo,333) 'down down ', &
               l,sop(l,2,2,1),sop(l,2,2,2),sop(l,2,2,3)
          write(stdo,333) 'up   down ', &
               l,sop(l,1,2,1),sop(l,1,2,2),sop(l,1,2,3)
          write(stdo,333) 'down up   ', &
               l,sop(l,2,1,1),sop(l,2,1,2),sop(l,2,1,3)
          write(stdo,333)
          if (lpzi(l) /= 0) then
             write(stdo,335) 'up   up   ', &
                  l,sopz(l,1,1,1),sopz(l,1,1,2),sopz(l,1,1,3)
             write(stdo,335) 'down down ', &
                  l,sopz(l,2,2,1),sopz(l,2,2,2),sopz(l,2,2,3)
             write(stdo,335) 'up   down ', &
                  l,sopz(l,1,2,1),sopz(l,1,2,2),sopz(l,1,2,3)
             write(stdo,335) 'down up   ', &
                  l,sopz(l,2,1,1),sopz(l,2,1,2),sopz(l,2,1,3)
             write(stdo,335)
          endif
       enddo
    endif
    !  333 format(1x,a,i3,3f20.15)
333 format(1x,a,i3,1x, 3f14.8)
335 format(1x,a,i3,'l',3f14.8)
  end subroutine soprm
end module m_potpus
