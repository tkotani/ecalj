module m_atwf !- Make properties related to core for one sphere
  use m_lgunit,only:stdo
  use m_ftox
  public atwf,rwftai,makrwf,phidx
  private
contains
  subroutine atwf(mode,a,lmxa,nr,nsp,pnu,pnz,rsml,ehl,rmt,z,v0,nphimx,ncore,konfig,ecore,gcore,gval,nmcore)
    use m_rhocor,only: getcor
    use m_lmfinit,only: n0
    !i   mode  :0 return ncore, and konfig, and nphimx only;
    !i         :  see description below for contents of nphimx
    !i         :1s digit
    !i         :1 return valence wave functions
    !i         :2 return core wave functions
    !i         :3 combination of 1+2
    !i       using large component only
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   lmxa  :augmentation l-cutoff
    !i   nr    :number of radial mesh points
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
    !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !i   pnz   :pnu por local orbitals
    !i   rmt   :MT boundary
    !i   z     :nuclear charge      (not used if mode=0)
    !i   v0    :spherical potential (not used if mode=0)
    !i   ehl   :energy of smoothed Hankel tail for extended local orbital
    !i   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
    ! o Inputs/Outputs
    ! o  nphimx:dimensions gval.  Must be at least as large as the
    ! o        :number of valence wave functions
    ! o        :For mode=0, nphimx is output and is assigned to
    !i         :maximum number radial wave functions for any l channel.
    !o Outputs
    !o   ncore :number of core levels
    !o   konfig:1s digit contains core configuration
    !o         :10s digit:
    !o         : 0 -> no local orbitals
    !o         : 1 -> local orbital with p.q.n. < pnu
    !o         : 2 -> local orbital with p.q.n. > pnu
    !o   ... The following are not used if mode=0
    !o   ecore :core eigenvalues
    !o   gcore :core wave functions
    !o   gval  :valence wave functions
    !o          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
    !o            i=0 : phi
    !o            i=1 : phidot
    !o            i=2 : local orbital
    implicit none
    integer :: mode,nr,nsp,lmxa,ncore,konfig(1+lmxa),nrmx,nphimx
    real(8):: rmt,z,a,v0(nr,nsp),pnu(n0,nsp),pnz(n0,nsp),gval(nr*2,0:lmxa,nphimx,nsp),ecore(*),gcore(nr,2,*),rsml(n0),ehl(n0)
    logical :: lpz
    integer :: l,isp,konf,konfz,k,mode0,  nmcore !,mode1
    real(8) :: sumtc,sumec,e,ez,xx
    real(8) :: rofi(nr),rwgt(nr),rhoc(nr,2),gp(2*nr*4)
    real(8) :: phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp
    logical:: l_dummy_isanrg
    mode0 = mod(mode,10)
    ! --- Count number of core states ---
    SetCoreIndex: block
      logical:: isanrg
      lpz = .false.
      ncore = 0
      do  l = 0, lmxa
         k = l+1
         konfig(k) = pnu(k,1)
         konfz = mod(pnz(k,1),10d0)
         if (konfz == 0) konfz = konfig(k)
         l_dummy_isanrg = isanrg(konfz,konfig(k)-1,konfig(k)+1,'atwf:','pnuz',.true.)
         do  konf = l+1, min(konfz,konfig(k))-1
            ncore = ncore+nsp
         enddo
         if (konfz < konfig(k)) then
            konfig(k) = konfz + 10
            lpz = .true.
         elseif (konfz > konfig(k)) then
            konfig(k) = konfig(k) + 20
            lpz = .true.
         endif
      enddo
    endblock SetCoreIndex
    if (mode0 == 0) then
       nphimx = 2
       if (lpz) nphimx = 3
       return
    endif
    call radmsh(rmt,a,nr,rofi)
    call radwgt(rmt,a,nr,rwgt)
    ValenceWaveFunctions: if (mod(mode0,2) == 1) then
       do  l = 0, lmxa
          k = l+1
          do  isp = 1, nsp
             konf = pnu(k,1)
             call makrwf(z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnu(1,isp),4, gval(1,l,1,isp),gp,e,phi,dphi,phip,dphip,p)
             gval(:,l,2,isp)=gp(1:2*nr) 
             if (konf /= konfig(k)) then
                call makrwf(z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnz(1,isp),2, gval(1,l,3,isp),gp,ez,phz,dphz,phzp,dphzp,p)
                call wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz, &
                     phzp,dphzp,pnz(1,isp),rsml,ehl, gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp))
             endif
          enddo
       enddo
    endif ValenceWaveFunctions
    if(mode0>=2) call getcor(1,z,a,pnu,pnz,nr,lmxa,rofi,v0,0,0,[0d0,0d0],sumec,sumtc,rhoc,ncore,ecore,gcore,nmcore) !nmcore jun2012 !Core eigenfunctions and eigenvalues --- 
  end subroutine atwf
  subroutine makrwf(z,rmax,l,v,a,nr,rofi,pnu,nptdif,g,gp, enu,phi,dphi,phip,dphip,p)!Radial wave functions and energy derivative
    !i   z     :nuclear charge
    !i   rmax  :augmentation radius, in a.u.
    !i   l     :l-quantum number
    !i   v     :spherical potential
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax, pnu = .5 - atan(Dl)/pi + (princ.quant.number).
    !o Outputs
    !o   g     :r * wave function corresponding to b.c. pnu
    !o   gp    :r * energy derivative of g; dimensioned 8*nr
    !o   enu   :eigenvalue of g
    !o   phi   :wave function at rmax, i.e. g/rmax
    !o   dphi  :slope of wave function at rmax, i.e. (d(g/r)/dr)_rmax
    !o   phip  :energy derivative of wave function at rmax
    !o   dphip :energy derivative of slope of wave function at rmax
    !o   p     :<gp**2> (potential parameter)
    !r Remarks
    !r   This routine makes r*phi and r*phidot, where phi and phidot are true radial wave function and energy derivatives.
    !r   phi is normalized, and p = <phidot**2>
    implicit none
    integer :: mode,l,nr,nptdif, konf,nn,nre,modep
    real(8):: a,rmax,z,rofi(1),v(nr,1),pnu(1:l+1),g(nr,2),gp(nr,2,4),phi,phip,dphi,dphip,p,b,dnu,eb1,eb2,enu,slo(5),sum,val(5)
    real(8),parameter:: pi = 4d0*datan(1d0), tol = 1d-12
    if(abs(rmax-rofi(nr))>1d-8) call rx('makrwf: rmax/=rofi(nr)')
    nn = mod(pnu(l+1),10d0)-l-1 !number of nodes
    val(1) = rmax
    slo(1) = dtan(pi*(.5d0-mod(pnu(l+1),10d0))) + 1d0
    call phidx(0,z,l,v,rofi,nr,nptdif,tol,enu,val,slo,nn,g,gp,phi,dphi,phip,dphip,p)!,0d0,[0d0],0d0,[0d0])
  end subroutine makrwf
  subroutine rwftai(rmt,a,nrmt,nrbig,ribig,phi,dphi,tphi,l, ehl,rsml,g)!Extend radial wave function outside MT boundary
    use m_hansr,only :hansr,hansmr,hansmronly
    ! Compute radial wave function on extended mesh using rsml,ehl for tail, scale gz so that value matches envelope h(rsm,eh)
    !i Inputs 
    !i   rmt   :augmentation radius, in a.u., by species
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nrmt  :number of radial mesh points between origin and rmt
    !i   l     :quantum number for this wave function
    !i   ehl   :energy of smoothed Hankel tail for local orbital
    !i   rsml  :smoothing radius for smoothed Hankel tail of local orbital
    !i  nrbig  :number of radial mesh points on extended mesh
    !i  ribig  :points for extended radial mesh
    !o Outputs
    !o   g      :radial wave function extended to points nrmt..nrbig
    !l Local variables
    !l   rbig   :rmax of extended mesh
    !l   nxtn   :number of points between rmt and rbig
    implicit none
    integer :: mode,nrmt,nrbig,l,idn,nxtn,ir,info,mode0,i
    integer,parameter:: nrx=1501,IPRTDB=10
    real(8):: a,rmt,phi,dphi,tphi,ribig(*),g(nrbig,2),ehl,rsml,rbig,rwgt(nrx),xi(nrx,0:l),r2(nrx),fac1,fac2,rsm,ekin,eh,alfa,beta
    rbig = ribig(nrbig)
    nxtn = nrbig-nrmt+1
    call radwgt(rbig,a,nrbig,rwgt)
    rwgt(nrmt) = rwgt(nrmt)/2
    r2(1:nxtn) = [(ribig(ir-1+nrmt)**2,ir=1,nxtn)]
    if(hansmronly) then
       do ir=1,nxtn
          call hansmr(r2(ir)**.5,ehl,1d0/rsml,xi(ir,0:l),l)
          xi(ir,0:l)=[(xi(ir,i)*(r2(ir)**.5)**i,i=0,l)]
       enddo
    else   
       call hansr(rsml,0,l,1,[l],[ehl],[r2],nrx,nxtn,[idn],11,xi)
    endif   
    fac1 = g(nrmt,1)/rmt/xi(1,l) ! Hankel function scaled to match g at nrmt
    fac2 = g(nrmt,2)/rmt/xi(1,l)
    g(1:nrmt,:)=g(1:nrmt,:)/fac1 
    fac2 = fac2 / fac1
    fac1 = 1
    g(nrmt:nrmt+nxtn-1,1) = [(xi(ir,l)*ribig(ir-1+nrmt) * fac1,ir=1,nxtn)]
    g(nrmt:nrmt+nxtn-1,2) = [(xi(ir,l)*ribig(ir-1+nrmt) * fac2,ir=1,nxtn)]
  end subroutine rwftai
  subroutine wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz,phzp,dphzp,pnz,rsml,ehl,g0,g1,gz)! Add a linear combination of two w.f. to a 3rd to make a local orbital
    !i Inputs
    !i   l     :l quantum number
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   rofi  :radial mesh points
    !i   rwgt  :radial mesh weights
    !i   phi   :1st wave function at rmax, i.e. r*g0
    !i   dphi  :radial derivative of r*g0
    !i   phip  :2nd wave function at rmax, i.e. r*g1
    !i   dphip :radial derivative of r*g1
    !i   g0    :1st radial w.f. to which to orthogonalize gz
    !i   g1    :2nd radial w.f. to which to orthogonalize gz
    !i   ehl   :energy of smoothed Hankel tail for extended local orbital
    !i   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
    ! o Inputs/Outputs
    ! o  phz   :3rd wave function at rmax, i.e. r*gz
    ! o        :Input for standard local orbital; scaled for extended loc. orb
    ! o  dphz  :radial derivative of r*gz
    ! o        :Input for standard local orbital; scaled for extended loc. orb
    ! o  phzp  :energy derivative of phz
    ! o        :scaled for extended local orbital
    ! o  dphzp :energy derivative of dphz
    ! o        :scaled for extended local orbital
    ! o  gz    :on input, radial w.f.
    ! o        :on output, gz is overwritten by:
    ! o        :(gz - alpha g0 - beta g1) so that value and slope
    ! o        :of the result are zero at rmt, standard local orbital
    ! o        :scaled to match sm. hankel at rmt, extended loc. orb.
    !u Updates
    !u   04 Sep 04 Adapted to extended local orbitals
    !u   06 Mar 02 New routine
    ! ----------------------------------------------------------------------
    implicit none
    integer :: k,lpzi,l,nr
    real(8) :: a,rofi(nr),rwgt(nr),rsml(*),ehl(*),phi,dphip,dphi,phip,phz,dphz,phzp,dphzp,pnz(*), g0(nr,2),g1(nr,2),gz(nr,2),&
         det,au,bu,as,bs,fac,x,xx,gzbig(nr,2)
    k = l+1
    lpzi = 0
    if (pnz(k) >  0) lpzi = 1
    if (pnz(k) >= 10) lpzi = 2
    if (pnz(k) >= 20) lpzi = 3
    if (lpzi == 0) then
       return
    elseif (lpzi == 1) then
       det = phi*dphip - dphi*phip
       au = dphip/det
       bu = -dphi/det
       as = -phip/det
       bs = phi/det
       gz(:,:)=gz(:,:)-(phz*au + dphz*as)*g0(:,:)
       gz(:,:)=gz(:,:)-(phz*bu + dphz*bs)*g1(:,:)
    elseif (lpzi == 2 .OR. lpzi == 3) then
       gzbig=gz
       call rwftai(rofi(nr),a,nr,nr,rofi,phz,dphz,xx,l, ehl(k),rsml(k),gzbig) 
       if(gzbig(nr,1) /= gz(nr,1)) then !If rwftai scales gzbig, rescale phz,gz
          xx = gzbig(nr,1)/gz(nr,1)
          phz   = phz*xx
          dphz  = dphz*xx
          phzp  = phzp*xx
          dphzp = dphzp*xx
          gz=gz*xx
          !call dscal(nr,xx,gz(1,1),1)
          !call dscal(nr,xx,gz(1,2),1)
       endif
    endif
  end subroutine wf2lo
  subroutine phidx(job,z,l,v,rofi,nr,nptdif,tol, e,val,slo,nn, g,gp,phi,dphi,phip,dphip,p)!Generate potential parameters for a prescribed energy or b.c.
    use m_rseq,only: rseq,rsq1
    !i Inputs:
    !i   job: 0 or 2
    !i      0, boundary conditions specified val,slo,nn 
    !i      2, boundary conditions specified by energy e 
    !ixxx          10s digit:    1, set dphip to satisfy Wronskian condition
    !i   z:     nuclear charge
    !i   l:     l quantum number for this g
    !i   v:     spherical potential on shifted logarithmic mesh
    !i   rofi:  list of points
    !i   nr:    number of mesh points
    !i   nptdif:2 or 4 for 3- or 5- point differentiation
    !i         :You may also set nptdif=0.  Then quantities related to
    !i         :energy differences are not calculated (phip,dphip,p)
    !i   tol:   precision to which wave function is integrated
    ! We supply e for job=2, while we supply val,slo,nn for job=0
    !io     e:   energy eigenvalue input for job=2, while output for job=0
    !io    val: val(1)=g(r)=r*phi(r) at rmax input for job=0, while output val(1)=value of normalized g(r) at rmax.
    !           Also on output, val(2..1+nptdif)=energy derivatives of val
    !io    slo: slo(1)=radial derivative of g(r) at rmax input for job=0, while output slo(1)=der. of normalized g(r) at rmax.
    !           Also on output, slo(2..1+nptdif) = energy derivatives of slo
    !io     nn:   number of nodes for input for job=0, while output for job=2.
    !o  Outputs:
    !o   g:     normalized wave function times r
    !o   gp:    first nptdif energy derivatives to G
    !o   phi:   wave function at rmax, i.e. g/rmax
    !o   dphi  :slope of wave function at rmax, i.e. (d(g/r)/dr)_rmax
    !o   phip:  energy derivative of wave function at rmax
    !o   dphip: energy derivative of slope of wave function at rmax
    !o   p:     <gp**2> (potential parameter)
    !r Remarks:
    !r   This version makes parameters related to wave function g(r)/r defined by potential v(r).
    !r
    !r   Boundary conditions are specified in one of two ways:
    !r     job=0   val,slo,nn are specified.  On output, val,slo are renormalized so that val=g(nr), with <gg>=1
    !r             Here energy eigenvalue e is calculated are assumed to correspond with g.
    !r     job=2   the energy eigenvalue e is specified. val,slo, and nn are calculated.
    !NOTE:
    !b      phi   = val(1)/rmax
    !b      dphi  = (slo(1) - phi)/rmax
    !b      phip  = vali(1)/rmax
    !b      dphip = (sloi(1) - phip)/rmax
    implicit none
    integer :: job,l,nr,nn,nptdif,nre,i,iprint
    real(8)::z,e,val(*),slo(*),phi,dphi,phip,dphip,p,tol,v(nr),rofi(nr),g(2*nr),gp(2*nr,4),& !vmtz,hcr,phia,phiap(*)dla,dlap(*)
         rmax,eb1,eb2,dele,ddde,sum,a,b,aold,dmach,tola,sloi(5),vali(5),phiai(4),ei(4),de1,de2,del1,del2,tole!dlax(1),,dlai(4)
    parameter (tole=1d-10)
    rmax = rofi(nr)
    dele = .002d0
    a = log(rofi(nr)/rofi(nr-1))
    tola = 8*dmach(1)
    do i = 1, 100
       aold = a
       a = log(rofi(nr)/rofi(nr-1)/(1-exp(-(nr-1)*a))*(1-exp(-(nr-2)*a)))
       if (i > 95) write(stdo,'(i4,1p,2e15.5)') i,a,a-aold
       if (abs(a-aold) <= tola) goto 1
    enddo
    call rx('phidx failed to determine ''a'' parameter')
1   continue
    b = rmax/(dexp(a*(nr-1)) - 1d0)
    if (job == 0) then ! --- Find energy e corresponding to val,slo ---
       eb1 = -20d0
       eb2 =  20d0
       e = (eb1+eb2)/2
       call rseq(eb1,eb2,e,tol,z,l,nn,val(1),slo(1),v,g,sum,a,b,rofi,nr,nre) !       This generates g, normalized to <g g> = 1
       val(1) = val(1)/dsqrt(sum) !Scale val, slo to correspond to normalization <g g> = 1
       slo(1) = slo(1)/dsqrt(sum)
    elseif (job == 2) then ! --- Find val,slo corresponding to energy e --- !       Initial estimate for val,slo
       call rsq1(0,e,l,z,v,nr,g,val(1),slo(1),nn,a,b,rofi,nr)  !Adjust slope iteratively until ei(slope)-e is within tole
       ei(1) = e
       eb1 = e-.1d0
       eb2 = e+.1d0
       do  i = 1, 5+5
          call rseq(eb1,eb2,ei(1),tol,z,l,nn,val(1),slo(1),v,g,sum,a,b,rofi,nr,nre)
          if (iprint()> 0.AND.i>5) &
               write(stdo,"('PHIDX Z=',f9.4,'  l nod=',i0,x,i0,' val slo=',2f9.4,' e(bc)=',f9.4,' e(bc)-e=',f9.4)") z,l,nn,val(1),&
               slo(1),ei(1),ei(1)-e
          if (abs(ei(1)-e) < tole) goto 2
          slo(1) = slo(1) + (ei(1)-e) * val(1) / g(nr)**2
       enddo
       call rx('phidx failed to converge')
2      continue
       val(1) = val(1)/dsqrt(sum)
       slo(1) = slo(1)/dsqrt(sum)
    else
       call rx('phidx: bad job')
    endif
!    if (hcr /= 0) call makdla(e-vmtz,l,hcr,slo(1),val(1),rmax,phia,dla)
    if (nptdif /= 0) then
       ddde = -rmax/g(nr)**2
       ei(1:4) = [1d0, -1d0, 1.5d0, -1.5d0]
       eb1 = e-.1d0
       eb2 = e+.1d0
       ! choice 2 takao. For deep semicore, to avoid error in rseq (Mg dimer in 10\A cubic cell).
       ei=ei/2d0 !In future, we have to use better solution. you may need to use ei/3.0 or better algolism
       !  I had a problem that it results in warning rseq, because two exponential solution makes huge changes due to slight difference of energy
       !  and node number can not be the same.
       !    ATOM= Mg Z= 12 R= 3.000
       !     RSMH=   1.500 1.500 1.500 1.500 EH=  -1.0 -1.0 -1.0 -1.0
       !     RSMH2=  1.500 1.500 1.500 1.500 EH2= -2.0 -2.0 -2.0 -2.0
       !     PZ=0,12.9 P=0,3.3     KMXA={kmxa}  LMX=3 LMXA=4
       ! cccccccccccccccccccccccccccccccccc
       do  10  i = 1, nptdif
          sloi(i) = slo(1) + dele*ei(i)*ddde*val(1)/rmax
          ei(i) = e + dele*ei(i)
          call rseq(eb1,eb2,ei(i),tol,z,l,nn,val(1),sloi(i),v,gp(1,i),   sum,a,b,rofi,nr,nre)
          vali(i) = val(1)/dsqrt(sum)
          sloi(i) = sloi(i)/dsqrt(sum)
!          if (hcr /= 0) call makdla(ei(i)-vmtz,l,hcr,sloi(i),vali(i), rmax,phiai(i),dlai(i))
10     enddo
       de1  = (ei(1) - ei(2))/2
       del1 = (ei(1) + ei(2))/2 - e
       de2  = (ei(3) - ei(4))/2
       del2 = (ei(3) + ei(4))/2 - e
       call dfphi(de1,del1,de2,del2,1,val,vali,nptdif.eq.4) !Energy derivatives of value and slope
       call dfphi(de1,del1,de2,del2,1,slo,sloi,nptdif.eq.4)
       call dfphi(de1,del1,de2,del2,2*nr,g,gp,nptdif.eq.4)!     Energy derivatives of g
       call gintsr(gp,gp,a,b,nr,z,e,l,v,rofi,p) !     p = integral <gp gp>
    endif
    phi = val(1)/rmax  !   phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
    dphi = (slo(1) - phi)/rmax
    if (nptdif /= 0) then
       phip  = vali(1)/rmax
       dphip = (sloi(1) - phip)/rmax
       call dcopy(nptdif,sloi,1,slo(2),1)!     Copy energy derivatives sloi to slo(2..)
       call dcopy(nptdif,vali,1,val(2),1)
    endif
    ! if(nptdif/=0.AND.job1/= 0) dphip=(dphi*phip-1d0/rmax**2)/phi !Set dphip from Wronskian condition:
    ! phi*dphip - dphi*phip = -1/rmax**2 => dphip = (dphi*phip - 1/rmax**2)/phi
    if (iprint() >= 111) write(stdo,"(' PHIDOT:  phi,phip,phip,dphip=',4f12.6)") phi,dphi,phip,dphip
  end subroutine phidx
  subroutine gintsr(g1,g2,a,b,nr,z,e,l,v,rofi,sum)
    use m_lmfinit,only: c=>cc
    !- Integrate inner product of two wave equations
    ! ----------------------------------------------------------------
    !i   g1,g2 :First and second radial wave functions
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nr    :number of radial mesh points
    !i   z     :nuclear charge
    !i   e     :energy
    !i   l     :l quantum number of g1,g2
    !i   v     :spherical potential
    !i   rofi  :radial mesh points
    !o Outputs:
    !o   sum   :inner product
    !r Remarks:
    !r   Uses Simpson's rule
    ! ----------------------------------------------------------------
    implicit none
    integer :: nr,l
    double precision :: a,b,z,e,sum,g1(nr,2),g2(nr,2),v(nr),rofi(nr)
    integer :: i,ir
    double precision :: fllp1,r,tmc,fi
    tmc(i,r) = c - (v(i) - 2d0*z/r - e)/c
    fi(i,r) = (r+b)*(g1(i,1)*g2(i,1)*(1 + fllp1/(tmc(i,r)*r)**2) + g1(i,2)*g2(i,2))
    fllp1 = l*(l+1)
    sum = 0d0
    do  10  ir = 2, nr-1, 2
       sum = sum + fi(ir,rofi(ir))
10  enddo
    sum = 2*sum
    do  11  ir = 3, nr-2, 2
       sum = sum + fi(ir,rofi(ir))
11  enddo
    sum = (2*sum + fi(nr,rofi(nr)))*a/3
  end subroutine gintsr
  subroutine dfphi(de1,del1,de2,del2,nr,g,gp,fivep)
    !- Numerically differentiates g using five point finite difference
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   de1   :see Remarks
    !i   del1  :see Remarks
    !i   de2   :see Remarks
    !i   del2  :see Remarks
    !i   nr    :number of mesh points
    !i   g     :normalized wave function times r
    !i   fivep :if true, five-point formula; else three-point formula
    ! o Inputs/Outputs
    ! o  gp    :On input,
    ! o        :contains g at de1+del1, -de1+del1, de2+del2, -de2+del2
    ! o        :On output, contains energy derivatives of g:
    ! o        :two derivatives or four derivatives, depending on fivep
    !r Remarks
    !r   A simple five-point estimate for the numerical differentiation for
    !r   the first four derivatives of phi is based on integration of phi
    !r   at four energies de1+del1, -de1+del1, de2+del2, -de2+del2.  The
    !r   offset del1 and del2 are present because (for consistency's sake)
    !r   RSEQ is called which takes evenly spaced increments in slopes
    !r   about the average.  But the deviations del1 and del2 are only
    !r   nonzero to second order in the increment de, and this can be
    !r   exploited to obtain accurate five point estimates without having
    !r   to solve for the four simultaneous equations.
    !r
    !r   A three point for energy differences about Enu of 0, ep = de+del,
    !r   em = -de+del gives
    !r
    !r   (1)  gp =  2 del/(de^2-del^2) g(0) +
    !r              (de-del)/(de+del)/(2 de) g(ep) -
    !r              (de+del)/(de-del)/(2 de) g(em) +
    !r
    !r              1/6 gppp (de^2-del^2) + 1/12 gpppp del (de^2-del^2) + ...
    !r
    !r
    !r   (2) gpp = -2/(de^2-del^2) g(0) +
    !r              1/(de+del)/de g(ep) +
    !r              1/(de-del)/de g(em) -
    !r
    !r              2/3 gppp del + 1/12 gpppp del (de^2 + 3 del^2) + ...
    !r
    !r
    !r   The gppp term in (1) can be knocked out by taking two three point
    !r   formulas in linear combination (de1^2-del1^2) and (de2^2-del2^2)
    !r   leaving only fourth and higher order terms.  Because del is of
    !r   order de^2, the fourth order term is of the same order as the
    !r   fifth order term and there is no advantage in eliminating it.
    !r   Also from the difference between this more accurate (5-point)
    !r   estimate for gp, an estimate for gppp can be made as
    !r
    !r   (3) gppp = -6 ( gp(five point) - gp(three point) ) /(de1^2-del1^2)
    !r
    !r             + 1/2 gpppp del1
    !r
    !r   which is again accurate to order de^2.  Once gppp is known to this
    !r   order the term proportional to gppp in three point estimate for
    !r   gpp can be subtracted out directly and the gpppp term can be
    !r   eliminated by taking three point formulas (with the gppp term
    !r   removed) in linear combinations (de1^2+3*del1^2) and
    !r   (de2^2+3*del2^2), and finally the fourth order derivative can be
    !r   estimated from the difference in the five-point estimate for gpp
    !r   and the three point estimate.
    !  ----------------------------------------------------------------
    implicit none
    integer :: nr,i
    logical :: fivep
    real(8) :: de1,del1,de2,del2,g(nr),gp(nr,4), gp5p,gpp5p,gppp,gpp3p,gpp32,xx1,xx2,xx3,xx4,&
         gp3p,w01,w11,w21,w02,w12,w22, w01d,w11d,w21d,w02d,w12d,w22d,wp1d,wp2d,gpppp
    ! --- Constants common to 3-point and 5-point formulas ---
    w01 = 2*del1/(de1**2-del1**2)
    w11 = (de1-del1)/(de1+del1)/(2*de1)
    w21 = (de1+del1)/(de1-del1)/(2*de1)
    w01d = -2/(de1**2-del1**2)
    w11d = 1/(de1+del1)/de1
    w21d = 1/(de1-del1)/de1
    if ( .NOT. fivep) goto 20
    if (dabs(del1/de1) > .1 .OR.  dabs(del1/de1) > .1) goto 20
    xx1 = de1**2 - del1**2 ! --- Extra constants for 5-point formula ---
    xx2 = de2**2 - del2**2
    xx3 = de1**2 + 3*del1**2
    xx4 = de2**2 + 3*del2**2
    w02 = 2*del2/(de2**2-del2**2)
    w12 = (de2-del2)/(de2+del2)/(2*de2)
    w22 = (de2+del2)/(de2-del2)/(2*de2)
    wp1d = 2d0*del1/3
    w02d = -2/(de2**2-del2**2)
    w12d = 1/(de2+del2)/de2
    w22d = 1/(de2-del2)/de2
    wp2d = 2d0*del2/3
    do  10  i = 1, nr
       gp3p = w01*g(i) + w11*gp(i,1) - w21*gp(i,2) ! Three point formula for gp; store in temporary gp3p
       gp5p = (xx2*gp3p - xx1*(w02*g(i) + w12*gp(i,3) - w22*gp(i,4))) /(xx2-xx1) ! Five point estimate for gp
       gppp = -6/xx1*(gp5p - gp3p) ! Difference between five point and three point gives estimate for gppp
       gpp3p = w01d * g(i) + w11d * gp(i,1) + w21d*gp(i,2) -wp1d * gppp ! Three point estimates for gpp with correction for gppp
       gpp32 = w02d * g(i) + w12d * gp(i,3) + w22d*gp(i,4) -wp2d * gppp
       gpp5p = (gpp3p*xx4 - gpp32*xx3) / (xx4 - xx3) ! Five point estimate for gpp with correction for gppp
       gpppp = -12/xx3*(gpp5p - gpp3p) ! Difference between five point and three point gives est for gpppp
       gp(i,1:4) = [gp5p, gpp5p, gppp, gpppp]
10  enddo
    return
20  continue
    do  30  i = 1, nr ! Three point formulae:  only gp, gpp calculated
       gp3p    = w01*g(i)  + w11*gp(i,1)  - w21*gp(i,2)
       gp(i,2) = w01d*g(i) + w11d*gp(i,1) + w21d*gp(i,2)
       gp(i,1) = gp3p
30  enddo
  end subroutine dfphi
end module m_atwf
