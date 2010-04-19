      subroutine hmltnc(mode,nbas,nl,indxsh,qspirl,eula,neul,pph,sod,
     .  lasa,ccor,lss,lnc,lso,ccd,vmtz,ldim,lihdim,sk,hk,ok,wk)
C- Three-center ASA noncollinear hamiltonian and overlap
C ---------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 hamiltonian-> hk and overlap -> ok (sk overwritten)
Ci         :1 small h -> hk
Ci         :2 1+oh -> hk
Ci         :3 small h -> hk, input sk is noncol srdel * S * srdel
Ci         :4 1+oh -> hk, input sk is noncol srdel * S * srdel
Ci         :  modes 1-4: 
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   qspirl:Parameters defining spin-spiral
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :pph(1..5,i,is): parms for ith RL and is(th) spin channel.
Ci         :pph(1) : enu
Ci         :pph(2) : calpha
Ci         :pph(3) : sqrdel
Ci         :pph(4) : palpha
Ci         :pph(5) : oalp
Ci   sod   :diagonal parms for S/O ham., LL block.
Ci         :sod(*,1..3,isp): 1-, 2-, 3- center terms, ++,-- blocks
Ci         :sod(*,1,   3):   1- center terms, L- S+ block
Ci         :sod(*,2..3,3):   2- center terms, L- S+ block
Ci         :sod(*,4,   3):   3- center terms, L- S+ block
Ci   lasa  :switch telling whether to add ASA H
Ci   ccor  :switch telling whether to add combined correction
Ci   lss   :switch for spin spiral
Ci   lnc   :switch turning on noncollinear hamiltonian
Ci   lso   :switch turning on S/O coupling
Ci   ccd:  :diagonal matrices for 1-, 2-, 3-center terms for the
Ci         :combined corr;they are terms in parentheses in Kanpur notes
Ci         :eq. 3.87 multiplied by w^2. (see makdia.f)
Ci   vmtz  :muffin-tin zero for combined correction (asamad.f)
Ci   ldim  :dimension of the hamiltonian
Ci   lihdim:dimensions ccd and pph
Ci   sk    :structure constants, s^beta
Ci   ok    :S^beta-dot
Ci   hk    :?i-wave 3-centre integrals (see i3cntre)
Ci   wk    :work array of length ldim*2
Co Outputs
Co   hk,ok,:(mode 0) hamiltonian-> hk and overlap -> ok
Co         :(mode 1) small h -> hk
Co         :(mode 2) 1+oh -> hk
Co   sk    :(mode 0) sk is changed to noncol sqrdel*sk*sqrdel
Co         :(mode>0) sk is changed to noncol (C-enu)+sqrdel*sk*sqrdel
Cr Remarks
Cr
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | psi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   The first form uses phidot = phidot^alpha; the second form uses
Cr     phidot^alpha = phidot^gamma + o phi; and calls 'dot'=phidot^gamma
Cr   Note that <phi|phi>=1 <phidot|phi>=o <phidot|phidot>=p
Cr   Considering the LL block only, the ASA part of the overlap is:
Cr    <psi|psi>_ASA = <phi|phi> + h<phidot|phi> + h.c.
Cr                    + h <phidot|phidot> h
Cr                  = 1 + ho + oh + hph
Cr
Cr   To work directly with  D = srdel S srdel, rather
Cr   that h = C-enu + D, the diagonal parts connecting C-enu
Cr   in the one-, two-, three-center terms are reshuffled.  Thus
Cr    <psi|psi>_ASA = 1 + ho + oh + hph
Cr                  = 1 + (C-e+D)o + o(C-e+D) + (C-e+D) p (C-e+D)
Cr                  = 1 + 2(C-e)o + (C-e)^2 p      (one center)
Cr                  + D(o+p(C-e)) + (o+p(C-e))D    (two center)
Cr                  + D p D                        (three center)
Cr
Cr   The hamiltonian <psi|H|psi>_ASA has a corresponding structure
Cr   with similar 1-, 2- and 3- center terms; but the diagonal parts
Cr   are calculated from <phi|H|phi>, <phidot|H|phi>, <phidot|H|phidot>.
Cr
Cr
Cr   In the noncollinear case the same formulas apply, with the
Cr   following extensions:
Cr     1.  potential parameters have a ++ and a -- part
Cr     2.  D is a 2x2 matrix in spin space
Cr     3.  There may be a spin-orbit hamiltonian which has a structure
Cr         similar to the noncollinear hamiltonian, i.e. terms one-,
Cr         two- and three-center in D, but the diagonal bits couple
Cr         (l,m) to (l,m+/-1).
Cr     4.  Contribution from the applied field.  This is done separately
Cr         in routine hmladb.f, which see.
Cr
Cr   ?Routine untested for ldim ne lihdim (downfolding not implemented)
Cr   <k|k> = <kappa|kappa>
Cu Updates
Cu   15 Nov 07 New mode, to generate noncollinear h or 1+oh
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nl,indxsh(*),lihdim,ldim,neul,mode
      logical lasa,ccor,lss,lnc,lso
      double precision ccd(lihdim,0:2),pph(5,lihdim,2),wk(ldim,2),
     .  eula(nbas,neul,3),sk(ldim,2,ldim,2*2),hk(ldim,2,ldim,2*2),
     .  ok(ldim,2,ldim,2*2),sod(ldim,4,3),vmtz,qspirl(4)
C Local parameters
      integer i,j,l2,i1,j1,nasa,iprint,ncsw
      double precision xx

      call tcn('hmltnc')
      l2 = ldim**2*4
      nasa = 0
      if (lasa) nasa = 1

C --- sk <- noncollinear d*S*d;  hk <- noncollinear d*Sdot*d ---
      if (mode .ne. 3 .and. mode .ne. 4) then
      if (lss .and. lso) call rx('S-O coupling incompatible with SS')
      if (lnc .or. lss .or. lso) then
        ncsw = 2000
        if (lss) ncsw = ncsw + 20000
        do  8  i1 = 1, 2
        do  8  i = 1, ldim
    8   wk(i,i1) = pph(3,i,i1)
        if (ccor) call rotspn(ncsw,nbas,nl,indxsh,eula,neul,qspirl(4),
     .    wk,wk,ldim,ldim,ldim,ldim,ldim,ok,hk)
        if (.not. lss) call dcopy(ldim**2*2,sk,1,ok,1)
        if (lss) call dcopy(ldim**2*4,sk,1,ok,1)
        call rotspn(ncsw,nbas,nl,indxsh,eula,neul,qspirl(4),wk,wk,
     .   ldim,ldim,ldim,ldim,ldim,ok,sk)
C       call yprm('d*S*d',2,sk,l2,ldim*2,ldim*2,ldim*2)
      endif
      endif

C --- mode 1-4: make h or 1+oh ---
      if (mode .ge. 1 .and. mode .le. 4) then
        do  i1 = 1, 2
          do  i = 1, ldim
            sk(i,i1,i,i1) = sk(i,i1,i,i1) + pph(2,i,i1) - pph(1,i,i1)
        enddo
        enddo
C       call yprm('h',12,sk,l2,ldim*2,ldim*2,ldim*2)
        if (mode .eq. 1 .or. mode .eq. 3) return
C       Make o
        do  i1 = 1, 2
        do  i = 1, ldim
          wk(i,i1) = pph(5,i,i1)
        enddo
        enddo
C       Make oh
        do  j1 = 1, 2
          do  j = 1, ldim
            do  i1 = 1, 2
              do  i = 1, ldim
                hk(i,i1,j,j1) =    wk(i,i1)*sk(i,i1,j,j1)
                hk(l2+i,i1,j,j1) = wk(i,i1)*sk(l2+i,i1,j,j1)
              enddo
            enddo
            hk(j,j1,j,j1) = 1 + hk(j,j1,j,j1)
          enddo
        enddo
C       call yprm('1+oh',2,hk,l2,ldim*2,ldim*2,ldim*2)
        return
      endif

C --- O += d*S*d*[o+p(c-e)] + h.c. + <k|k>_linear ---
      if (.not. lasa) goto 26
      do  20  i1 = 1, 2
      do  20  i = 1, ldim
   20 wk(i,i1) = pph(5,i,i1) + pph(4,i,i1)*(pph(2,i,i1)-pph(1,i,i1))
C     call prmx('hmltnc diag. matrix for 2C O',wk,ldim,ldim,2)
      if (ccor) call daxpy(ldim,1d0,ccd(1,1),1,wk,1)
      if (ccor) call daxpy(ldim,1d0,ccd(1,1),1,wk(1,2),1)
      do  22  i1 = 1, 2
      do  22  j1 = 1, 2
      if (ccor) then
        do  23  j = 1, ldim
        do  23  i = 1, ldim
          ok(i,i1,j,j1)    = sk(i,i1,j,j1)*(wk(i,i1)+wk(j,j1)) +
     .                       hk(i,i1,j,j1)
          ok(l2+i,i1,j,j1) = sk(l2+i,i1,j,j1)*(wk(i,i1)+wk(j,j1)) +
     .                       hk(l2+i,i1,j,j1)
   23   continue
      else
        do  24  j = 1, ldim
        do  24  i = 1, ldim
          ok(i,i1,j,j1) =    sk(i,i1,j,j1)*(wk(i,i1)+wk(j,j1))
          ok(l2+i,i1,j,j1) = sk(l2+i,i1,j,j1)*(wk(i,i1)+wk(j,j1))
   24   continue
      endif
   22 continue
   26 continue

C --- H += d*S*d*[1/2+oe+(c-e)(o+pe)] + h.c. + vmtz*<k|k>_linear ---
      do  120  i1 = 1, 2
        do  121  i = 1, ldim
          wk(i,i1) = (.5d0 + pph(1,i,i1)*pph(5,i,i1) +
     .      (pph(2,i,i1) - pph(1,i,i1)) *
     .      (pph(5,i,i1) + pph(4,i,i1)*pph(1,i,i1)))*nasa
  121   continue
        if (lso)  call daxpy(ldim,1d0,sod(1,2,i1),1,wk(1,i1),1)
        if (ccor) call daxpy(ldim,vmtz,ccd(1,1),1,wk(1,i1),1)
  120 continue
      do  122  i1 = 1, 2
      do  122  j1 = 1, 2
      if (ccor) then
        do  123  j = 1, ldim
        do  123  i = 1, ldim
          hk(i,i1,j,j1)    = sk(i,i1,j,j1)*(wk(i,i1)+wk(j,j1)) +
     .                       hk(i,i1,j,j1)*vmtz
          hk(l2+i,i1,j,j1) = sk(l2+i,i1,j,j1)*(wk(i,i1)+wk(j,j1)) +
     .                       hk(l2+i,i1,j,j1)*vmtz
  123   continue
      else
        do  124  j = 1, ldim
        do  124  i = 1, ldim
          hk(i,i1,j,j1)    = sk(i,i1,j,j1)*(wk(i,i1)+wk(j,j1))
          hk(l2+i,i1,j,j1) = sk(l2+i,i1,j,j1)*(wk(i,i1)+wk(j,j1))
  124   continue
      endif
  122 continue

C --- O += (d*S*d) p^alpha (d*S*d) + <k|k>_bilinear ---
      if (.not. lasa) goto 32
      do  30  i1 = 1, 2
      do  30  i = 1, ldim
        wk(i,i1) = pph(4,i,i1)
        if (ccor) wk(i,i1) = wk(i,i1) + ccd(i,2)/pph(3,i,i1)**2
   30 continue
C     call prmx('hmltnc diag. matrix for 3C O',wk,ldim,ldim,2)
      call yyhmul(ldim*2,ldim*2,ldim*2,ldim*2,sk,l2,wk,ok,l2)
   32 continue

C --- H += d*S*d*(ep+o)*d*S*d + vmtz*<k|k>_bilinear ---
      do  130  i1 = 1, 2
      do  132  i = 1, ldim
        wk(i,i1) = (pph(1,i,i1)*pph(4,i,i1) + pph(5,i,i1))*nasa
        if (ccor) wk(i,i1) = wk(i,i1) + ccd(i,2)*vmtz/pph(3,i,i1)**2
  132 continue
      if (lso) call daxpy(ldim,1d0,sod(1,3,i1),1,wk(1,i1),1)
  130 continue
      call yyhmul(ldim*2,ldim*2,ldim*2,ldim*2,sk,l2,wk,hk,l2)

C --- O += 1 + 2(c-e)o + p(c-e)^2 + <k|k>_constant ---
      if (.not. lasa) goto 14
      do  10  i1 = 1, 2
      do  12  i = 1, ldim
        wk(i,i1) = 1d0 + 2d0*(pph(2,i,i1) - pph(1,i,i1))*pph(5,i,i1)
     .        + pph(4,i,i1)*(pph(2,i,i1) - pph(1,i,i1))**2
   12 continue
      do  10  i = 1, ldim
        ok(i,i1,i,i1) = ok(i,i1,i,i1) + wk(i,i1)
        if (ccor) ok(i,i1,i,i1) = ok(i,i1,i,i1) +ccd(i,0)*pph(3,i,i1)**2
   10 continue
   14 continue

C --- H += C + (C-e)[2oe+(C-e)(o+pe)] + vmtz*<k|k>_constant ---
      do  110  i1 = 1, 2
      do  112  i = 1, ldim
        wk(i,i1) = (pph(2,i,i1)  + (pph(2,i,i1)-pph(1,i,i1))*
     .    (2d0*pph(5,i,i1)*pph(1,i,i1) + (pph(2,i,i1)-pph(1,i,i1))*
     .    (pph(5,i,i1) + pph(4,i,i1)*pph(1,i,i1))))*nasa
        if (ccor) wk(i,i1) = wk(i,i1) + vmtz*ccd(i,0)*pph(3,i,i1)**2
  112 continue
      if (lso) call daxpy(ldim,1d0,sod(1,1,i1),1,wk(1,i1),1)
      do  110  i = 1, ldim
  110 hk(i,i1,i,i1) = hk(i,i1,i,i1) + wk(i,i1)

C ------ (1-2) block of S-O hamiltonian ------
      if (.not. lso) goto 99

C --- ho<phi H phi>oh + ho<phi H dot>h + h<dot H phi>oh + h<dot H dot>h
      call yhmlso(ldim,sk,sod(1,4,3),hk)

C --- H += <phi_lm L-S+ phi_lm+1> oh + <phi_lm L-S+ dot_lm+1> h ---
C This loop multipies matrices with structure
C          ( s11 s12 )   ( 0 sod )   ( 0  s11*sod)
C          (         ) * (       ) = (           )
C          ( s21 s22 )   ( 0  0  )   ( 0  s21*sod)
C Hermitian conjugate of the s21*sod block added for L+.S- term
      do  42  i1 = 1, 2
        do  43  j = 1, ldim-1
        xx = sod(j,2,3)
        do  43  i = 1, ldim
          hk(i,i1,j+1,2)    = hk(i,i1,j+1,2)    + sk(i,i1,j,1)*xx
          hk(l2+i,i1,j+1,2) = hk(l2+i,i1,j+1,2) + sk(l2+i,i1,j,1)*xx
          hk(j+1,2,i,i1)    = hk(j+1,2,i,i1)    + sk(i,i1,j,1)*xx
          hk(l2+j+1,2,i,i1) = hk(l2+j+1,2,i,i1) - sk(l2+i,i1,j,1)*xx
   43 continue
   42 continue

C --- H +=  ho <phi_lm-1 L-S+ phi_lm> + h <dot_lm-1 L-S+ phi_lm> ---
C Loop 42  multipies matrices that have the form in spin space
C         ( 0 sod )   ( s11 s12 )   ( s21*sod  s22*sod)
C         (       ) * (         ) = (                 )
C         ( 0  0  )   ( s21 s22 )   ( 0        0      )
C Hermitian conjugate of the s21*sod block added for L+.S- term
      do  44  j1 = 1, 2
      do  45  i = 2, ldim
      xx = sod(i,3,3)
      do  45  j = 1, ldim
        hk(i-1,1,j,j1)    = hk(i-1,1,j,j1)    + xx*sk(i,2,j,j1)
        hk(l2+i-1,1,j,j1) = hk(l2+i-1,1,j,j1) + xx*sk(l2+i,2,j,j1)
        hk(j,j1,i-1,1)    = hk(j,j1,i-1,1)    + xx*sk(i,2,j,j1)
        hk(l2+j,j1,i-1,1) = hk(l2+j,j1,i-1,1) - xx*sk(l2+i,2,j,j1)
   45 continue
   44 continue

C     call yprm('s.sod.3',2,hk,(ldim*2)**2,ldim*2,ldim*2)

C --- H += F/2 <phi_lm | L-S+ | phi_lm+1> ---
      do  46  i = 1, ldim-1
      hk(i+1,2,i,1) = hk(i+1,2,i,1) + sod(i,1,3)
   46 hk(i,1,i+1,2) = hk(i,1,i+1,2) + sod(i,1,3)

C --- Cleanup ---
   99 continue
      if (iprint() .gt. 110) then
        call yprm('h in hmltnc',12,hk,(ldim*2)**2,ldim*2,ldim*2,ldim*2)
        call yprm('o in hmltnc',12,ok,(ldim*2)**2,ldim*2,ldim*2,ldim*2)
      endif

      call tcx('hmltnc')
      end
C      subroutine mkdsdc(ladd,lzero,nbas,nl,ipc,lmx,indxsh,theta,
C     .  eula,neul,u,pph,ldim,sll,sk)
CC- srdel*s*srdel for combined spin-spiral and noncollinear hamiltonian
C      implicit none
C      logical ladd,lzero
C      integer ldim,nbas,neul,nl,lmx(*),ipc(nbas),indxsh(*)
C      double complex u(2,2,nl*nl,nbas)
C      double precision pph(5,ldim,2),sll(ldim,ldim,2,2),
C     .  sk(ldim,2,ldim,2,2),eula(nbas,neul,3),theta
C      double precision xxc,xxs,alpha,gamma,s,cc,ss
C      double complex uij(2,2),s1,s2,s11,s12,s22
C      integer lmi,ic,ib,i,il,im,lmj,jc,jb,j,jl,jm,i1,j1,ilm,jlm
C
CC      call prmx('sll in mkdsdc',sll,ldim,ldim,ldim)
C
CC --- Rotation matrix for SS ---
C      cc = dcos(theta/2)**2
C      ss = dsin(theta/2)**2
C      s  = dsin(theta)/2
C      if (lzero) call dpzero(sk,ldim*2*ldim*2*2)
C
CC --- Spinor rotation matrices for all sites ---
C      do  10  ib = 1, nbas
CC ... Assume initially euler angles are not l-dependent
C      xxc = dcos(eula(ib,1,2)/2)
C      xxs = dsin(eula(ib,1,2)/2)
C      alpha = eula(ib,1,1)
C      gamma = eula(ib,1,3)
C      ilm = 0
C      do  10  il = 1, nl
CC ... If euler angles are l dependent
C      if (neul .eq. nl) then
C        xxc = dcos(eula(ib,il,2)/2)
C        xxs = dsin(eula(ib,il,2)/2)
C        alpha = eula(ib,il,1)
C        gamma = eula(ib,il,3)
C      endif
C      do  10  im = -il+1, il-1
C      ilm = ilm+1
CC ... If euler angles are lm dependent
C      if (neul .eq. nl*nl) then
C        xxc = dcos(eula(ib,ilm,2)/2)
C        xxs = dsin(eula(ib,ilm,2)/2)
C        alpha = eula(ib,ilm,1)
C        gamma = eula(ib,ilm,3)
C      endif
CC     u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
CC     u(1,2,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
CC     u(2,1,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
CC     u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(-alpha-gamma)/2))
C      u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
C      u(2,1,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
C      u(1,2,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
C      u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
CC     call zprm('u in mkdsdc',2,u(1,1,ilm,ib),2,2,2)
C   10 continue
C
CC --- srdel*S*srdel ---
C      lmi = 0
C      do  20  ib = 1, nbas
C        ic = ipc(ib)
C        ilm = 0
C        do  25  il = 0, nl-1
C        do  25  im = -il, il
C        lmi = lmi+1
C        ilm = ilm+1
C        if (il .gt. lmx(ic)) goto 25
C        i = indxsh(lmi)
C        lmj = 0
C        do  30  jb = 1, nbas
C          jc = ipc(jb)
C          jlm = 0
C          do  35  jl = 0, nl-1
C          do  35  jm = -jl, jl
C            jlm = jlm+1
C            lmj = lmj+1
C            if (jl .gt. lmx(jc)) goto 35
C            j = indxsh(lmj)
C            s1 = dcmplx(sll(i,j,1,1),sll(i,j,2,1))
C            s2 = dcmplx(sll(i,j,1,2),sll(i,j,2,2))
C            s12 = s*s2-s*s1
C            s11 = ss*s2+cc*s1
C            s22 = ss*s1+cc*s2
C            uij(1,1) =
C     .      dconjg(u(2,1,ilm,ib))*(s22*u(2,1,jlm,jb)+s12*u(1,1,jlm,jb))+
C     .      dconjg(u(1,1,ilm,ib))*(s12*u(2,1,jlm,jb)+s11*u(1,1,jlm,jb))
C            uij(1,2) =
C     .      dconjg(u(2,1,ilm,ib))*(s22*u(2,2,jlm,jb)+s12*u(1,2,jlm,jb))+
C     .      dconjg(u(1,1,ilm,ib))*(s12*u(2,2,jlm,jb)+s11*u(1,2,jlm,jb))
C            uij(2,1) =
C     .      dconjg(u(2,2,ilm,ib))*(s22*u(2,1,jlm,jb)+s12*u(1,1,jlm,jb))+
C     .      dconjg(u(1,2,ilm,ib))*(s12*u(2,1,jlm,jb)+s11*u(1,1,jlm,jb))
C            uij(2,2) =
C     .      dconjg(u(2,2,ilm,ib))*(s22*u(2,2,jlm,jb)+s12*u(1,2,jlm,jb))+
C     .      dconjg(u(1,2,ilm,ib))*(s12*u(2,2,jlm,jb)+s11*u(1,2,jlm,jb))
C            do  36  i1 = 1, 2
C            do  36  j1 = 1, 2
C              sk(i,i1,j,j1,1) = sk(i,i1,j,j1,1) +
C     .          dble(uij(i1,j1))*pph(3,i,i1)*pph(3,j,j1)
C              sk(i,i1,j,j1,2) = sk(i,i1,j,j1,2) +
C     .          dimag(uij(i1,j1))*pph(3,i,i1)*pph(3,j,j1)
C   36       continue
C   35     continue
C   30   continue
C   25 continue
C   20 continue
C
CC     call prmx('sk in mkdsdc',sk,ldim*2,ldim*2,ldim*2)
C
C      if (ladd) then
C        do  40  i1 = 1, 2
C        do  40  i = 1, ldim
C   40   sk(i,i1,i,i1,1) = sk(i,i1,i,i1,1) + pph(2,i,i1)-pph(1,i,i1)
C      endif
C
C      end
C      subroutine mkdsds(ladd,lzero,theta,pph,ldim,sll,sk)
CC- srdel*s*srdel for spin spiral hamiltonian
C      implicit none
C      logical ladd,lzero
C      integer ldim
C      double precision theta,pph(5,ldim,2),sll(ldim,ldim,2,2),
C     .  sk(ldim,2,ldim,2,2)
C      double precision s,cc,ss,u1(2,2),u2(2,2)
C      integer i,j,i1,j1
C
CC --- Rotation matrix for SS ---
C      cc = dcos(theta/2)**2
C      ss = dsin(theta/2)**2
C      s  = dsin(theta)/2
C      if (lzero) call dpzero(sk,ldim*2*ldim*2*2)
C      u1(1,1) = cc
C      u1(1,2) = -s
C      u1(2,1) = -s
C      u1(2,2) = ss
C      u2(1,1) = ss
C      u2(1,2) =  s
C      u2(2,1) =  s
C      u2(2,2) = cc
C
CC --- srdel*S*srdel ---
C      do  25  i1 = 1, 2
C      do  25  j1 = 1, 2
C      do  25  j = 1, ldim
C      do  25  i = 1, ldim
C        sk(i,i1,j,j1,1) = sk(i,i1,j,j1,1) + pph(3,i,i1)*pph(3,j,j1)*
C     .    (sll(i,j,1,1)*u1(i1,j1) + sll(i,j,1,2)*u2(i1,j1))
C        sk(i,i1,j,j1,2) = sk(i,i1,j,j1,2) + pph(3,i,i1)*pph(3,j,j1)*
C     .    (sll(i,j,2,1)*u1(i1,j1) + sll(i,j,2,2)*u2(i1,j1))
C   26   continue
C   25 continue
C
CC     call yprm('sk in mkdsds',02,sk,(ldim*2)**2,ldim*2,ldim*2,ldim*2)
C
CC --- Add C-enu ---
C      if (ladd) then
C        do  40  i1 = 1, 2
C        do  40  i = 1, ldim
C   40   sk(i,i1,i,i1,1) = sk(i,i1,i,i1,1) + pph(2,i,i1)-pph(1,i,i1)
C      endif
C
C      end
C      subroutine mkdsdn(ladd,lzero,nbas,nl,ipc,lmx,indxsh,eula,neul,u,
C     .  pph,ldim,sll,sk)
CC- srdel*s*srdel for noncollinear hamiltonian
C      implicit none
C      logical ladd,lzero
C      integer ldim,nbas,neul,nl,lmx(*),ipc(nbas),indxsh(*)
C      double complex u(2,2,nl*nl,nbas)
C      double precision pph(5,ldim,2),sll(ldim,ldim,2),
C     .  sk(ldim,2,ldim,2,2),eula(nbas,neul,3)
C      double precision xxc,xxs,alpha,gamma
C      double complex xx,uij
C      integer lmi,ic,ib,i,il,im,lmj,jc,jb,j,jl,jm,i1,j1,ilm,jlm
C
CC     call prmx('sll in mkdsdn',sll,ldim,ldim,ldim)
C      if (lzero) call dpzero(sk,ldim*2*ldim*2*2)
C
CC --- Spinor rotation matrices for all sites ---
C      do  10  ib = 1, nbas
CC ... Assume initially euler angles are not l-dependent
C      xxc = dcos(eula(ib,1,2)/2)
C      xxs = dsin(eula(ib,1,2)/2)
C      alpha = eula(ib,1,1)
C      gamma = eula(ib,1,3)
C      ilm = 0
C      do  10  il = 1, nl
CC   ... If euler angles are l dependent
C        if (neul .eq. nl) then
C          xxc = dcos(eula(ib,il,2)/2)
C          xxs = dsin(eula(ib,il,2)/2)
C          alpha = eula(ib,il,1)
C          gamma = eula(ib,il,3)
C        endif
C        do  10  im = -il+1, il-1
C        ilm = ilm+1
CC   ... If euler angles are lm dependent
C        if (neul .eq. nl*nl) then
C          xxc = dcos(eula(ib,ilm,2)/2)
C          xxs = dsin(eula(ib,ilm,2)/2)
C          alpha = eula(ib,ilm,1)
C          gamma = eula(ib,ilm,3)
C        endif
CC       u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
CC       u(1,2,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
CC       u(2,1,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
CC       u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
C        u(1,1,ilm,ib) =  xxc*cdexp(dcmplx(0d0,-(alpha+gamma)/2))
C        u(2,1,ilm,ib) =  xxs*cdexp(dcmplx(0d0,(alpha-gamma)/2))
C        u(1,2,ilm,ib) = -xxs*cdexp(dcmplx(0d0,(-alpha+gamma)/2))
C        u(2,2,ilm,ib) =  xxc*cdexp(dcmplx(0d0,(alpha+gamma)/2))
CC       call zprm('u in mkdsdn',2,u(1,1,ilm,ib),2,2,2)
C   10 continue
C
CC --- srdel*S*srdel ---
C      lmi = 0
C      do  20  ib = 1, nbas
C        ic = ipc(ib)
C        ilm = 0
C        do  25  il = 0, nl-1
C        do  25  im = -il, il
C        lmi = lmi+1
C        ilm = ilm+1
C        if (il .gt. lmx(ic)) goto 25
C        i = indxsh(lmi)
C        lmj = 0
C        do  30  jb = 1, nbas
C          jc = ipc(jb)
C          jlm = 0
C          do  35  jl = 0, nl-1
C          do  36  jm = -jl, jl
C            lmj = lmj+1
C            jlm = jlm+1
C            if (jl .gt. lmx(jc)) goto 36
C            j = indxsh(lmj)
C            do  32  i1 = 1, 2
C            do  32  j1 = 1, 2
C              uij = dconjg(u(1,i1,ilm,ib))*u(1,j1,jlm,jb) +
C     .              dconjg(u(2,i1,ilm,ib))*u(2,j1,jlm,jb)
C              xx  = dcmplx(sll(i,j,1),sll(i,j,2))*uij*
C     .              pph(3,i,i1)*pph(3,j,j1)
C              sk(i,i1,j,j1,1) = sk(i,i1,j,j1,1) + dble(xx)
C              sk(i,i1,j,j1,2) = sk(i,i1,j,j1,2) + dimag(xx)
C   32       continue
C   37       continue
C   36     continue
C   35     continue
C   30   continue
C   25 continue
C   20 continue
C
CC     call yprm('sk in mkdsdn',02,sk,(ldim*2)**2,ldim*2,ldim*2,ldim*2)
C
C      if (ladd) then
C        do  40  i1 = 1, 2
C        do  40  i = 1, ldim
C   40   sk(i,i1,i,i1,1) = sk(i,i1,i,i1,1) + pph(2,i,i1)-pph(1,i,i1)
C      endif
C
C      end
C      subroutine mkdsdo(ladd,lzero,nbas,nl,ipc,lmx,indxsh,
C     .  pph,ldim,sll,sk)
CC- srdel*s*srdel for S-O hamiltonian (diagonal in +,-)
C      implicit none
C      logical ladd,lzero
C      integer ldim,nbas,nl,lmx(*),ipc(nbas),indxsh(*)
C      double precision pph(5,ldim,2),sll(ldim,ldim,2),
C     .  sk(ldim,2,ldim,2,2)
C      double complex xx
C      integer lmi,ic,ib,i,il,im,lmj,jc,jb,j,jl,jm,i1,j1
C
CC     call prmx('sll in mkdsdn',sll,ldim,ldim,ldim)
CC     call yprm('sll in mkdsdn',02,sll,ldim*ldim,ldim,ldim,ldim)
C      if (lzero) call dpzero(sk,ldim*2*ldim*2*2)
C
CC --- srdel*S*srdel ---
C      lmi = 0
C      do  20  ib = 1, nbas
C        ic = ipc(ib)
C        do  25  il = 0, nl-1
C        do  25  im = -il, il
C        lmi = lmi+1
C        if (il .gt. lmx(ic)) goto 25
C        i = indxsh(lmi)
C        lmj = 0
C        do  30  jb = 1, nbas
C          jc = ipc(jb)
C          do  35  jl = 0, nl-1
C          do  36  jm = -jl, jl
C            lmj = lmj+1
C            if (jl .gt. lmx(jc)) goto 36
C            j = indxsh(lmj)
C            do  37  i1 = 1, 2
C            j1 = i1
CC           do  37  j1 = 1, 2
C              xx = dcmplx(sll(i,j,1),sll(i,j,2))*pph(3,i,i1)*pph(3,j,j1)
C              sk(i,i1,j,j1,1) = sk(i,i1,j,j1,1) + dble(xx)
C              sk(i,i1,j,j1,2) = sk(i,i1,j,j1,2) + dimag(xx)
C   37       continue
C   36     continue
C   35     continue
C   30   continue
C   25 continue
C   20 continue
C
Cc     call prmx('sk in mkdsdo',sk,ldim*2,ldim*2,ldim*2)
CC     call yprm('sk in mkdsdo',12,sk,(ldim*2)**2,ldim*2,ldim*2,ldim*2)
C
C      if (ladd) then
C        do  40  i1 = 1, 2
C        do  40  i = 1, ldim
C   40   sk(i,i1,i,i1,1) = sk(i,i1,i,i1,1) + pph(2,i,i1)-pph(1,i,i1)
C      endif
C
C      end
C

