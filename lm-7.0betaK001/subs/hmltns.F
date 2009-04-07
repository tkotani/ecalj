      subroutine hmltns(mode,ccor,ccd,vmtz,ldim,lihdim,pph,sk,hk,ok,wk)
C- Generate ASA and CCOR part of Hamiltonian and Overlap
C ---------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 hamiltonian-> hk and overlap -> ok (sk overwritten)
Ci         :1 small h -> hk
Ci         :2 1+oh -> hk
Ci         :3 small h -> hk, input sk is srdel * S * srdel
Ci         :4 1+oh -> hk, input sk is srdel * S * srdel
Ci   ccor  :switch telling whether to add combined correction
Ci   ccd:  :diagonal matrices for 1-, 2-, 3-center terms in the strux S;
Ci         :they are terms in parentheses in Kanpur notes
Ci         :eq. 3.87 multiplied by w^2. (see makdia.f)
Ci   vmtz  :muffin-tin zero (asamad.f)
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lihdim:dimensions ccd and pph
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :pph(1..5,i,is): parms for ith RL and is(th) spin channel.
Ci         :pph(1) : enu
Ci         :pph(2) : calpha
Ci         :pph(3) : sqrdel
Ci         :pph(4) : palpha
Ci         :pph(5) : oalp
Ci   sk    :structure constants, s^beta
Ci   ok    :S^beta-dot
Ci   hk    :i-wave 3-centre integrals (see i3cntre)
Ci   wk    :work array of length ldim
Co Outputs
Co   hk,ok,:(mode 0) hamiltonian-> hk and overlap -> ok
Co   sk    :         Also sk is changed to sqrdel*sk*sqrdel
Co         :(mode 1) small h -> hk
Co         :(mode 2) 1+oh -> hk
Cr Remarks
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
Cr   are calculated from <phi|H|phi>, <phidot|H|phi>, <phidot|H|phidot>
Cr   and are passed in array ccd.
Cr
Cr   Also, In the notation below, <k|k> = <kappa|kappa>
Cu Updates
Cu   15 Nov 07 New mode, to generate h or 1+oh
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer mode,lihdim,ldim
      logical ccor
      double precision pph(5,lihdim),wk(ldim),ccd(lihdim,0:2),vmtz
      double precision sk(ldim,ldim*2),hk(ldim,ldim*2),ok(ldim,ldim*2)
C Local parameters
      integer i,j,l2
C External calls
      external yyhmul,daxpy

      call tcn('hmltns')

      l2 = ldim**2
      call isanrg(mode,0,4,'hmltns:','mode',.true.)

C --- mode 1-4: make h or 1+oh ---
      if (mode .ge. 1 .and. mode .le. 4) then
        i = 1
        if (mode .eq. 3 .or. mode .eq. 4) i = 31
        call makdsd(i,ldim,ldim,ldim,ldim,0,0,pph,sk,hk)
C       call yprm('h',12,hk,ldim*ldim,ldim,ldim,ldim)
        if (mode .eq. 1 .or. mode .eq. 3) goto 999
        do   i = 1, ldim
          wk(i) = pph(5,i)
        enddo
        do  j = 1, ldim
          do  i = 1, ldim
            hk(i,j)    = wk(i)*hk(i,j)
            hk(l2+i,j) = wk(i)*hk(l2+i,j)
          enddo
          hk(j,j) = 1 + hk(j,j)
        enddo
C       call yprm('1+oh',2,hk,ldim*ldim,ldim,ldim,ldim)
        goto 999
      endif

C --- Make d*Sdot*d, Add vmtz*d*Sdot*d into H ---
      call makdsd(0,ldim,ldim,ldim,ldim,0,0,pph,ok,ok)
      call daxpy(2*l2,vmtz,ok,1,hk,1)

C --- Make d*S*d with d = sqrt(delta) ---
      call makdsd(0,ldim,ldim,ldim,ldim,0,0,pph,sk,sk)
C     call yprm('dsd',12,sk,ldim*ldim,ldim,ldim,ldim)

C --- O += d*S*d*p^alpha*d*S*d + <k|k>_bilinear ---
      do  1  i = 1, ldim
        wk(i) = pph(4,i)
        if (ccor) wk(i) = wk(i) + ccd(i,2)/pph(3,i)**2
    1 continue
      call yyhmul(ldim,ldim,ldim,ldim,sk,l2,wk,ok,l2)
C     call yprm('O += 3c',12,ok,ldim*ldim,ldim,ldim,ldim)

C --- O += d*S*d*[o+p(c-e)] + h.c. + <k|k>_linear ---
      do  2  i = 1, ldim
    2 wk(i) = pph(5,i) + pph(4,i)*(pph(2,i) - pph(1,i))
      if (ccor) call daxpy(ldim,1d0,ccd(1,1),1,wk,1)
#if SGI_PARALLEL
C$DOACROSS local(j,i)
#endif
      do  3  j = 1, ldim
      do  3  i = 1, ldim
        ok(i,j) = ok(i,j) + sk(i,j)*(wk(i) + wk(j))
        ok(l2+i,j) = ok(l2+i,j) + sk(l2+i,j)*(wk(i) + wk(j))
    3 continue
C     call yprm('O+= 2c',12,ok,ldim*ldim,ldim,ldim,ldim)

C --- O += 1 + 2(c-e)o + p(c-e)^2 + <k|k>_constant ---
      do  4  i = 1, ldim
        wk(i) = 1d0 + 2d0*(pph(2,i) - pph(1,i))*pph(5,i)
     .        + pph(4,i)*(pph(2,i) - pph(1,i))**2
    4 continue
      do  5  i = 1, ldim
        ok(i,i) = ok(i,i) + wk(i)
        if (ccor) ok(i,i) = ok(i,i) + ccd(i,0)*pph(3,i)**2
    5 continue
C     call yprm('O+= 1c',12,ok,ldim*ldim,ldim,ldim,ldim)

C --- H += d*S*d*(ep+o)*d*S*d + vmtz*<k|k>_bilinear ---
      do  6  i = 1, ldim
        wk(i) = pph(1,i)*pph(4,i) + pph(5,i)
        if (ccor) wk(i) = wk(i) + ccd(i,2)*vmtz/pph(3,i)**2
    6 continue
      call yyhmul(ldim,ldim,ldim,ldim,sk,l2,wk,hk,l2)
C     call yprm('H += 3c',12,hk,ldim*ldim,ldim,ldim,ldim)

C --- H += d*S*d*[1/2+oe+(c-e)(o+pe)] + h.c. + vmtz*<k|k>_linear ---
      do  7  i = 1, ldim
        wk(i) = .5d0 + pph(1,i)*pph(5,i) + (pph(2,i) - pph(1,i))
     .        * (pph(5,i) + pph(4,i)*pph(1,i))
    7 continue
      if (ccor) call daxpy(ldim,vmtz,ccd(1,1),1,wk,1)

#if SGI_PARALLEL
C$DOACROSS local(j,i)
#endif
      do  8  j = 1, ldim
      do  8  i = 1, ldim
        hk(i,j) = hk(i,j) + sk(i,j)*( wk(i) + wk(j) )
        hk(l2+i,j) = hk(l2+i,j) + sk(l2+i,j)*( wk(i) + wk(j) )
    8 continue
C     call yprm('H += 2c',12,hk,ldim*ldim,ldim,ldim,ldim)

C --- H += c + (c-e)[2oe+(c-e)(o+pe)] + vmtz*<k|k>_constant ---
      do  9  i = 1, ldim
        wk(i) = pph(2,i) + (pph(2,i) - pph(1,i))*(2d0*pph(5,i)*pph(1,i)
     .        + (pph(2,i) - pph(1,i))*(pph(5,i) + pph(4,i)*pph(1,i)))
    9 continue
      do  10  i = 1, ldim
        hk(i,i) = hk(i,i) + wk(i)
        if (ccor) hk(i,i) = hk(i,i) + vmtz*ccd(i,0)*pph(3,i)**2
   10 continue
C     call yprm('H += 1c',12,hk,ldim*ldim,ldim,ldim,ldim)

C ... Exit
  999 continue
      call tcx('hmltns')
      end

