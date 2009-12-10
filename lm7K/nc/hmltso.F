      subroutine hmltso(lasa,ccor,lso,ldim,idim,pph,
     .  sod,diags,ndim,vmtz,sk,hk,ok,wk)
C- ASA + spin-orbit hamiltonian, and overlap
C ---------------------------------------------------------------------
Ci Inputs
Ci   lasa,ccor,lso: switches telling whether to add (1) asa H,
Ci                  (2) combined correction and (3) spin-orbit H
Ci   ldim,idim, dimension of lower and intermediate waves
Ci   pph: vector of pot pars, make by makpph
Ci   sod: diagonal arrays for SO coupling, made by mksod:
Co   sod(*,1..3,isp): arrays for 1-, 2-, 3- center terms, ++,-- blocks
Co   sod(*,1,3): for 1- center terms,    L- S+ block
Co   sod(*,2..3,3): for 2- center terms, L- S+ block
Co   sod(*,4,3):    for 3- center terms, L- S+ block
Cr Remarks
Cr   (H_so)_(i,j) = <psi_i| I L . S | psi_j>.
Cr   where L.S = Lz.Sz + .5(L+.S- + L-.S+) and I = (2/(c^2) dV/dr*(1/r).
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | psi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   Subroutine mksod combines radial matrix elements I between phi's
Cr   and dot's (phidot^gamma) and the L.S operator into diagonal
Cr   matrices associated with the 1-, 2- and 3-center terms comprising
Cr   matrix elements (H_so)_(i,j).
Cr
Cr   The up-up and down-down (Lz.Sz) terms calculated first
Cr   the up-down (L-S+) calculated next; (L+S-) obtained from (L-S+).
Cr
Cr   The ASA hamiltonian works directly with D = srdel S srdel, rather
Cr   that h = C-e+D.  To combine SO with the ASA H, the diagonal parts
Cr   of the one-, two-, three-center terms are reshuffled.  Thus
Cr   w1 + w2 (C-e+D) + h.c. + (C-e+D) w3 (C-e+D)  becomes
Cr   w1+2*w2*(C-e)+w3(C-e)**2 + (w2+(C-e)w3) D + h.c. + D w3 D.
C ---------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lasa,ccor,lso
      integer ldim,idim,ndim
      double precision vmtz,
     .  sk(ldim,2,ldim,2),pph(5,ldim,2),diags(ndim,0:2),
     .  ok(ldim,2,ldim,2),hk(ldim,2,ldim,2),sod(ldim,4,3),wk(ldim,2)
C Local parameters
      integer i,j,i1,j1,l2,iprint,nasa

      l2 = ldim**2*4
      nasa = 0
      if (lasa) nasa = 1
C     call yprm('sk on entry',12,sk,ldim*2*ldim*2,ldim*2,ldim*2,ldim*2)

C --- O += d*S*d*[o+p(c-e)] + h.c. + <k|k>_linear ---
      if (.not. lasa) goto 26
      do  20  i1 = 1, 2
      do  20  i = 1, ldim
   20 wk(i,i1) = pph(5,i,i1) + pph(4,i,i1)*(pph(2,i,i1)-pph(1,i,i1))
      if (ccor) call daxpy(ldim,1d0,diags(1,1),1,wk,1)
      if (ccor) call daxpy(ldim,1d0,diags(1,1),1,wk(1,2),1)
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
C       print *, 'zap so(2)'
        if (lso)  call daxpy(ldim,1d0,sod(1,2,i1),1,wk(1,i1),1)
        if (ccor) call daxpy(ldim,vmtz,diags(1,1),1,wk(1,i1),1)
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

C --- O += d*S*d*p^alpha*d*S*d + <k|k>_bilinear ---
      if (.not. lasa) goto 32
      do  30  i1 = 1, 2
      do  30  i = 1, ldim
        wk(i,i1) = pph(4,i,i1)
        if (ccor) wk(i,i1) = wk(i,i1) + diags(i,2)/pph(3,i,i1)**2
   30 continue
      call yyhmul(ldim*2,ldim*2,ldim*2,ldim*2,sk,l2,wk,ok,l2)
   32 continue

C --- H += d*S*d*(ep+o)*d*S*d + vmtz*<k|k>_bilinear ---
      do  130  i1 = 1, 2
      do  132  i = 1, ldim
        wk(i,i1) = (pph(1,i,i1)*pph(4,i,i1) + pph(5,i,i1))*nasa
        if (ccor) wk(i,i1) = wk(i,i1) + diags(i,2)*vmtz/pph(3,i,i1)**2
  132 continue
C     print *, 'zap so(3)'
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
        if (ccor) ok(i,i1,i,i1) = ok(i,i1,i,i1) +
     .    diags(i,0)*pph(3,i,i1)**2
   10 continue
   14 continue

C      print *, 'zap hk'
C      call dpzero(hk,l2*2)

C --- H += C + (C-e)[2oe+(C-e)(o+pe)] + vmtz*<k|k>_constant ---
      do  110  i1 = 1, 2
      do  112  i = 1, ldim
        wk(i,i1) = (pph(2,i,i1)  + (pph(2,i,i1)-pph(1,i,i1))*
     .    (2d0*pph(5,i,i1)*pph(1,i,i1) + (pph(2,i,i1)-pph(1,i,i1))*
     .    (pph(5,i,i1) + pph(4,i,i1)*pph(1,i,i1))))*nasa
        if (ccor) wk(i,i1) = wk(i,i1) + vmtz*diags(i,0)*pph(3,i,i1)**2
  112 continue
c      print *, 'zap so(1)'
C      print *, 'zap wk(1) to get soham only'
C      call dpzero(wk(1,i1),ldim)
      if (lso) call daxpy(ldim,1d0,sod(1,1,i1),1,wk(1,i1),1)
      do  110  i = 1, ldim
  110 hk(i,i1,i,i1) = hk(i,i1,i,i1) + wk(i,i1)

      if (.not. lso) goto 99
C      print *, 'goto 99'
C      goto 99

C ------ Make up-down, down-up blocks of S-O hamiltonian ------

C ---    H += h_IL+ <phi_I|H_so|phi_I> h_IL
C        do  22  k = ldim+1, ldim+idim-1
C          elkp1 = soph(4,1,1,k+1)
C          emkp1 = soph(4,2,1,k+1)
C          elmkp1 = dsqrt((elkp1+emkp1)*(elkp1-emkp1+1))
C
C          hk(i,1,j,2) = hk(i,1,j,2) +
C     .      hil(k-ldim,i,1,1)*soph(1,1,2,k)*sp*elmkp1*
C     .      hil(k+1-ldim,j,1,2) +
C     .      hil(k-ldim,i,2,1)*soph(1,1,2,k)*sp*elmkp1*
C     .      hil(k+1-ldim,j,2,2)
C
C          hk(l2+i,1,j,2) = hk(l2+i,1,j,2) +
C     .      hil(k-ldim,i,1,1)*soph(1,1,2,k)*sp*elmkp1*
C     .      hil(k+1-ldim,j,2,2) -
C     .      hil(k-ldim,i,2,1)*soph(1,1,2,k)*sp*elmkp1*
C     .      hil(k+1-ldim,j,1,2)
C
C   22   continue

C --- ho<phi H phi>oh + ho<phi H dot>h + h<dot H phi>oh + h<dot H dot>h
C     Make sk(i,k) sod(k,4,3). NB: this uses h(i,2,j,1) as a workspace;
C     also initial contents of h(i,1,j,2) are destroyed.
      do  41  j = 1, ldim
      do  41  i = 1, ldim
      hk(i,2,j,1)    = sk(i,1,j,1)*sod(j,4,3)
   41 hk(l2+i,2,j,1) = sk(l2+i,1,j,1)*sod(j,4,3)
      call yygemm('N','N',ldim,ldim,ldim-1,1d0,hk(1,2,1,1),
     .  hk(l2+1,2,1,1),2*ldim,sk(2,2,1,2),sk(l2+2,2,1,2),2*ldim,
     .  0d0,hk(1,1,1,2),hk(l2+1,1,1,2),2*ldim)

C --- H += <phi_lm L-S+ phi_lm+1> oh + <phi_lm L-S+ dot_lm+1> h ---
      do  43  j = 1, ldim-1
      do  43  i = 1, ldim
        hk(i,1,j+1,2)    = hk(i,1,j+1,2)    + sk(i,1,j,1)*sod(j,2,3)
        hk(l2+i,1,j+1,2) = hk(l2+i,1,j+1,2) + sk(l2+i,1,j,1)*sod(j,2,3)
   43 continue

C --- H +=  ho <phi_lm-1 L-S+ phi_lm> + h <dot_lm-1 L-S+ phi_lm> ---
      do  45  j = 1, ldim
      do  45  i = 2, ldim
        hk(i-1,1,j,2)    = hk(i-1,1,j,2)    + sod(i,3,3)*sk(i,2,j,2)
        hk(l2+i-1,1,j,2) = hk(l2+i-1,1,j,2) + sod(i,3,3)*sk(l2+i,2,j,2)
   45 continue

C --- H += F/2 <phi_lm | L-S+ | phi_lm+1> ---
      do  46  i = 1, ldim-1
   46 hk(i,1,i+1,2) = hk(i,1,i+1,2) + sod(i,1,3)

C --- Get L+ . S-  from L- . S+ ---
      do  50  j = 1, ldim
      do  50  i = 1, ldim
        hk(i,2,j,1)    =  hk(j,1,i,2)
        hk(l2+i,2,j,1) = -hk(l2+j,1,i,2)
   50 continue

C     call yprm('hk at end',2,hk,ldim*2*ldim*2,ldim*2,ldim*2,ldim*2)

C --- Cleanup ---
   99 continue
      if (iprint() .gt. 110) then
      call yprm('hk at end of hmltso',12,hk,
     .    (ldim*2)**2,ldim*2,ldim*2,ldim*2)
      call yprm('ok at end of hmltso',12,ok,
     .  (ldim*2)**2,ldim*2,ldim*2,ldim*2)
      endif
      end


