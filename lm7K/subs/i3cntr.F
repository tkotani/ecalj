      subroutine i3cntr(ldim,idim,lnc,sil,pph,vmtz,h,wk)
C- make 3-centre integrals in H over the i-waves
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   lnc   :T, suppress scaling sqrdel
Ci   sil   :lower-intermediate block, rotated to 'beta' repsn
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci   vmtz  :muffin-tin zero (asamad.f)
Ci    wk   :work array
Co  Output:
Co    h    :i-wave contribution to hamiltonian is added; see remarks
Co    sil:  sil is OVERWRITTEN
Cr  Remarks:
Cr    These are the Omega^dagger*e_nu*p^alpha*Omega terms over the
Cr    i-waves in the hamiltonian (eq. 2.40, Kanpur notes).
Cr    Beta^dot in the i-set is chosen so that corresponding terms in
Cr    the overlap matrix vanish. In this case the i-wave bilinear term
Cr    in H becomes (sqrdel*S)_li * D_i * (S*sqrdel)_il
Cr    where the diagonal matrix D is (e_nu - v_mtz) DELTA / (C - e_nu)^2
Cr    which is representation - independent.
Cr
Cr    In the noncollinear case, the decoration by sqrdel ( ) sqrdel
Cr    is suppressed, pending rotation by Euler angles.
Cu Updates
Cu   29 Jan 00 added lnc switch
C-----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical lnc
      integer ldim,idim
      double precision sil(idim,ldim,2),h(ldim,ldim,2),
     .                 pph(5,ldim),vmtz,wk(idim)
C ... Local parameters
      integer i,j,iprint
C ... External calls
      external tcn,tcx,yprm,yyhmul

      if (idim .eq. 0) return

      call tcn('i-waves H')

      do  1  i = 1, idim
        j = i + ldim
        wk(i) = ( (pph(1,j) - vmtz)*(pph(3,j)/(pph(2,j)-pph(1,j)))**2 )
        if (.not. lnc) then
        do  2  j = 1, ldim
          sil(i,j,1) = sil(i,j,1)*pph(3,j)
          sil(i,j,2) = sil(i,j,2)*pph(3,j)
    2   continue
        endif
    1 continue
C     call yprm('wk',1,wk,ldim*ldim,ldim,ldim,1)

      call yyhmul(idim,idim,ldim,ldim,sil,idim*ldim,wk,h,ldim**2)
      if (iprint() .ge. 110 .and. .not. lnc)
     .  call yprm('iwaves H',12,h,ldim*ldim,ldim,ldim,ldim)

      call tcx('i-waves H')
      end

