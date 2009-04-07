      subroutine trS(ldim,idim,alpha,pti,sll,sil,sii,a,wk)
C- Transform the structure constant matrix for downfolding
C-----------------------------------------------------------------------
Ci Inputs
Ci    ldim :dimension of lower block of hamiltonian
Ci    idim :dimension of intermediate block of hamiltonian
Ci   alpha, pti, screening constants and inverse potential function
Ci               vectors in downfolding order
Ci   wk (work array length nl**2*nbas)
Ci   sll, sil, sii: respectively l-l, i-l and i-i blocks of the
Ci   Bloch transformed structures constants
Co Output
Co   s (see remarks), alpha - the lower-wave screening constants
Co   are unchanged; on exit, the i-wave alpha'a are returned as
Co   (beta - alpha)_i^-1; a: the matrix
Co   A_il = [(beta - alpha)_i^-1 - S^alpha_ii]^-1 S_il^alpha
Cr Remarks
Cr   s is Bloch transformed and of dimension lidim X lidim. The upper
Cr   left block is the l-l block of lower waves, lower right is the
Cr   i-i block of intermediate (down-folded) waves. On exit, the l-l
Cr   block and the i-l block are s in the new representation, beta, in
Cr   which the i-i can eventually be discarded. The i-i block is
Cr   overwritten
Cu Updates
Cu    6 Jun 00 revised to distingish int from higher blocks
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ldim,idim
      double precision alpha(*),pti(*),
     .                 wk(2*idim),sll(0:ldim-1,0:2*ldim-1),
     .                 sil(0:idim-1,0:ldim-1,2),
     .                 sii(0:idim-1,0:2*idim-1),a(0:idim-1,0:ldim-1,2)
C Local parameters
      integer ierr,i2,li,i,j,oiwk,ldw
      double precision alp

C heap
      integer w(1)
      common  / w / w

C External procedures
      external invalp,dscal,zmpy,yyhifa,yyhidi

      call tcn('trS')

      i2 = idim**2
      li = ldim * idim

C --- xsi := beta - alpha for intermediate waves ---
      call daxpy(idim,-1d0,pti(1+ldim),1,alpha(1+ldim),1)
      call dscal(idim,-1d0,alpha(1+ldim),1)
      call invalp(idim,1,alpha(1+ldim))

C --- Invert i-i block of -(xsi^-1 - S) ---
      call daxpy(idim,-1d0,alpha(1+ldim),1,sii,idim+1)
      call defask(i)
#ifndef YQINV
      i = 0
#endif
C     call tcn('Make Ail')
      if (i .gt. idim*(idim+1)) then
        ldw = max(ldim,idim)
        call defrr(oiwk,ldw*(idim+1))
C       call yqinv('h',sii,i2,idim,2,idim,w(oiwk),idim,ierr)
        call dpcopy(sil,a,1,2*li,-1d0)
C   ... A_il := (xsi^-1 - S)_ii^-1 S_il
        call yqinvb('h',sii,i2,idim,idim,ldim,w(oiwk),ldw,
     .    w(oiwk),a,idim*ldim,idim,ierr)
        call rxx(ierr.ne.0,'trS: (xsi^-1 - S) was singular')
        call rlse(oiwk)
      else
        call defi(oiwk,idim)
        call yyhifa(sii,sii(i2,0),idim,idim,w(oiwk),ierr)
        call rxx(ierr.ne.0,'trS: (xsi^-1 - S) was singular')
        call yyhidi(sii,sii(i2,0),idim,idim,w(oiwk),0,0,wk,wk(1+idim),1)
        call rlse(oiwk)
        do  1  j = 0, idim-1
        do  1  i = j, idim-1
          sii(i,j) = sii(j,i)
          sii(i2+i,j) = - sii(i2+j,i)
    1   continue
        call dscal(i2*2,-1d0,sii,1)
C       call yprm('-(xsi^1-S)^-1',2,sii,idim*idim,idim,idim,idim)
C   ... A_il := (xsi^-1 - S)_ii^-1 S_il
        call zmpy(sii,idim,1,i2,sil,idim,1,li,a,idim,1,li,idim,ldim,
     .    idim)
      endif
C     call yprm('A_il',2,a,idim*ldim,idim,idim,ldim)
C     call tcx('Make Ail')

C --- S_ll --> S_ll + S_il^dagger A_il ---
#ifndef CRAY
      call yyhmpy('c','n',ldim,idim,sil,a,.true.,sll)
#else
      call dscal(li,-1d0,sil(0,0,2),1)
      call zhmpy(sil,1,idim,li,a,idim,1,li,sll,ldim,.true.,idim)
#endif

C --- S_il --> xsi^-1 A_il ---
      do  3  i = 0, idim-1
        alp = alpha(1+ldim+i)
        do  2  j = 0, ldim-1
          sil(i,j,1) = alp*a(i,j,1)
          sil(i,j,2) = alp*a(i,j,2)
    2   continue
    3 continue

C     call yprm('S_il(beta)',2,sil,idim*ldim,idim,idim,ldim)
C     call yprm('S_ll(beta)',2,sll,ldim*ldim,ldim,ldim,ldim)

      call tcx('trS')
      end

