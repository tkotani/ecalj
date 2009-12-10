      subroutine trSdot(ldim,idim,alpha,adot,a,sll,sil,sii)
C- transform s^dot for downfolding
C-----------------------------------------------------------------------
Ci  Inputs:
Ci    ldim :dimension of lower block of hamiltonian
Ci    idim :dimension of intermediate block of hamiltonian
Ci    alpha:(only i-wave entries used)  (beta-alpha)^-1
Ci         :in downfolding order 
Ci    adot :(only i-wave entries used) w^-2(beta^dot - alpha^dot)
Ci         :in downfolding order 
Ci    a    :A_il (see trS); sll,sil,sii (S^alpha^dot on entry)
Co  Outputs:
Co    sll  :input S^alpha^dot_ll transformed for to S^beta^dot_ll
Co    adot :(i-set only) is changed to w^-2(beta-alpha)^-1^dot
Co    A_il :is returned as A_li^dagger
Cu Updates
Cu    6 Jun 00 revised to distingish int from higher blocks
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ldim,idim
      double precision alpha(*),adot(*),
     .                 a(0:idim-1,0:ldim-1,2),sil(0:idim-1,0:ldim-1,2),
     .                 sll(0:ldim-1,0:ldim*2-1),sii(0:idim-1,0:idim*2-1)
C Local parameters
      integer i2,li

C External procedures
      external scalpv,dscal,daxpy,zmpy

      call tcn('trSdot')
          
      i2 = idim**2
      li = ldim * idim

C --- S^dot_ll --> S^dot_ll + A_li S^dot_il + S^dot_li A_il ---
C     call tcn('call yyhmpy')
#ifndef CRAY
      call yyhmpy('c','n',ldim,idim,a,sil,.false.,sll)
      call yyhmpy('c','n',ldim,idim,sil,a,.false.,sll)
#else
c      call dscal(li,-1d0,a(0,0,2),1)
c      call zhmpy(a,1,idim,li,sil,idim,1,li,sll,ldim,.false.,idim)
c      call dscal(li,-1d0,a(0,0,2),1)
c      call dscal(li,-1d0,sil(0,0,2),1)
c      call zhmpy(sil,1,idim,li,a,idim,1,li,sll,ldim,.false.,idim)
#endif
C     call tcx('call yyhmpy')

C --- xsi^dot := -(beta-alpha)^-1(beta^dot-alpha^dot)(beta-alpha)^-1 ---
      call scalpv(adot(1+ldim),idim,1,1,alpha(1+ldim),0)
      call scalpv(adot(1+ldim),idim,1,1,alpha(1+ldim),0)
      call dscal(idim,-1d0,adot(1+ldim),1)

C --- S^dot_il -->  - (xsi^dot - S^dot)_ii A_il ---
C     call tcn('call zmpy')
      call daxpy(idim,-1d0,adot(1+ldim),1,sii,idim+1)
      call zmpy(sii,idim,1,i2,a,idim,1,li,sil,idim,1,li,idim,ldim,idim)
C     call tcx('call zmpy')

C --- S^dot_ll --> S^dot_ll - A_li (xsi^dot - S^dot)_ii A_il ---
C     call tcn('call yyhmpy')
#ifndef CRAY
      call yyhmpy('c','n',ldim,idim,a,sil,.true.,sll)
#else
c      call dscal(li,-1d0,a(0,0,2),1)
c      call zhmpy(a,1,idim,li,sil,idim,1,li,sll,ldim,.true.,idim)
#endif
C     call tcx('call yyhmpy')
      call tcx('trSdot')
      end

