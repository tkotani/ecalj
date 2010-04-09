      subroutine salph2(nbas,nenv,lqinv,ldot,loka,siinv,ndimW,lmaxW,
     .  nenvx,npr,plat,bas,alphv,adotv,iax,cy,cg,indxcg,jcg,el,
     .  s0a,s0b,s01,s2,wk,kpvt,salphx,sdotx)

C- Calculate one-center expansion coefficients for one cluster, 2 energies
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nenv  :number of envelope function on head sphere; see Bugs
Ci   lqinv :switch governing how matrix inversion is accomplished
Ci          0 uses lapack dsytrf
Ci          1 uses dqinvb (not implemented)
Ci   ldot  :T, make energy derivatives of bdot also
Ci   loka: :Old Andersen conventions for Hankel and Bessels
Ci   siinv: parameters for iterative inversion of strux (not used now)
Ci          1: number of sites for cutoff
Ci          2: number of iterations
Ci          3: allowed tolerance
Ci   ndimW :dimensions of strux, including space for Watson sphere
Ci   lmaxW :maximum l for Watson-sphere (not used now)
Ci   nenvx :max number of envelope function on head sphere (eg nl**2)
Ci   npr   :number of pairs in cluster about this site
Ci   plat  :primitive lattice vectors, in atomic units
Ci   bas   :basis vectors, in atomic units
Ci   alphv :vector of inverse of alphas, of length ndimW
Ci   adotv :vector of adot/alpha**2, of length ndimW
Ci         :(ldot=F) not used
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   el    :energy of Hankels, atomic units
Ci   Work arrays:
Ci   s0a   :matrix s0 for first energy, dimension ndimW*ndimW
Ci   s0b   :matrix s0 for second energy, dimension ndimW*ndimW
Ci   s01   :for iterative solutions (not used now)
Ci   s2    :for iterative solutions (not used now)
Ci   wk    :work array, of length ndimW*max(2+4*nenv,64)
Ci          Used in the following contexts:
Ci          work space for LU decomposition of s0a,s0b
Ci          wk(1,ndimW) used to hold s12
Ci          wk(ndimW+1,2*ndimW) used to hold s21
Ci          The last ndimS*4*nenv used by dsolv2
Co Outputs
Co  salphx :structure constants, same as salph but with fixed dimensions
Co  sdotx  :energy derivative of salphx 
Cl Local variables
Cl  salph  :structure constants as described in Remarks
Cl  sdot   :energy derivative of salph 
Cb Bugs
Cb   This routine should have the following:
Cb  ndimW  :The dimension of strux matrix, including the Watson sphere
Cb  ndimS  :Number of functions to be determined by direct inversion
Cb  nenvx  :max number of envelope functions on head sphere 
Cb  nenv   :actual number of envelope functions on head sphere,
Cb          used as a dimension of local arrays salph and sdot       
Cb
Cb ldot=T is not implemented
Cr Remarks
Cr   Salph is treated as an array of dimension (ndimW,2,nenv,2) or
Cr   (for sim eqns) as an array of dimension  (2*ndimW,2*nenv):
Cr       nenv         nenv
Cr      +-----+     +-----+
Cr      |     |     |     |
Cr      | E=0 |     | E=0 |   ndimW
Cr      | E'=0|     | E'=1|
Cr      |     |     |     |
Cr      +- - -+     +- - -+
Cr      |     |     |     |
Cr      | E=1 |     | E=1 |
Cr      | E'=1|     | E'=1|   ndimW
Cr      |     |     |     |
Cr      +-----+     +-----+
Cr   The 2X2 block in energy for a given 0L,R'L':
Cr     columns are separated by 2*ndimW*nenv elements; rows by ndimW
Cr   The nenv X nenv block in for a given E0,E'R':
Cr     columns are separated by 2*ndimW elements; rows by 1
Cr
Cr   for debugging
Cr   mc -f6f17.12 out.gas -p -v2dia -v2dia -- -a b1 tmp/out.gas -p -v2dia -v2dia -- -a b2 -1:nr -s0 -a z b1 z -ccat z b2 -ccat -rcat >b0
Cr   mc -f6f17.12 alp -i -a ialp ialp ialp b0 -+ -i ialp -x -x  ialp -tog -- >b
Cr   mc -f6f17.12 alp -i -a ialp ialp b0 -+ -i ialp -x >s
Cu Updates
Cu   23 Jul 08 Adjusted to work with species-dependent lmax
Cu   03 Aug 06 Adapted from ASA v5
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lqinv,loka,nbas,ndimW,nenv,npr,kpvt(1),niax,lmaxW,nenvx
      logical ldot
      parameter (niax=10)
      integer iax(niax,npr)
      double precision alphv(ndimW,2,2),adotv(ndimW,2,2)
C Do not specifiy dimensions of s0a,s0b, for future iterative solutions
      double precision bas(3,nbas),plat(3,3),wk(*),el(2),siinv(5),
     .  s0a(*),s0b(*),s01(ndimW,*),s2(*),
     .  salphx(ndimW,2,nenvx,2),sdotx(ndimW,2,nenvx,2)
C For Clebsch Gordon coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local parameters
      integer j,i,lm1,lm2
      double precision tmp(2,2)
      double precision salph(ndimW,2,nenv,2),sdot(ndimW,2,nenv,2)

      call tcn('salph2')
      
      if (ldot)  
     .  call rx('salph2: sdot is not implemented for double-kappa SSSW')

C     call xxalp2(ndimW,alphv)
C     call invalp(nenv*nbas,2,alpha)

C      print *, '!!'
C      alpha = 0

      i = 11
      if (ldot) i = 30
      if (loka .eq. 1) i = i+1200

      call maks0(npr,plat,bas,indxcg,jcg,cg,cy,
     .  alphv,adotv,iax,el(1),i,lmaxw,s0a,s0a,ndimW)

      call maks0(npr,plat,bas,indxcg,jcg,cg,cy,
     .  alphv(1,2,2),adotv(1,2,2),iax,el(2),i,lmaxw,s0b,s0b,ndimW)

C     call prmx('s0a',s0a,ndimW,ndimW,ndimW)
C     call prmx('s0b',s0b,ndimW,ndimW,ndimW)

C     call xxalp2(nenv,nbas,ndimW,npr,iax,alpha)
C     call xxalp3(ndimW,s0a,s0b)

      call decmp2(s0a,alphv(1,1,2),alphv(1,2,1),s0b,ndimW,kpvt,wk)

C --- Add term delta_RL-,R'L' to rhs ---
      call dpzero(salph, 4*ndimW*nenv)
      if (ldot) call dpzero(sdot, 4*ndimW*nenv)
      do  i = 1, nenv
        salph(i,1,i,1) = salph(i,1,i,1) + 1
        salph(i,2,i,2) = salph(i,2,i,2) + 1
      enddo

C --- (alp^-1 + s0)^-1 ---
      call dsolv2(s0a,alphv(1,1,2),alphv(1,2,1),s0b,ndimW,kpvt,
     .            wk,salph,2*nenv)

C     call prmx('(alpha^-1 + s0)^-1',salph,2*ndimW,2*ndimW,2*nenv)

C --- alp^-1 * (alp^-1 + s0)^-1 * alp^-1 ---
      do  i = 1, ndimW
        do  lm2 = 1, nenv
C        (alp^-1 + s0)^-1 alp^-1
          call dmx22(salph(i,1,lm2,1),ndimW*2*nenv,ndimW,
     .              alphv(lm2,1,1),ndimW*2,ndimW,tmp,2,1)
C         call dmpy(salph(i,1,lm2,1),ndimW*2*nenv,ndimW,
C                   alphv(lm2,1,1),ndimW*2,ndimW,tmp,2,1,2,2,2)
C         alp^-1 (alp^-1 + s0)^-1 alp^-1
          call dmx22(alphv(i,1,1),ndimW*2,ndimW,tmp,2,1,
     .              salph(i,1,lm2,1),2*ndimW*nenv,ndimW)
C         call dmpy(alphv(i,1,1),ndimW*2,ndimW,tmp,2,1,
C    .              salph(i,1,lm2,1),2*ndimW*nenv,ndimW,2,2,2)
        enddo
      enddo

C --- alp^-1 -  alp^-1 (alp^-1 + s0)^-1 alp^-1 ---
      call dscal(ndimW*2*nenv*2,-1d0,salph,1)
      do  i = 1, 2
      do  j = 1, 2
      do  lm2 = 1, nenv
        salph(lm2,i,lm2,j) = salph(lm2,i,lm2,j) + alphv(lm2,i,j)
      enddo
      enddo
      enddo

C     call prmx('alpha^-1  -  alpha^-1 (alpha^-1 + s0)^-1 alpha^-1',
C    .  salph,2*ndimW,2*ndimW,2*nenv)

C --- copy salph and sdot into fixed dimension arrays
C     call dpzero(salphx, 4*ndimW*nenvx)
C     call dpzero(sdotx, 4*ndimW*nenvx)
      do  j = 1, 2
        call dcopy(2*ndimW*nenv,salph(1,1,1,j),1,salphx(1,1,1,j),1)
        if (ldot) 
     .    call dcopy(2*ndimW*nenv,sdot(1,1,1,j),1,sdotx(1,1,1,j),1)
      enddo

C --- Restore alpha before exit ---
C     call invalp(nenv*nbas,2,alpha)

      call tcx('salph2')
Cgetarg       end
       end subroutine salph2 


C      subroutine xxalp2(ndimW,alphv)
CC- print 2-kap alpha in matrix form
C      implicit none
C      integer ndimW
C      double precision alphv(ndimW,2,2)
C      double precision b(ndimW,2,ndimW,2)
C      integer j,k,iRL
C
C      b = 0d0
C      do  j = 1, 2
C      do  k = 1, 2
C      do  iRL = 1, ndimW
C        b(iRL,j,iRL,k) = alphv(iRL,j,k)
C      enddo
C      enddo
C      enddo
C
C      call prmx('alpha',b,2*ndimW,2*ndimW,2*ndimW)
C      end

C      subroutine xxalp3(ndimW,b1,b2)
CC- Copy B1, B2 into 2-kap form for testing
C      implicit none
C      integer ndimW
C      double precision b1(ndimW,ndimW),b2(ndimW,ndimW)
C      double precision b0(ndimW,2,ndimW,2)
C      double precision wk(ndimW,ndimW)
C
C      b0 = 0
C
C      call dcopy(ndimW**2,b1,1,wk,1)
C      call dupack(wk,ndimW)
C      b0(:,1,:,1) = wk(:,:)
C      call dcopy(ndimW**2,b2,1,wk,1)
C      call dupack(wk,ndimW)
C      b0(:,2,:,2) = wk(:,:)
C
C      call prmx('b0',b0,2*ndimW,2*ndimW,2*ndimW)
C      end


