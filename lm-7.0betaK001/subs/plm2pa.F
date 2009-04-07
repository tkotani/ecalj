      subroutine plm2pa(p,nbas,lpdim,ldim,indxsh,iopt,nf,pa)
C- Contract l- or lm- dependent components of real p into real pa
C ----------------------------------------------------------------------
Ci Inputs
Ci   iopt   1s digit
Ci          0  one-dimensional matrix
Ci          1  two-dimensional matrix
Ci         10s digit
Ci          0  dimension nbas
Ci          1  dimension nbas*nl
Ci          2  dimension nbas*nlm
Ci        100s digit
Ci          1  accumulate p/nl
Ci          2  accumulate p/nlm
Ci       1000s digit
Ci          1  do not initialize pa to zero before starting
Ci   P     :ASA Green function like object G_ij (or slice),
Ci         :depending on opt1 = 10s digit iopt
Ci         :P can be one of:
Ci          :  opt1 = 0 => G contracted over l and m 
Ci          :  opt1 = 1 => G contracted over m
Ci          :  opt1 = 2 => full G_ij = (not contracted)
Ci   nbas  :size of basis
Ci   lpdim :dimensions P; 
Ci          :  opt1 = 0 => lpdim = nbas
Ci          :  opt1 = 1 => lpdim = nnrl(1) ~ nbas*nl (see function nnrl)
Ci          :  opt1 = 2 => lpdim = nnrl(0) ~ nbas*nlm (see function nnrl)
Co Outputs
Co   pa: all l- or lm-dependent parts of p are contracted.
C ----------------------------------------------------------------------
C     implicit none
      integer lpdim,nbas,ldim,indxsh(1),iopt,nf
      double precision pa(nbas,nbas,nf),p(lpdim,lpdim,nf)
C Local variables
      integer ib,ii,ilm,jb,jlm,nd,nli,nlj,n0,nkap0,
     .  norbi,norbj,ntorbi,ntorbj,offi,offj,opt1,opt2,opt3
      parameter (n0=10,nkap0=3)
      integer ltabi(n0*nkap0),ltabj(n0*nkap0),ktab(n0*nkap0),
     .        offl(n0*nkap0)

      double precision scale

      nd = 1
      if (mod(iopt,10) .eq. 1) nd = 2
      opt1 = mod(iopt/10,10)
      opt2 = mod(iopt/100,10)
      opt3 = mod(iopt/1000,10)

      if (opt3 .eq. 0) call dpzero(pa, nbas**nd*nf)

C ... Simple copy if no indices to contract over
C     BUG: assumes p is dimensioned (nbas,nbas)
      if (opt1 .eq. 0) then
        call daxpy(nbas**nd*nf,1d0,p,1,pa,1)
        return
      endif

      nlj = 0
      do  20  jb = 1, nbas
C     Uses norbj,ntorbj
      call orbl(jb,0,ldim,indxsh,norbj,ltabj,ktab,offj,offl,ntorbj)
      scale = 1
      if (opt2 .eq. 1) scale = 1/dble(norbj)
      if (opt1 .eq. 1) ntorbj = norbj
      do  20  jlm = 1, ntorbj
        nlj = nlj+1
C       Sanity check
C       if (nlj .gt. lpdim) call rx('bad indexing in plm2pa')

        if (nd .eq. 1) then
          do  22  ii = 1, nf
   22     pa(jb,ii,1) = pa(jb,ii,1) + p(nlj,ii,1)*scale
        else
          nli = 0
          do  30  ib = 1, nbas
C         Uses norbi,ntorbi
          call orbl(ib,0,ldim,indxsh,norbi,ltabi,ktab,offi,offl,ntorbi)
          scale = 1
          if (opt2 .eq. 1) scale = 1/dble(norbi)
          if (opt1 .eq. 1) ntorbi = norbi
          do  30  ii = 1, nf
          do  30  ilm = 1, ntorbi
            nli = nli+1
            pa(ib,jb,ii) = pa(ib,jb,ii) + p(nli,nlj,ii)*scale
   30     continue
        endif

   20 continue

      end

