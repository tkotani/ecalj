      subroutine pa2plm(pa,nbas,lpdim,ldim,indxsh,iopt,nf,p)
C- Expand l- independent matrix into matrix l- or lm- components
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpdim  dimension of p
Ci   nbas
Ci   iopt   1s digit
Ci          0  one-dimensional matrix
Ci          1  two-dimensional matrix
Ci          10s digit
Ci          0  p has dimension nbas
Ci          1  p has dimension nbas*nl
Ci          2  p has dimension nbas*nlm
Ci       1000s digit
Ci          1  do not initialize pa to zero before starting
Co Outputs
Co   p:  the l-indepent parts of pa copied into l- (lm-) dependent p
Cr Remarks
Cr   pa2plm performs the inverse function to plm2pa 
C ----------------------------------------------------------------------
C     implicit none
      integer lpdim,nbas,ldim,indxsh(1),iopt,nf
      double precision pa(nbas,nbas),p(lpdim,lpdim)
C Local variables
      integer ib,ilm,jb,jlm,nli,nlj,norbi,norbj,ntorbi,ntorbj,nd,
     .  opt1,n0,offj,offi,opt3,nkap0
      parameter (n0=10,nkap0=3)
      integer ltabi(n0*nkap0),ltabj(n0*nkap0),ktab(n0*nkap0),
     .        offl(n0*nkap0)
C     double precision scale

      nd = 1
      if (mod(iopt,10) .eq. 1) nd = 2
      opt1 = mod(iopt/10,10)
      opt3 = mod(iopt/1000,10)

C ... Simple copy if no indices to contract over
      if (opt1 .eq. 0) then
        call dcopy(nbas**nd,pa,1,p,1)
        return
      endif

      if (opt3 .eq. 0) call dpzero(p, lpdim**nd)

      nlj = 0
      do  20  jb = 1, nbas
C     Uses norbj,ntorbj
      call orbl(jb,0,ldim,indxsh,norbj,ltabj,ktab,offj,offl,ntorbj)
      if (opt1 .eq. 1) ntorbj = norbj
      do  20  jlm = 1, ntorbj
        nlj = nlj+1

        if (nd .eq. 1) then
          p(nlj,1) = pa(jb,1)
        else
          nli = 0
          do  30  ib = 1, nbas
C         Uses norbi,ntorbi
          call orbl(ib,0,ldim,indxsh,norbi,ltabi,ktab,offi,offl,ntorbi)
          if (opt1 .eq. 1) ntorbi = norbi
          do  30  ilm = 1, ntorbi
              nli = nli+1
              p(nli,nlj) = pa(ib,jb)
   30     continue
        endif

   20 continue

      end

