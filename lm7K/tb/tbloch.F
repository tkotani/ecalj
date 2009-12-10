      subroutine tbloch(lgamma,bk,nl,nspc,nspu,ispu,nbas,plat,lmx,ipc,
     .  indxsh,nsite,iax,npr,s,vso,hso,lso1,ldim,sk,indxH)
C- Bloch transform of real-space matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   lgamma: gamma point only, sk real
Ci   bk:    k point / 2 pi
Ci   nl,nbas,plat,lmx,ipc;
Ci   nspc = 2 for coupled spins (empirical S-O)
Ci   nspu = 2 if TB+U, ispu is then current spin
Ci   indxsh: pointers for constructing permutated H and S from makidx
Ci   nsite: total number of neighbors in all clusters;
Ci   iax: neighbor lists; npr: see tbham;
Ci   s:     real-space matrix that is to be summed;
Ci   vso,hso: table of spin-orbit parameters and the hamiltonian;
Ci   lso1:  if T add spin-orbit terms
Ci   ldim: dimension of sk
Ci   indxH:  work array of length nbas
Co Outputs
Co   sk: Bloch sum of real-space matrix S
Cr Remarks
Cr   s(k;r1,l1,r2,l2) = sum_T s(r1,l1,T+r2,l2) * exp(-i k . T)
Cr   where r1 and r2 are basis vectors and T = t2-t1 is the difference
Cr   in primitive lattice translation vectors.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nspc,nspu,ispu,nbas,nsite,ldim,niax
      parameter (niax=10)
      integer lmx(*),ipc(nbas),indxsh(*),iax(0:niax-1,nsite),
     .        npr(0:1,nbas),indxH(nbas)
      double precision bk(3),plat(3,3)
      double precision s(nl**2,nl**2,nsite*nspc**2,nspu),vso(0:nl-1,1),
     .                 hso(nl**2,nl**2,4,2),sk(ldim,ldim*2)
      logical lgamma,lso1

C Local parameters
      integer iat,isite,j,k,ix0,ixR,lm0,lmR,nlm0,nlmR,i,l2,lsp,icl0
      integer nlm,ll,iprint
      double precision twopi,TdotK,cosT,sinT,v0
      logical lso

C External procedures
      external zinit

C Intrinsic functions
      intrinsic datan, dsin, dcos, dabs
      nlm(i) = (1+lmx(i))**2

      l2 = ldim**2
      lsp = ldim / nspc
      if (lgamma) then
        call dpzero(sk,l2)
      else
        call zinit(sk,l2)
      endif
      twopi = 8*datan(1d0)

C --- Form table of indices that mark offset of Rth block in S(k) ---
      indxH(1) = 0
      do  iat = 2, nbas
        indxH(iat) = indxH(iat-1) + nlm(ipc(iat-1))
      enddo

C --- For each RR' pair add contribution to Bloch sum ---
      if (nspc .eq. 2) goto 1
      do  isite = 1, nsite
        if (.not. lgamma) then
          TdotK = 0
          do  j = 1,3
            do  k = 1,3
              TdotK = TdotK + twopi*bk(j)*plat(j,k)*iax(1+k,isite)
            enddo
          enddo
          cosT = dcos(TdotK)
          sinT = dsin(TdotK)
        endif
        ix0 = indxH(iax(0,isite))
        ixR = indxH(iax(1,isite))
        nlm0 = nlm(ipc(iax(0,isite)))
        nlmR = nlm(ipc(iax(1,isite)))
        do  lm0 = 1, nlm0
          do  lmR = 1, nlmR
            j = indxsh(ixR+lmR)
            i = indxsh(ix0+lm0)
            if (i .le. lsp .and. j .le. lsp) then
              if (lgamma) then
                sk(i,j) = sk(i,j) + s(lm0,lmR,isite,ispu)
              else
                sk(i,j) = sk(i,j) + s(lm0,lmR,isite,ispu)*cosT
                sk(l2+i,j) = sk(l2+i,j) + s(lm0,lmR,isite,ispu)*sinT
              endif
            endif
          enddo
        enddo
      enddo
        if (lgamma .and. ((iprint() .gt. 30 .and. ldim .lt. 19)
     .       .or. iprint() .gt. 60)) then
        print *, 'TBLOCH real hamiltonian :'
        do  j = 1, ldim
          write (*,10) (sk(j,k),k=1,ldim)
        enddo
      endif
   10 format (18f10.6)
      return

C --- NSPIN=2 and spin-orbit branch ---
    1 continue 
      if (lgamma) call rx('tbloch: nspc not set up for gamma only')
      lso = lso1
      do  isite = 1, nsite
        TdotK = 0
        do  j = 1,3
          do  k = 1,3
            TdotK = TdotK + twopi*bk(j)*plat(j,k)*iax(1+k,isite)
          enddo
        enddo
        cosT = dcos(TdotK)
        sinT = dsin(TdotK)

C --- Determine whether to add spin-orbit interaction terms ---
        if (lso1) then
          lso = nl .ge. 2 .and. isite .eq. npr(1,iax(0,isite))+1
          icl0 = ipc(iax(0,isite))
        endif

        ix0 = indxH(iax(0,isite))
        ixR = indxH(iax(1,isite))
        nlm0 = nlm(ipc(iax(0,isite)))
        nlmR = nlm(ipc(iax(1,isite)))
        do  lm0 = 1, nlm0
          do  lmR = 1, nlmR
            j = indxsh(ixR+lmR)
            i = indxsh(ix0+lm0)
            if (i .le. lsp .and. j .le. lsp) then
              sk(i,j) = sk(i,j) + s(lm0,lmR,isite,ispu)*cosT
              sk(l2+i,j) = sk(l2+i,j) + s(lm0,lmR,isite,ispu)*sinT
              sk(i,j+lsp) = sk(i,j+lsp)
     .          + s(lm0,lmR,isite+2*nsite,ispu)*cosT
              sk(i+lsp,j) = sk(i+lsp,j)
     .          + s(lm0,lmR,isite+1*nsite,ispu)*cosT
              sk(i+lsp,j+lsp) = sk(i+lsp,j+lsp)
     .          + s(lm0,lmR,isite+3*nsite,ispu)*cosT
              sk(l2+i,j+lsp) = sk(l2+i,j+lsp)
     .          + s(lm0,lmR,isite+2*nsite,ispu)*sinT
              sk(l2+i+lsp,j) = sk(l2+i+lsp,j)
     .          + s(lm0,lmR,isite+1*nsite,ispu)*sinT
              sk(l2+i+lsp,j+lsp) = sk(l2+i+lsp,j+lsp)
     .          + s(lm0,lmR,isite+3*nsite,ispu)*sinT
              if (lso) then
C --- Spin-orbit interaction terms ---
                v0 = vso(ll(lm0),icl0)
                sk(i,j) = sk(i,j) + v0*hso(lm0,lmR,1,1)
                sk(i,j+lsp) = sk(i,j+lsp) + v0*hso(lm0,lmR,2,1)
                sk(i+lsp,j) = sk(i+lsp,j) + v0*hso(lm0,lmR,3,1)
                sk(i+lsp,j+lsp) = sk(i+lsp,j+lsp)
     .            + v0*hso(lm0,lmR,4,1)
                sk(l2+i,j) = sk(l2+i,j)
     .            + v0*hso(lm0,lmR,1,2)
                sk(l2+i,j+lsp) = sk(l2+i,j+lsp)
     .            + v0*hso(lm0,lmR,2,2)
                sk(l2+i+lsp,j) = sk(l2+i+lsp,j)
     .            + v0*hso(lm0,lmR,3,2)
                sk(l2+i+lsp,j+lsp) = sk(l2+i+lsp,j+lsp)
     .            + v0*hso(lm0,lmR,4,2)
              endif
            endif
          enddo
        enddo
      enddo

      end
