      subroutine invbt1(lov,nbas,nl,nsp,ldim,plat,nsites,npr,iax,lmx,
     .                  indxsh,ipc,ik,nk,bk,wtkp,hk,ok,indxH,hrs,ors)
C- Inverse Bloch transform to extract real space hamiltonian and overlap
C ----------------------------------------------------------------------
Ci Inputs:  lov,nbas,nl,ldim,alat,pos,npr,iax,nk,hk,ok
Ci          indxH:  work array of length nbas
Co Outputs: none, matrix elements tabulated to disc
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer niax
      parameter (niax=10)
      integer nbas,nl,ldim,nsp,nsites,npr(0:nbas-1),indxH(nbas),lmx(1),
     .        ik,nk,iax(0:niax-1,nsites),indxsh(1),ipc(nbas)
      logical lov
      double precision bk(3),plat(3,3),hk(ldim,ldim,2),ok(ldim,ldim,2),
     .                 wtkp(nk),
     .                 hrs(nl**2,nl**2,nsites),ors(nl**2,nl**2,nsites)
C Local Variables
      integer ifi,i,j,k,iat,isite,nlm,ix0,ixR,lm0,
     .        lmR,nlm0,nlmR
      integer fopn,i1mach,iprint
      double precision TdotK,twopi,cosT,sinT,wt
      double precision dcos

C Intrinsic functions
      intrinsic datan, dcos, dsin
      nlm(i) = (1+lmx(i))**2

      call rxx(nsp.ne.1,'INVBTB: not set up for spin')
      
      wt = wtkp(ik) / 2d0
      ifi = i1mach(2)
      if (bk(1) .eq. 0d0 .and. bk(2) .eq. 0d0 .and. bk(3) .eq. 0d0
     .    .and. iprint() .ge. 30) then
        call awrit0('Hamiltonian at Gamma:',' ',120,ifi)
        do  i = 1, ldim
          write (ifi,1) (hk(i,j,1),j=1,ldim)
        enddo
        if (lov) then
          call awrit0('Overlap at Gamma:',' ',120,ifi)
          do  i = 1, ldim
            write (ifi,1) (ok(i,j,1),j=1,ldim)
          enddo
        endif
      endif
    1 format (128f10.6)

      twopi = 8*datan(1d0)
C --- Form table of indices that mark offsets to unpermuted h(k) ---
      indxH(1) = 0
      do  iat = 2, nbas
        indxH(iat) = indxH(iat-1) + nlm(ipc(iat-1))
      enddo

C      write(*,3)
C    3 format('isite iax0 iax1 ix0 ixR',
C     .       3x'lm0',2x,'lmR',3x,'i',4x,'j',7x,'hr',8x,'cosT',10x,'hi',
C     .       8x,'sinT',8x'hrs')

C --- loop over all pairs ---
      do  isite = 1, nsites
        TdotK = 0
        do  j = 1, 3
          do  k = 1, 3
            TdotK = TdotK + twopi*bk(j)*plat(j,k)*iax(1+k,isite)
          enddo
        enddo
        cosT = dcos(TdotK)
        sinT = dsin(TdotK)
        ix0 = indxH(iax(0,isite))
        ixR = indxH(iax(1,isite))
        nlm0 = nlm(ipc(iax(0,isite)))
        nlmR = nlm(ipc(iax(1,isite)))
        do  lm0 = 1, nlm0
          do  lmR = 1, nlmR
            j = indxsh(ixR+lmR)
            i = indxsh(ix0+lm0)

       
       ifi = fopn('RSH')
       call awrit7('isite=%i, ia=%i, ib=%i, j=%i, i=%i, lmR=%i'//
     .             ' lm0=%i',' ',256,ifi,isite,iax(1,isite),
     .             iax(0,isite),j,i,lmR,lm0)

            if (i .le. ldim .and. j .le. ldim) then
              hrs(lm0,lmR,isite) = hrs(lm0,lmR,isite) +
     .                          (hk(i,j,1)*cosT + hk(i,j,2)*sinT) * wt

C        write(*,2) isite,iax(0,isite),iax(1,isite),ix0,ixR,
C     .             lm0,lmR,i,j,hk(i,j,1),cosT,hk(i,j,2),sinT,
C     .             hrs(lm0,lmR,isite)
C    2   format (1x,3i4,1x,2i4,4i5,5f12.6)

              if (lov) then
                ors(lm0,lmR,isite) = ors(lm0,lmR,isite) +
     .                            (ok(i,j,1)*cosT + ok(i,j,2)*sinT) * wt
              endif
            endif
          enddo
        enddo
      enddo
      end

