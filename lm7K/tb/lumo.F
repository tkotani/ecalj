      subroutine lumo(lmol,zval,ldim,nsp,nspc,nbmax,eband,wtkb)
C- Try and make a HOMO-LUMO excitation by manipulating band weights
C ----------------------------------------------------------------------
Ci Inputs:
Ci   qval, nsp, wtkb
Co Outputs:
Co   wtkb
Cr Remarks
Cr   Try to make a HOMO-LUMO excitation by exchanging the weights of the
Cr   HOMO and LUMO.
Cb Bugs
Cb   If the HOMO or LUMO are degenerate then this routine wil fail
Cb   with an error exit. In that case the weights should have to be
Cb   remade as in evlwsm in mol; this is not yet implemented in tbe.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nsp,nspc,ldim,nbmax
      double precision zval,wtkb(ldim,2),eband(nbmax,2)
      logical lmol
C Local Variables
      integer obot,otop,obmap,owk
      integer i,isp,ibh,isph,ibl,ispl,nbpw,nkp,nevxx,nbmxx
      integer i1mach,iprint
      double precision elo,ehi,ehomo,elumo,wt
      logical lfill,efrng2,fndh,fndl,T,F
C Dynamic memory
      integer w(1)
      common /w/ w
      data T /.true./ F /.false./

C --- Sanity checks ---
      call rxx(.not.lmol, 'Use LUMO for molecules only')
      call rxx(nsp.ne.2,  'LUMO needs spin pol')
      call rxx(nspc.eq.2, 'LUMO needs npsc=1')

      if (iprint() .gt. 40) then
        print *,'  Bands and weights ..'
        print *,'  Eigenvalue            weight'
        write (*,10) ((eband(i,isp),wtkb(i,isp),i=1,ldim),isp=1,nsp)
      endif

      nkp=1
      nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
      call defi(obmap,-(ldim*nsp*nkp/nbpw+1))
      call defdr(owk, ldim*nsp)
      call ebcpl(0,nbmax,ldim,nsp,nspc,nkp,nbpw,w(obmap),w(owk),eband)
      call defdr(obot,ldim*nsp)
      call defdr(otop,ldim*nsp)
      nevxx = ldim*nsp
      nbmxx = nbmax*nsp
      lfill = efrng2(1,nkp,nbmxx,nevxx,zval*2,eband,w(obot),w(otop),
     .               elo,ehi,ehomo,elumo)
      if (iprint() .ge. 40) then
        call awrit2(' LUMO efrng2 found eHOMO=%d, eLUMO=%d',' ',120,
     .              i1mach(2),ehomo,elumo)
      endif
      call ebcpl(1,nbmax,ldim,nsp,nspc,nkp,nbpw,w(obmap),w(owk),eband)
      call rlse(obmap)

      isph = 0
      fndh = F
      do  i = 1, ldim
        do  isp = 1, nsp
          if (dabs(eband(i,isp)-ehomo) .lt. 1d-12) then
            if (fndh) then
              if (isp .eq. isph) then
                call rx('LUMO: HOMO is degenerate')
              endif
            endif
            ibh = i
            isph = isp
            fndh = T
          endif
        enddo
      enddo
      fndl = F
      do  i = 1, ldim
C ... flip the spin: Hund's rule 
        do  isp = nsp, 1, -1
          if (dabs(eband(i,isp)-elumo) .lt. 1d-12) then
            if (fndl) then
              if (isp .eq. ispl) then
                call rx('LUMO: LUMO is degenerate')
              endif
            endif
            ibl = i
            ispl = isp
            fndl = T
          endif
        enddo
      enddo

      if (iprint() .ge. 30) then
        call awrit4(' LUMO found eHOMO=%d (weight=%d), eLUMO=%d'//
     .              ' (weight=%d)',' ',120,i1mach(2),
     .              eband(ibh,isph),wtkb(ibh,isph),
     .              eband(ibl,ispl),wtkb(ibl,ispl))
      endif

C --- exchange weights ---
      wt = wtkb(ibh,isph)
      wtkb(ibh,isph) = wtkb(ibl,ispl)
      wtkb(ibl,ispl) = wt

      if (iprint() .gt. 31) then
        print *,'  Modified bands and weights ..'
        print *,'  Eigenvalue            weight'
        write (*,10) ((eband(i,isp),wtkb(i,isp),i=1,ldim),isp=1,nsp)
      endif
   10 format (f10.6,4x,f10.6)
      
      end
