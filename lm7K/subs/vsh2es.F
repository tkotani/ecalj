      subroutine vsh2es(mode,nbas,ipc,nrc,vshft,ves)
C- Adds/subtracts vshft from ves to produce true effective potentials
C ----------------------------------------------------------------
C  mode   1s digit   0 add vshft to ves   1 subtract from ves
C        10s digit   1 initialize input ves to zero
C       100s digit   1 reverse role of ves and vshft
C ----------------------------------------------------------------
C     implicit none
      integer mode,nbas,ipc(nbas),nrc(nbas)
      double precision vshft(-7:nbas),ves(nbas),scal
      logical mode2
      integer ib,ic,iprint,lgunit,nclass

C ... Initialize
      mode2 = mod(mode/100,10) .eq. 1
      if (mod(mode/10,10) .eq. 1)  then
        nclass = 0
        do  10  ib = 1, nbas
        ic = ipc(ib)
        nclass = max(nclass,ic)
        if (.not. mode2) ves(ic) = 0
        if (mode2) vshft(ib) = 0
   10 continue
      endif

C ... Scale +/- 1
      scal = 1
      if (mod(mode/1,10) .eq. 1) scal = -1

C ... Do the copy
      do  20  ib = 1, nbas
        ic = ipc(ib)
        if (.not. mode2) ves(ic) = ves(ic) + scal*vshft(ib)/nrc(ic)
        if (mode2) vshft(ib) = vshft(ib) + scal*ves(ic)
   20 continue

      if (iprint() .ge. 30) then

        if (.not. mode2) write(lgunit(1),500) 'Class'
        if (mode2) write(lgunit(1),500) ' Site'
  500   format(/1x,a,3x,'Shifted V')
  501   format(i4,f14.6)
        if (mode2) then
          do  30  ic = 1, nclass
   30     write(lgunit(1),501) ic, ves(ic)
        else
          do  32  ib = 1, nbas
   32     write(lgunit(1),501) ib, vshft(ib)
        endif
      endif

      end

