      subroutine vmix2(nclass,dclabl,betv,vold,vnew)
C- Extra Madelung potential mixing
C ----------------------------------------------------------------
Ci Inputs
Ci   nclas,clabl
Ci   betv,vold,vnew (see Outputs)
Co Outputs
Co   vnew is replaced by beta*vnew + (1-beta)*vold.
Cr Remarks
Cr   Purpose is to mix potentials more slowly than the charges owing
Cr   to the amplification in potential changes by the dielectric
Cr   constant.  vold should be starting potential; vnew the potential
Cr   arising from the mixed charge densities.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nclass
      double precision betv,vold(nclass),vnew(nclass),dclabl(1)
C Local variables
      integer i,ic,lgunit,ipr
      double precision vmix,rms
      character*8 clabl*8, outs*80

      call getpr(ipr)
      call awrit1(' vmix2: independent potential mixing with betv = %d',
     .  outs,80,0,betv)
      if (ipr .gt. 30) then
        do  5  i = 1, 2
          call awrit0('%a',outs,80,-lgunit(i))
          write(lgunit(i),500)
    5   continue
  500   format(' class        Vold       Vnew       diff      Vmix')
      endif

      rms = 0
      do 10  ic = 1, nclass
        vmix = betv*vnew(ic) + (1-betv)*vold(ic)
        call r8tos8(dclabl(ic),clabl)
        do  12  i = 1, 2
   12   if (ipr .gt. 30) write(lgunit(i),783)
     .    clabl,vold(ic),vnew(ic),vnew(ic)-vold(ic),vmix
  783   format(2x,a8,4f11.6)
        rms = rms + (vnew(ic)-vold(ic))**2
        vnew(ic) = vmix
   10 continue
      rms = dsqrt(rms/nclass)

      if (ipr .ge. 10 .and. ipr .le. 30)
     .  call awrit1('%a rms delta V = %d',outs,80,-lgunit(1),rms)
      call awrit1(' vmix2: rms delta V = %d',outs,80,-lgunit(2),rms)
      if (ipr .gt. 30) call awrit0('%a',outs,80,-lgunit(1))
      if (ipr .ge. 20) print '(1x)'

      end

