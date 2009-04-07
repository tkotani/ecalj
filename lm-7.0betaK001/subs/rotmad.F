      subroutine rotmad(nbas,nclass,ipc,rmax,dmad,ves1,ves2,
     .  nrclas,emad0,modcst,dmad2,igroup,evmad,u,v1til,v2til,nvmix)
C- Rotates Madelung matrix, calculate rotated Madelung potential
C ----------------------------------------------------------------
Ci Inputs
Ci   same as madpot
Ci   dmad2: temporary work array of length nbas**2
Ci   modcst: 0, uses ordinary Madelung matrix
Ci           1, averages potentials according to igroup
Ci          -1, same as 1, but just two groups, depending on sgn(igroup)
Ci           2, rotates by eigenvector of Madelung matrix, with idea
Ci              to weight low eval differently from high.
C    emad0:  freeze rotated Mad to vtil for the emad(j) .lt. emad0
Co Outputs
Co   evmad: eigenvalues of madelung matrix  (modcst=2 only)
Co          ratio of M2/M1 (modcst=-1 only)
Co   u:     eigenvectors of madelung matrix (modcst=2 only)
Co   v1til:  ves2 * u
Co   v2til:  ves1 * u
Cr Remarks                                (M2  -M1)
Cr   In the case modcst=-1, the matrix M= (       )  is constructed
Cr                                        (M1   M2)
Cr
Cr   and vtil is replaced with M * vtil.  In vmix the change in
Cr   this rotated vtil(2) is constrained to be zero.
C ----------------------------------------------------------------
C     implicit none
      integer nbas,nclass,modcst,igroup(7),nvmix,nrclas(1)
      double precision ves1(1),rmax(1),evmad(2),emad0
      double precision dmad(nbas,nbas),ves2(1),
     . dmad2(nbas,nbas),u(nbas,nbas),v1til(2),v2til(2)
      integer ipc(1)
C Heap
      integer w(1)
      common /w/ w
C Local variables
      character*1 modchr
      integer ibas,jbas,iprint,lgunit,ic,owk2,owk3,i
      logical keepld
      double precision e1mad, e2mad, e3mad, xx(2)
      integer ngrp,now

      goto (1,2), iabs(modcst)

C --- Case modcst = 0 ---
      nvmix = nclass
      call dcopy(nclass,ves1,1,v1til,1)
      call dcopy(nclass,ves2,1,v2til,1)
      return

C --- Case modcst = +/- 1 ---
    1 continue
      if (modcst .eq. -1) then
        do  60  ic = 1, nclass
          if (igroup(ic) .ge. 0) then
            igroup(ic) = 1
          else
            igroup(ic) =-1
          endif
   60   continue
        call vtoq(nbas,1,1,nclass,ipc,nrclas,0d0,4,
     .    0,'dum',0,0,rmax,ves1,dmad,0,0,u,u,0)
        evmad(1) = 0
        evmad(2) = 0
      endif
      call ivshel(1,nclass,igroup,igroup(nclass),.true.)

C --- Accumulate V by group ---
      if (iprint() .ge. 30) print
     .'(''Rotmad:  ic  group  nvmix    Ves1        Ves2        Diff'')'
      nvmix = 0
      ngrp = 0
      now = igroup(igroup(nclass))-1
      call dpzero(v1til,nclass)
      call dpzero(v2til,nclass)
      do  40  ic = 0, nclass-1
        if (now .ne. igroup(igroup(nclass+ic))) then
          now = igroup(igroup(nclass+ic))
          ngrp = 0
          nvmix = nvmix+1
        endif
        ngrp = ngrp+1
        i = igroup(nclass+ic)
        v1til(nvmix) = v1til(nvmix) + ves1(i)
        v2til(nvmix) = v2til(nvmix) + ves2(i)
        if (modcst .eq. -1) evmad(nvmix) = evmad(nvmix) + u(i+1,1)
C -- Printout --
        if (iprint() .ge. 30) then
          if (iprint() .ge. 40)
     .    print 333, i,now,nvmix,ves1(i),ves2(i),ves2(i)-ves1(i)
          if (ic.eq.nclass-1 .or. now.ne.igroup(igroup(nclass+ic+1)))
     .      then
            v1til(nvmix) = v1til(nvmix)/ngrp
            v2til(nvmix) = v2til(nvmix)/ngrp
            print 334, ngrp, nvmix, v1til(nvmix), v2til(nvmix),
     .      v2til(nvmix)-v1til(nvmix)
          endif
  333     format(i11,i5,i7,3f12.6)
  334     format(i6,' elts in group',i3,3f12.6,'#')
        endif
   40 continue
C --- Case modcst = -1 ---
      if (modcst .eq. -1) then
        xx(1) = dsqrt(evmad(1)**2 + evmad(2)**2)
        evmad(1) = evmad(1)/xx(1)
        evmad(2) = evmad(2)/xx(1)
        xx(1) =    evmad(1)*v1til(1) + evmad(2)*v1til(2)
        v1til(1) = evmad(2)*v1til(1) - evmad(1)*v1til(2)
        v1til(2) = xx(1)
        xx(1) =    evmad(1)*v2til(1) + evmad(2)*v2til(2)
        v2til(1) = evmad(2)*v2til(1) - evmad(1)*v2til(2)
        v2til(2) = xx(1)
        nvmix = 1
      endif
      return

C --- Case modcst = 2 ---
    2 continue
      nvmix = nbas
C --- Find eigenvectors of madelung matrix; store in mad2 ---
      call dcopy(nbas**2,dmad,1,dmad2,1)
      do  10  ibas = 1, nbas
        ic = ipc(ibas)
        dmad2(ibas,ibas) = dmad(ibas,ibas) + 1/rmax(ic)
   10 continue
      call defdr(owk2,nbas)
      call defdr(owk3,nbas)
      call rs(nbas,nbas,dmad2,evmad,1,u,w(owk2),w(owk3),i)
      if (i .ne. 0) call fexit(-1,119,'rotmad: bad madelung matrix',0)
      call rlse(owk2)

C --- Make v2tilda and v1tilda ---
      call dpzero(v2til,nbas)
      call dpzero(v1til,nbas)
      do  20  ibas = 1, nbas
        ic = ipc(ibas)
        do  20  jbas = 1, nbas
        v1til(jbas) = v1til(jbas) + ves1(ic)*u(ibas,jbas)
        v2til(jbas) = v2til(jbas) + ves2(ic)*u(ibas,jbas)
   20 continue

C --- Printout ---
      if (iprint() .gt. 100) then
        print 100
  100   format(/' Eigenvectors of madelung matrix:')
        do  32  ibas = 1, nbas
          print 110, (u(ibas,jbas),jbas=1,nbas)
  110     format(5f11.7)
   32   continue
        print 105
  105   format(/' V by bas:')
        do  34  ibas = 1, nbas
          ic = ipc(ibas)
          print 110, ves2(ic)
   34   continue

      endif

      if (iprint() .ge. 30) then
        e1mad = 0
        e2mad = 0
        e3mad = 0
        do  22  i = 1, 2
          write(lgunit(i),369)
  369     format('Rotmad:'/
     .      '  i       Emad       V1til       V2til        V2-V1')
        do  22  ibas = 1, nbas
          modchr = ' '
          keepld = (modcst.eq.2 .and. evmad(ibas).gt.emad0)
          if (keepld) modchr = '*'
          write(lgunit(i),440) ibas, evmad(ibas),
     .      v1til(ibas), v2til(ibas), v2til(ibas)-v1til(ibas), modchr
          if (i .eq. 1) e1mad = e1mad + v1til(ibas)**2/evmad(ibas)/4
          if (i .eq. 1) e2mad = e2mad + v2til(ibas)**2/evmad(ibas)/4
          if (keepld) v2til(ibas) = v1til(ibas)
          if (i .eq. 1) e3mad = e3mad + v2til(ibas)**2/evmad(ibas)/4
  440     format(i3, 4f12.5, a1)
   22   continue
        print 442, e1mad, e2mad, e3mad
  442   format(' Mads:',9x,3f12.6,'*')
      endif
      end

