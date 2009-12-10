      subroutine iosv(ifi,iter,dmxprm,etot,amgm,ioflg)
C- saves info for convergence (true for write, false for read)
C     implicit none
      double precision dmxprm(20),etot,amgm
      logical ioflg
      integer ifi,iter,lbroy,isw,ineg
      character cflg*1,mflg*1,outs*80
      integer nmix

      if (ifi .le. 0) then
        ineg = 0
        cflg = ' '
        if (dmxprm(11) .lt. 0) cflg = '*'
        lbroy = dmxprm(14)
        nmix = nint(dmxprm(13))
        mflg = 'A'
        if (nmix .eq. 0) mflg = 'L'
        if (lbroy .eq. 1 .and. nmix .gt. 0) then
          nmix = 1
          mflg = 'B'
        endif
        if (lbroy .eq. 2) then
          nmix = 1
          mflg = 'C'
          ineg = isw(dmxprm(16).lt.0)
          if (ineg .ne. 0) dmxprm(16) = -dmxprm(16)
        endif
        write(outs,334) 'SV:',iter,cflg,dabs(dmxprm(11)),
     .    dmxprm(15),dmxprm(12),etot,amgm,mflg
        if (nmix.gt.0) call awrit3('%a%?#n#*##%n:1;3,3g',outs,
     .    len(outs),0,ineg,nmix,dmxprm(16))
        call awrit0('%a',outs,-len(outs),ifi)
        if (ineg .ne. 0) dmxprm(16) = -dmxprm(16)
      else
        call rx('iosv: not ready for this branch')
        ioflg = .false.
C       read(ifi,333,err=10,end=10) dum,iter,dmxprm(11),dmxprm(12),etot,amgm
C       beta = dmxprm(11)
        ioflg = .true.
C  10   continue
      endif
C  333 format(1x,a4,i3,1pe10.3,0p,f7.4,1pe10.3,0p,f17.8,f9.6:
C     .       ' |',10f8.4)
  334 format(1x,a3,i4,a,1pd9.3,0p,f7.4,1pe10.3,0p,f17.8,f9.6,1x,a,1x,
     .  10f8.4)
      end

