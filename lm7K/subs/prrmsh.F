      subroutine prrmsh(strn,rofi,s,ns,nr,nc)
C     implicit none
      integer nr,nc,ns,ifi
      double precision rofi(nr),s(ns,nc)
      character*(20) fmt, strn*(*), outs*80
      integer i,j,fopna,i1mach,i1
C     fmt = '(f10.6,10f15.10)'
      fmt = '(f11.8,10g18.10)'
      ifi = fopna('out',-1,0)
      i1 = 1
      if (s(1,1)  .eq. 0) i1 = 2
      call awrit2('%% rows %i cols %i',' ',80,ifi,nr-i1+1,nc+1)
      do  10  i = i1, nr
        write(ifi,fmt) rofi(i), (s(i,j), j=1,nc)
   10 continue
      call fclose(ifi)
      outs = 'prrmsh: done writing to file out data '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx0('quit in prrmsh')
      end

