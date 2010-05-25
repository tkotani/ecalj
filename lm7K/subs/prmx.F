      subroutine prm(sfmt,msga,s,ns,nr,nc)
C- writes matrix into file with message (for debugging)
C     implicit none
      integer nr,nc,ns,ifi,msga
      double precision s(ns,nc,2)
      character*(10) fmt, sfmt*(*), outs*80
      integer i,j,fopna,i1mach
      fmt = '(9f22.17)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
      do  10  i = 1, nr
        write(ifi,fmt) (s(i,j,1), j=1,nc)
   10 continue
      write(ifi,*)
C      do  12  i = 1, nr
C   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
      call fclose(ifi)
      call awrit1(sfmt,outs,80,0,msga)
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx0('quit in prmx')
      end
      subroutine prmx(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
C     implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2)
      character*(14) fmt, fmt0, strn*(*), outs*80
      integer i,j,fopna,i1mach
      save fmt
      data fmt /'(9f20.15)'/
C      fmt = '(1p9e20.10)'
      ifi = fopna('out',29,0)
      rewind ifi
#if NOAWRITE
      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc
#else
      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
#endif
      do  10  i = 1, nr
        write(ifi,fmt) (s(i,j,1), j=1,nc)
   10 continue
      write(ifi,*)
C      do  12  i = 1, nr
C   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
      call fclose(ifi)

#if NOAWRITE
      outs = ' prm: wrote '//strn//' continue?'
      print *, outs
#else
      outs = ' prm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
#endif
      read(*,'(a80)') outs

#if NOAWRITE
      if (outs .eq. 'q') stop 'quit in prmx'
#else
      if (outs .eq. 'q') call rx0('quit in prmx')
#endif
      return

      entry prmx0(fmt0)
      fmt = fmt0
      end
      subroutine prmr(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
C     implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2),xx
      character*(10) fmt, strn*(*)
      integer i,j,fopna
      fmt = '(9f15.10)'
      ifi = fopna('out',29,0)
      write(ifi,*) nr, nc
      do  10  i = 1, nr
        write(ifi,fmt) (s(i,j,1), j=1,nc)
   10 continue
      write(ifi,*)
      call fclose(ifi)
      print *, 'prm: pausing after writing data ',strn
      read(5,*) xx !pause
      end
      subroutine yprm(strn,icast,s,ofi,ns,nr,nc)
C ofi used only for kcplx=0
C ns,nr,nc are formal dimensions, not real ones
C     implicit none
      integer icast,ofi,ns,nr,nc,ifi
*      double precision s(ns,nsc,2)
      integer s(ns,nc)
      character*(20) fmt, fmt0, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      save fmt
      data fmt /'(9f15.10)'/
*     fmt = '(9f15.10)'
*     fmt = '(%9;10,10d)'
      ifi = fopna('out',29,0)

      if (icast .ne. 0) then
        call ywrm(0,' ',icast,ifi,fmt,s,ofi,ns,nr,nc)
      else
        call awrit2('%% rows %i cols %i',' ',80,ifi,nr,nc)
        do  10  i = 1, nr
          write(ifi,'(22i7)') (s(i,j), j=1,nc)
   10   continue
      endif

      call fclose(ifi)
      outs = ' prm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx0('quit in prmx')
      return
      entry yprm0(fmt0)
      fmt = fmt0
      end
      subroutine zprm(strn,icast,s,ns,nr,nc)
C     implicit none
      integer icast,nr,nc,ns,ifi
      double precision s(2,ns,nc)
      character*(10) fmt, outs*80, strn*(*)
      integer i,j,fopna,i1mach
      fmt = '(9f15.10)'
      fmt = '(9f18.11)'
C     fmt = '(5f20.15)'
      outs = ' '
      if (icast .eq. 1)  outs = ' real'
      if (icast .eq. 11) outs = ' symm'
      if (icast .eq. 2)  outs = ' complex'
      if (icast .eq. 12) outs = ' herm'
      ifi = fopna('out',29,0)
      rewind ifi
#if NOAWRITE
      write(ifi,'(''% rows'',i5,'' cols'',i5,a)') nr,nc,outs(1:10)
#else
      call awrit2('%% rows %i cols %i'//outs,' ',80,ifi,nr,nc)
#endif
      do  10  i = 1, nr
        write(ifi,fmt) (s(1,i,j), j=1,nc)
   10 continue
      if (mod(icast,10) .gt. 1) then
        write(ifi,*)
        do  20  i = 1, nr
          write(ifi,fmt) (s(2,i,j), j=1,nc)
   20   continue
      endif
      call fclose(ifi)
#if NOAWRITE
      outs = ' prm: wrote '//strn//' continue?'
      print *, outs
#else
      outs = ' zprm: wrote '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
#endif
      read(*,'(a80)') outs

#if NOAWRITE
      if (outs .eq. 'q') stop 'quit in zprm'
#else
      if (outs .eq. 'q') call rx0('quit in zprm')
#endif
      end
      subroutine zprm3(sfmt,msga,s,n1,n2,n3)
C     implicit none
      character *(*) sfmt
      integer n1,n2,n3,msga
      double precision s(2,n1,n2,n3)
      character*(20) fmt
      character*80 outs
      integer i,j,k,i1mach,ifi,fopna
      fmt = '(9f15.10)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i complex',' ',80,ifi,n1*n2,n3)
      do  10  i = 1, n1
        do  10  j = 1, n2
          write(ifi,fmt) (s(1,i,j,k), k=1,n3)
   10 continue
      do  20  i = 1, n1
        do  20  j = 1, n2
          write(ifi,fmt) (s(2,i,j,k), k=1,n3)
   20 continue
      call fclose(ifi)
      call awrit1(sfmt,outs,80,0,msga)
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx0('quit in zprm3')
      end
      subroutine prm3(sfmt,msga,s,n1,n2,n3)
C     implicit none
      character *(*) sfmt
      integer n1,n2,n3,msga
      double precision s(n1,n2,n3)
      character*(20) fmt
      character*80 outs
      integer i,j,k,i1mach,ifi,fopna
      fmt = '(9f15.10)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i real',' ',80,ifi,n1*n2,n3)
      do  10  i = 1, n1
        do  10  j = 1, n2
          write(ifi,fmt) (s(i,j,k), k=1,n3)
   10 continue
      call fclose(ifi)
      call awrit1(sfmt,outs,80,0,msga)
      call awrit0('%a.  Continue?',outs,80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx0('quit in zprm3')
      end
      subroutine rdm3(s,n1,n2,n3)
C     implicit none
      integer n1,n2,n3
      double precision s(n1,n2,n3)
      integer rdm,fopna,ifi,i,j,k
      double precision xx

      ifi = fopna('in',29,0)
      if (rdm(ifi,0,0,' ',xx,n1*n2,n3) .ne. 1) call
     .rx('rdm3:  file mismatch')
      do  10  i = 1, n1
        do  10  j = 1, n2
          read(ifi,*) (s(i,j,k), k=1,n3)
   10 continue

      call fclose(ifi)

      end

