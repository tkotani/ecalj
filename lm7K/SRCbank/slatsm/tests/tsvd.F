      subroutine fmain
C     Driver for routine SVBKSB, which calls routine SVDCMP
      implicit none
      integer nm,np,j,k,l,m,n,ierr
      parameter(nm=5,np=5)
      double precision wmax,wmin
      double precision a(nm,nm),b(nm,nm),u(nm,np),w(np),wk(nm,11)
      double precision v(nm,np),x(nm,np),rv1(np)
      integer ww(10000)
      common /w/ ww

C      data a /
C     .-2.3572015481d0, 1.7851881585d0, 0.4416605027d0, 0.1303528870d0,
C     . 1.7851881585d0,-3.0510328841d0, 0.4391160749d0, 0.8267286507d0,
C     . 0.4416605027d0, 0.4391160749d0,-0.9103618572d0, 0.0295852797d0,
C     . 0.1303528870d0, 0.8267286507d0, 0.0295852797d0,-0.9866668174d0/
C      data b /1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,
C     .        0d0,0d0,1d0,0d0,0d0,0d0,0d0,1d0/

      data a /
     .-2.3572015481d0, 1.7851881585d0, 0.4416605027d0, 0.1303528870d0,
     .  0d0,
     . 1.7851881585d0,-3.0510328841d0, 0.4391160749d0, 0.8267286507d0,
     .  0d0,
     . 0.4416605027d0, 0.4391160749d0,-0.9103618572d0, 0.0295852797d0,
     .  0d0,
     . 0.1303528870d0, 0.8267286507d0, 0.0295852797d0,-0.9866668174d0,
     .  0d0,
     . 0d0,0d0,0d0,0d0/
      data b /1d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,
     .        0d0,0d0,1d0,0d0,0d0,0d0,0d0,0d0,1d0,0d0/


      call wkinit(10000)
      n = 4
      m = 4

      call dcopy(nm*np,a,1,u,1)
C ... this for svd
      call svd(nm,n,n,u,w,.true.,u,.true.,v,ierr,rv1)
C ... or this for symmetric eigenvalues
C      call dtridx(nm,n,a,wk,wk(1,4),wk(1,5),wk(1,2))
C      call imtqlv(n,wk,wk(1,4),wk(1,5),w,wk(1,11),ierr,wk(1,6))
C      call rxx(ierr.ne.0,'imtqlv cannot find all evals') 
C      call tinvit(nm,n,wk(1,1),wk(1,4),wk(1,5),n,w,wk(1,11),v,
C     .  ierr,wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10))
C      call rxx(ierr.ne.0,'tinvit cannot find all evecs') 
C      call dtribx(nm,n,a,wk(1,2),n,v)
C      print *, 'ierr is',ierr
C ... end

C     Find maximum singular value
      wmax=0.0d0
      do  13  k = 1, n
   13 if (dabs(w(k)) .gt. wmax) wmax=dabs(w(k))
C     Define "small"
      wmin=wmax*(1d-10)
C     Zero the "small" singular values
      do  14  k = 1, n
   14 if (dabs(w(k)) .lt. wmin) w(k)=0.0d0
C     Backsubstitute for each right-hand side vector
* ... for svd
      call svbksb(nm,n,n,m,w,u,v,b,x,rv1)
* ... for eigenvector routine
*     call svbksb(nm,n,n,m,w,v,v,b,x,rv1)
*     call dcopy(nm*np,u,1,a,1)
C ... end
      do  18  l = 1, m
        write(*,'(1x,a,i2)') 'vector number ',l
        write(*,*) '    solution vector is:'
        write(*,'(1x,6f18.12)') (x(k,l), k=1,n)
        write(*,*) '    original right-hand side vector:'
        write(*,'(1x,6f18.12)') (b(k,l), k=1,n)
        write(*,*) '    result of (matrix)*(sol''n vector):'
        do  17  k = 1, n
          b(k,l)=0.0d0
          do  16  j = 1, n
   16     b(k,l) = b(k,l) + a(k,j)*x(j,l)
   17   continue
        write(*,'(1x,6f18.12)') (b(k,l), k=1,n)
   18 continue
      end
      subroutine prmx(strn,s,ns,nr,nc)
C- writes matrix into out file (for debugging)
C     implicit none
      integer nr,nc,ns,ifi
      double precision s(ns,nc,2)
      character*(10) fmt, strn*(*), outs*80
      integer i,j,fopna,i1mach
      fmt = '(9f22.17)'
      ifi = fopna('out',29,0)
      call awrit2('%% rows %i cols %i real',' ',80,ifi,nr,nc)
      do  10  i = 1, nr
   10 write(ifi,fmt) (s(i,j,1), j=1,nc)
      write(ifi,*)
C      do  12  i = 1, nr
C   12 write(ifi,fmt) (s(i,j,2), j=1,nc)
      call fclose(ifi)
      outs = 'prm: done writing data '//strn
      call awrit0('%a.  Continue?',outs,-80,-i1mach(2))
      read(*,'(a80)') outs
      if (outs .eq. 'q') call rx0('quit in prmx')
      end
