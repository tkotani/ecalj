      subroutine fmain
      implicit none
      integer i,j,n,ndim,ier,nevmx,nev,itime,i1mach
      logical lov,lnv,lx
      parameter (ndim=4)
      complex*16 a(ndim,ndim),z(ndim,ndim),wk(32),c(ndim,ndim),
     .  s(ndim,ndim),zz(ndim,ndim),scale
      double precision e(ndim),eevmx
      
      eevmx=100
      nevmx=100
      lov = .true.
      lx  = .true.
      lnv = .false.
      itime = i1mach(2)
      call sa(a)
      print *, 'h: (see imsl)'
      print 333, a

C --- Check cplx2r ---
C      call zprm('(5f12.6)',a,ndim,ndim)
C      call cplx2r(ndim**2,0,a,z)
C      call zprm2('(5f12.6)',a,ndim,ndim)
C      call cplx2r(ndim**2,1,a,z)
C      call zprm('(5f12.6)',a,ndim,ndim)
C      stop

      call ss(s)
      call zhev(ndim,a,s,lov,lx,nevmx,eevmx,nev,wk,lnv,itime,e,z)

      print *, 'eigenvalues: (see imsl)'
      print 333, e
      print *, 'eigenvectors:'
      print 333, z
  333 format(4f14.5)

      call sa(a)
      n = ndim
      call zmpy(a,2*n,2,1,z,2*n,2,1,c,2*n,2,1,n,n,n)
      call ss(s)
      call zmpy(s,2*n,2,1,z,2*n,2,1,zz,2*n,2,1,n,n,n)

      print *, 'H Z  -  E O Z'
      do  10  i = 1, ndim
      do  10  j = 1, ndim
   10 print 333, c(i,j) - e(j)*zz(i,j)

      end

      subroutine sa(a)
      complex*16 a(4,4)
      a(1,1) = 3+1
      a(2,2) = 3+1
      a(3,3) = 1+1
      a(4,4) = 1+1
      a(2,1) = 1
      a(3,1) = 0
      a(4,1) = (0d0,-2d0)
      a(3,2) = (0d0,+2d0)
      a(4,2) = 0
      a(4,3) = 1
      a(1,2) = 1
      a(1,3) = 0
      a(1,4) = (0d0,+2d0)
      a(2,3) = (0d0,-2d0)
      a(2,4) = 0
      a(3,4) = 1

      end
      subroutine ss(s)
      complex*16 s(4,4)

      call zinit(s,16)
      s(1,1) = 1
      s(1,3) = (.1d0,.2d0)
      s(3,1) = (.1d0,-.2d0)
      s(2,2) = 2
      s(3,3) = 3
      s(4,4) = 4
      end

      subroutine zprm(fmt,s,nr,nc)
      double precision s(2,nr,nc)
C#ifdef APOLLO_BUG
      character*(20) fmt
C#elseC
C      character*(*) fmt
C#endif
      print *, nr, nc
      do  10  i = 1, nr
   10 print fmt, (s(1,i,j), j=1,nc)
      do  20  i = 1, nr
   20 print fmt, (s(2,i,j), j=1,nc)
      end
      subroutine zprm2(fmt,s,nr,nc)
      double precision s(nr,nc,2)
C#ifdef APOLLO_BUG
      character*(20) fmt
C#elseC
C      character*(*) fmt
C#endif
      print *, nr, nc
      do  10  i = 1, nr
   10 print fmt, (s(i,j,1), j=1,nc)
      do  20  i = 1, nr
   20 print fmt, (s(i,j,2), j=1,nc)
      end
