      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
c
      integer i,j,k,l,m,n,nm
      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
      double precision h,s,si
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  htridi.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  htridi  in their
c          full lower triangles except for the diagonal of ar.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c     Aug 1990 MvS altered into daxpy-style loops
c     ------------------------------------------------------------------
c
      if (m .eq. 0) return
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........
      do  k = 1, n
        do j = 1, m
          zi(k,j) = -zr(k,j) * tau(2,k)
          zr(k,j) = zr(k,j) * tau(1,k)
       enddo
      enddo

      if (n .eq. 1) return
C     .......... recover and apply the householder matrices ..........
      do  140  i = 2, n
        l = i - 1
        h = ai(i,i)
        if (h .eq. 0.0d0) go to 140
c        do 100 j = 1, m
c          tau(1,j) = 0.0d0
c          tau(2,j) = 0.0d0
c     100   continue
        tau=0d0
        do  110  k = 1, l
c          call yyaxpy(m,ar(i,k),ai(i,k),zr(k,1),zi(k,1),nm,
c     .          tau,tau(2,1),2,.true.)
          tau(1,:)= ar(i,k)*zr(k,:) - ai(i,k)*zi(k,:) + tau(1,:)
          tau(2,:)= ar(i,k)*zi(k,:) + ai(i,k)*zr(k,:) + tau(2,:) 
  110   continue
C     ... double divisions avoid possible underflow ...
        tau=tau/h**2
c        call dscal(m,1/h,tau,2)
c        call dscal(m,1/h,tau,2)
c        call dscal(m,1/h,tau(2,1),2)
c        call dscal(m,1/h,tau(2,1),2)
        do  120  j = 1, m
c          call yyxcpy(l,-tau(1,j),-tau(2,j),ar(i,1),ai(i,1),nm,
c     .    zr(1,j),zi(1,j),1,.true.)
           zr(:,j)= -tau(1,j)*ar(i,:) + tau(2,j)*ai(i,:) + zr(:,j)
           zi(:,j)= +tau(1,j)*ai(i,:) - tau(2,j)*ar(i,:) + zi(:,j)
 120    continue
  140 continue
      end

