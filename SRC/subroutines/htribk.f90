subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)

  integer :: i,j,k,l,m,n,nm
  double precision :: ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
  double precision :: h,s,si

  !     this subroutine is a translation of a complex analogue of
  !     the algol procedure trbak1, num. math. 11, 181-195(1968)
  !     by martin, reinsch, and wilkinson.
  !     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

  !     this subroutine forms the eigenvectors of a complex hermitian
  !     matrix by back transforming those of the corresponding
  !     real symmetric tridiagonal matrix determined by  htridi.

  !     on input

  !        nm must be set to the row dimension of two-dimensional
  !          array parameters as declared in the calling program
  !          dimension statement.

  !        n is the order of the matrix.

  !        ar and ai contain information about the unitary trans-
  !          formations used in the reduction by  htridi  in their
  !          full lower triangles except for the diagonal of ar.

  !        tau contains further information about the transformations.

  !        m is the number of eigenvectors to be back transformed.

  !        zr contains the eigenvectors to be back transformed
  !          in its first m columns.

  !     on output

  !        zr and zi contain the real and imaginary parts,
  !          respectively, of the transformed eigenvectors
  !          in their first m columns.

  !     note that the last component of each returned vector
  !     is real and that vector euclidean norms are preserved.

  !     questions and comments should be directed to burton s. garbow,
  !     mathematics and computer science div, argonne national laboratory

  !     this version dated august 1983.
  !     Aug 1990 MvS altered into daxpy-style loops
  !     ------------------------------------------------------------------

  if (m == 0) return
  !     .......... transform the eigenvectors of the real symmetric
  !                tridiagonal matrix to those of the hermitian
  !                tridiagonal matrix. ..........
  do  k = 1, n
     do j = 1, m
        zi(k,j) = -zr(k,j) * tau(2,k)
        zr(k,j) = zr(k,j) * tau(1,k)
     enddo
  enddo

  if (n == 1) return
  !     .......... recover and apply the householder matrices ..........
  do  140  i = 2, n
     l = i - 1
     h = ai(i,i)
     if (h == 0.0d0) go to 140
     !        do 100 j = 1, m
     !          tau(1,j) = 0.0d0
     !          tau(2,j) = 0.0d0
     !     100   continue
     tau=0d0
     do  110  k = 1, l
        !          call yyaxpy(m,ar(i,k),ai(i,k),zr(k,1),zi(k,1),nm,
        !     .          tau,tau(2,1),2,.true.)
        tau(1,:)= ar(i,k)*zr(k,:) - ai(i,k)*zi(k,:) + tau(1,:)
        tau(2,:)= ar(i,k)*zi(k,:) + ai(i,k)*zr(k,:) + tau(2,:)
110  enddo
     !     ... double divisions avoid possible underflow ...
     tau=tau/h**2
     !        call dscal(m,1/h,tau,2)
     !        call dscal(m,1/h,tau,2)
     !        call dscal(m,1/h,tau(2,1),2)
     !        call dscal(m,1/h,tau(2,1),2)
     do  120  j = 1, m
        !          call yyxcpy(l,-tau(1,j),-tau(2,j),ar(i,1),ai(i,1),nm,
        !     .    zr(1,j),zi(1,j),1,.true.)
        zr(:,j)= -tau(1,j)*ar(i,:) + tau(2,j)*ai(i,:) + zr(:,j)
        zi(:,j)= +tau(1,j)*ai(i,:) - tau(2,j)*ar(i,:) + zi(:,j)
120  enddo
140 enddo
end subroutine htribk

