!> Wigner-Seitz supercell of the k mesh (lattice vectors R and their degeneracy) for the real-space matrix elements of
!> v and W in hwmatK_MPI. Moved here 2026-10-02 from m_maxloc0 (the Wannier functions, removed; tag last-wannier).
module m_wigner_seitz
  implicit none
  public:: wigner_seitz, sortvec2
  private
contains
  subroutine  sortvec2(ndat,vec,dist,idat)
    implicit integer (i-n)
    implicit real*8(a-h,o-z)
    real(8) :: vec(3,ndat),vtmp(3,ndat),dist(ndat)
    integer :: idat(ndat)

    vtmp = vec
    do i = 1,ndat
      dist(i) = dsqrt(sum(vtmp(:,i)**2))
      idat(i) = i
    enddo

    do j = 2,ndat
      d = dist(j)
      do i = j-1,1,-1
        if (dist(i)<=d) goto 999
        dist(i+1) = dist(i)
        idat(i+1) = idat(i)
      enddo
      i = 0
999   dist(i+1) = d
      idat(i+1) = j
    enddo

    do i = 1,ndat
      vec(1:3,i) = vtmp(1:3,idat(i))
    enddo

    do i = 1,ndat-1
      d1 = dsqrt(sum(vec(:,i)**2))
      d2 = dsqrt(sum(vec(:,i+1)**2))
      if (d1 > d2) stop 'sortvec: sorting error!'
      if (abs(d1-dist(i)) > 1.d-4) &
           stop 'sortvec: sorting error in d!'
    enddo

    return
  end subroutine sortvec2
  subroutine wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
    implicit real*8(a-h,o-z)
    implicit integer (i-n)
    integer :: n1,n2,n3,nrws
    real(8) :: alat,plat(3,3)
    integer :: irws(n1*n2*n3*8)
    real(8) :: rws(3,n1*n2*n3*8),drws(n1*n2*n3*8)
    integer(4):: ii0(3,8),isort(8), &
         iwork1(n1*n2*n3*8),iwork2(n1*n2*n3*8)
    real(8):: rr(3,8),dd(8)
    parameter (tol=1.d-6)
    nrws = 0
    do i1=0,n1-1
      do i2=0,n2-1
        do i3=0,n3-1
          n = 0
          do j1=0,1
            do j2=0,1
              do j3=0,1
                n = n+1
                ii0(1,n) = i1 - j1*n1
                ii0(2,n) = i2 - j2*n2
                ii0(3,n) = i3 - j3*n3
              enddo ! j3
            enddo ! j2
          enddo ! j1
          do n=1,8
            rr(1:3,n) =  ( plat(1:3,1)*dble(ii0(1,n)) &
                 +   plat(1:3,2)*dble(ii0(2,n)) &
                 +   plat(1:3,3)*dble(ii0(3,n)) )
          enddo
          call sortvec2(8,rr,dd,isort)
          ndegen = 1
          do n=2,8
            if ((dd(n)-dd(1)) <= tol) ndegen = n
          enddo
          do n=1,ndegen
            nrws = nrws + 1
            rws(1:3,nrws) = rr(1:3,n)
            drws(nrws) = dd(n)
            irws(nrws) = ndegen
          enddo
        enddo ! i3
      enddo ! i2
    enddo ! i1
    call sortvec2(nrws,rws,drws,iwork1)
    iwork2(1:nrws) = irws(1:nrws)
    do n=1,nrws
      irws(n) = iwork2(iwork1(n))
    enddo
  end subroutine wigner_seitz
end module m_wigner_seitz
