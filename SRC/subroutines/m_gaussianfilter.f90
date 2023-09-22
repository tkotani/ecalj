!>Gaussian filtering for rcxq, imaginary part of x0
module m_GaussianFilter
  use m_freq,only: frhis,nwhis,npm
  implicit none
  public:: GaussianFilter
  private
  logical:: eginit=.true.
  real(8),allocatable,save:: gfmat(:,:)
  complex(8),allocatable:: rcxqin(:)
contains
  subroutine GaussianFilter(rcxq,nmbas1,nmbas2, egauss,iprint)
    intent(in):: nmbas1,nmbas2,egauss,iprint
    complex(8):: rcxq(nmbas1,nmbas2,nwhis,npm)
    logical:: iprint
    integer:: nmbas1,nmbas2,ipm,imbas1,imbas2
    real(8):: egauss
    if(eginit) then
       if(iprint) write(6,'("GaussianFilterX0= ",d13.6)') egauss
       allocate(gfmat(nwhis,nwhis))
       gfmat=gaussianfilterhis(egauss,frhis,nwhis)
       allocate(rcxqin(1:nwhis))
       eginit=.false.
    endif
    do ipm=1,npm
       do imbas1=1,nmbas1
          do imbas2=1,nmbas2
             rcxqin = rcxq(imbas1,imbas2,1:nwhis,ipm)
             rcxq(imbas1,imbas2,1:nwhis,ipm) = matmul(gfmat,rcxqin)
          enddo
       enddo
    enddo
    write(6,"(' End of Gaussian Filter egauss=',f9.4)") egauss
  end subroutine GaussianFilter
  pure function gaussianfilterhis(egauss, frhis,nwhis) result(gfmat)
    implicit none
    integer,intent(in):: nwhis
    real(8),intent(in):: egauss,frhis(nwhis+1)
    real(8):: gfmat(nwhis,nwhis)
    real(8),allocatable:: frc(:),gfm(:)
    real(8):: ggg
    integer:: i,j
    allocate(frc(nwhis),gfm(nwhis))
    do i=1,nwhis
       frc(i)=(frhis(i)+frhis(i+1))/2d0
    enddo
    do i=1,nwhis
       do j=1,nwhis
          gfm(j)= exp( -(frc(i)-frc(j))**2/(2d0*egauss))
       enddo
       ggg = sum(gfm(:))
       do j=1,nwhis
          gfmat(j,i)= gfm(j)/ggg
       enddo
    enddo
    deallocate(frc,gfm)
  end function gaussianfilterhis
end module m_GaussianFilter
