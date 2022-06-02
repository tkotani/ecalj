module m_readqgcou
  !! this is somehow duplicated with reqdqg.F---> I hope they are unified... jun2012takao
  !! ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc,ikp) is added on jun2012takao
  use NaNum,only: NaN
  implicit none
  integer,protected:: nqnum=NaN, nqbz=NaN, imxc=NaN, ngcmx=NaN, nqbz_=NaN
  integer,allocatable,protected:: ngc(:)
  real(8),allocatable,protected:: qtt_(:,:)
  integer,allocatable,protected:: ngvecc(:,:,:), ngveccrev(:,:,:,:),itimermap(:,:)
contains

  !!--- Readin QGcou ---
  subroutine readqgcou()
    integer:: ifiqg,nqi,ikp,nnnn! igc,ikpm
    real(8):: QpGcut_cou,qmm(3)
    open(newunit=ifiqg,file='QGcou',form='unformatted')
    read(ifiqg) nqnum , ngcmx, QpGcut_cou, nqbz_, nqi,imxc
    write(6,"('read QGcou file',2i6,f8.3,10i6)") nqnum , ngcmx, QpGcut_cou, nqbz_,nqi,imxc
    allocate( qtt_(3,nqnum),ngc(nqnum) )
    allocate( ngvecc(3,ngcmx,nqnum))
    allocate( ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc,nqnum) )
    allocate( itimermap(ngcmx,nqnum))
    do ikp = 1,nqnum
       read (ifiqg) qtt_(:,ikp), ngc(ikp)
       read (ifiqg) ngvecc(1:3, 1:ngc(ikp),ikp),ngveccrev(-imxc:imxc,-imxc:imxc,-imxc:imxc,ikp)
       !! time reversal mapping. only needed for time reversal case.
    enddo
    ! takao is trying to a mapping for time-reversal. But not yet. it can be complicated
    ! when we used current q-point mesh where -q do not exist in the list of q vector.
    !      do ikp=1,nqnum
    !        call iqindx2(-qtt_(:,ikp),ginv,qtt_,nqnum, ikpm,qmm) !qinv is true q
    !        do igc=1,ngc(ikp)
!!!  ikp ---> ikpm
!!!  q+G (igc,ikp) is mapped to  qmm+G (itimermap,ikpm)
!!!
    ! xxxxxxxxxxxxxxxxxxx qmm= -q + delta G xxxxxxxx
    ! xxxxxxxx this should be taken into account xxxxxxxx
    ! xxxxxxxxx in addition, we have phase problem when there is delta G vector for inversion xxxxxx.
    !          itimermap(igc,ikp)= ngveccrev(-ngvecc(1,igc,ikp),-ngvecc(2,igc,ikp),-ngvecc(3,igc,ikp),ikpm)
    !          !ikp is mapped to itimermap
    !        enddo
    !     enddo
    close(ifiqg)
  end subroutine readqgcou
end module m_readqgcou
