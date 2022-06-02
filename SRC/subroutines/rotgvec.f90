subroutine rotgvec(symops, nqibz, &
     ngcmx,ngcn, qbas,ngveccB, &
     ngveccBr)
  !> Rotated ngveccB
  !! G' = R(G), where R denotes rotation and
  !! We determine ngveccBr so that
  !!      matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq)))
  !!      =  matmul(qbas, ngveccBr(1:3, igc, iq)).
  !! See the variable sumchck.
  implicit none
  real(8),intent(in)    :: symops(3,3), qbas(3,3)
  integer(4),intent(in) :: nqibz, ngcmx, ngcn(nqibz)
  integer(4), intent(in) :: ngveccB(3,ngcmx,nqibz)
  integer(4), intent(out) :: ngveccBr(3,ngcmx,nqibz)

  integer(4) :: iq,igc
  real(8) :: qbasinv(3,3),rotnvg(3,3),vec(3),det,sumchk
  integer:: verbose
  logical:: debug=.false.

  if(verbose()>=100) debug= .TRUE. 
  if(debug) write(6,*)' rotgvec: '
  call minv33(qbas,qbasinv)
  !      symops = grp(:,:,irot)
  if(debug) write(6,'(3f12.6)') symops(1,1:3)
  if(debug) write(6,'(3f12.6)') symops(2,1:3)
  if(debug) write(6,'(3f12.6)') symops(3,1:3)

  rotnvg = matmul(qbasinv,matmul(symops,qbas))
  if(debug) write(6,*) 'sum rotnvg=',sum(rotnvg),ngcmx,ngcn(1:nqibz)
  sumchk  =0d0
  !      sumchk2 =0d0
  do iq  = 1, nqibz
     if(debug) write(6,*)' iq=',iq, ' sum ngveccB=', sum(abs(ngveccB(1:3,1:ngcn(iq), iq)))
     do igc = 1, ngcn(iq)
        vec  = matmul( rotnvg, ngveccB(1:3,igc, iq))
        ! vec should be the almost integer and ngveccBr = vec.
        ! But we need this procedure in order to get correct integer value.
        ngveccBr(1:3, igc, iq) = idint( vec + dsign(.5d0,vec))
     enddo
     do igc= 1, ngcn(iq)
        sumchk = sumchk + &
             sum(abs( &
             matmul(qbas, ngveccBr(1:3, igc, iq)) &
             - matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq)))))
     enddo
     ! cccccccccccccc
     !        do igc= 1, ngcn(iq)
     !          write(6,*) " igc=",igc
     !          write(6,'(3f13.5)')
     !     &       matmul(qbas, ngveccBr(1:3, igc, iq))
     !          write(6,'(3f13.5)')
     !     &      matmul(symops, matmul(qbas, ngveccB(1:3, igc, iq)))
     !          sumchk2 = sumchk2 +
     !     &    sum(abs(
     !     &     qrot + matmul(qbas, ngveccBr(1:3, igc, iq)) ))
     !        enddo
     ! cccccccccccccc
     if(debug) write(6,*)" rotgvec: nmin nmax=" &
          ,minval(ngveccBr(1:3, 1:ngcn(iq), iq)) &
          ,maxval(ngveccBr(1:3, 1:ngcn(iq), iq))
  enddo
  if(abs(sumchk)/nqibz/minval(ngcn)>1d-4) then
     write(6,*)" rotgvec: sum chk error sumchk=",sumchk
     call rx( "rotgvec: sum chk error >1d-4")
  endif
end subroutine rotgvec
