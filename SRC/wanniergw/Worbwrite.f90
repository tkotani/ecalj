  ! ccc identify d-orb:
  use m_read_Worb,only: s_read_Worb, s_cal_Worb, &
       nwf, nclass_mlwf, cbas_mlwf, nbasclass_mlwf, &
       classname_mlwf, iclassin, &
       iphi, iphidot, nphi, nphix

  integer :: iclass2

  ! input parameters specific to MAXLOC
  call s_read_Worb()

  do iclass2=1,nclass_mlwf
     write(*,*)'output:',iclassin(iclass2), nwf &
          ,trim(classname_mlwf(iclass2)),cbas_mlwf(1:nbasclass_mlwf(iclass2),iclass2)
  enddo

  !      call s_cal_Worb()



END PROGRAM
