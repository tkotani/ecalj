  !! Determine Fermi energy ef for given valn (legas case), or corresponding charge given by z and konf.==
  !!    When esmr is negative, esmr is geven automatically by efsimplef.
  !      write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
  legas=.false. ! if legas=T, homogenius electron gas test case.
  call efsimplef2a(nspin,wibz,qibz,ginv, &
       nband,nqibz &
       ,konf,z,nl,natom,iclass,nclass &
       ,valn, legas, esmref,    ! & !! valn is input for legas=T, output otherwise.
  qbz,nqbz                 ! &  index_qbz, n_index_qbz,
  ,efnew)
  if(ixc/=3) ef = efnew
  eftrue = efnew
