subroutine writeboltztrap(evlall,eferm) !write input file for boltztrap !test by gomi at year2020 around
  use m_lgunit,only:stdo
  use m_lmfinit,only: nlmax,nsp,nbas,nlmax,nspc,qbg=>zbak,alat=>lat_alat
  use m_suham,only: ndhamx=>ham_ndhamx,ndham=>ham_ndham
  use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize,procid
  use m_qplist,only: nkp,xdatt,qplist
  use m_suham,only: nspx=>ham_nspx
  use m_mkpot,only:  qval
  use m_ext,only: sname
  use m_hamindex, only: ngrp,symops !,norbmto,ibastab,ltab,ktab,offl, symops_af
  use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, plat=>lat_plat,pos=>rv_a_opos
  real(8):: evlall(:,:,:)
!  use m_bandcal,only: evlall
  character strn*120,strn2*120
  integer:: iqread,iqindex,job,ist,ip,ni,ix,ifi,jsp,ncount,iq,iband,i,j,nbandx,ig
  real(8):: eferm,qvec(3),symxx(3,3)
  character(10):: i2char
  !   qbg = !homogenious background charge.  qcore + qval-qbg = \sum_i Z_i
  iqread=0
  iqindex = index(strn(12:),'nb=')+2
  nbandx = ndhamx
  if(iqindex/=2) then
     read(strn(12+iqindex:),*) nbandx
     nbandx=min(nbandx,ndhamx)
  endif
  !      open(newunit=ifi,file='efermi.lmf') !readin fermi energy from efermi.lmf
  !      read(ifi,*)  eferm
  !      close(ifi)
  open(newunit=ifi,file=trim(sname)//'.intrans_template.boltztrap')
  write(ifi,"(a)")'GENE          # format '
  write(ifi,"(a)")'0 0 0 0       # iskip (not presently used) idebug setgap shiftgap'
  write(ifi,"(f20.16,a,f9.4,a)") eferm,' 0.0005 0.4 ', qval-qbg, &
       '       # efermi.lmf (Ry), energygrid, energy span around Fermilevel, number of electrons'
  write(ifi,"(a)")'CALC          # CALC (calculate expansion coeff), NOCALC read from file'
  write(ifi,"(a)")'5             # lpfac, number of latt-points per k-point'
  write(ifi,"(a)")'BOLTZ         # run mode (only BOLTZ is supported)'
  write(ifi,"(a)")'0.15          # (efcut) energy range of chemical potential'
  write(ifi,"(a)")'800.0 50.0  # Tmax, temperature grid'
  write(ifi,"(a)")'-1.0  # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)'
  write(ifi,"(a)")'TETRA'
  close(ifi)
  open(newunit=ifi,file=trim(sname)//'.struct.boltztrap')
  write(ifi,"(a)") trim(sname)
  write(ifi,"(3d16.8, ' # plat1  ')") plat(:,1)*alat
  write(ifi,"(3d16.8, ' # plat2  ')") plat(:,2)*alat
  write(ifi,"(3d16.8, ' # plat3  ')") plat(:,3)*alat
  write(ifi,*) ngrp
  do ig=1,ngrp
     symxx = matmul(transpose(qlat), matmul(symops(:,:,ig),plat))
     if(abs(sum( nint(symxx(:,:))-symxx(:,:) )) >1d-6) then !sanity check
        call rx('boltztrap: rotation-matrix elements are not integers --- probably strange')
     endif
     write(ifi,"(9i3)") ((nint(symxx(i,j)),j=1,3),i=1,3)
  enddo
  close(ifi)
  do jsp=1,nspx             !=isp.  nspx=1 for so=1
     open(newunit=ifi,file=trim(sname)//'.energy.isp'//trim(i2char(jsp))//'.boltztrap')
     write(ifi,"(a)") trim(sname)
     write(ifi,"(i10)") nkp
     do iq=1,nkp
        !! true q vector = 2pi/alat *qplist(:,iq) in cartesian
        qvec= matmul(transpose(plat),qplist(:,iq)) !qvec in qlat unit
        ncount=0
        do iband=1,nbandx
           if(evlall(iband,jsp,iq)>1d98) cycle
           ncount=ncount+1
        enddo
        write(ifi,"(3f15.8,i10)") qvec,ncount
        do iband=1,nbandx
           if(evlall(iband,jsp,iq)>1d98) cycle !evlall=1d99 is dummy
           write(ifi,"(d23.16)") evlall(iband,jsp,iq)
        enddo
     enddo
     close(ifi)
  enddo
end subroutine writeboltztrap
