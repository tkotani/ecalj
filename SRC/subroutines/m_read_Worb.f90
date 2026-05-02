!--------------------------------------
!! This module is for readling informations of <Worb> section in GWinput.
!! Then this is transformed to the form for the calculation of MLWF.
!!
!! nclass_mlwf    : the number of atoms for Wannier
!! classname_mlwf : the name of atoms
!! nbasclass_mlwf : the number of Wannier orbitals for each line.
!! cbas_mlwf : each Wannier orbital
!!             (e.g) s = 1, py = 2,... (see texts in the script job_band)
!-------------------------------------
module m_read_Worb  !! all output
  integer,protected:: nwf, nclass_mlwf
  integer,protected,allocatable:: cbas_mlwf(:,:),nbasclass_mlwf(:)
  character(20),protected,allocatable:: classname_mlwf(:)
  integer,protected,allocatable:: iclassin(:)
  integer,parameter,private::maxdat=1024
  ! t is convenient if you use this type of number.
  integer(4),protected,allocatable::  iphi(:,:),iphidot(:,:), nphi(:)
  integer(4),protected:: nphix
  integer(4):: natom
contains
  subroutine s_read_Worb()
    use m_keyvalue,only: getkeyvalue
    use m_GWinput, only: gwinput_init, gwinput_loaded, &
                         tg_n_worb => n_worb, tg_worb_iatom => worb_iatom, &
                         tg_worb_label => worb_label, tg_worb_lm => worb_lm, &
                         tg_worb_nlm => worb_nlm
    implicit none
    integer(4):: ret
    integer:: i, ifile_handle,il,ix,iclass,iline,nbasclassMax,nline
    character*256:: a,aaa
    integer::ifmloc
    integer,allocatable:: cbastemp(:,:)
    call gwinput_init()
    if (gwinput_loaded) then
       nclass_mlwf = tg_n_worb
       if (nclass_mlwf == 0) call rx('s_read_Worb: empty Worb in GWinput.toml')
       allocate(iclassin(nclass_mlwf), nbasclass_mlwf(nclass_mlwf), classname_mlwf(nclass_mlwf))
       nbasclassMax = 0
       do iclass = 1, nclass_mlwf
          iclassin(iclass)       = tg_worb_iatom(iclass)
          classname_mlwf(iclass) = tg_worb_label(iclass)
          nbasclass_mlwf(iclass) = tg_worb_nlm(iclass)
          if (tg_worb_nlm(iclass) > nbasclassMax) nbasclassMax = tg_worb_nlm(iclass)
       enddo
       allocate(cbas_mlwf(nbasclassMax, nclass_mlwf))
       cbas_mlwf = -999
       do iclass = 1, nclass_mlwf
          cbas_mlwf(1:nbasclass_mlwf(iclass), iclass) = &
               tg_worb_lm(1:nbasclass_mlwf(iclass), iclass)
       enddo
    else
       call getkeyvalue("GWinput","<Worb>",unit=ifmloc,status=ret)
       iline = 0
       nline = 0
       nclass_mlwf = 0
       do
          iline = iline + 1
          read(ifmloc,"(a)") aaa
          if (aaa(1:1) == '!') cycle
          if(aaa(1:7) == "</Worb>") then
             exit
          else
             nclass_mlwf = nclass_mlwf+1
          end if
       end do
       close(ifmloc)
       nline = iline-1
       allocate(iclassin(nclass_mlwf),cbastemp(maxdat,nclass_mlwf), &
            nbasclass_mlwf(nclass_mlwf),classname_mlwf(nclass_mlwf))
       call getkeyvalue("GWinput","<Worb>",unit=ifmloc,status=ret)
       cbastemp=-999
       iclass = 0
       do 1001 iline=1,nline
          read(ifmloc,"(a)") aaa
          if (aaa(1:1) == '!') then
             read(aaa,*)
             cycle
          end if
          iclass = iclass + 1
          read(aaa,*,end=1201) iclassin(iclass),a,(cbastemp(i,iclass),i=1,maxdat)
1201      continue
          classname_mlwf(iclass) = trim(a)
          do i=1,maxdat
             if(cbastemp(i,iclass)==-999) then
                nbasclass_mlwf(iclass)=i-1
                exit
             endif
          enddo
1001   enddo
       nbasclassMax = maxval(nbasclass_mlwf(1:nclass_mlwf))
       allocate(cbas_mlwf(nbasclassMax,nclass_mlwf))
       cbas_mlwf = cbastemp(1:nbasclassMax,1:nclass_mlwf)
       deallocate(cbastemp)
       close(ifmloc)
    endif
    nwf = 0
    do iclass=1,nclass_mlwf
       nwf = nwf + nbasclass_mlwf(iclass)
    end do
  end subroutine s_read_Worb
  subroutine s_cal_Worb()
    use m_ll,only:ll
    use m_keyvalue,only: getkeyvalue
    use m_GWinput, only: gwinput_init, gwinput_loaded, &
                         tg_pb_nlx => pb_nlx, tg_pb_n_nlx => pb_n_nlx
    use m_struct_from_lmf, only: natom !    use m_HamPMT,  only: natom=>nbas !NOT 2023-8-4
    implicit none
    integer:: iclass, iclass2, iphidot_plus, ifmloc, iphi_tmp
    integer :: i, j, l_number, correction
    integer :: tmp_atom, tmp_l, iatom, il, iwf, ret,ixatom
    integer :: nnvv,ix,ioffset(0:10,natom),ioffadd,mm
    integer,allocatable:: l_numbermx(:)
    !! Read index for cphi from GWinput: it should be essentially the same as @MNLA_CPHI
    call gwinput_init()
    if (gwinput_loaded) then
       ! pb_nlx(1:4, 1:pb_n_nlx) = [iatom, l, nnvv, nnc]
       ioffadd = 0
       do ix = 1, tg_pb_n_nlx
          ixatom = tg_pb_nlx(1, ix)
          il     = tg_pb_nlx(2, ix)
          nnvv   = tg_pb_nlx(3, ix)
          ioffset(il, ixatom) = ioffadd
          ioffadd = ioffadd + nnvv * (2*il + 1)
       enddo
    else
       call getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifmloc,status=ret)
       read(ifmloc,*)
       read(ifmloc,*)
       read(ifmloc,*)
       read(ifmloc,*)
       read(ifmloc,*)
       ioffadd = 0
       do
          read(ifmloc,*,err=888)ixatom,il,nnvv
          ioffset(il,ixatom) = ioffadd
          ioffadd = ioffadd + nnvv * (2*il+1)
       end do
888    continue
       close(ifmloc)
    endif
    !! real harmonics case
    allocate (nphi(nwf)) ! number of radial waves for each iwf.
    nphi  = 1  ! We use a simple setting.
    nphix = 1
    allocate (iphi(nphix,nwf),iphidot(nphix,nwf))
    iphi = 0
    iphidot = 0
    iwf = 1
    do iclass=1,nclass_mlwf !atom
       do i=1,nbasclass_mlwf(iclass) !wannier index for atom from GWinput Worb.
          l_number = ll(cbas_mlwf(i,iclass)) !l for Wannier
          mm = cbas_mlwf(i,iclass) - l_number**2
          iphi  (nphix,iwf) = ioffset(l_number,iclassin(iclass)) + mm !bugfix 2022-8-7 (based on Suzuki report.)
          iphidot(nphix,iwf) = iphi(nphix,iwf) + l_number*2+1
          iwf = iwf +1
       enddo
    enddo
  end subroutine s_cal_Worb
end module m_read_Worb

