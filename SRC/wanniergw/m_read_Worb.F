c--------------------------------------
!! This module is for readling informations of <Worb> section in GWinput.
!! Then this is transformed to the form for the calculation of MLWF.
!!
!! nclass_mlwf    : the number of atoms for Wannier
!! classname_mlwf : the name of atoms
!! nbasclass_mlwf : the number of Wannier orbitals for each line.
!! cbas_mlwf : each Wannier orbital 
!!             (e.g) s = 1, py = 2,... (see texts in the script job_band)
c-------------------------------------
      module m_read_Worb
!! all output
      integer,protected:: nwf, nclass_mlwf
  ! cbas_mlwf(nbasclassMax,nclass_mlwf) !We may need nbasclassMax to pass array to subroutine.
      integer,protected,allocatable:: cbas_mlwf(:,:),nbasclass_mlwf(:)
      character(20),protected,allocatable:: classname_mlwf(:)
      integer,protected,allocatable:: iclassin(:)
      integer,parameter,private::maxdat=1024 
                   !it is convenient if you use this type of number.
      integer(4),protected,allocatable::  iphi(:,:),iphidot(:,:), nphi(:)
      integer(4),protected:: nphix
      integer(4):: natom

      contains

      subroutine s_read_Worb()
      use m_keyvalue,only: getkeyvalue
      implicit none
      integer(4):: ret
      integer:: i, ifile_handle,il,ix,iclass,iline,nbasclassMax,nline
      character*256:: a,aaa
      integer::ifmloc
      integer,allocatable:: cbastemp(:,:)
      
      call getkeyvalue("GWinput","<Worb>",unit=ifmloc,status=ret) 
c      write(6,*)' ifmloc ret=',ifmloc,ret

      iline = 0
      nline = 0
      nclass_mlwf = 0
      do 
        iline = iline + 1
        read(ifmloc,"(a)") aaa
c        print *, aaa
        if (aaa(1:1) == '!') cycle
        if(aaa(1:7) == "</Worb>") then
          exit
        else 
          nclass_mlwf = nclass_mlwf+1
        end if
      end do
      close(ifmloc)

      nline = iline-1
      write(6,'("nline nclass mlwf=",3i5)') nline,nclass_mlwf
      allocate(iclassin(nclass_mlwf),cbastemp(maxdat,nclass_mlwf),
     &        nbasclass_mlwf(nclass_mlwf),classname_mlwf(nclass_mlwf))

      call getkeyvalue("GWinput","<Worb>",unit=ifmloc,status=ret) 
      cbastemp=-999
      iclass = 0
      do 1001 iline=1,nline
        read(ifmloc,"(a)") aaa
        if (aaa(1:1) == '!') then
          read(aaa,*)
          cycle
        end if
c        print *,' read line===',trim(aaa),'==='
        iclass = iclass + 1
        read(aaa,*,end=1201) iclassin(iclass),a,(cbastemp(i,iclass),i=1,maxdat)
 1201   continue
        classname_mlwf(iclass) = trim(a)
c       write(*,*) iclassin(iclass),trim(a),(cbastemp(i,iclass),i=1,10)
        do i=1,maxdat
          if(cbastemp(i,iclass)==-999) then
            nbasclass_mlwf(iclass)=i-1
            exit
          endif
        enddo  
      !write(*,"(i5,a,i5)") iclassin(iline),trim(classname(iline)) ,nbasclass_mlwf(iline) (cbas_mlwf(i,iline),i=1,nbasclass_mlwf(iline))
 1001 continue
      nbasclassMax = maxval(nbasclass_mlwf(1:nclass_mlwf))
      print *,'nbasclass_mlwf=',nbasclass_mlwf,nbasclassMaX
      allocate(cbas_mlwf(nbasclassMax,nclass_mlwf))
      cbas_mlwf = cbastemp(1:nbasclassMax,1:nclass_mlwf)
      deallocate(cbastemp)
      nwf = 0
      do iclass=1,nclass_mlwf
        nwf = nwf + nbasclass_mlwf(iclass)
      end do
      close(ifmloc)
      end subroutine s_read_Worb

!!------------------------------------------------
      subroutine s_cal_Worb()
      use m_keyvalue,only: getkeyvalue
      use m_genallcf_v3, only : natom
	implicit none
      integer:: iclass, iclass2, iphidot_plus, ifmloc, iphi_tmp
      integer :: i, j, l_number, correction
      integer :: tmp_atom, tmp_l, iatom, il, iwf, ret,ixatom
      integer :: nnvv(0:10,natom),ix,ioffset(0:10,natom),ll,ioffadd,mm
      integer,allocatable:: l_numbermx(:)

!! Read index for cphi from GWinput: it should be essentially the same as @MNLA_CPHI
      call getkeyvalue("GWinput","<PRODUCT_BASIS>",unit=ifmloc,status=ret) 
      read(ifmloc,*)
      read(ifmloc,*)
      read(ifmloc,*)
      read(ifmloc,*)
      read(ifmloc,*)
      ioffadd = 0
      do 
         read(ifmloc,*,err=888) ixatom, il ,nnvv(il,ixatom)
         ioffset(il,ixatom) = ioffadd
c         write(6,"(' iatom l nnvv=',4i5)") ixatom, il, nnvv(il,ixatom),ioffset(il,ixatom)
         ioffadd = ioffadd + nnvv(il,ixatom) * (2*il+1)
      end do
 888  continue
      close(ifmloc)
      
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
c          print *, "l_number, iclassin(iclass)", l_number , iclassin(iclass)
          mm = cbas_mlwf(i,iclass) - l_number**2
          iphi   (nphix,iwf) = ioffset(l_number,iclass) + mm
          iphidot(nphix,iwf) = iphi(nphix,iwf) + l_number*2+1
          iwf = iwf +1
      enddo
      enddo

      
c      write(6,"('iii:iphi     = ',10i5)") iphi
c      write(6,"('iii:iphidot  = ',10i5)") iphidot
      end subroutine s_cal_Worb

      end module m_read_Worb

