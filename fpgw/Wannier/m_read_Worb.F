c--------------------------------------
c This module is for readling informations of <Worb> section in GWinput.
c Then this is transformed to the form for the calculation of MLWF.

c nclass_mlwf : the number of material
c classname_mlwf : the name of materials
c nbasclass_mlwf : the number of Wannier orbitals for each line.
c cbas_mlwf : each Wannier orbital (e.g) s = 1, py = 2,... 
c-------------------------------------

      module m_read_Worb
!! all output
      integer,protected:: nwf, nclass_mlwf
        ! cbas_mlwf(nbasclassMax,nclass_mlwf) !We may need nbasclassMax to pass array to subroutine.
      integer,protected,allocatable:: cbas_mlwf(:,:),nbasclass_mlwf(:)
      character(20),protected,allocatable:: classname_mlwf(:)
      integer,protected,allocatable:: iclassin(:)
      integer,parameter,private::maxdat=1024 !it is convenient if you use this type of number.
      integer(4),protected,allocatable::  iphi(:,:),iphidot(:,:), nphi(:)
      integer(4),protected:: nphix

      contains

      subroutine s_read_Worb()
      use keyvalue,only: getkeyvalue
      implicit none
      integer(4):: ret
      integer:: i, ifile_handle,il,ix,iclass,iline,nbasclassMax
      character*256:: a,aaa
      integer::ifmloc
      integer,allocatable:: cbastemp(:,:)

      
      call getkeyvalue("GWinput","<Worb>",unit=ifmloc,status=ret) 
      write(6,*)' ifmloc ret=',ifmloc,ret

      iline=0
      do 
        read(ifmloc,"(a)") aaa
        print *, aaa
        if(aaa(1:7) == "</Worb>") then
          goto 999
        else 
          iline=iline+1
        end if
      end do
 999  continue


      nclass_mlwf=iline
      print *, nclass_mlwf

      allocate(iclassin(nclass_mlwf),cbastemp(maxdat,nclass_mlwf),nbasclass_mlwf(nclass_mlwf),classname_mlwf(nclass_mlwf))

      close(ifmloc)

      call getkeyvalue("GWinput","<Worb>",unit=ifmloc,status=ret) 

      cbastemp=-999
      do 1001 iclass=1,nclass_mlwf
        read(ifmloc,"(a)") aaa
c        print *,' read line===',trim(aaa),'==='
        read(aaa,*,end=1201) iclassin(iclass),a,(cbastemp(i,iclass),i=1,maxdat)
 1201   continue
        !print *, iclassin(iclass), a, cbastemp(i,iclass)
        !if(iclassin(iclass)/=iclass) call rx('iclass is not i')
        classname_mlwf(iclass)=trim(a)
c        write(*,*) iclassin(iclass),trim(a),(cbastemp(i,iclass),i=1,10)
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





      subroutine s_cal_Worb()
	implicit none
      integer:: iclass, l_number, correction, i, j, iwf

      allocate (nphi(nwf))
      nphix = 0
      do i =1,nwf
        nphi(i) = 1
        if(nphi(i).gt.nphix)nphix = nphi(i)

      end do

      do i = 1,nwf
         if(nphi(i).gt.nphix)nphix = nphi(i)
      enddo

      allocate (iphi(nphix,nwf),iphidot(nphix,nwf))

      iphi = 0
      print *, iphi
      iphidot = 0
      print *, iphidot
      !print *, iphi, iphidot


      iwf = 1
      do iclass=1,nclass_mlwf
        do i=1,nbasclass_mlwf(iclass)
          !print *, iclassin(iclass)
          if (cbas_mlwf(i,iclass) == 1) then
            l_number = 0
          else if (int(cbas_mlwf(i,iclass))/4.0 <= 1) then
            l_number = 1
          else if (int(cbas_mlwf(i,iclass))/9.0 <= 1) then
            l_number = 2
          else if (int(cbas_mlwf(i,iclass))/16.0 <= 1) then
            l_number = 3
          else if (int(cbas_mlwf(i,iclass))/25.0 <= 1) then
            l_number = 4
          end if

	    print *, l_number
	   
	    correction = 0

	    do j=0, l_number-1
	      correction = correction + 2*(2*j+1)
	    end do		

	    print *, iclassin(iclass)
          iphi(nphix,iwf) = 2 * (iclassin(iclass)-1) * (1+3+5+7+9) +
     &                     correction  + cbas_mlwf(i,iclass)
          iphidot(nphix,iwf) = iphi(nphix,iwf) + l_number*2+1
          iwf = iwf +1
        enddo
      enddo
      print *, "iphi = ", iphi
      print *, "iphidot = ", iphidot

      end subroutine s_cal_Worb

      end module m_read_Worb

!------------------------------------------------------------------
