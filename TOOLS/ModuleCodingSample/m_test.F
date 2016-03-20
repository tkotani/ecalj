c---------------------------
c To read a following file.
c 1 Ndup 3 4 8
c 2 Nddn 5 6 7 8
c 3 Cu 2 3 4 123 556 45
c----------------------------
c You may need -132 option
c----------------------------------------
!! ----------------------------------------------------------


!!----------------------------------------------------------------------------
      module m_readline
      integer,protected:: nclass,nbasclassMax 
        ! cbas(nbasclassMax,nclass) !We may need nbasclassMax to pass array to subroutine.
      integer,protected,allocatable:: cbas(:,:),nbasclass(:)
      character(20),protected,allocatable:: classname(:)

      integer,parameter,private::maxdat=1024 !it is convenient if you use this type of number.


      contains

      subroutine s_readclass()
      implicit none
      integer:: i,ifix,ifile_handle,il,ix,iclass,iline
      integer,allocatable:: iclassin(:)
      character*256:: a,aaa
      integer,allocatable:: cbastemp(:,:)

      ifix=ifile_handle()
      open(ifix,file='sample.dat')
      iline=0
      do 
        read(ifix,"(a)",end=999) aaa
        iline=iline+1
      enddo
 999  continue
      nclass=iline
      allocate(iclassin(nclass),cbastemp(maxdat,nclass),nbasclass(nclass),classname(nclass))

      rewind(ifix)
      cbastemp=-999
      do 1001, iclass=1,nclass
        read(ifix,"(a)") aaa
c        print *,' read line===',trim(aaa),'==='
        read(aaa,*,end=1201) iclassin(iclass),a,(cbastemp(i,iclass),i=1,maxdat)
 1201   continue
        if(iclassin(iclass)/=iclass) call rx('iclass is not i')
        classname(iclass)=trim(a)
c        write(*,*) iclassin(iclass),trim(a),(cbastemp(i,iclass),i=1,10)
        do i=1,maxdat
          if(cbastemp(i,iclass)==-999) then
            nbasclass(iclass)=i-1
            exit
          endif
        enddo  
c    write(*,"(i5,a,i5)") iclassin(iline),trim(classname(iline)) ,nbasclass(iline) (cbas(i,iline),i=1,nbasclass(iline))
 1001 continue
      nbasclassMax = maxval(nbasclass(1:nclass))
      print *,'nbasclass=',nbasclass,nbasclassMaX
      allocate(cbas(nbasclassMax,nclass))
      cbas = cbastemp(1:nbasclassMax,1:nclass)
      deallocate(cbastemp)
      end subroutine s_readclass
      end module m_readline

!----------------------------------------------------------------
!! These are already in fpgw/ library
      integer function ifile_handle()
!! find open file handle
      implicit none
      integer:: i
      logical:: nexist
      do i=5001,9999
         inquire(unit=i,opened=nexist)
         if(.not.nexist) then
            ifile_handle=i
            return
         endif
      enddo
c      call rx('ifile_hundle: we did not find open file hundle')
      end
      subroutine rx(aaa)
      character*(*):: aaa
      print *,aaa
      stop 'error exit'
      end

!------------------------------------------------------------------
      program test
      use m_readline,only: s_readclass, 
     &  nbasclass, nclass, cbas, classname, nbasclassmax
      integer:: i,ix,iclass
!!
      call s_readclass()

      write(*,*) '=== Read lines nclass=',nclass
      do iclass=1,nclass
        write(*,*)'output:',iclass
     & ,trim(classname(iclass)),cbas(1:nbasclass(iclass),iclass)
      enddo
      end
