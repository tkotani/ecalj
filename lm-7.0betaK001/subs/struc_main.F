c---------------------------------------------------
c     internal subroutines
c---------------------------------------------------
      subroutine struc_checkclass(sname,iout)
      implicit none
      character(*),intent(in):: sname
      integer,intent(out):: iout
      character(17):: thisname='struc_checkclass'

      iout=0
      select case(trim(adjustl(sname)))
      case('array')
      case('bz')
      case('ctrl')
      case('gw')
      case('ham')
      case('lat')
      case('mix')
      case('move')
      case('optic')
      case('ordn')
      case('pot')
      case('site')
         iout=1
      case('spec')
         iout=1
      case('str')
      case('strn')
      case('tb')
      case default
         write(*,*) thisname,': error unknown sname=',sname 
         stop
      end select

      end subroutine struc_checkclass


      subroutine pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x5,rw)
      use m_struc_def
      use m_struc_func
      implicit none
      character(*),intent(in):: sname
c     integer(8):: struc(:)
c      integer(8):: struc(*)
      real(8):: struc(*)
      integer:: x1,x2,x3,x4,x5
      character,intent(in):: rw
      type (s_ri):: v_ri 

      integer,parameter:: mlist=20
      integer:: nlist,ilist(mlist,2)

      integer:: ia1,ia2
      integer:: ib1,ib2

      integer:: itype,ith,i,nsize

      character:: arr='a'

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      character(22):: thisname='pack5upack5_realbody'

#if STRUC_DEBUG_WRITE
      write(*,*) thisname,': address=',%loc(struc),len(sname)
#endif



      call salias(sname,trunam,alias,switch,range,mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      if ( nlist <2 ) then
         write(*,*) thisname,': parameter error: nlist=',nlist
         write(*,*) 'arg1=',sname
         stop
      endif

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)

      if (itype.eq.0) then

#if STRUC_DEBUG_WRITE
         write(*,*) thisname,': pass x1,...,x5'
#endif

         do i=2,nlist
            ib1=ilist(i,1)
            ib2=ilist(i,2)
            select case(i)
         case(2)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
     .           arr,rw,range(1,i-1),range(2,i-1), x1,v_ri) 
         case(3)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
     .           arr,rw,range(1,i-1),range(2,i-1), x2,v_ri) 
         case(4)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
     .           arr,rw,range(1,i-1),range(2,i-1), x3,v_ri) 
         case(5)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
     .           arr,rw,range(1,i-1),range(2,i-1), x4,v_ri) 
         case(6)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
     .           arr,rw,range(1,i-1),range(2,i-1), x5,v_ri) 
         case default
            write(*,*) thisname,': error, nlist>max, nlist=',nlist
            write(*,*) thisname,': sname=<',sname,'>'
            stop
         end select
      enddo

      else

#if STRUC_DEBUG_WRITE
         write(*,*) thisname,': x1 is id, x1=',x1
         write(*,*) thisname,': pass x2,,,.x5'
#endif

         nsize= nint(struc(1))
         if (nsize<=0) then
            write(*,*) thisname,': error, size of the structure is ',nsize
            write(*,*) thisname,': sname=',sname
            stop
         endif
         ith=1+(x1-1)*nsize

         do i=2,nlist
            ib1=ilist(i,1)
            ib2=ilist(i,2)
            select case(i)
         case(2)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x2,v_ri)
         case(3)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x3,v_ri)
         case(4)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x4,v_ri)
         case(5)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x5,v_ri)
         case default
            write(*,*) thisname,': error ith-mode, nlist>max, nlist=',nlist
            write(*,*) thisname,': sname=<',sname,'>'
            stop
         end select
      enddo

      endif

      end subroutine pack5upack5_realbody

      subroutine  struc_spackv_iv(lpack,sname,struc,is1,is2,n,x1)
      use m_struc_def
      use m_struc_func
      implicit none
      integer,intent(in):: lpack 
      character(*):: sname
      real(8):: struc(*)
      integer,intent(in):: is1,is2,n
      integer::x1(*)

      character:: rw,arr
      integer,parameter:: mlist=20
      integer:: nlist,ilist(mlist,2)

      integer:: ia1,ia2
      integer:: ib1,ib2
      integer:: i,itype
      type (s_ri) :: v_ri
      integer:: size_struc

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      character(15):: thisname='struc_spackv_iv'

      arr='a'
      select case ( mod(lpack,10) )
      case (0) 
         rw='u'
      case (1)
         rw='p'
      end select


      call salias(sname,trunam,alias,switch,range,mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)

      size_struc=nint(struc(1))
      if ( size_struc<=0) then
         write(*,*) thisname,': error size =', size_struc
         write(*,*) thisname,': sname=',sname
      endif

      do i=is1,is2
         select case ( lpack/10 )
         case(0)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1,v_ri)
         case(1)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1(1+(i-is1)*n),v_ri)
         end select
      enddo 

      end subroutine  struc_spackv_iv


      subroutine  struc_spackv_r8v(lpack,sname,struc,is1,is2,n,x1)
      use m_struc_def
      use m_struc_func
      implicit none
      integer,intent(in):: lpack 
      character(*):: sname
      real(8):: struc(*)
      integer,intent(in):: is1,is2,n
      real(8) :: x1(*)

      character:: rw,arr
      integer,parameter:: mlist=20
      integer:: nlist,ilist(mlist,2)

      integer:: ia1,ia2
      integer:: ib1,ib2
      integer:: i,itype
      type (s_ri) :: v_ri
      integer:: size_struc

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      character(16):: thisname='struc_spackv_r8v'

      arr='a'
      select case ( mod(lpack,10) )
      case (0) 
         rw='u'
      case (1)
         rw='p'
      end select


      call salias(sname,trunam,alias,switch,range,mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)

      size_struc=nint(struc(1))
      if (size_struc<=0) then
         write(*,*) thisname,': size=',size_struc
         write(*,*) thisname,': sname=',sname
         stop
      endif

      do i=is1,is2
         select case ( lpack/10 )
         case(0)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1,v_ri)
         case(1)
            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1(1+(i-is1)*n),v_ri)
         end select
      enddo 

      end subroutine  struc_spackv_r8v


c----------------------------------------------------------
c     public 
c----------------------------------------------------------

      subroutine pack1(sname,struc,x1)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1
      character :: rw
      rw='p'
c     write(*,*) 'pack2: ',sname
      call pack5upack5_realbody(sname,struc,x1,x1,x1,x1,x1,rw)
      end subroutine pack1


      subroutine pack2(sname,struc,x1,x2)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2
      character :: rw
      rw='p'
c     write(*,*) 'pack2: ',sname
      call pack5upack5_realbody(sname,struc,x1,x2,x2,x2,x2,rw)
      end subroutine pack2

      subroutine pack3(sname,struc,x1,x2,x3)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3
      character :: rw
      rw='p'
      call pack5upack5_realbody(sname,struc,x1,x2,x3,x3,x3,rw)
      end subroutine pack3

      subroutine pack4(sname,struc,x1,x2,x3,x4)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3,x4
      character :: rw
      rw='p'
      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x4,rw)
      end subroutine pack4

      subroutine pack5(sname,struc,x1,x2,x3,x4,x5)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3,x4,x5
      character :: rw
      rw='p'
c     write(*,*) 'pack5: address=',%loc(struc)
      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x5,rw)
      end subroutine pack5

      subroutine pack(sname,struc,x1,x2,x3,x4,x5)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3,x4,x5
      call pack5(sname,struc,x1,x2,x3,x4,x5)
      end subroutine pack


      subroutine upack1(sname,struc,x1)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1
      character :: rw
      rw='u'
      call pack5upack5_realbody(sname,struc,x1,x1,x1,x1,x1,rw)
      end subroutine upack1

      
      subroutine upack2(sname,struc,x1,x2)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2
      character :: rw
      rw='u'
      call pack5upack5_realbody(sname,struc,x1,x2,x2,x2,x2,rw)
      end subroutine upack2

      subroutine upack3(sname,struc,x1,x2,x3)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3
      character :: rw
      rw='u'
      call pack5upack5_realbody(sname,struc,x1,x2,x3,x3,x3,rw)
      end subroutine upack3

      subroutine upack4(sname,struc,x1,x2,x3,x4)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3,x4
      character :: rw
      rw='u'
      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x4,rw)
      end subroutine upack4

      subroutine upack5(sname,struc,x1,x2,x3,x4,x5)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3,x4,x5
      character :: rw
      rw='u'
      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x5,rw)
      end subroutine upack5

      subroutine upack(sname,struc,x1,x2,x3,x4,x5)
      implicit none
      character(*),intent(in):: sname
      real(8):: struc
      integer:: x1,x2,x3,x4,x5
      call upack5(sname,struc,x1,x2,x3,x4,x5)
      end subroutine upack


      real(8) function dgets(sname,struc)
C     - Returns one double precision entry in a structure, referred by name
      use m_struc_def
      use m_struc_func
      implicit none
      character(*),intent(in):: sname
      real(8),intent(in):: struc
      real(8):: ret
      integer,parameter:: mlist=10
      integer:: nlist,ilist(mlist,2)
      integer:: i,ia1,ia2,itype,ib1,ib2
      type (s_ri):: v_ri 

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      call salias(sname,trunam,alias,switch,range,mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)

      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
     .     '1','u',range(1,1),range(2,1),ret,v_ri) 

      dgets=v_ri%r

      end function dgets


      real(8) function dgetss(sname,is,struc)
C     - Returns one double precision entry in a structure, referred by name
      use m_struc_def
      use m_struc_func
      implicit none
      character(*),intent(in):: sname
      real(8),intent(in):: struc(*)
      integer,intent(in):: is
      real(8):: ret
      integer,parameter:: mlist=10
      integer:: nlist,ilist(mlist,2)
      integer:: i,ia1,ia2,itype,ib1,ib2

      integer:: ith
      type (s_ri):: v_ri

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      character(6):: thisname='dgetss'

      call salias(sname,trunam,alias,switch,range,mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)
      if (itype.eq.0) then
         write(*,*) thisname,' is called, with ',sname
         write(*,*)'but ',sname,' is scalar structure.'
         stop
      endif

      if (struc(1)<=0) then
         write(*,*) thisname,': size=',struc(1)
         write(*,*) thisname,': sname=',sname
         stop
      endif
      ith = struc(1)*(is-1)+1

      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .     struc(ith),'1','u',range(1,1),range(2,1),ret,v_ri) 

      dgetss=v_ri%r
      
      end function dgetss



      integer function igets(sname,struc)
C     - Returns one integer entry in a structure, referred by name
      use m_struc_def
      use m_struc_func
      implicit none
      character(*),intent(in):: sname
      real(8),intent(in):: struc
      real(8):: ret
      integer,parameter:: mlist=10
      integer:: nlist,ilist(mlist,2)
      integer:: i,ia1,ia2,itype,ib1,ib2
      type (s_ri):: v_ri 

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      integer:: bitand

      call salias(sname,trunam,alias,switch,range,mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)

      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .     struc,'1','u',range(1,1),range(2,1),ret,v_ri) 

      igets=v_ri%i

      if (mask(1) .gt. 0) igets = bitand(mask,v_ri%i)


      end function igets


      integer function igetss(sname,is,struc)
C     - Returns one integer entry in one species of a structure, referred by name
      use m_struc_def
      use m_struc_func
      implicit none
      character(*),intent(in):: sname
      real(8),intent(in):: struc(*)
      integer,intent(in):: is
      real(8):: ret
      integer,parameter:: mlist=10
      integer:: nlist,ilist(mlist,2)
      integer:: i,ia1,ia2,itype,ib1,ib2

      integer:: ith
      type (s_ri):: v_ri

      integer mask(5),range(2,5),struc_size
      character(80) trunam,alias,switch
      character(6):: thisname='igetss'

      call salias(sname,trunam,alias,switch,range,mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)
      if (itype.eq.0) then
         write(*,*) thisname,' is called, with ',sname
         write(*,*)'but ',sname,' is scalar structure.'
         stop
      endif

      struc_size= nint(struc(1))
      if (struc_size <=0) then
         write(*,*) thisname,': error, size of the structure is ',struc(1),struc_size
         write(*,*) thisname,': sname=',sname
         stop
      endif

      ith = struc_size*(is-1)+1

      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
     .     struc(ith),'1','u',range(1,1),range(2,1),ret,v_ri) 

      igetss=v_ri%i
      
      end function igetss



      subroutine lsets(sname,struc,lval,mask)
C     - Set a logical bits to a structure entry
C     ----------------------------------------------------------------
C     i Inputs
C     i   lval is a logical T or F, of length of sname entry
C     i   mask should be an integer multiple of 2.
C     i        Only the lowest bit of mask is used.
C     o Outputs
C     o  struc element corresponding to label 'name' is modified.
C     o        The mask bit of that entry is set to lval.
C     ----------------------------------------------------------------
      use m_struc_def 
      use m_struc_func
      implicit none
      character(*),intent(in):: sname
      real(8),intent(in):: struc
      logical,intent(in):: lval
      integer,intent(in):: mask

      real(8):: ret
      integer,parameter:: mlist=10
      integer:: nlist,ilist(mlist,2)
      integer:: i,ia1,ia2,ib1,ib2,itype
      type (s_ri):: v_ri 
      integer:: newval,lmask,n

      integer v_mask(5),range(2,5)
      character(80) trunam,alias,switch


      integer:: bitor, bitand ,bitlow

      integer,allocatable:: v_iv(:)

      character(5):: thisname='lsets'

      call salias(sname,trunam,alias,switch,range,v_mask)

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)

      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,'s','u',range(1,1),range(2,1),ret,v_ri) 
      n= v_ri%i
      allocate( v_iv( n ) )

#if STRUC_DEBUG_WRITE
      write(*,*) thisname,' note:  ',struc,' size=',n 
#endif

      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,'a','u',range(1,1),range(2,1),v_iv,v_ri) 

      if ( v_ri%t.ne.2 ) then
         write(*,*) thisname,' warning: ',sname,' is not integer, OK?'
         stop
      endif

      lmask=bitlow(mask)

      do i=1,n

         if (lval) then
            newval = bitor(v_iv(i),lmask)
         else
            newval = v_iv(i) - bitand(v_iv(i),lmask)
         endif
         
         v_iv(i) = newval 

      enddo

      call struc_packupack_val1(sname(ia1:ia2), sname(ib1:ib2),struc,'a','p',range(1,1),range(2,1),v_iv,v_ri) 
      deallocate( v_iv )

      end subroutine lsets

#if 1
C     not tested, probabaly works

      logical function lgors(sname,struc)
C     - Logical OR of bits in one (masked) entry in a structure

      implicit none
      character*(*),intent(in):: sname
      real(8),intent(in):: struc(*)
      integer mask(5),bitand,range(2,5),ival,i,k
      double precision x1(20),x2(20),x3(20),x4(20),x5(20)
      character*40 trunam,alias,switch

      integer:: igets
      character(5)::thisname='lgors'

      call salias(sname,trunam,alias,switch,range,mask)

      ival= igets(trunam,struc)
#if STRUC_DEBUG_WRITE
      write(*,*) thisname,' ival=',ival
#endif
      if (mask(1)<0) then
         lgors = ival .ne. 0
      else 
         lgors = bitand(mask(1),ival) .ne. 0
      endif

      end 

      logical function lgands(sname,struc)
C     - Logical AND of all bits corresponding to mask in a structure entry

      implicit none
      character*(*),intent(in):: sname
      real(8),intent(in):: struc(*)
      integer mask(5),bitand,range(2,5),ival,i,k
      double precision x1(20),x2(20),x3(20),x4(20),x5(20)
      character*40 trunam,alias,switch

      integer:: igets

      call salias(sname,trunam,alias,switch,range,mask)

      ival= igets(trunam,struc)

      lgands = bitand(mask(1),ival) .eq. mask(1) 

      end
#endif


      subroutine spackv(lpack0,sname,struc,is1,is2,x1)
C     - Pack/unpack a vector of one element type from a struc containing spec
C     ----------------------------------------------------------------
C     i Inputs
C     i   lpack  1s digit
C     i           0 unpacks, 1 pack
C     i         10s digit
C     i           0 x1 is same for all species (packing only)
C     i           1 x1 is different for each species
C     i           2 x1 and its size is different for each species
C     i             (not implemented)
C     i is1,is2
C     i    x1    packed into is1..is2.
C     r Remarks
C     r   Structure types must be consistent with those in lstra, above
      use m_struc_def
      use m_struc_func
      implicit none
      integer,intent(in):: lpack0
      character(*),intent(in):: sname
      real(8):: struc(*)
      integer,intent(in):: is1,is2 
      integer::x1
      
      integer:: lpack  
      character:: rw
      integer:: n,itype

      integer,parameter:: mlist=20
      integer:: nlist,ilist(mlist,2)
      integer:: ia1,ia2
      integer:: ib1,ib2
      integer:: i
      real(8):: ret
      type(s_ri):: v_ri

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      character(6):: thisname='spackv'

      lpack=lpack0

      call salias(sname,trunam,alias,switch,range,mask)


      select case (lpack)
c      case (00)  ! unpack, same for all,  not supported
      case (01)                 ! pack,  same for all
      case (10)                 ! unpack, different
      case (11)                 ! pack, different
      case default
         if (lpack.eq.0 .and. is1.eq.is2) then
            lpack=10
c            write(*,*) thisname,': change lpack=',lpack,' and continue'
         else

         write(*,*) thisname,': unsupported lpack, lpack=',lpack
         write(*,*) thisname,': sname=<',sname,'>'
         write(*,*) thisname,': is1,is2=',is1,is2
         stop

         endif
      end select

      select case ( mod(lpack,10) )
      case (0)
         rw='u'
      case (1)
         rw='p'
      end select

      call struc_strtok(trunam, mlist,nlist,ilist)

      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)

      call struc_checkclass(trunam(ia1:ia2), itype)

      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
     .     's','u',range(1,1),range(1,2),ret,v_ri) 
      n= v_ri%i 
      itype= v_ri%t


c it is ok to call using sname, not trunam
      select case (itype)
      case (2)
         call struc_spackv_iv(lpack,sname,struc,is1,is2,n,x1)
      case (1)
         call struc_spackv_r8v(lpack,sname,struc,is1,is2,n,x1)
      case (4)
         call struc_spackv_r8v(lpack,sname,struc,is1,is2,n,x1)
      case default
         write(*,*) thisname,': wrong itype=',itype
         write(*,*) thisname,': sname=<',sname,'>'
         write(*,*) thisname,': is1,is2=',is1,is2
         stop
      end select

      end subroutine spackv
      


      subroutine spacks(lpack,sname,struc,name,is1,is2)
C     - Pack/unpack string element in a struc, element referenced by name
C     ----------------------------------------------------------------
C     i   lpack   0 for unpack, 1 for pack, 2 to return parameters
C     i   is1..is2 range of species is for which to pack
C     For structure without species, set both to 0
C     ----------------------------------------------------------------

      implicit none
      integer,intent(in):: lpack
      character(*),intent(in):: sname
      real(8):: struc(*)
      character(*):: name(is1:is2)
      integer,intent(in):: is1,is2

      real(8),allocatable:: v_r(:)
      integer:: i
      character(6):: thisname='spacks'

      allocate( v_r( is1:is2 ) )


      select case (lpack )
      case (1)
         do i=is1,is2
            call s8tor8(name(i),v_r(i))
         enddo
         call spackv(11,sname, struc,is1,is2,v_r)
      case (0)
         call spackv(10,sname, struc,is1,is2,v_r)
         do i=is1,is2
            call r8tos8(v_r(i),name(i))
         enddo
      case default
         write(*,*) thisname,': unknown lpack, lpack=',lpack
         stop
      end select 

      deallocate(v_r)

      end subroutine spacks


      subroutine sp2cls(sname,sspec,ics,nelt,is1,is2,oarr)
C     - Allocate and extract one element of a structure for a range of species
C     ----------------------------------------------------------------------
C     i Inputs
C     i   sname :a species-based structure and the element to extract, e.g.
C     i         :'spec idmod'
C     i   sspec :array holding structure sname
C     i   ics   :species table: class ic belongs to species ics(ic)
C     i   nelt  :number of row entries this element has in structure
C     i   is1   :extract elements starting at species is1
C     i   is2   :extract elements ending at species is2
C     o Outputs
C     o   oarr  :Array oarr is allocated and loaded,
C     o         :corresponding to arr(1:nelt,is1:is2)

      use m_struc_def
      use m_struc_func
      implicit none
      character(*),intent(in):: sname
      real(8):: sspec(*)
      integer,intent(in):: ics(*)
      integer,intent(in):: nelt, is1,is2
      integer,intent(out):: oarr
      
      integer:: is
      integer w(1)
      common /w/ w

      integer:: i,ioarr,factor

      type (s_ri):: v_ri 

      integer,parameter:: mlist=20
      integer:: nlist,ilist(mlist,2)

      integer:: ia1,ia2
      integer:: ib1,ib2

      integer:: v_i,cast
      real(8):: v_r

      integer mask(5),range(2,5)
      character(80) trunam,alias,switch

      character(6):: thisname='sp2cls'

      call salias(sname,trunam,alias,switch,range,mask)

      
c     must query type
      call struc_strtok(trunam, mlist,nlist,ilist)
      i=1
      ia1=ilist(i,1)
      ia2=ilist(i,2)
      i=2
      ib1=ilist(i,1)
      ib2=ilist(i,2)
      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),sspec,'1','r', range(1,1),range(2,1),v_ri,v_ri)

      cast= v_ri%t

      select case ( cast ) 
      case (2)
         factor=1
         call defi(oarr,(is2-is1+1)*nelt)
      case (4)
         factor=kind(v_r)/kind(v_i)
         call defrr(oarr,(is2-is1+1)*nelt)
      case default
         write(*,*) thisname,': error, unknown type=',v_ri%t
         stop
      end select

      do i=is1,is2
         if (ics(i).eq.0) then; is=i;
         else; is=ics(i); endif
         ioarr=oarr+ (i-is1+1)*nelt*factor
         call spackv(10,sname,sspec,is1,is2,w(ioarr))
      enddo

      end subroutine sp2cls


      subroutine shstru
      end subroutine shstru

      subroutine s_mix_init(v)
      use m_struc_def
      implicit none
      type (s_mix) :: v
      real(8) :: v_r
      integer :: v_i

        v_i = sizeof (v) / kind(v_r) +2
        if (mod(v_i,2).eq.1) v_i=v_i+1
        v%size = v_i 

        v%b = 0.0d0
        v%bl =0.0d0
        v%bv =0.0d0
        v%elind =0.0d0
        v%fn =0.0d0 !string
         v%kill=0
         v%lxpot=0
         v%mmix=0
         v%mode=0
         v%model=0
         v%n=0
         v%nitu=0
         v%nmix=0
         v%nsave=0
        v%r(:)=0.0d0 !string
        v%rms1=0.0d0
        v%rms2=0.0d0
        v%tj(:)=0.0d0
        v%tolu=0.0d0
        v%umix=0.0d0
        v%w(:)=0.0d0
        v%wc=0.0d0

      end subroutine s_mix_init 

