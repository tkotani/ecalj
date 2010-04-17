C$$$c---------------------------------------------------
C$$$c     internal subroutines
C$$$c---------------------------------------------------
C$$$      subroutine struc_checkclass(sname,iout)
C$$$      implicit none
C$$$      character(*),intent(in):: sname
C$$$      integer,intent(out):: iout
C$$$      character(17):: thisname='struc_checkclass'

C$$$      iout=0
C$$$      select case(trim(adjustl(sname)))
C$$$      case('array')
C$$$      case('bz')
C$$$      case('ctrl')
C$$$      case('gw')
C$$$      case('ham')
C$$$      case('lat')
C$$$      case('mix')
C$$$      case('move')
C$$$      case('optic')
C$$$      case('ordn')
C$$$      case('pot')
C$$$      case('site')
C$$$         iout=1
C$$$      case('spec')
C$$$         iout=1
C$$$      case('str')
C$$$      case('strn')
C$$$      case('tb')
C$$$      case default
C$$$         write(*,*) thisname,': error unknown sname=',sname 
C$$$         stop
C$$$      end select

C$$$      end subroutine struc_checkclass


c$$$      subroutine pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x5,rw)
c$$$      use m_struc_def
c$$$      use m_struc_func
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$c     integer(8):: struc(:)
c$$$c      integer(8):: struc(*)
c$$$      real(8):: struc(*)
c$$$      integer:: x1,x2,x3,x4,x5
c$$$      character,intent(in):: rw
c$$$      type (s_ri):: v_ri 
c$$$
c$$$      integer,parameter:: mlist=20
c$$$      integer:: nlist,ilist(mlist,2)
c$$$
c$$$      integer:: ia1,ia2
c$$$      integer:: ib1,ib2
c$$$
c$$$      integer:: itype,ith,i,nsize
c$$$
c$$$      character:: arr='a'
c$$$
c$$$      integer mask(5),range(2,5)
c$$$      character(80) trunam,alias,switch
c$$$
c$$$      character(22):: thisname='pack5upack5_realbody'
c$$$
c$$$#if STRUC_DEBUG_WRITE
c$$$      write(*,*) thisname,': address=',%loc(struc),len(sname)
c$$$#endif
c$$$
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,mask)
c$$$
c$$$      call struc_strtok(trunam, mlist,nlist,ilist)
c$$$
c$$$      if ( nlist <2 ) then
c$$$         write(*,*) thisname,': parameter error: nlist=',nlist
c$$$         write(*,*) 'arg1=',sname
c$$$         stop
c$$$      endif
c$$$
c$$$      i=1
c$$$      ia1=ilist(i,1)
c$$$      ia2=ilist(i,2)
c$$$
c$$$      call struc_checkclass(trunam(ia1:ia2), itype)
c$$$
c$$$      if (itype.eq.0) then
c$$$
c$$$#if STRUC_DEBUG_WRITE
c$$$         write(*,*) thisname,': pass x1,...,x5'
c$$$#endif
c$$$
c$$$         do i=2,nlist
c$$$            ib1=ilist(i,1)
c$$$            ib2=ilist(i,2)
c$$$            select case(i)
c$$$         case(2)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
c$$$     .           arr,rw,range(1,i-1),range(2,i-1), x1,v_ri) 
c$$$         case(3)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
c$$$     .           arr,rw,range(1,i-1),range(2,i-1), x2,v_ri) 
c$$$         case(4)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
c$$$     .           arr,rw,range(1,i-1),range(2,i-1), x3,v_ri) 
c$$$         case(5)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
c$$$     .           arr,rw,range(1,i-1),range(2,i-1), x4,v_ri) 
c$$$         case(6)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
c$$$     .           arr,rw,range(1,i-1),range(2,i-1), x5,v_ri) 
c$$$         case default
c$$$            write(*,*) thisname,': error, nlist>max, nlist=',nlist
c$$$            write(*,*) thisname,': sname=<',sname,'>'
c$$$            stop
c$$$         end select
c$$$      enddo
c$$$
c$$$      else
c$$$
c$$$#if STRUC_DEBUG_WRITE
c$$$         write(*,*) thisname,': x1 is id, x1=',x1
c$$$         write(*,*) thisname,': pass x2,,,.x5'
c$$$#endif
c$$$
c$$$         nsize= nint(struc(1))
c$$$         if (nsize<=0) then
c$$$            write(*,*) thisname,': error, size of the structure is ',nsize
c$$$            write(*,*) thisname,': sname=',sname
c$$$            stop
c$$$         endif
c$$$         ith=1+(x1-1)*nsize
c$$$
c$$$         do i=2,nlist
c$$$            ib1=ilist(i,1)
c$$$            ib2=ilist(i,2)
c$$$            select case(i)
c$$$         case(2)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
c$$$     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x2,v_ri)
c$$$         case(3)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
c$$$     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x3,v_ri)
c$$$         case(4)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
c$$$     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x4,v_ri)
c$$$         case(5)
c$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
c$$$     .           struc(ith),arr,rw,range(1,i-1),range(2,i-1),x5,v_ri)
c$$$         case default
c$$$            write(*,*) thisname,': error ith-mode, nlist>max, nlist=',nlist
c$$$            write(*,*) thisname,': sname=<',sname,'>'
c$$$            stop
c$$$         end select
c$$$      enddo
c$$$
c$$$      endif
c$$$
c$$$      end subroutine pack5upack5_realbody

C$$$      subroutine  struc_spackv_iv(lpack,sname,struc,is1,is2,n,x1)
C$$$      use m_struc_def
C$$$      use m_struc_func
C$$$      implicit none
C$$$      integer,intent(in):: lpack 
C$$$      character(*):: sname
C$$$      real(8):: struc(*)
C$$$      integer,intent(in):: is1,is2,n
C$$$      integer::x1(*)

C$$$      character:: rw,arr
C$$$      integer,parameter:: mlist=20
C$$$      integer:: nlist,ilist(mlist,2)

C$$$      integer:: ia1,ia2
C$$$      integer:: ib1,ib2
C$$$      integer:: i,itype
C$$$      type (s_ri) :: v_ri
C$$$      integer:: size_struc

C$$$      integer mask(5),range(2,5)
C$$$      character(80) trunam,alias,switch

C$$$      character(15):: thisname='struc_spackv_iv'

C$$$      arr='a'
C$$$      select case ( mod(lpack,10) )
C$$$      case (0) 
C$$$         rw='u'
C$$$      case (1)
C$$$         rw='p'
C$$$      end select


C$$$      call salias(sname,trunam,alias,switch,range,mask)

C$$$      call struc_strtok(trunam, mlist,nlist,ilist)

C$$$      i=1
C$$$      ia1=ilist(i,1)
C$$$      ia2=ilist(i,2)
C$$$      i=2
C$$$      ib1=ilist(i,1)
C$$$      ib2=ilist(i,2)

C$$$      call struc_checkclass(trunam(ia1:ia2), itype)

C$$$      size_struc=nint(struc(1))
C$$$      if ( size_struc<=0) then
C$$$         write(*,*) thisname,': error size =', size_struc
C$$$         write(*,*) thisname,': sname=',sname
C$$$      endif

C$$$      do i=is1,is2
C$$$         select case ( lpack/10 )
C$$$         case(0)
C$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
C$$$     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1,v_ri)
C$$$         case(1)
C$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
C$$$     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1(1+(i-is1)*n),v_ri)
C$$$         end select
C$$$      enddo 

C$$$      end subroutine  struc_spackv_iv


C$$$      subroutine  struc_spackv_r8v(lpack,sname,struc,is1,is2,n,x1)
C$$$      use m_struc_def
C$$$      use m_struc_func
C$$$      implicit none
C$$$      integer,intent(in):: lpack 
C$$$      character(*):: sname
C$$$      real(8):: struc(*)
C$$$      integer,intent(in):: is1,is2,n
C$$$      real(8) :: x1(*)

C$$$      character:: rw,arr
C$$$      integer,parameter:: mlist=20
C$$$      integer:: nlist,ilist(mlist,2)

C$$$      integer:: ia1,ia2
C$$$      integer:: ib1,ib2
C$$$      integer:: i,itype
C$$$      type (s_ri) :: v_ri
C$$$      integer:: size_struc

C$$$      integer mask(5),range(2,5)
C$$$      character(80) trunam,alias,switch

C$$$      character(16):: thisname='struc_spackv_r8v'

C$$$      arr='a'
C$$$      select case ( mod(lpack,10) )
C$$$      case (0) 
C$$$         rw='u'
C$$$      case (1)
C$$$         rw='p'
C$$$      end select


C$$$      call salias(sname,trunam,alias,switch,range,mask)

C$$$      call struc_strtok(trunam, mlist,nlist,ilist)

C$$$      i=1
C$$$      ia1=ilist(i,1)
C$$$      ia2=ilist(i,2)
C$$$      i=2
C$$$      ib1=ilist(i,1)
C$$$      ib2=ilist(i,2)

C$$$      call struc_checkclass(trunam(ia1:ia2), itype)

C$$$      size_struc=nint(struc(1))
C$$$      if (size_struc<=0) then
C$$$         write(*,*) thisname,': size=',size_struc
C$$$         write(*,*) thisname,': sname=',sname
C$$$         write(*,*) thisname,': arr,rw=',arr,rw
C$$$         if ( trunam(ia1:ia2).eq.'gw' ) call ugw_show(struc)
C$$$         stop
C$$$      endif

C$$$      do i=is1,is2
C$$$         select case ( lpack/10 )
C$$$         case(0)
C$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
C$$$     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1,v_ri)
C$$$         case(1)
C$$$            call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
C$$$     .           struc(size_struc*(i-1)+1),arr,rw,range(1,1),range(2,1),x1(1+(i-is1)*n),v_ri)
C$$$         end select
C$$$      enddo 

C$$$      end subroutine  struc_spackv_r8v


c----------------------------------------------------------
c     public 
c----------------------------------------------------------
c$$$
c$$$      subroutine pack1(sname,struc,x1)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1
c$$$      character :: rw
c$$$      rw='p'
c$$$c     write(*,*) 'pack2: ',sname
c$$$      call pack5upack5_realbody(sname,struc,x1,x1,x1,x1,x1,rw)
c$$$      end subroutine pack1
c$$$
c$$$
c$$$      subroutine pack2(sname,struc,x1,x2)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2
c$$$      character :: rw
c$$$      rw='p'
c$$$c     write(*,*) 'pack2: ',sname
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x2,x2,x2,rw)
c$$$      end subroutine pack2
c$$$
c$$$      subroutine pack3(sname,struc,x1,x2,x3)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3
c$$$      character :: rw
c$$$      rw='p'
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x3,x3,x3,rw)
c$$$      end subroutine pack3
c$$$
c$$$      subroutine pack4(sname,struc,x1,x2,x3,x4)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3,x4
c$$$      character :: rw
c$$$      rw='p'
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x4,rw)
c$$$      end subroutine pack4
c$$$
c$$$      subroutine pack5(sname,struc,x1,x2,x3,x4,x5)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3,x4,x5
c$$$      character :: rw
c$$$      rw='p'
c$$$c     write(*,*) 'pack5: address=',%loc(struc)
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x5,rw)
c$$$      end subroutine pack5
c$$$
c$$$      subroutine pack(sname,struc,x1,x2,x3,x4,x5)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3,x4,x5
c$$$      call pack5(sname,struc,x1,x2,x3,x4,x5)
c$$$      end subroutine pack
c$$$
c$$$
c$$$      subroutine upack1(sname,struc,x1)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1
c$$$      character :: rw
c$$$      rw='u'
c$$$      call pack5upack5_realbody(sname,struc,x1,x1,x1,x1,x1,rw)
c$$$      end subroutine upack1
c$$$
c$$$      
c$$$      subroutine upack2(sname,struc,x1,x2)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2
c$$$      character :: rw
c$$$      rw='u'
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x2,x2,x2,rw)
c$$$      end subroutine upack2
c$$$
c$$$      subroutine upack3(sname,struc,x1,x2,x3)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3
c$$$      character :: rw
c$$$      rw='u'
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x3,x3,x3,rw)
c$$$      end subroutine upack3
c$$$
c$$$      subroutine upack4(sname,struc,x1,x2,x3,x4)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3,x4
c$$$      character :: rw
c$$$      rw='u'
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x4,rw)
c$$$      end subroutine upack4
c$$$
c$$$      subroutine upack5(sname,struc,x1,x2,x3,x4,x5)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3,x4,x5
c$$$      character :: rw
c$$$      rw='u'
c$$$      call pack5upack5_realbody(sname,struc,x1,x2,x3,x4,x5,rw)
c$$$      end subroutine upack5
c$$$
c$$$      subroutine upack(sname,struc,x1,x2,x3,x4,x5)
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8):: struc
c$$$      integer:: x1,x2,x3,x4,x5
c$$$      call upack5(sname,struc,x1,x2,x3,x4,x5)
c$$$      end subroutine upack
c$$$
c$$$
c$$$      real(8) function dgets(sname,struc)
c$$$C     - Returns one double precision entry in a structure, referred by name
c$$$      use m_struc_def
c$$$      use m_struc_func
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8),intent(in):: struc
c$$$      real(8):: ret
c$$$      integer,parameter:: mlist=10
c$$$      integer:: nlist,ilist(mlist,2)
c$$$      integer:: i,ia1,ia2,itype,ib1,ib2
c$$$      type (s_ri):: v_ri 
c$$$
c$$$      integer mask(5),range(2,5)
c$$$      character(80) trunam,alias,switch
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,mask)
c$$$
c$$$      call struc_strtok(trunam, mlist,nlist,ilist)
c$$$
c$$$      i=1
c$$$      ia1=ilist(i,1)
c$$$      ia2=ilist(i,2)
c$$$      i=2
c$$$      ib1=ilist(i,1)
c$$$      ib2=ilist(i,2)
c$$$
c$$$      call struc_checkclass(trunam(ia1:ia2), itype)
c$$$
c$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
c$$$     .     '1','u',range(1,1),range(2,1),ret,v_ri) 
c$$$
c$$$      dgets=v_ri%r
c$$$
c$$$      end function dgets
c$$$
c$$$
c$$$      real(8) function dgetss(sname,is,struc)
c$$$C     - Returns one double precision entry in a structure, referred by name
c$$$      use m_struc_def
c$$$      use m_struc_func
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8),intent(in):: struc(*)
c$$$      integer,intent(in):: is
c$$$      real(8):: ret
c$$$      integer,parameter:: mlist=10
c$$$      integer:: nlist,ilist(mlist,2)
c$$$      integer:: i,ia1,ia2,itype,ib1,ib2
c$$$
c$$$      integer:: ith
c$$$      type (s_ri):: v_ri
c$$$
c$$$      integer mask(5),range(2,5)
c$$$      character(80) trunam,alias,switch
c$$$
c$$$      character(6):: thisname='dgetss'
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,mask)
c$$$
c$$$      call struc_strtok(trunam, mlist,nlist,ilist)
c$$$
c$$$      i=1
c$$$      ia1=ilist(i,1)
c$$$      ia2=ilist(i,2)
c$$$      i=2
c$$$      ib1=ilist(i,1)
c$$$      ib2=ilist(i,2)
c$$$
c$$$      call struc_checkclass(trunam(ia1:ia2), itype)
c$$$      if (itype.eq.0) then
c$$$         write(*,*) thisname,' is called, with ',sname
c$$$         write(*,*)'but ',sname,' is scalar structure.'
c$$$         stop
c$$$      endif
c$$$
c$$$      if (struc(1)<=0) then
c$$$         write(*,*) thisname,': size=',struc(1)
c$$$         write(*,*) thisname,': sname=',sname
c$$$         stop
c$$$      endif
c$$$      ith = struc(1)*(is-1)+1
c$$$
c$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
c$$$     .     struc(ith),'1','u',range(1,1),range(2,1),ret,v_ri) 
c$$$
c$$$      dgetss=v_ri%r
c$$$      
c$$$      end function dgetss
c$$$
c$$$
c$$$
c$$$      integer function igets(sname,struc)
c$$$C     - Returns one integer entry in a structure, referred by name
c$$$      use m_struc_def
c$$$      use m_struc_func
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8),intent(in):: struc
c$$$      real(8):: ret
c$$$      integer,parameter:: mlist=10
c$$$      integer:: nlist,ilist(mlist,2)
c$$$      integer:: i,ia1,ia2,itype,ib1,ib2
c$$$      type (s_ri):: v_ri 
c$$$
c$$$      integer mask(5),range(2,5)
c$$$      character(80) trunam,alias,switch
c$$$
c$$$      integer:: bitand
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,mask)
c$$$
c$$$      call struc_strtok(trunam, mlist,nlist,ilist)
c$$$
c$$$      i=1
c$$$      ia1=ilist(i,1)
c$$$      ia2=ilist(i,2)
c$$$      i=2
c$$$      ib1=ilist(i,1)
c$$$      ib2=ilist(i,2)
c$$$
c$$$      call struc_checkclass(trunam(ia1:ia2), itype)
c$$$
c$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
c$$$     .     struc,'1','u',range(1,1),range(2,1),ret,v_ri) 
c$$$
c$$$      igets=v_ri%i
c$$$
c$$$      if (mask(1) .gt. 0) igets = bitand(mask,v_ri%i)
c$$$
c$$$
c$$$      end function igets
c$$$
c$$$
c$$$      integer function igetss(sname,is,struc)
c$$$C     - Returns one integer entry in one species of a structure, referred by name
c$$$      use m_struc_def
c$$$      use m_struc_func
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8),intent(in):: struc(*)
c$$$      integer,intent(in):: is
c$$$      real(8):: ret
c$$$      integer,parameter:: mlist=10
c$$$      integer:: nlist,ilist(mlist,2)
c$$$      integer:: i,ia1,ia2,itype,ib1,ib2
c$$$
c$$$      integer:: ith
c$$$      type (s_ri):: v_ri
c$$$
c$$$      integer mask(5),range(2,5),struc_size
c$$$      character(80) trunam,alias,switch
c$$$      character(6):: thisname='igetss'
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,mask)
c$$$
c$$$      call struc_strtok(trunam, mlist,nlist,ilist)
c$$$
c$$$      i=1
c$$$      ia1=ilist(i,1)
c$$$      ia2=ilist(i,2)
c$$$      i=2
c$$$      ib1=ilist(i,1)
c$$$      ib2=ilist(i,2)
c$$$
c$$$      call struc_checkclass(trunam(ia1:ia2), itype)
c$$$      if (itype.eq.0) then
c$$$         write(*,*) thisname,' is called, with ',sname
c$$$         write(*,*)'but ',sname,' is scalar structure.'
c$$$         stop
c$$$      endif
c$$$
c$$$      struc_size= nint(struc(1))
c$$$      if (struc_size <=0) then
c$$$         write(*,*) thisname,': error, size of the structure is ',struc(1),struc_size
c$$$         write(*,*) thisname,': sname=',sname
c$$$         stop
c$$$      endif
c$$$
c$$$      ith = struc_size*(is-1)+1
c$$$
c$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),
c$$$     .     struc(ith),'1','u',range(1,1),range(2,1),ret,v_ri) 
c$$$
c$$$      igetss=v_ri%i
c$$$      
c$$$      end function igetss
c$$$
c$$$
c$$$
c$$$      subroutine lsets(sname,struc,lval,mask)
c$$$C     - Set a logical bits to a structure entry
c$$$C     ----------------------------------------------------------------
c$$$C     i Inputs
c$$$C     i   lval is a logical T or F, of length of sname entry
c$$$C     i   mask should be an integer multiple of 2.
c$$$C     i        Only the lowest bit of mask is used.
c$$$C     o Outputs
c$$$C     o  struc element corresponding to label 'name' is modified.
c$$$C     o        The mask bit of that entry is set to lval.
c$$$C     ----------------------------------------------------------------
c$$$      use m_struc_def 
c$$$      use m_struc_func
c$$$      implicit none
c$$$      character(*),intent(in):: sname
c$$$      real(8),intent(in):: struc
c$$$      logical,intent(in):: lval
c$$$      integer,intent(in):: mask
c$$$
c$$$      real(8):: ret
c$$$      integer,parameter:: mlist=10
c$$$      integer:: nlist,ilist(mlist,2)
c$$$      integer:: i,ia1,ia2,ib1,ib2,itype
c$$$      type (s_ri):: v_ri 
c$$$      integer:: newval,lmask,n
c$$$
c$$$      integer v_mask(5),range(2,5)
c$$$      character(80) trunam,alias,switch
c$$$
c$$$
c$$$      integer:: bitor, bitand ,bitlow
c$$$
c$$$      integer,allocatable:: v_iv(:)
c$$$
c$$$      character(5):: thisname='lsets'
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,v_mask)
c$$$
c$$$      call struc_strtok(trunam, mlist,nlist,ilist)
c$$$
c$$$      i=1
c$$$      ia1=ilist(i,1)
c$$$      ia2=ilist(i,2)
c$$$      i=2
c$$$      ib1=ilist(i,1)
c$$$      ib2=ilist(i,2)
c$$$
c$$$      call struc_checkclass(trunam(ia1:ia2), itype)
c$$$
c$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,'s','u',range(1,1),range(2,1),ret,v_ri) 
c$$$      n= v_ri%i
c$$$      allocate( v_iv( n ) )
c$$$
c$$$#if STRUC_DEBUG_WRITE
c$$$      write(*,*) thisname,' note:  ',struc,' size=',n 
c$$$#endif
c$$$
c$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,'a','u',range(1,1),range(2,1),v_iv,v_ri) 
c$$$
c$$$      if ( v_ri%t.ne.2 ) then
c$$$         write(*,*) thisname,' warning: ',sname,' is not integer, OK?'
c$$$         stop
c$$$      endif
c$$$
c$$$      lmask=bitlow(mask)
c$$$
c$$$      do i=1,n
c$$$
c$$$         if (lval) then
c$$$            newval = bitor(v_iv(i),lmask)
c$$$         else
c$$$            newval = v_iv(i) - bitand(v_iv(i),lmask)
c$$$         endif
c$$$         
c$$$         v_iv(i) = newval 
c$$$
c$$$      enddo
c$$$
c$$$      call struc_packupack_val1(sname(ia1:ia2), sname(ib1:ib2),struc,'a','p',range(1,1),range(2,1),v_iv,v_ri) 
c$$$      deallocate( v_iv )
c$$$
c$$$      end subroutine lsets
c$$$
c$$$#if 1
c$$$C     not tested, probabaly works
c$$$
c$$$      logical function lgors(sname,struc)
c$$$C     - Logical OR of bits in one (masked) entry in a structure
c$$$
c$$$      implicit none
c$$$      character*(*),intent(in):: sname
c$$$      real(8),intent(in):: struc(*)
c$$$      integer mask(5),bitand,range(2,5),ival,i,k
c$$$      double precision x1(20),x2(20),x3(20),x4(20),x5(20)
c$$$      character*40 trunam,alias,switch
c$$$
c$$$      integer:: igets
c$$$      character(5)::thisname='lgors'
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,mask)
c$$$
c$$$      ival= igets(trunam,struc)
c$$$#if STRUC_DEBUG_WRITE
c$$$      write(*,*) thisname,' ival=',ival
c$$$#endif
c$$$      if (mask(1)<0) then
c$$$         lgors = ival .ne. 0
c$$$      else 
c$$$         lgors = bitand(mask(1),ival) .ne. 0
c$$$      endif
c$$$
c$$$      end 
c$$$
c$$$      logical function lgands(sname,struc)
c$$$C     - Logical AND of all bits corresponding to mask in a structure entry
c$$$
c$$$      implicit none
c$$$      character*(*),intent(in):: sname
c$$$      real(8),intent(in):: struc(*)
c$$$      integer mask(5),bitand,range(2,5),ival,i,k
c$$$      double precision x1(20),x2(20),x3(20),x4(20),x5(20)
c$$$      character*40 trunam,alias,switch
c$$$
c$$$      integer:: igets
c$$$
c$$$      call salias(sname,trunam,alias,switch,range,mask)
c$$$
c$$$      ival= igets(trunam,struc)
c$$$
c$$$      lgands = bitand(mask(1),ival) .eq. mask(1) 
c$$$
c$$$      end
c$$$#endif
c$$$

C$$$      subroutine spackv(lpack0,sname,struc,is1,is2,x1)
C$$$C     - Pack/unpack a vector of one element type from a struc containing spec
C$$$C     ----------------------------------------------------------------
C$$$C     i Inputs
C$$$C     i   lpack  1s digit
C$$$C     i           0 unpacks, 1 pack
C$$$C     i         10s digit
C$$$C     i           0 x1 is same for all species (packing only)
C$$$C     i           1 x1 is different for each species
C$$$C     i           2 x1 and its size is different for each species
C$$$C     i             (not implemented)
C$$$C     i is1,is2
C$$$C     i    x1    packed into is1..is2.
C$$$C     r Remarks
C$$$C     r   Structure types must be consistent with those in lstra, above
C$$$      use m_struc_def
C$$$      use m_struc_func
C$$$      implicit none
C$$$      integer,intent(in):: lpack0
C$$$      character(*),intent(in):: sname
C$$$      real(8):: struc(*)
C$$$      integer,intent(in):: is1,is2 
C$$$      integer::x1
      
C$$$      integer:: lpack  
C$$$      character:: rw
C$$$      integer:: n,itype

C$$$      integer,parameter:: mlist=20
C$$$      integer:: nlist,ilist(mlist,2)
C$$$      integer:: ia1,ia2
C$$$      integer:: ib1,ib2
C$$$      integer:: i
C$$$      real(8):: ret
C$$$      type(s_ri):: v_ri

C$$$      integer mask(5),range(2,5)
C$$$      character(80) trunam,alias,switch

C$$$      character(6):: thisname='spackv'

C$$$      lpack=lpack0

C$$$      call salias(sname,trunam,alias,switch,range,mask)


C$$$      select case (lpack)
C$$$c      case (00)  ! unpack, same for all,  not supported
C$$$      case (01)                 ! pack,  same for all
C$$$      case (10)                 ! unpack, different
C$$$      case (11)                 ! pack, different
C$$$      case default
C$$$         if (lpack.eq.0 .and. is1.eq.is2) then
C$$$            lpack=10
C$$$c            write(*,*) thisname,': change lpack=',lpack,' and continue'
C$$$         else

C$$$         write(*,*) thisname,': unsupported lpack, lpack=',lpack
C$$$         write(*,*) thisname,': sname=<',sname,'>'
C$$$         write(*,*) thisname,': is1,is2=',is1,is2
C$$$         stop

C$$$         endif
C$$$      end select

C$$$      select case ( mod(lpack,10) )
C$$$      case (0)
C$$$         rw='u'
C$$$      case (1)
C$$$         rw='p'
C$$$      end select

C$$$      call struc_strtok(trunam, mlist,nlist,ilist)

C$$$      i=1
C$$$      ia1=ilist(i,1)
C$$$      ia2=ilist(i,2)
C$$$      i=2
C$$$      ib1=ilist(i,1)
C$$$      ib2=ilist(i,2)

C$$$      call struc_checkclass(trunam(ia1:ia2), itype)

C$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),struc,
C$$$     .     's','u',range(1,1),range(1,2),ret,v_ri) 
C$$$      n= v_ri%i 
C$$$      itype= v_ri%t


C$$$c it is ok to call using sname, not trunam
C$$$      select case (itype)
C$$$      case (2)
C$$$         call struc_spackv_iv(lpack,sname,struc,is1,is2,n,x1)
C$$$      case (1)
C$$$         call struc_spackv_r8v(lpack,sname,struc,is1,is2,n,x1)
C$$$      case (4)
C$$$         call struc_spackv_r8v(lpack,sname,struc,is1,is2,n,x1)
C$$$      case default
C$$$         write(*,*) thisname,': wrong itype=',itype
C$$$         write(*,*) thisname,': sname=<',sname,'>'
C$$$         write(*,*) thisname,': is1,is2=',is1,is2
C$$$         stop
C$$$      end select

C$$$      end subroutine spackv
      


C$$$      subroutine spacks(lpack,sname,struc,name,is1,is2)
C$$$C     - Pack/unpack string element in a struc, element referenced by name
C$$$C     ----------------------------------------------------------------
C$$$C     i   lpack   0 for unpack, 1 for pack, 2 to return parameters
C$$$C     i   is1..is2 range of species is for which to pack
C$$$C     For structure without species, set both to 0
C$$$C     ----------------------------------------------------------------

C$$$      implicit none
C$$$      integer,intent(in):: lpack
C$$$      character(*),intent(in):: sname
C$$$      real(8):: struc(*)
C$$$      character(*):: name(is1:is2)
C$$$      integer,intent(in):: is1,is2

C$$$      real(8),allocatable:: v_r(:)
C$$$      integer:: i
C$$$      character(6):: thisname='spacks'

C$$$      allocate( v_r( is1:is2 ) )


C$$$      select case (lpack )
C$$$      case (1)
C$$$         do i=is1,is2
C$$$            call s8tor8(name(i),v_r(i))
C$$$         enddo
C$$$         call spackv(11,sname, struc,is1,is2,v_r)
C$$$      case (0)
C$$$         call spackv(10,sname, struc,is1,is2,v_r)
C$$$         do i=is1,is2
C$$$            call r8tos8(v_r(i),name(i))
C$$$         enddo
C$$$      case default
C$$$         write(*,*) thisname,': unknown lpack, lpack=',lpack
C$$$         stop
C$$$      end select 

C$$$      deallocate(v_r)

C$$$      end suroutine spacks


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
c      use m_struc_func
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
      end


C$$$      subroutine sp2cls_org(sname,sspec,ics,nelt,is1,is2,oarr)
C$$$C     - Allocate and extract one element of a structure for a range of species
C$$$C     ----------------------------------------------------------------------
C$$$C     i Inputs
C$$$C     i   sname :a species-based structure and the element to extract, e.g.
C$$$C     i         :'spec idmod'
C$$$C     i   sspec :array holding structure sname
C$$$C     i   ics   :species table: class ic belongs to species ics(ic)
C$$$C     i   nelt  :number of row entries this element has in structure
C$$$C     i   is1   :extract elements starting at species is1
C$$$C     i   is2   :extract elements ending at species is2
C$$$C     o Outputs
C$$$C     o   oarr  :Array oarr is allocated and loaded,
C$$$C     o         :corresponding to arr(1:nelt,is1:is2)

C$$$      use m_struc_def
C$$$      use m_struc_func
C$$$      implicit none
C$$$      character(*),intent(in):: sname
C$$$      real(8):: sspec(*)
C$$$      integer,intent(in):: ics(*)
C$$$      integer,intent(in):: nelt, is1,is2
C$$$      integer,intent(out):: oarr
      
C$$$      integer:: is
C$$$      integer w(1)
C$$$      common /w/ w

C$$$      integer:: i,ioarr,factor

C$$$      type (s_ri):: v_ri 

C$$$      integer,parameter:: mlist=20
C$$$      integer:: nlist,ilist(mlist,2)

C$$$      integer:: ia1,ia2
C$$$      integer:: ib1,ib2

C$$$      integer:: v_i,cast
C$$$      real(8):: v_r

C$$$      integer mask(5),range(2,5)
C$$$      character(80) trunam,alias,switch

C$$$      character(6):: thisname='sp2cls'




C$$$      call salias(sname,trunam,alias,switch,range,mask)

      
C$$$c     must query type
C$$$      call struc_strtok(trunam, mlist,nlist,ilist)
C$$$      i=1
C$$$      ia1=ilist(i,1)
C$$$      ia2=ilist(i,2)
C$$$      i=2
C$$$      ib1=ilist(i,1)
C$$$      ib2=ilist(i,2)
C$$$      call struc_packupack_val1(trunam(ia1:ia2), trunam(ib1:ib2),sspec,'1','r', range(1,1),range(2,1),v_ri,v_ri)

C$$$      cast= v_ri%t

C$$$      select case ( cast ) 
C$$$      case (2)
C$$$         factor=1
C$$$         call defi(oarr,(is2-is1+1)*nelt)
C$$$      case (4)
C$$$         factor=kind(v_r)/kind(v_i)
C$$$         call defrr(oarr,(is2-is1+1)*nelt)
C$$$      case default
C$$$         write(*,*) thisname,': error, unknown type=',v_ri%t
C$$$         stop
C$$$      end select

C$$$      do i=is1,is2
C$$$         if (ics(i).eq.0) then; is=i;
C$$$         else; is=ics(i); endif
C$$$         ioarr=oarr+ (i-is1+1)*nelt*factor
C$$$         call spackv(10,sname,sspec,is1,is2,w(ioarr))
C$$$      enddo

C$$$      end subroutine


