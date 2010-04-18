      module m_struc_func

      end module


c$$$
c$$$      interface struc_eval_io
c$$$      module procedure struc_eval_io_r8, struc_eval_io_i8
c$$$     .     , struc_eval_io_r8v, struc_eval_io_i8v
c$$$      end interface struc_eval_io 
c$$$
c$$$      contains
c$$$
c$$$
c$$$        integer function uarray_size()
c$$$        integer:: n
c$$$        n= 32
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        uarray_size=n
c$$$        end function uarray_size
c$$$       
c$$$        integer function ubz_size()
c$$$        integer:: n
c$$$        n= 47
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        ubz_size=n
c$$$        end function ubz_size
c$$$       
c$$$        integer function uctrl_size()
c$$$        integer:: n
c$$$        n= 88
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        uctrl_size=n
c$$$        end function uctrl_size
c$$$       
c$$$        integer function ugw_size()
c$$$        integer:: n
c$$$        n= 22
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        ugw_size=n
c$$$        end function ugw_size
c$$$       
c$$$        integer function uham_size()
c$$$        integer:: n
c$$$        n= 129
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        uham_size=n
c$$$        end function uham_size
c$$$       
c$$$        integer function ulat_size()
c$$$        integer:: n
c$$$        n= 120
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        ulat_size=n
c$$$        end function ulat_size
c$$$       
c$$$        integer function umix_size()
c$$$        integer:: n
c$$$        n= 36
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        umix_size=n
c$$$        end function umix_size
c$$$       
c$$$        integer function umove_size()
c$$$        integer:: n
c$$$        n= 34
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        umove_size=n
c$$$        end function umove_size
c$$$       
c$$$        integer function uoptic_size()
c$$$        integer:: n
c$$$        n= 33
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        uoptic_size=n
c$$$        end function uoptic_size
c$$$       
c$$$        integer function uordn_size()
c$$$        integer:: n
c$$$        n= 13
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        uordn_size=n
c$$$        end function uordn_size
c$$$       
c$$$        integer function upot_size()
c$$$        integer:: n
c$$$        n= 63
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        upot_size=n
c$$$        end function upot_size
c$$$       
c$$$        integer function usite_size()
c$$$        integer:: n
c$$$        n= 86
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        usite_size=n
c$$$        end function usite_size
c$$$       
c$$$        integer function uspec_size()
c$$$        integer:: n
c$$$        n= 352
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        uspec_size=n
c$$$        end function uspec_size
c$$$       
c$$$        integer function ustr_size()
c$$$        integer:: n
c$$$        n= 46
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        ustr_size=n
c$$$        end function ustr_size
c$$$       
c$$$        integer function utb_size()
c$$$        integer:: n
c$$$        n= 10
c$$$C  +2 is margin
c$$$C         n=n+2
c$$$#ifdef ALIGNMENT128
c$$$        if (mod(n,2).eq.1) n=n+1
c$$$#endif
c$$$        utb_size=n
c$$$        end function utb_size
c$$$
c$$$      subroutine struc_eval_io_r8( v,rw, n,irange1,irange2,val)
c$$$      implicit none
c$$$      real(8) ::  v
c$$$      integer::val
c$$$      integer,intent(in):: n,irange1,irange2
c$$$      character,intent(in):: rw
c$$$      call struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
c$$$      end subroutine struc_eval_io_r8
c$$$
c$$$      subroutine struc_eval_io_r8v( v, rw,n,irange1,irange2,val)
c$$$      implicit none
c$$$      real(8) ::  v(:)
c$$$      integer::val
c$$$      integer,intent(in):: n,irange1,irange2
c$$$      character,intent(in):: rw
c$$$      call struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
c$$$      end subroutine struc_eval_io_r8v
c$$$
c$$$      subroutine struc_eval_io_i8( v, rw,n,irange1,irange2,val)
c$$$      implicit none
c$$$      integer(8) ::  v
c$$$      integer ::val
c$$$      integer,intent(in):: n,irange1,irange2
c$$$      character,intent(in):: rw
c$$$      call struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
c$$$      end subroutine struc_eval_io_i8
c$$$
c$$$      subroutine struc_eval_io_i8v( v, rw,n,irange1,irange2,val)
c$$$      implicit none
c$$$      integer(8) ::  v(:)
c$$$      integer ::val
c$$$      integer,intent(in):: n,irange1,irange2
c$$$      character,intent(in):: rw
c$$$      call struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
c$$$      end subroutine struc_eval_io_i8v
c$$$
c$$$      subroutine struc_strtok(str0, mlist,  nlist,ilist )
c$$$      character(*),intent(in):: str0
c$$$      integer,intent(in):: mlist
c$$$      integer:: nlist
c$$$      integer:: ilist(mlist,2)
c$$$
c$$$      integer:: i1,i2,n
c$$$      logical:: flag
c$$$
c$$$      character(12):: thisname='struc_strtok'
c$$$
c$$$      ilist=0
c$$$      nlist=0
c$$$
c$$$      n=len_trim(str0) 
c$$$#if STRUC_DEBUG_WRITE
c$$$      write(*,*)thisname,' <',str0(:n),'>'
c$$$#endif
c$$$c     accept 0-9,A-Z,a-z,' '
c$$$      do i1=1,n
c$$$         flag = .false.
c$$$         flag = flag .or.
c$$$     .        ( ichar('0') <= ichar(str0(i1:i1))  .and. ichar(str0(i1:i1)) <= ichar('9') )
c$$$         flag = flag .or.
c$$$     .        ( ichar('A') <= ichar(str0(i1:i1))  .and. ichar(str0(i1:i1)) <= ichar('Z') )
c$$$         flag = flag .or.
c$$$     .        ( ichar('a') <= ichar(str0(i1:i1))  .and. ichar(str0(i1:i1)) <= ichar('z') )
c$$$         flag = flag .or. str0(i1:i1).eq.' ' 
c$$$         if ( .not. flag ) then
c$$$            write(*,*) thisname,' unknown character in <',str0(1:n),'>'
c$$$            stop
c$$$         endif 
c$$$      enddo
c$$$
c$$$      i1=1
c$$$      do while (i1<=n)
c$$$         
c$$$                                ! find non' '
c$$$         do while ( str0(i1:i1) ==' ' .and. i1<=n) ; i1=i1+1; enddo
c$$$
c$$$            i2=i1
c$$$            do while ( str0(i2:i2) .ne.' ' .and. i2<=n ); i2=i2+1; enddo
c$$$               if (i2.eq.n) i2=i2+1
c$$$
c$$$               nlist=nlist+1
c$$$               ilist(nlist,1)=i1 
c$$$               ilist(nlist,2)=i2-1
c$$$
c$$$               if ( nlist>= mlist ) exit;
c$$$
c$$$               i1=i2+1
c$$$
c$$$            enddo 
c$$$
c$$$            end subroutine struc_strtok
c$$$
c$$$      end module m_struc_func
c$$$
c$$$
c$$$
c$$$      subroutine struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
c$$$      implicit none
c$$$      real(8):: v(n), val(n)
c$$$      integer,intent(in):: n,irange1,irange2
c$$$      character,intent(in):: rw
c$$$      integer:: ir1,ir2
c$$$      integer:: i
c$$$
c$$$      ir1=irange1
c$$$      ir2=irange2
c$$$
c$$$      if (ir1.eq.-1) then
c$$$         if (rw.eq.'u') then
c$$$            val=v
c$$$         else if (rw.eq.'p') then
c$$$            v=val
c$$$         endif
c$$$      else 
c$$$         if (ir2.eq.-1)  ir2=n 
c$$$         do i=ir1,ir2
c$$$            if (rw.eq.'u') then
c$$$               val(i-ir1+1)=v(i)
c$$$            else if (rw.eq.'p') then
c$$$               v(i)=val(i-ir1+1)
c$$$            endif
c$$$         enddo
c$$$      endif
c$$$      end subroutine struc_eval_io_r8_realbody
c$$$
c$$$
c$$$      subroutine struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
c$$$      implicit none
c$$$      integer(8):: v(n)
c$$$      integer:: val(n)
c$$$      integer,intent(in):: n,irange1,irange2
c$$$      character,intent(in):: rw
c$$$      integer:: ir1,ir2
c$$$      integer:: i
c$$$
c$$$      ir1=irange1
c$$$      ir2=irange2
c$$$
c$$$      if (ir1.eq.-1) then
c$$$         if (rw.eq.'u') then
c$$$            val=v
c$$$         else if (rw.eq.'p') then
c$$$            v=val
c$$$         endif
c$$$      else
c$$$         if (ir2.eq.-1)  ir2=n
c$$$         do i=ir1,ir2
c$$$            if (rw.eq.'u') then
c$$$               val(i-ir1+1)=v(i)
c$$$            else if (rw.eq.'p') then
c$$$               v(i)=val(i-ir1+1)
c$$$            endif
c$$$         enddo
c$$$      endif
c$$$
c$$$      end subroutine struc_eval_io_i8_realbody
c$$$
c$$$      
c$$$
