
      module m_struc_func

      interface struc_eval_io
      module procedure struc_eval_io_r8, struc_eval_io_i8
     .     , struc_eval_io_r8v, struc_eval_io_i8v
      end interface struc_eval_io 

      contains


        integer function uarray_size()
        integer:: n
        n= 32
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        uarray_size=n
        end function uarray_size
       
        integer function ubz_size()
        integer:: n
        n= 47
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        ubz_size=n
        end function ubz_size
       
        integer function uctrl_size()
        integer:: n
        n= 88
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        uctrl_size=n
        end function uctrl_size
       
        integer function ugw_size()
        integer:: n
        n= 22
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        ugw_size=n
        end function ugw_size
       
        integer function uham_size()
        integer:: n
        n= 129
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        uham_size=n
        end function uham_size
       
        integer function ulat_size()
        integer:: n
        n= 120
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        ulat_size=n
        end function ulat_size
       
        integer function umix_size()
        integer:: n
        n= 36
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        umix_size=n
        end function umix_size
       
        integer function umove_size()
        integer:: n
        n= 34
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        umove_size=n
        end function umove_size
       
        integer function uoptic_size()
        integer:: n
        n= 33
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        uoptic_size=n
        end function uoptic_size
       
        integer function uordn_size()
        integer:: n
        n= 13
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        uordn_size=n
        end function uordn_size
       
        integer function upot_size()
        integer:: n
        n= 63
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        upot_size=n
        end function upot_size
       
        integer function usite_size()
        integer:: n
        n= 86
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        usite_size=n
        end function usite_size
       
        integer function uspec_size()
        integer:: n
        n= 352
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        uspec_size=n
        end function uspec_size
       
        integer function ustr_size()
        integer:: n
        n= 46
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        ustr_size=n
        end function ustr_size
       
        integer function utb_size()
        integer:: n
        n= 10
C  +2 is margin
C         n=n+2
#ifdef ALIGNMENT128
        if (mod(n,2).eq.1) n=n+1
#endif
        utb_size=n
        end function utb_size
       

      subroutine struc_eval_io_r8( v,rw, n,irange1,irange2,val)
      implicit none
      real(8) ::  v
      integer::val
      integer,intent(in):: n,irange1,irange2
      character,intent(in):: rw
      call struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_r8

      subroutine struc_eval_io_r8v( v, rw,n,irange1,irange2,val)
      implicit none
      real(8) ::  v(:)
      integer::val
      integer,intent(in):: n,irange1,irange2
      character,intent(in):: rw
      call struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_r8v

      subroutine struc_eval_io_i8( v, rw,n,irange1,irange2,val)
      implicit none
      integer(8) ::  v
      integer ::val
      integer,intent(in):: n,irange1,irange2
      character,intent(in):: rw
      call struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_i8

      subroutine struc_eval_io_i8v( v, rw,n,irange1,irange2,val)
      implicit none
      integer(8) ::  v(:)
      integer ::val
      integer,intent(in):: n,irange1,irange2
      character,intent(in):: rw
      call struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_i8v

      subroutine struc_strtok(str0, mlist,  nlist,ilist )
      character(*),intent(in):: str0
      integer,intent(in):: mlist
      integer:: nlist
      integer:: ilist(mlist,2)

      integer:: i1,i2,n
      logical:: flag

      character(12):: thisname='struc_strtok'

      ilist=0
      nlist=0

      n=len_trim(str0) 
#if STRUC_DEBUG_WRITE
      write(*,*)thisname,' <',str0(:n),'>'
#endif
c     accept 0-9,A-Z,a-z,' '
      do i1=1,n
         flag = .false.
         flag = flag .or.
     .        ( ichar('0') <= ichar(str0(i1:i1))  .and. ichar(str0(i1:i1)) <= ichar('9') )
         flag = flag .or.
     .        ( ichar('A') <= ichar(str0(i1:i1))  .and. ichar(str0(i1:i1)) <= ichar('Z') )
         flag = flag .or.
     .        ( ichar('a') <= ichar(str0(i1:i1))  .and. ichar(str0(i1:i1)) <= ichar('z') )
         flag = flag .or. str0(i1:i1).eq.' ' 
         if ( .not. flag ) then
            write(*,*) thisname,' unknown character in <',str0(1:n),'>'
            stop
         endif 
      enddo

      i1=1
      do while (i1<=n)
         
                                ! find non' '
         do while ( str0(i1:i1) ==' ' .and. i1<=n) ; i1=i1+1; enddo

            i2=i1
            do while ( str0(i2:i2) .ne.' ' .and. i2<=n ); i2=i2+1; enddo
               if (i2.eq.n) i2=i2+1

               nlist=nlist+1
               ilist(nlist,1)=i1 
               ilist(nlist,2)=i2-1

               if ( nlist>= mlist ) exit;

               i1=i2+1

            enddo 

            end subroutine struc_strtok
      
      end module m_struc_func




      subroutine struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
      implicit none
      real(8):: v(n), val(n)
      integer,intent(in):: n,irange1,irange2
      character,intent(in):: rw
      integer:: ir1,ir2
      integer:: i

      ir1=irange1
      ir2=irange2

      if (ir1.eq.-1) then
         if (rw.eq.'u') then
            val=v
         else if (rw.eq.'p') then
            v=val
         endif
      else 
         if (ir2.eq.-1)  ir2=n 
         do i=ir1,ir2
            if (rw.eq.'u') then
               val(i-ir1+1)=v(i)
            else if (rw.eq.'p') then
               v(i)=val(i-ir1+1)
            endif
         enddo
      endif
      end subroutine struc_eval_io_r8_realbody


      subroutine struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
      implicit none
      integer(8):: v(n)
      integer:: val(n)
      integer,intent(in):: n,irange1,irange2
      character,intent(in):: rw
      integer:: ir1,ir2
      integer:: i

      ir1=irange1
      ir2=irange2

      if (ir1.eq.-1) then
         if (rw.eq.'u') then
            val=v
         else if (rw.eq.'p') then
            v=val
         endif
      else
         if (ir2.eq.-1)  ir2=n
         do i=ir1,ir2
            if (rw.eq.'u') then
               val(i-ir1+1)=v(i)
            else if (rw.eq.'p') then
               v(i)=val(i-ir1+1)
            endif
         enddo
      endif

      end subroutine struc_eval_io_i8_realbody

      

