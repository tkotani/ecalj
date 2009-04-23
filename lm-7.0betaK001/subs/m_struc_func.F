
      module m_struc_func

      interface struc_eval_io
      module procedure struc_eval_io_r8, struc_eval_io_i8
     .     , struc_eval_io_r8v, struc_eval_io_i8v
      end interface struc_eval_io 

      contains

      subroutine struc_eval_io_r8( v,rw, n,irange1,irange2,val)
      implicit none
      real(8) ::  v
      integer::val
      integer,intent(in):: n,irange1,irange2
      character:: rw
      call struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_r8

      subroutine struc_eval_io_r8v( v, rw,n,irange1,irange2,val)
      implicit none
      real(8) ::  v(:)
      integer::val
      integer,intent(in):: n,irange1,irange2
      character:: rw
      call struc_eval_io_r8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_r8v

      subroutine struc_eval_io_i8( v, rw,n,irange1,irange2,val)
      implicit none
      integer(8) ::  v
      integer ::val
      integer,intent(in):: n,irange1,irange2
      character:: rw
      call struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_i8

      subroutine struc_eval_io_i8v( v, rw,n,irange1,irange2,val)
      implicit none
      integer(8) ::  v(:)
      integer ::val
      integer,intent(in):: n,irange1,irange2
      character:: rw
      call struc_eval_io_i8_realbody(v,rw,n,irange1,irange2,val)
      end subroutine struc_eval_io_i8v

      subroutine struc_strtok(str0, mlist,  nlist,ilist )
      character(*):: str0
      integer:: mlist
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
      character:: rw
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
      character:: rw
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

      

