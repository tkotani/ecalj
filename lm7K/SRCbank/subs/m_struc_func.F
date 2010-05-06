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

      

