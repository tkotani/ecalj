!-----------------------------------------------
      module m_alloclist
      implicit none
      integer,parameter :: m_var_total=5000
      character(30):: var_name(m_var_total)
      integer(8):: var_size(m_var_total)=m_var_total*0 
      integer:: var_total=0 
      integer(8):: var_unit=1024*1024 
      end module m_alloclist
!-----------------------------------------------
      subroutine add_alloclist(name,ik,isi)
      use m_mpi, only: mpi__rank
      use m_alloclist
      implicit none
      character(*),intent(in)::name
      integer,intent(in):: ik 
      integer(8),intent(in)::isi

      integer:: i,iadd
      integer(8):: total,isize
      logical:: idone

      isize=int8(ik)*isi
      iadd=0

      idone=.false.
      do i=1,var_total
         if (name.eq.var_name(i) .and. var_size(i).eq.0) then 
            var_size(i)=isize
            idone=.true.
            iadd=i
            exit 
         endif
      enddo 

      if (.not.idone) then 
      var_total=var_total+1
      if (var_total>m_var_total) then
        write(6,*)'add_alloclist can not add a list, increase m_var_total=' , m_var_total 
        return 
      else 
        var_name(var_total)=name
        var_size(var_total)=isize 
        iadd=var_total
      endif
      endif

      total=0
      do i=1, var_total
         total=total+var_size(i)  
      enddo
      write(6,'(a,i4,1x,a,a20,a,i10,a,i10,a,i5)')
     &  'mpi=',mpi__rank,  
     & 'add_alloc_list name=',trim(name),
     & ' size=',isize/var_unit,'(MB) total=',total/var_unit ,'(MB)',iadd
 
      end subroutine add_alloclist
!-----------------------------------------------
      subroutine del_alloclist(name) 
      use m_mpi, only: mpi__rank
      use m_alloclist
      implicit none
      character(*),intent(in):: name
      integer:: i,idel
      integer(8):: total,isize
 
      idel=0
      do i=var_total,1,-1
          if (name.eq.var_name(i)) then
             var_size(i)=0
             idel=i
             exit 
          endif
      enddo 
 
      total=0
      do i= 1,var_total
         total=total+var_size(i)  
      enddo
      write(6,'(a,i4,1x,a,a20,a,i10,a,i1)')
     & 'mpi=',mpi__rank,
     & 'delete_alloc_list name=', trim(name),' total=',total/var_unit ,'(MB)' !,idel

      end subroutine del_alloclist


