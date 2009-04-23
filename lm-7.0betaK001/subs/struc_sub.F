C ubz 30
       
C offe=-1 is only supported
        subroutine ubz(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_bz):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'ubz offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ubz
       
        subroutine struc_ubz_io(v_bz,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_bz):: v_bz
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call ubz(v_bz,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_bz%size,rw,1,range1,range2,val)
          endif
        case ('def')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%def,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%def)/kind(v_bz%def)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%def)/kind(v_bz%def)
            call struc_eval_io( v_bz%def,rw,n,range1,range2,val)
          endif
        case ('dosw')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%dosw,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%dosw)/kind(v_bz%dosw)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%dosw)/kind(v_bz%dosw)
            call struc_eval_io( v_bz%dosw,rw,n,range1,range2,val)
          endif
        case ('ef')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%ef,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%ef)/kind(v_bz%ef)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%ef)/kind(v_bz%ef)
            call struc_eval_io( v_bz%ef,rw,n,range1,range2,val)
          endif
        case ('efmax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%efmax,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%efmax)/kind(v_bz%efmax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%efmax)/kind(v_bz%efmax)
            call struc_eval_io( v_bz%efmax,rw,n,range1,range2,val)
          endif
        case ('fsmom')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%fsmom,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%fsmom)/kind(v_bz%fsmom)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%fsmom)/kind(v_bz%fsmom)
            call struc_eval_io( v_bz%fsmom,rw,n,range1,range2,val)
          endif
        case ('lcond')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%lcond,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%lcond)/kind(v_bz%lcond)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%lcond)/kind(v_bz%lcond)
            call struc_eval_io( v_bz%lcond,rw,n,range1,range2,val)
          endif
        case ('lio')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%lio,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%lio)/kind(v_bz%lio)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%lio)/kind(v_bz%lio)
            call struc_eval_io( v_bz%lio,rw,n,range1,range2,val)
          endif
        case ('lmet')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%lmet,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%lmet)/kind(v_bz%lmet)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%lmet)/kind(v_bz%lmet)
            call struc_eval_io( v_bz%lmet,rw,n,range1,range2,val)
          endif
        case ('lmull')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%lmull,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%lmull)/kind(v_bz%lmull)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%lmull)/kind(v_bz%lmull)
            call struc_eval_io( v_bz%lmull,rw,n,range1,range2,val)
          endif
        case ('lopt')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%lopt,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%lopt)/kind(v_bz%lopt)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%lopt)/kind(v_bz%lopt)
            call struc_eval_io( v_bz%lopt,rw,n,range1,range2,val)
          endif
        case ('lshft')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%lshft,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%lshft)/kind(v_bz%lshft)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%lshft)/kind(v_bz%lshft)
            call struc_eval_io( v_bz%lshft,rw,n,range1,range2,val)
          endif
        case ('n')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%n,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%n)/kind(v_bz%n)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%n)/kind(v_bz%n)
            call struc_eval_io( v_bz%n,rw,n,range1,range2,val)
          endif
        case ('ndos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%ndos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%ndos)/kind(v_bz%ndos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%ndos)/kind(v_bz%ndos)
            call struc_eval_io( v_bz%ndos,rw,n,range1,range2,val)
          endif
        case ('nevmx')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%nevmx,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%nevmx)/kind(v_bz%nevmx)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%nevmx)/kind(v_bz%nevmx)
            call struc_eval_io( v_bz%nevmx,rw,n,range1,range2,val)
          endif
        case ('nkabc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%nkabc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%nkabc)/kind(v_bz%nkabc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%nkabc)/kind(v_bz%nkabc)
            call struc_eval_io( v_bz%nkabc,rw,n,range1,range2,val)
          endif
        case ('nkp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%nkp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%nkp)/kind(v_bz%nkp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%nkp)/kind(v_bz%nkp)
            call struc_eval_io( v_bz%nkp,rw,n,range1,range2,val)
          endif
        case ('ntet')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%ntet,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%ntet)/kind(v_bz%ntet)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%ntet)/kind(v_bz%ntet)
            call struc_eval_io( v_bz%ntet,rw,n,range1,range2,val)
          endif
        case ('odos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%odos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%odos)/kind(v_bz%odos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%odos)/kind(v_bz%odos)
            call struc_eval_io( v_bz%odos,rw,n,range1,range2,val)
          endif
        case ('oidtet')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%oidtet,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%oidtet)/kind(v_bz%oidtet)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%oidtet)/kind(v_bz%oidtet)
            call struc_eval_io( v_bz%oidtet,rw,n,range1,range2,val)
          endif
        case ('oipq')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%oipq,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%oipq)/kind(v_bz%oipq)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%oipq)/kind(v_bz%oipq)
            call struc_eval_io( v_bz%oipq,rw,n,range1,range2,val)
          endif
        case ('opdos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%opdos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%opdos)/kind(v_bz%opdos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%opdos)/kind(v_bz%opdos)
            call struc_eval_io( v_bz%opdos,rw,n,range1,range2,val)
          endif
        case ('oqp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%oqp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%oqp)/kind(v_bz%oqp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%oqp)/kind(v_bz%oqp)
            call struc_eval_io( v_bz%oqp,rw,n,range1,range2,val)
          endif
        case ('ostar')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%ostar,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%ostar)/kind(v_bz%ostar)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%ostar)/kind(v_bz%ostar)
            call struc_eval_io( v_bz%ostar,rw,n,range1,range2,val)
          endif
        case ('owtkp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_bz%owtkp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_bz%owtkp)/kind(v_bz%owtkp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%owtkp)/kind(v_bz%owtkp)
            call struc_eval_io( v_bz%owtkp,rw,n,range1,range2,val)
          endif
        case ('range')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%range,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%range)/kind(v_bz%range)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%range)/kind(v_bz%range)
            call struc_eval_io( v_bz%range,rw,n,range1,range2,val)
          endif
        case ('semsh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%semsh,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%semsh)/kind(v_bz%semsh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%semsh)/kind(v_bz%semsh)
            call struc_eval_io( v_bz%semsh,rw,n,range1,range2,val)
          endif
        case ('stnr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%stnr,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%stnr)/kind(v_bz%stnr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%stnr)/kind(v_bz%stnr)
            call struc_eval_io( v_bz%stnr,rw,n,range1,range2,val)
          endif
        case ('w')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%w,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%w)/kind(v_bz%w)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%w)/kind(v_bz%w)
            call struc_eval_io( v_bz%w,rw,n,range1,range2,val)
          endif
        case ('zval')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_bz%zval,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_bz%zval)/kind(v_bz%zval)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_bz%zval)/kind(v_bz%zval)
            call struc_eval_io( v_bz%zval,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ubz_io: failed to find ',name
          stop
        end select
        end subroutine struc_ubz_io
       
C uarray 32
       
C offe=-1 is only supported
        subroutine uarray(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_array):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'uarray offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uarray
       
        subroutine struc_uarray_io(v_array,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_array):: v_array
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call uarray(v_array,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_array%size,rw,1,range1,range2,val)
          endif
        case ('ncl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ncl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ncl)/kind(v_array%ncl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ncl)/kind(v_array%ncl)
            call struc_eval_io( v_array%ncl,rw,n,range1,range2,val)
          endif
        case ('nclasp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%nclasp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%nclasp)/kind(v_array%nclasp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%nclasp)/kind(v_array%nclasp)
            call struc_eval_io( v_array%nclasp,rw,n,range1,range2,val)
          endif
        case ('nofgl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%nofgl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%nofgl)/kind(v_array%nofgl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%nofgl)/kind(v_array%nofgl)
            call struc_eval_io( v_array%nofgl,rw,n,range1,range2,val)
          endif
        case ('nofgr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%nofgr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%nofgr)/kind(v_array%nofgr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%nofgr)/kind(v_array%nofgr)
            call struc_eval_io( v_array%nofgr,rw,n,range1,range2,val)
          endif
        case ('npadl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%npadl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%npadl)/kind(v_array%npadl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%npadl)/kind(v_array%npadl)
            call struc_eval_io( v_array%npadl,rw,n,range1,range2,val)
          endif
        case ('npadr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%npadr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%npadr)/kind(v_array%npadr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%npadr)/kind(v_array%npadr)
            call struc_eval_io( v_array%npadr,rw,n,range1,range2,val)
          endif
        case ('npl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%npl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%npl)/kind(v_array%npl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%npl)/kind(v_array%npl)
            call struc_eval_io( v_array%npl,rw,n,range1,range2,val)
          endif
        case ('oclabl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oclabl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oclabl)/kind(v_array%oclabl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oclabl)/kind(v_array%oclabl)
            call struc_eval_io( v_array%oclabl,rw,n,range1,range2,val)
          endif
        case ('ocllst')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ocllst,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ocllst)/kind(v_array%ocllst)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ocllst)/kind(v_array%ocllst)
            call struc_eval_io( v_array%ocllst,rw,n,range1,range2,val)
          endif
        case ('oclp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oclp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oclp)/kind(v_array%oclp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oclp)/kind(v_array%oclp)
            call struc_eval_io( v_array%oclp,rw,n,range1,range2,val)
          endif
        case ('oclssl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oclssl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oclssl)/kind(v_array%oclssl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oclssl)/kind(v_array%oclssl)
            call struc_eval_io( v_array%oclssl,rw,n,range1,range2,val)
          endif
        case ('ogroup')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ogroup,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ogroup)/kind(v_array%ogroup)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ogroup)/kind(v_array%ogroup)
            call struc_eval_io( v_array%ogroup,rw,n,range1,range2,val)
          endif
        case ('ogrp2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ogrp2,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ogrp2)/kind(v_array%ogrp2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ogrp2)/kind(v_array%ogrp2)
            call struc_eval_io( v_array%ogrp2,rw,n,range1,range2,val)
          endif
        case ('ohave')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ohave,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ohave)/kind(v_array%ohave)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ohave)/kind(v_array%ohave)
            call struc_eval_io( v_array%ohave,rw,n,range1,range2,val)
          endif
        case ('oics')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oics,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oics)/kind(v_array%oics)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oics)/kind(v_array%oics)
            call struc_eval_io( v_array%oics,rw,n,range1,range2,val)
          endif
        case ('oifrlx')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oifrlx,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oifrlx)/kind(v_array%oifrlx)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oifrlx)/kind(v_array%oifrlx)
            call struc_eval_io( v_array%oifrlx,rw,n,range1,range2,val)
          endif
        case ('oipc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oipc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oipc)/kind(v_array%oipc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oipc)/kind(v_array%oipc)
            call struc_eval_io( v_array%oipc,rw,n,range1,range2,val)
          endif
        case ('oipcp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oipcp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oipcp)/kind(v_array%oipcp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oipcp)/kind(v_array%oipcp)
            call struc_eval_io( v_array%oipcp,rw,n,range1,range2,val)
          endif
        case ('oips')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%oips,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%oips)/kind(v_array%oips)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%oips)/kind(v_array%oips)
            call struc_eval_io( v_array%oips,rw,n,range1,range2,val)
          endif
        case ('omxcst')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%omxcst,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%omxcst)/kind(v_array%omxcst)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%omxcst)/kind(v_array%omxcst)
            call struc_eval_io( v_array%omxcst,rw,n,range1,range2,val)
          endif
        case ('onrc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%onrc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%onrc)/kind(v_array%onrc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%onrc)/kind(v_array%onrc)
            call struc_eval_io( v_array%onrc,rw,n,range1,range2,val)
          endif
        case ('onrcp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%onrcp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%onrcp)/kind(v_array%onrcp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%onrcp)/kind(v_array%onrcp)
            call struc_eval_io( v_array%onrcp,rw,n,range1,range2,val)
          endif
        case ('opgfsl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%opgfsl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%opgfsl)/kind(v_array%opgfsl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%opgfsl)/kind(v_array%opgfsl)
            call struc_eval_io( v_array%opgfsl,rw,n,range1,range2,val)
          endif
        case ('opgfvl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%opgfvl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%opgfvl)/kind(v_array%opgfvl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%opgfvl)/kind(v_array%opgfvl)
            call struc_eval_io( v_array%opgfvl,rw,n,range1,range2,val)
          endif
        case ('opgord')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%opgord,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%opgord)/kind(v_array%opgord)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%opgord)/kind(v_array%opgord)
            call struc_eval_io( v_array%opgord,rw,n,range1,range2,val)
          endif
        case ('opgplp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%opgplp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%opgplp)/kind(v_array%opgplp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%opgplp)/kind(v_array%opgplp)
            call struc_eval_io( v_array%opgplp,rw,n,range1,range2,val)
          endif
        case ('opos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%opos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%opos)/kind(v_array%opos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%opos)/kind(v_array%opos)
            call struc_eval_io( v_array%opos,rw,n,range1,range2,val)
          endif
        case ('ormax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ormax,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ormax)/kind(v_array%ormax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ormax)/kind(v_array%ormax)
            call struc_eval_io( v_array%ormax,rw,n,range1,range2,val)
          endif
        case ('osordn')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%osordn,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%osordn)/kind(v_array%osordn)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%osordn)/kind(v_array%osordn)
            call struc_eval_io( v_array%osordn,rw,n,range1,range2,val)
          endif
        case ('ostni')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ostni,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ostni)/kind(v_array%ostni)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ostni)/kind(v_array%ostni)
            call struc_eval_io( v_array%ostni,rw,n,range1,range2,val)
          endif
        case ('ozos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_array%ozos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_array%ozos)/kind(v_array%ozos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_array%ozos)/kind(v_array%ozos)
            call struc_eval_io( v_array%ozos,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uarray_io: failed to find ',name
          stop
        end select
        end subroutine struc_uarray_io
       
C uctrl 56
       
C offe=-1 is only supported
        subroutine uctrl(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_ctrl):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'uctrl offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uctrl
       
        subroutine struc_uctrl_io(v_ctrl,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_ctrl):: v_ctrl
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call uctrl(v_ctrl,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_ctrl%size,rw,1,range1,range2,val)
          endif
        case ('defm')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%defm,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%defm)/kind(v_ctrl%defm)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%defm)/kind(v_ctrl%defm)
            call struc_eval_io( v_ctrl%defm,rw,n,range1,range2,val)
          endif
        case ('elin')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%elin,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%elin)/kind(v_ctrl%elin)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%elin)/kind(v_ctrl%elin)
            call struc_eval_io( v_ctrl%elin,rw,n,range1,range2,val)
          endif
        case ('lasa')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lasa,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lasa)/kind(v_ctrl%lasa)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lasa)/kind(v_ctrl%lasa)
            call struc_eval_io( v_ctrl%lasa,rw,n,range1,range2,val)
          endif
        case ('lbas')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lbas,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lbas)/kind(v_ctrl%lbas)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lbas)/kind(v_ctrl%lbas)
            call struc_eval_io( v_ctrl%lbas,rw,n,range1,range2,val)
          endif
        case ('lcd')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lcd,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lcd)/kind(v_ctrl%lcd)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lcd)/kind(v_ctrl%lcd)
            call struc_eval_io( v_ctrl%lcd,rw,n,range1,range2,val)
          endif
        case ('lcgf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lcgf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lcgf)/kind(v_ctrl%lcgf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lcgf)/kind(v_ctrl%lcgf)
            call struc_eval_io( v_ctrl%lcgf,rw,n,range1,range2,val)
          endif
        case ('ldos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%ldos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%ldos)/kind(v_ctrl%ldos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%ldos)/kind(v_ctrl%ldos)
            call struc_eval_io( v_ctrl%ldos,rw,n,range1,range2,val)
          endif
        case ('lfp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lfp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lfp)/kind(v_ctrl%lfp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lfp)/kind(v_ctrl%lfp)
            call struc_eval_io( v_ctrl%lfp,rw,n,range1,range2,val)
          endif
        case ('lfrce')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lfrce,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lfrce)/kind(v_ctrl%lfrce)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lfrce)/kind(v_ctrl%lfrce)
            call struc_eval_io( v_ctrl%lfrce,rw,n,range1,range2,val)
          endif
        case ('lgen3')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lgen3,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lgen3)/kind(v_ctrl%lgen3)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lgen3)/kind(v_ctrl%lgen3)
            call struc_eval_io( v_ctrl%lgen3,rw,n,range1,range2,val)
          endif
        case ('lham')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lham,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lham)/kind(v_ctrl%lham)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lham)/kind(v_ctrl%lham)
            call struc_eval_io( v_ctrl%lham,rw,n,range1,range2,val)
          endif
        case ('lmet')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lmet,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lmet)/kind(v_ctrl%lmet)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lmet)/kind(v_ctrl%lmet)
            call struc_eval_io( v_ctrl%lmet,rw,n,range1,range2,val)
          endif
        case ('lncol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lncol,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lncol)/kind(v_ctrl%lncol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lncol)/kind(v_ctrl%lncol)
            call struc_eval_io( v_ctrl%lncol,rw,n,range1,range2,val)
          endif
        case ('loptc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%loptc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%loptc)/kind(v_ctrl%loptc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%loptc)/kind(v_ctrl%loptc)
            call struc_eval_io( v_ctrl%loptc,rw,n,range1,range2,val)
          endif
        case ('lordn')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lordn,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lordn)/kind(v_ctrl%lordn)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lordn)/kind(v_ctrl%lordn)
            call struc_eval_io( v_ctrl%lordn,rw,n,range1,range2,val)
          endif
        case ('lpgf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lpgf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lpgf)/kind(v_ctrl%lpgf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lpgf)/kind(v_ctrl%lpgf)
            call struc_eval_io( v_ctrl%lpgf,rw,n,range1,range2,val)
          endif
        case ('lqp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lqp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lqp)/kind(v_ctrl%lqp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lqp)/kind(v_ctrl%lqp)
            call struc_eval_io( v_ctrl%lqp,rw,n,range1,range2,val)
          endif
        case ('lrel')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lrel,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lrel)/kind(v_ctrl%lrel)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lrel)/kind(v_ctrl%lrel)
            call struc_eval_io( v_ctrl%lrel,rw,n,range1,range2,val)
          endif
        case ('lrs')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lrs,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lrs)/kind(v_ctrl%lrs)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lrs)/kind(v_ctrl%lrs)
            call struc_eval_io( v_ctrl%lrs,rw,n,range1,range2,val)
          endif
        case ('lscr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lscr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lscr)/kind(v_ctrl%lscr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lscr)/kind(v_ctrl%lscr)
            call struc_eval_io( v_ctrl%lscr,rw,n,range1,range2,val)
          endif
        case ('lstonr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lstonr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lstonr)/kind(v_ctrl%lstonr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lstonr)/kind(v_ctrl%lstonr)
            call struc_eval_io( v_ctrl%lstonr,rw,n,range1,range2,val)
          endif
        case ('lstr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lstr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lstr)/kind(v_ctrl%lstr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lstr)/kind(v_ctrl%lstr)
            call struc_eval_io( v_ctrl%lstr,rw,n,range1,range2,val)
          endif
        case ('lsx')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lsx,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lsx)/kind(v_ctrl%lsx)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lsx)/kind(v_ctrl%lsx)
            call struc_eval_io( v_ctrl%lsx,rw,n,range1,range2,val)
          endif
        case ('ltb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%ltb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%ltb)/kind(v_ctrl%ltb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%ltb)/kind(v_ctrl%ltb)
            call struc_eval_io( v_ctrl%ltb,rw,n,range1,range2,val)
          endif
        case ('lves')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lves,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lves)/kind(v_ctrl%lves)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lves)/kind(v_ctrl%lves)
            call struc_eval_io( v_ctrl%lves,rw,n,range1,range2,val)
          endif
        case ('lxcf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%lxcf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%lxcf)/kind(v_ctrl%lxcf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%lxcf)/kind(v_ctrl%lxcf)
            call struc_eval_io( v_ctrl%lxcf,rw,n,range1,range2,val)
          endif
        case ('maxit')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%maxit,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%maxit)/kind(v_ctrl%maxit)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%maxit)/kind(v_ctrl%maxit)
            call struc_eval_io( v_ctrl%maxit,rw,n,range1,range2,val)
          endif
        case ('mdprm')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%mdprm,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%mdprm)/kind(v_ctrl%mdprm)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%mdprm)/kind(v_ctrl%mdprm)
            call struc_eval_io( v_ctrl%mdprm,rw,n,range1,range2,val)
          endif
        case ('modep')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%modep,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%modep)/kind(v_ctrl%modep)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%modep)/kind(v_ctrl%modep)
            call struc_eval_io( v_ctrl%modep,rw,n,range1,range2,val)
          endif
        case ('nbas')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nbas,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nbas)/kind(v_ctrl%nbas)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nbas)/kind(v_ctrl%nbas)
            call struc_eval_io( v_ctrl%nbas,rw,n,range1,range2,val)
          endif
        case ('nbasp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nbasp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nbasp)/kind(v_ctrl%nbasp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nbasp)/kind(v_ctrl%nbasp)
            call struc_eval_io( v_ctrl%nbasp,rw,n,range1,range2,val)
          endif
        case ('nclass')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nclass,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nclass)/kind(v_ctrl%nclass)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nclass)/kind(v_ctrl%nclass)
            call struc_eval_io( v_ctrl%nclass,rw,n,range1,range2,val)
          endif
        case ('nesabc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nesabc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nesabc)/kind(v_ctrl%nesabc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nesabc)/kind(v_ctrl%nesabc)
            call struc_eval_io( v_ctrl%nesabc,rw,n,range1,range2,val)
          endif
        case ('nitmv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nitmv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nitmv)/kind(v_ctrl%nitmv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nitmv)/kind(v_ctrl%nitmv)
            call struc_eval_io( v_ctrl%nitmv,rw,n,range1,range2,val)
          endif
        case ('nl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nl)/kind(v_ctrl%nl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nl)/kind(v_ctrl%nl)
            call struc_eval_io( v_ctrl%nl,rw,n,range1,range2,val)
          endif
        case ('nmap')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nmap,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nmap)/kind(v_ctrl%nmap)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nmap)/kind(v_ctrl%nmap)
            call struc_eval_io( v_ctrl%nmap,rw,n,range1,range2,val)
          endif
        case ('npl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%npl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%npl)/kind(v_ctrl%npl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%npl)/kind(v_ctrl%npl)
            call struc_eval_io( v_ctrl%npl,rw,n,range1,range2,val)
          endif
        case ('nsite')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nsite,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nsite)/kind(v_ctrl%nsite)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nsite)/kind(v_ctrl%nsite)
            call struc_eval_io( v_ctrl%nsite,rw,n,range1,range2,val)
          endif
        case ('nspec')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nspec,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nspec)/kind(v_ctrl%nspec)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nspec)/kind(v_ctrl%nspec)
            call struc_eval_io( v_ctrl%nspec,rw,n,range1,range2,val)
          endif
        case ('nspin')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nspin,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nspin)/kind(v_ctrl%nspin)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nspin)/kind(v_ctrl%nspin)
            call struc_eval_io( v_ctrl%nspin,rw,n,range1,range2,val)
          endif
        case ('nvario')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%nvario,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%nvario)/kind(v_ctrl%nvario)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%nvario)/kind(v_ctrl%nvario)
            call struc_eval_io( v_ctrl%nvario,rw,n,range1,range2,val)
          endif
        case ('omax1')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%omax1,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%omax1)/kind(v_ctrl%omax1)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%omax1)/kind(v_ctrl%omax1)
            call struc_eval_io( v_ctrl%omax1,rw,n,range1,range2,val)
          endif
        case ('omax2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%omax2,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%omax2)/kind(v_ctrl%omax2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%omax2)/kind(v_ctrl%omax2)
            call struc_eval_io( v_ctrl%omax2,rw,n,range1,range2,val)
          endif
        case ('pfloat')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%pfloat,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%pfloat)/kind(v_ctrl%pfloat)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%pfloat)/kind(v_ctrl%pfloat)
            call struc_eval_io( v_ctrl%pfloat,rw,n,range1,range2,val)
          endif
        case ('quit')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%quit,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%quit)/kind(v_ctrl%quit)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%quit)/kind(v_ctrl%quit)
            call struc_eval_io( v_ctrl%quit,rw,n,range1,range2,val)
          endif
        case ('rmaxes')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%rmaxes,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%rmaxes)/kind(v_ctrl%rmaxes)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%rmaxes)/kind(v_ctrl%rmaxes)
            call struc_eval_io( v_ctrl%rmaxes,rw,n,range1,range2,val)
          endif
        case ('rmines')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%rmines,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%rmines)/kind(v_ctrl%rmines)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%rmines)/kind(v_ctrl%rmines)
            call struc_eval_io( v_ctrl%rmines,rw,n,range1,range2,val)
          endif
        case ('sclwsr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%sclwsr,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%sclwsr)/kind(v_ctrl%sclwsr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%sclwsr)/kind(v_ctrl%sclwsr)
            call struc_eval_io( v_ctrl%sclwsr,rw,n,range1,range2,val)
          endif
        case ('sdmod')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%sdmod,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%sdmod)/kind(v_ctrl%sdmod)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%sdmod)/kind(v_ctrl%sdmod)
            call struc_eval_io( v_ctrl%sdmod,rw,n,range1,range2,val)
          endif
        case ('sdprm')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%sdprm,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%sdprm)/kind(v_ctrl%sdprm)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%sdprm)/kind(v_ctrl%sdprm)
            call struc_eval_io( v_ctrl%sdprm,rw,n,range1,range2,val)
          endif
        case ('sdxsi')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%sdxsi,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%sdxsi)/kind(v_ctrl%sdxsi)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%sdxsi)/kind(v_ctrl%sdxsi)
            call struc_eval_io( v_ctrl%sdxsi,rw,n,range1,range2,val)
          endif
        case ('smalit')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ctrl%smalit,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ctrl%smalit)/kind(v_ctrl%smalit)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%smalit)/kind(v_ctrl%smalit)
            call struc_eval_io( v_ctrl%smalit,rw,n,range1,range2,val)
          endif
        case ('tol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%tol,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%tol)/kind(v_ctrl%tol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%tol)/kind(v_ctrl%tol)
            call struc_eval_io( v_ctrl%tol,rw,n,range1,range2,val)
          endif
        case ('wsrmax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%wsrmax,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%wsrmax)/kind(v_ctrl%wsrmax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%wsrmax)/kind(v_ctrl%wsrmax)
            call struc_eval_io( v_ctrl%wsrmax,rw,n,range1,range2,val)
          endif
        case ('zbak')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ctrl%zbak,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ctrl%zbak)/kind(v_ctrl%zbak)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ctrl%zbak)/kind(v_ctrl%zbak)
            call struc_eval_io( v_ctrl%zbak,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uctrl_io: failed to find ',name
          stop
        end select
        end subroutine struc_uctrl_io
       
C ugw 19
       
C offe=-1 is only supported
        subroutine ugw(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_gw):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'ugw offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ugw
       
        subroutine struc_ugw_io(v_gw,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_gw):: v_gw
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call ugw(v_gw,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_gw%size,rw,1,range1,range2,val)
          endif
        case ('delre')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%delre,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%delre)/kind(v_gw%delre)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%delre)/kind(v_gw%delre)
            call struc_eval_io( v_gw%delre,rw,n,range1,range2,val)
          endif
        case ('deltaw')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%deltaw,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%deltaw)/kind(v_gw%deltaw)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%deltaw)/kind(v_gw%deltaw)
            call struc_eval_io( v_gw%deltaw,rw,n,range1,range2,val)
          endif
        case ('deltax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%deltax,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%deltax)/kind(v_gw%deltax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%deltax)/kind(v_gw%deltax)
            call struc_eval_io( v_gw%deltax,rw,n,range1,range2,val)
          endif
        case ('ecuts')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%ecuts,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%ecuts)/kind(v_gw%ecuts)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%ecuts)/kind(v_gw%ecuts)
            call struc_eval_io( v_gw%ecuts,rw,n,range1,range2,val)
          endif
        case ('gcutb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%gcutb,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%gcutb)/kind(v_gw%gcutb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%gcutb)/kind(v_gw%gcutb)
            call struc_eval_io( v_gw%gcutb,rw,n,range1,range2,val)
          endif
        case ('gcutx')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%gcutx,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%gcutx)/kind(v_gw%gcutx)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%gcutx)/kind(v_gw%gcutx)
            call struc_eval_io( v_gw%gcutx,rw,n,range1,range2,val)
          endif
        case ('gsmear')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%gsmear,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%gsmear)/kind(v_gw%gsmear)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%gsmear)/kind(v_gw%gsmear)
            call struc_eval_io( v_gw%gsmear,rw,n,range1,range2,val)
          endif
        case ('lgw')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_gw%lgw,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_gw%lgw)/kind(v_gw%lgw)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%lgw)/kind(v_gw%lgw)
            call struc_eval_io( v_gw%lgw,rw,n,range1,range2,val)
          endif
        case ('lmxpb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_gw%lmxpb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_gw%lmxpb)/kind(v_gw%lmxpb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%lmxpb)/kind(v_gw%lmxpb)
            call struc_eval_io( v_gw%lmxpb,rw,n,range1,range2,val)
          endif
        case ('mksig')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_gw%mksig,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_gw%mksig)/kind(v_gw%mksig)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%mksig)/kind(v_gw%mksig)
            call struc_eval_io( v_gw%mksig,rw,n,range1,range2,val)
          endif
        case ('nband')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_gw%nband,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_gw%nband)/kind(v_gw%nband)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%nband)/kind(v_gw%nband)
            call struc_eval_io( v_gw%nband,rw,n,range1,range2,val)
          endif
        case ('nime')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_gw%nime,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_gw%nime)/kind(v_gw%nime)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%nime)/kind(v_gw%nime)
            call struc_eval_io( v_gw%nime,rw,n,range1,range2,val)
          endif
        case ('nkabc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_gw%nkabc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_gw%nkabc)/kind(v_gw%nkabc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%nkabc)/kind(v_gw%nkabc)
            call struc_eval_io( v_gw%nkabc,rw,n,range1,range2,val)
          endif
        case ('pb1')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%pb1,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_gw%pb1)/kind(v_gw%pb1)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%pb1)/kind(v_gw%pb1)
            call struc_eval_io( v_gw%pb1,rw,n,range1,range2,val)
          endif
        case ('pb2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%pb2,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_gw%pb2)/kind(v_gw%pb2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%pb2)/kind(v_gw%pb2)
            call struc_eval_io( v_gw%pb2,rw,n,range1,range2,val)
          endif
        case ('pbtol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%pbtol,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%pbtol)/kind(v_gw%pbtol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%pbtol)/kind(v_gw%pbtol)
            call struc_eval_io( v_gw%pbtol,rw,n,range1,range2,val)
          endif
        case ('qoffp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_gw%qoffp,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_gw%qoffp)/kind(v_gw%qoffp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%qoffp)/kind(v_gw%qoffp)
            call struc_eval_io( v_gw%qoffp,rw,n,range1,range2,val)
          endif
        case ('rdsig')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_gw%rdsig,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_gw%rdsig)/kind(v_gw%rdsig)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_gw%rdsig)/kind(v_gw%rdsig)
            call struc_eval_io( v_gw%rdsig,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ugw_io: failed to find ',name
          stop
        end select
        end subroutine struc_ugw_io
       
C uham 56
       
C offe=-1 is only supported
        subroutine uham(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_ham):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'uham offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uham
       
        subroutine struc_uham_io(v_ham,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_ham):: v_ham
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call uham(v_ham,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_ham%size,rw,1,range1,range2,val)
          endif
        case ('alfsi')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%alfsi,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%alfsi)/kind(v_ham%alfsi)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%alfsi)/kind(v_ham%alfsi)
            call struc_eval_io( v_ham%alfsi,rw,n,range1,range2,val)
          endif
        case ('amgm')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%amgm,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%amgm)/kind(v_ham%amgm)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%amgm)/kind(v_ham%amgm)
            call struc_eval_io( v_ham%amgm,rw,n,range1,range2,val)
          endif
        case ('bandw')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%bandw,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%bandw)/kind(v_ham%bandw)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%bandw)/kind(v_ham%bandw)
            call struc_eval_io( v_ham%bandw,rw,n,range1,range2,val)
          endif
        case ('dabc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%dabc,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%dabc)/kind(v_ham%dabc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%dabc)/kind(v_ham%dabc)
            call struc_eval_io( v_ham%dabc,rw,n,range1,range2,val)
          endif
        case ('ehf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%ehf,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%ehf)/kind(v_ham%ehf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ehf)/kind(v_ham%ehf)
            call struc_eval_io( v_ham%ehf,rw,n,range1,range2,val)
          endif
        case ('ehk')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%ehk,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%ehk)/kind(v_ham%ehk)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ehk)/kind(v_ham%ehk)
            call struc_eval_io( v_ham%ehk,rw,n,range1,range2,val)
          endif
        case ('elind')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%elind,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%elind)/kind(v_ham%elind)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%elind)/kind(v_ham%elind)
            call struc_eval_io( v_ham%elind,rw,n,range1,range2,val)
          endif
        case ('eterms')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%eterms,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%eterms)/kind(v_ham%eterms)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%eterms)/kind(v_ham%eterms)
            call struc_eval_io( v_ham%eterms,rw,n,range1,range2,val)
          endif
        case ('hord')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%hord,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%hord)/kind(v_ham%hord)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%hord)/kind(v_ham%hord)
            call struc_eval_io( v_ham%hord,rw,n,range1,range2,val)
          endif
        case ('kmto')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%kmto,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%kmto)/kind(v_ham%kmto)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%kmto)/kind(v_ham%kmto)
            call struc_eval_io( v_ham%kmto,rw,n,range1,range2,val)
          endif
        case ('lasa')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lasa,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lasa)/kind(v_ham%lasa)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lasa)/kind(v_ham%lasa)
            call struc_eval_io( v_ham%lasa,rw,n,range1,range2,val)
          endif
        case ('ldham')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ldham,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ldham)/kind(v_ham%ldham)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ldham)/kind(v_ham%ldham)
            call struc_eval_io( v_ham%ldham,rw,n,range1,range2,val)
          endif
        case ('lgen3')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lgen3,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lgen3)/kind(v_ham%lgen3)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lgen3)/kind(v_ham%lgen3)
            call struc_eval_io( v_ham%lgen3,rw,n,range1,range2,val)
          endif
        case ('lham')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lham,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lham)/kind(v_ham%lham)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lham)/kind(v_ham%lham)
            call struc_eval_io( v_ham%lham,rw,n,range1,range2,val)
          endif
        case ('lmaxu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lmaxu,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lmaxu)/kind(v_ham%lmaxu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lmaxu)/kind(v_ham%lmaxu)
            call struc_eval_io( v_ham%lmaxu,rw,n,range1,range2,val)
          endif
        case ('lmxax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lmxax,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lmxax)/kind(v_ham%lmxax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lmxax)/kind(v_ham%lmxax)
            call struc_eval_io( v_ham%lmxax,rw,n,range1,range2,val)
          endif
        case ('lncol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lncol,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lncol)/kind(v_ham%lncol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lncol)/kind(v_ham%lncol)
            call struc_eval_io( v_ham%lncol,rw,n,range1,range2,val)
          endif
        case ('lsig')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lsig,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lsig)/kind(v_ham%lsig)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lsig)/kind(v_ham%lsig)
            call struc_eval_io( v_ham%lsig,rw,n,range1,range2,val)
          endif
        case ('ltb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ltb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ltb)/kind(v_ham%ltb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ltb)/kind(v_ham%ltb)
            call struc_eval_io( v_ham%ltb,rw,n,range1,range2,val)
          endif
        case ('lxcf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%lxcf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%lxcf)/kind(v_ham%lxcf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%lxcf)/kind(v_ham%lxcf)
            call struc_eval_io( v_ham%lxcf,rw,n,range1,range2,val)
          endif
        case ('nbf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%nbf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%nbf)/kind(v_ham%nbf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%nbf)/kind(v_ham%nbf)
            call struc_eval_io( v_ham%nbf,rw,n,range1,range2,val)
          endif
        case ('ndham')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ndham,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ndham)/kind(v_ham%ndham)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ndham)/kind(v_ham%ndham)
            call struc_eval_io( v_ham%ndham,rw,n,range1,range2,val)
          endif
        case ('ndhrs')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ndhrs,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ndhrs)/kind(v_ham%ndhrs)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ndhrs)/kind(v_ham%ndhrs)
            call struc_eval_io( v_ham%ndhrs,rw,n,range1,range2,val)
          endif
        case ('ndofH')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ndofH,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ndofH)/kind(v_ham%ndofH)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ndofH)/kind(v_ham%ndofH)
            call struc_eval_io( v_ham%ndofH,rw,n,range1,range2,val)
          endif
        case ('neula')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%neula,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%neula)/kind(v_ham%neula)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%neula)/kind(v_ham%neula)
            call struc_eval_io( v_ham%neula,rw,n,range1,range2,val)
          endif
        case ('nkaph')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%nkaph,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%nkaph)/kind(v_ham%nkaph)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%nkaph)/kind(v_ham%nkaph)
            call struc_eval_io( v_ham%nkaph,rw,n,range1,range2,val)
          endif
        case ('nlibu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%nlibu,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%nlibu)/kind(v_ham%nlibu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%nlibu)/kind(v_ham%nlibu)
            call struc_eval_io( v_ham%nlibu,rw,n,range1,range2,val)
          endif
        case ('nmto')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%nmto,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%nmto)/kind(v_ham%nmto)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%nmto)/kind(v_ham%nmto)
            call struc_eval_io( v_ham%nmto,rw,n,range1,range2,val)
          endif
        case ('npwmin')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%npwmin,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%npwmin)/kind(v_ham%npwmin)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%npwmin)/kind(v_ham%npwmin)
            call struc_eval_io( v_ham%npwmin,rw,n,range1,range2,val)
          endif
        case ('npwpad')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%npwpad,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%npwpad)/kind(v_ham%npwpad)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%npwpad)/kind(v_ham%npwpad)
            call struc_eval_io( v_ham%npwpad,rw,n,range1,range2,val)
          endif
        case ('nqsig')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%nqsig,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%nqsig)/kind(v_ham%nqsig)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%nqsig)/kind(v_ham%nqsig)
            call struc_eval_io( v_ham%nqsig,rw,n,range1,range2,val)
          endif
        case ('obdots')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%obdots,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%obdots)/kind(v_ham%obdots)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%obdots)/kind(v_ham%obdots)
            call struc_eval_io( v_ham%obdots,rw,n,range1,range2,val)
          endif
        case ('oeula')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%oeula,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%oeula)/kind(v_ham%oeula)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%oeula)/kind(v_ham%oeula)
            call struc_eval_io( v_ham%oeula,rw,n,range1,range2,val)
          endif
        case ('ohrs')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ohrs,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ohrs)/kind(v_ham%ohrs)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ohrs)/kind(v_ham%ohrs)
            call struc_eval_io( v_ham%ohrs,rw,n,range1,range2,val)
          endif
        case ('oiaxs')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%oiaxs,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%oiaxs)/kind(v_ham%oiaxs)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%oiaxs)/kind(v_ham%oiaxs)
            call struc_eval_io( v_ham%oiaxs,rw,n,range1,range2,val)
          endif
        case ('oindxo')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%oindxo,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%oindxo)/kind(v_ham%oindxo)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%oindxo)/kind(v_ham%oindxo)
            call struc_eval_io( v_ham%oindxo,rw,n,range1,range2,val)
          endif
        case ('olmxa')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%olmxa,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%olmxa)/kind(v_ham%olmxa)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%olmxa)/kind(v_ham%olmxa)
            call struc_eval_io( v_ham%olmxa,rw,n,range1,range2,val)
          endif
        case ('omagf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%omagf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%omagf)/kind(v_ham%omagf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%omagf)/kind(v_ham%omagf)
            call struc_eval_io( v_ham%omagf,rw,n,range1,range2,val)
          endif
        case ('ontabs')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ontabs,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ontabs)/kind(v_ham%ontabs)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ontabs)/kind(v_ham%ontabs)
            call struc_eval_io( v_ham%ontabs,rw,n,range1,range2,val)
          endif
        case ('ooffH')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%ooffH,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%ooffH)/kind(v_ham%ooffH)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%ooffH)/kind(v_ham%ooffH)
            call struc_eval_io( v_ham%ooffH,rw,n,range1,range2,val)
          endif
        case ('oqsig')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%oqsig,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%oqsig)/kind(v_ham%oqsig)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%oqsig)/kind(v_ham%oqsig)
            call struc_eval_io( v_ham%oqsig,rw,n,range1,range2,val)
          endif
        case ('oveps')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%oveps,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%oveps)/kind(v_ham%oveps)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%oveps)/kind(v_ham%oveps)
            call struc_eval_io( v_ham%oveps,rw,n,range1,range2,val)
          endif
        case ('pmax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%pmax,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%pmax)/kind(v_ham%pmax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%pmax)/kind(v_ham%pmax)
            call struc_eval_io( v_ham%pmax,rw,n,range1,range2,val)
          endif
        case ('pmin')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%pmin,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%pmin)/kind(v_ham%pmin)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%pmin)/kind(v_ham%pmin)
            call struc_eval_io( v_ham%pmin,rw,n,range1,range2,val)
          endif
        case ('pwemax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%pwemax,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%pwemax)/kind(v_ham%pwemax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%pwemax)/kind(v_ham%pwemax)
            call struc_eval_io( v_ham%pwemax,rw,n,range1,range2,val)
          endif
        case ('pwemin')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%pwemin,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%pwemin)/kind(v_ham%pwemin)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%pwemin)/kind(v_ham%pwemin)
            call struc_eval_io( v_ham%pwemin,rw,n,range1,range2,val)
          endif
        case ('pwmode')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%pwmode,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%pwmode)/kind(v_ham%pwmode)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%pwmode)/kind(v_ham%pwmode)
            call struc_eval_io( v_ham%pwmode,rw,n,range1,range2,val)
          endif
        case ('qpoff')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%qpoff,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%qpoff)/kind(v_ham%qpoff)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%qpoff)/kind(v_ham%qpoff)
            call struc_eval_io( v_ham%qpoff,rw,n,range1,range2,val)
          endif
        case ('qss')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%qss,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%qss)/kind(v_ham%qss)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%qss)/kind(v_ham%qss)
            call struc_eval_io( v_ham%qss,rw,n,range1,range2,val)
          endif
        case ('rsrnge')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%rsrnge,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%rsrnge)/kind(v_ham%rsrnge)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%rsrnge)/kind(v_ham%rsrnge)
            call struc_eval_io( v_ham%rsrnge,rw,n,range1,range2,val)
          endif
        case ('rsstol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%rsstol,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%rsstol)/kind(v_ham%rsstol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%rsstol)/kind(v_ham%rsstol)
            call struc_eval_io( v_ham%rsstol,rw,n,range1,range2,val)
          endif
        case ('seref')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%seref,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%seref)/kind(v_ham%seref)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%seref)/kind(v_ham%seref)
            call struc_eval_io( v_ham%seref,rw,n,range1,range2,val)
          endif
        case ('sigp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%sigp,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%sigp)/kind(v_ham%sigp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%sigp)/kind(v_ham%sigp)
            call struc_eval_io( v_ham%sigp,rw,n,range1,range2,val)
          endif
        case ('thrpv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ham%thrpv,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ham%thrpv)/kind(v_ham%thrpv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%thrpv)/kind(v_ham%thrpv)
            call struc_eval_io( v_ham%thrpv,rw,n,range1,range2,val)
          endif
        case ('udiag')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ham%udiag,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ham%udiag)/kind(v_ham%udiag)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ham%udiag)/kind(v_ham%udiag)
            call struc_eval_io( v_ham%udiag,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uham_io: failed to find ',name
          stop
        end select
        end subroutine struc_uham_io
       
C ulat 43
       
C offe=-1 is only supported
        subroutine ulat(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_lat):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'ulat offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ulat
       
        subroutine struc_ulat_io(v_lat,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_lat):: v_lat
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call ulat(v_lat,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_lat%size,rw,1,range1,range2,val)
          endif
        case ('alat')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%alat,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%alat)/kind(v_lat%alat)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%alat)/kind(v_lat%alat)
            call struc_eval_io( v_lat%alat,rw,n,range1,range2,val)
          endif
        case ('as')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%as,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%as)/kind(v_lat%as)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%as)/kind(v_lat%as)
            call struc_eval_io( v_lat%as,rw,n,range1,range2,val)
          endif
        case ('avw')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%avw,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%avw)/kind(v_lat%avw)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%avw)/kind(v_lat%avw)
            call struc_eval_io( v_lat%avw,rw,n,range1,range2,val)
          endif
        case ('awald')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%awald,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%awald)/kind(v_lat%awald)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%awald)/kind(v_lat%awald)
            call struc_eval_io( v_lat%awald,rw,n,range1,range2,val)
          endif
        case ('dist')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%dist,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%dist)/kind(v_lat%dist)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%dist)/kind(v_lat%dist)
            call struc_eval_io( v_lat%dist,rw,n,range1,range2,val)
          endif
        case ('gam')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%gam,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%gam)/kind(v_lat%gam)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%gam)/kind(v_lat%gam)
            call struc_eval_io( v_lat%gam,rw,n,range1,range2,val)
          endif
        case ('gmax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%gmax,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%gmax)/kind(v_lat%gmax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%gmax)/kind(v_lat%gmax)
            call struc_eval_io( v_lat%gmax,rw,n,range1,range2,val)
          endif
        case ('ldist')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%ldist,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%ldist)/kind(v_lat%ldist)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%ldist)/kind(v_lat%ldist)
            call struc_eval_io( v_lat%ldist,rw,n,range1,range2,val)
          endif
        case ('nabc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%nabc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%nabc)/kind(v_lat%nabc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%nabc)/kind(v_lat%nabc)
            call struc_eval_io( v_lat%nabc,rw,n,range1,range2,val)
          endif
        case ('ng')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%ng,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%ng)/kind(v_lat%ng)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%ng)/kind(v_lat%ng)
            call struc_eval_io( v_lat%ng,rw,n,range1,range2,val)
          endif
        case ('nkd')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%nkd,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%nkd)/kind(v_lat%nkd)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%nkd)/kind(v_lat%nkd)
            call struc_eval_io( v_lat%nkd,rw,n,range1,range2,val)
          endif
        case ('nkdmx')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%nkdmx,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%nkdmx)/kind(v_lat%nkdmx)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%nkdmx)/kind(v_lat%nkdmx)
            call struc_eval_io( v_lat%nkdmx,rw,n,range1,range2,val)
          endif
        case ('nkq')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%nkq,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%nkq)/kind(v_lat%nkq)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%nkq)/kind(v_lat%nkq)
            call struc_eval_io( v_lat%nkq,rw,n,range1,range2,val)
          endif
        case ('nkqmx')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%nkqmx,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%nkqmx)/kind(v_lat%nkqmx)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%nkqmx)/kind(v_lat%nkqmx)
            call struc_eval_io( v_lat%nkqmx,rw,n,range1,range2,val)
          endif
        case ('npgrp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%npgrp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%npgrp)/kind(v_lat%npgrp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%npgrp)/kind(v_lat%npgrp)
            call struc_eval_io( v_lat%npgrp,rw,n,range1,range2,val)
          endif
        case ('nsgrp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%nsgrp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%nsgrp)/kind(v_lat%nsgrp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%nsgrp)/kind(v_lat%nsgrp)
            call struc_eval_io( v_lat%nsgrp,rw,n,range1,range2,val)
          endif
        case ('oag')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%oag,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%oag)/kind(v_lat%oag)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%oag)/kind(v_lat%oag)
            call struc_eval_io( v_lat%oag,rw,n,range1,range2,val)
          endif
        case ('obgv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%obgv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%obgv)/kind(v_lat%obgv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%obgv)/kind(v_lat%obgv)
            call struc_eval_io( v_lat%obgv,rw,n,range1,range2,val)
          endif
        case ('ocg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%ocg,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%ocg)/kind(v_lat%ocg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%ocg)/kind(v_lat%ocg)
            call struc_eval_io( v_lat%ocg,rw,n,range1,range2,val)
          endif
        case ('ocy')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%ocy,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%ocy)/kind(v_lat%ocy)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%ocy)/kind(v_lat%ocy)
            call struc_eval_io( v_lat%ocy,rw,n,range1,range2,val)
          endif
        case ('odlv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%odlv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%odlv)/kind(v_lat%odlv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%odlv)/kind(v_lat%odlv)
            call struc_eval_io( v_lat%odlv,rw,n,range1,range2,val)
          endif
        case ('ogv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%ogv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%ogv)/kind(v_lat%ogv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%ogv)/kind(v_lat%ogv)
            call struc_eval_io( v_lat%ogv,rw,n,range1,range2,val)
          endif
        case ('oidxcg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%oidxcg,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%oidxcg)/kind(v_lat%oidxcg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%oidxcg)/kind(v_lat%oidxcg)
            call struc_eval_io( v_lat%oidxcg,rw,n,range1,range2,val)
          endif
        case ('oips0')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%oips0,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%oips0)/kind(v_lat%oips0)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%oips0)/kind(v_lat%oips0)
            call struc_eval_io( v_lat%oips0,rw,n,range1,range2,val)
          endif
        case ('oistab')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%oistab,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%oistab)/kind(v_lat%oistab)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%oistab)/kind(v_lat%oistab)
            call struc_eval_io( v_lat%oistab,rw,n,range1,range2,val)
          endif
        case ('ojcg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%ojcg,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%ojcg)/kind(v_lat%ojcg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%ojcg)/kind(v_lat%ojcg)
            call struc_eval_io( v_lat%ojcg,rw,n,range1,range2,val)
          endif
        case ('okv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%okv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%okv)/kind(v_lat%okv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%okv)/kind(v_lat%okv)
            call struc_eval_io( v_lat%okv,rw,n,range1,range2,val)
          endif
        case ('opos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%opos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%opos)/kind(v_lat%opos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%opos)/kind(v_lat%opos)
            call struc_eval_io( v_lat%opos,rw,n,range1,range2,val)
          endif
        case ('oqlv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%oqlv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%oqlv)/kind(v_lat%oqlv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%oqlv)/kind(v_lat%oqlv)
            call struc_eval_io( v_lat%oqlv,rw,n,range1,range2,val)
          endif
        case ('osymgr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_lat%osymgr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_lat%osymgr)/kind(v_lat%osymgr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%osymgr)/kind(v_lat%osymgr)
            call struc_eval_io( v_lat%osymgr,rw,n,range1,range2,val)
          endif
        case ('plat')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%plat,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%plat)/kind(v_lat%plat)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%plat)/kind(v_lat%plat)
            call struc_eval_io( v_lat%plat,rw,n,range1,range2,val)
          endif
        case ('plat0')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%plat0,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%plat0)/kind(v_lat%plat0)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%plat0)/kind(v_lat%plat0)
            call struc_eval_io( v_lat%plat0,rw,n,range1,range2,val)
          endif
        case ('plat2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%plat2,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%plat2)/kind(v_lat%plat2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%plat2)/kind(v_lat%plat2)
            call struc_eval_io( v_lat%plat2,rw,n,range1,range2,val)
          endif
        case ('plate')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%plate,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%plate)/kind(v_lat%plate)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%plate)/kind(v_lat%plate)
            call struc_eval_io( v_lat%plate,rw,n,range1,range2,val)
          endif
        case ('platl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%platl,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%platl)/kind(v_lat%platl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%platl)/kind(v_lat%platl)
            call struc_eval_io( v_lat%platl,rw,n,range1,range2,val)
          endif
        case ('platr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%platr,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%platr)/kind(v_lat%platr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%platr)/kind(v_lat%platr)
            call struc_eval_io( v_lat%platr,rw,n,range1,range2,val)
          endif
        case ('qlat')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%qlat,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%qlat)/kind(v_lat%qlat)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%qlat)/kind(v_lat%qlat)
            call struc_eval_io( v_lat%qlat,rw,n,range1,range2,val)
          endif
        case ('rpad')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%rpad,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%rpad)/kind(v_lat%rpad)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%rpad)/kind(v_lat%rpad)
            call struc_eval_io( v_lat%rpad,rw,n,range1,range2,val)
          endif
        case ('slat')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%slat,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%slat)/kind(v_lat%slat)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%slat)/kind(v_lat%slat)
            call struc_eval_io( v_lat%slat,rw,n,range1,range2,val)
          endif
        case ('tol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%tol,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%tol)/kind(v_lat%tol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%tol)/kind(v_lat%tol)
            call struc_eval_io( v_lat%tol,rw,n,range1,range2,val)
          endif
        case ('tolft')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%tolft,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%tolft)/kind(v_lat%tolft)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%tolft)/kind(v_lat%tolft)
            call struc_eval_io( v_lat%tolft,rw,n,range1,range2,val)
          endif
        case ('vol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_lat%vol,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_lat%vol)/kind(v_lat%vol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_lat%vol)/kind(v_lat%vol)
            call struc_eval_io( v_lat%vol,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ulat_io: failed to find ',name
          stop
        end select
        end subroutine struc_ulat_io
       
C umix 23
       
C offe=-1 is only supported
        subroutine umix(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_mix):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'umix offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine umix
       
        subroutine struc_umix_io(v_mix,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_mix):: v_mix
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call umix(v_mix,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_mix%size,rw,1,range1,range2,val)
          endif
        case ('b')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%b,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%b)/kind(v_mix%b)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%b)/kind(v_mix%b)
            call struc_eval_io( v_mix%b,rw,n,range1,range2,val)
          endif
        case ('bl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%bl,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%bl)/kind(v_mix%bl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%bl)/kind(v_mix%bl)
            call struc_eval_io( v_mix%bl,rw,n,range1,range2,val)
          endif
        case ('bv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%bv,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%bv)/kind(v_mix%bv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%bv)/kind(v_mix%bv)
            call struc_eval_io( v_mix%bv,rw,n,range1,range2,val)
          endif
        case ('elind')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%elind,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%elind)/kind(v_mix%elind)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%elind)/kind(v_mix%elind)
            call struc_eval_io( v_mix%elind,rw,n,range1,range2,val)
          endif
        case ('fn')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%fn,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_mix%fn)/kind(v_mix%fn)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%fn)/kind(v_mix%fn)
            call struc_eval_io( v_mix%fn,rw,n,range1,range2,val)
          endif
        case ('kill')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%kill,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%kill)/kind(v_mix%kill)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%kill)/kind(v_mix%kill)
            call struc_eval_io( v_mix%kill,rw,n,range1,range2,val)
          endif
        case ('lxpot')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%lxpot,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%lxpot)/kind(v_mix%lxpot)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%lxpot)/kind(v_mix%lxpot)
            call struc_eval_io( v_mix%lxpot,rw,n,range1,range2,val)
          endif
        case ('mmix')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%mmix,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%mmix)/kind(v_mix%mmix)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%mmix)/kind(v_mix%mmix)
            call struc_eval_io( v_mix%mmix,rw,n,range1,range2,val)
          endif
        case ('mode')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%mode,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%mode)/kind(v_mix%mode)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%mode)/kind(v_mix%mode)
            call struc_eval_io( v_mix%mode,rw,n,range1,range2,val)
          endif
        case ('model')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%model,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%model)/kind(v_mix%model)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%model)/kind(v_mix%model)
            call struc_eval_io( v_mix%model,rw,n,range1,range2,val)
          endif
        case ('n')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%n,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%n)/kind(v_mix%n)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%n)/kind(v_mix%n)
            call struc_eval_io( v_mix%n,rw,n,range1,range2,val)
          endif
        case ('nitu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%nitu,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%nitu)/kind(v_mix%nitu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%nitu)/kind(v_mix%nitu)
            call struc_eval_io( v_mix%nitu,rw,n,range1,range2,val)
          endif
        case ('nmix')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%nmix,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%nmix)/kind(v_mix%nmix)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%nmix)/kind(v_mix%nmix)
            call struc_eval_io( v_mix%nmix,rw,n,range1,range2,val)
          endif
        case ('nsave')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_mix%nsave,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_mix%nsave)/kind(v_mix%nsave)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%nsave)/kind(v_mix%nsave)
            call struc_eval_io( v_mix%nsave,rw,n,range1,range2,val)
          endif
        case ('r')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%r,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_mix%r)/kind(v_mix%r)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%r)/kind(v_mix%r)
            call struc_eval_io( v_mix%r,rw,n,range1,range2,val)
          endif
        case ('rms1')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%rms1,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%rms1)/kind(v_mix%rms1)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%rms1)/kind(v_mix%rms1)
            call struc_eval_io( v_mix%rms1,rw,n,range1,range2,val)
          endif
        case ('rms2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%rms2,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%rms2)/kind(v_mix%rms2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%rms2)/kind(v_mix%rms2)
            call struc_eval_io( v_mix%rms2,rw,n,range1,range2,val)
          endif
        case ('tj')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%tj,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%tj)/kind(v_mix%tj)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%tj)/kind(v_mix%tj)
            call struc_eval_io( v_mix%tj,rw,n,range1,range2,val)
          endif
        case ('tolu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%tolu,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%tolu)/kind(v_mix%tolu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%tolu)/kind(v_mix%tolu)
            call struc_eval_io( v_mix%tolu,rw,n,range1,range2,val)
          endif
        case ('umix')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%umix,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%umix)/kind(v_mix%umix)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%umix)/kind(v_mix%umix)
            call struc_eval_io( v_mix%umix,rw,n,range1,range2,val)
          endif
        case ('w')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%w,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%w)/kind(v_mix%w)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%w)/kind(v_mix%w)
            call struc_eval_io( v_mix%w,rw,n,range1,range2,val)
          endif
        case ('wc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_mix%wc,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_mix%wc)/kind(v_mix%wc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_mix%wc)/kind(v_mix%wc)
            call struc_eval_io( v_mix%wc,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_umix_io: failed to find ',name
          stop
        end select
        end subroutine struc_umix_io
       
C umove 11
       
C offe=-1 is only supported
        subroutine umove(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_move):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'umove offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine umove
       
        subroutine struc_umove_io(v_move,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_move):: v_move
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call umove(v_move,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_move%size,rw,1,range1,range2,val)
          endif
        case ('ct')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%ct,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%ct)/kind(v_move%ct)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%ct)/kind(v_move%ct)
            call struc_eval_io( v_move%ct,rw,n,range1,range2,val)
          endif
        case ('gyro')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%gyro,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%gyro)/kind(v_move%gyro)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%gyro)/kind(v_move%gyro)
            call struc_eval_io( v_move%gyro,rw,n,range1,range2,val)
          endif
        case ('kt')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%kt,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%kt)/kind(v_move%kt)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%kt)/kind(v_move%kt)
            call struc_eval_io( v_move%kt,rw,n,range1,range2,val)
          endif
        case ('modt')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_move%modt,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_move%modt)/kind(v_move%modt)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%modt)/kind(v_move%modt)
            call struc_eval_io( v_move%modt,rw,n,range1,range2,val)
          endif
        case ('nmodt')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_move%nmodt,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_move%nmodt)/kind(v_move%nmodt)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%nmodt)/kind(v_move%nmodt)
            call struc_eval_io( v_move%nmodt,rw,n,range1,range2,val)
          endif
        case ('prmint')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%prmint,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%prmint)/kind(v_move%prmint)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%prmint)/kind(v_move%prmint)
            call struc_eval_io( v_move%prmint,rw,n,range1,range2,val)
          endif
        case ('tnow')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%tnow,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%tnow)/kind(v_move%tnow)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%tnow)/kind(v_move%tnow)
            call struc_eval_io( v_move%tnow,rw,n,range1,range2,val)
          endif
        case ('ts')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%ts,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%ts)/kind(v_move%ts)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%ts)/kind(v_move%ts)
            call struc_eval_io( v_move%ts,rw,n,range1,range2,val)
          endif
        case ('tsequ')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%tsequ,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%tsequ)/kind(v_move%tsequ)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%tsequ)/kind(v_move%tsequ)
            call struc_eval_io( v_move%tsequ,rw,n,range1,range2,val)
          endif
        case ('tstot')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_move%tstot,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_move%tstot)/kind(v_move%tstot)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_move%tstot)/kind(v_move%tstot)
            call struc_eval_io( v_move%tstot,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_umove_io: failed to find ',name
          stop
        end select
        end subroutine struc_umove_io
       
C uoptic 13
       
C offe=-1 is only supported
        subroutine uoptic(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_optic):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'uoptic offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uoptic
       
        subroutine struc_uoptic_io(v_optic,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_optic):: v_optic
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call uoptic(v_optic,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_optic%size,rw,1,range1,range2,val)
          endif
        case ('axes')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%axes,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%axes)/kind(v_optic%axes)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%axes)/kind(v_optic%axes)
            call struc_eval_io( v_optic%axes,rw,n,range1,range2,val)
          endif
        case ('cll')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%cll,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%cll)/kind(v_optic%cll)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%cll)/kind(v_optic%cll)
            call struc_eval_io( v_optic%cll,rw,n,range1,range2,val)
          endif
        case ('cln')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%cln,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%cln)/kind(v_optic%cln)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%cln)/kind(v_optic%cln)
            call struc_eval_io( v_optic%cln,rw,n,range1,range2,val)
          endif
        case ('cls')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%cls,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%cls)/kind(v_optic%cls)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%cls)/kind(v_optic%cls)
            call struc_eval_io( v_optic%cls,rw,n,range1,range2,val)
          endif
        case ('clsite')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%clsite,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%clsite)/kind(v_optic%clsite)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%clsite)/kind(v_optic%clsite)
            call struc_eval_io( v_optic%clsite,rw,n,range1,range2,val)
          endif
        case ('esciss')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_optic%esciss,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_optic%esciss)/kind(v_optic%esciss)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%esciss)/kind(v_optic%esciss)
            call struc_eval_io( v_optic%esciss,rw,n,range1,range2,val)
          endif
        case ('mode')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%mode,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%mode)/kind(v_optic%mode)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%mode)/kind(v_optic%mode)
            call struc_eval_io( v_optic%mode,rw,n,range1,range2,val)
          endif
        case ('nchi2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%nchi2,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%nchi2)/kind(v_optic%nchi2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%nchi2)/kind(v_optic%nchi2)
            call struc_eval_io( v_optic%nchi2,rw,n,range1,range2,val)
          endif
        case ('ne')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%ne,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%ne)/kind(v_optic%ne)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%ne)/kind(v_optic%ne)
            call struc_eval_io( v_optic%ne,rw,n,range1,range2,val)
          endif
        case ('ocrng')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%ocrng,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%ocrng)/kind(v_optic%ocrng)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%ocrng)/kind(v_optic%ocrng)
            call struc_eval_io( v_optic%ocrng,rw,n,range1,range2,val)
          endif
        case ('unrng')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_optic%unrng,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_optic%unrng)/kind(v_optic%unrng)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%unrng)/kind(v_optic%unrng)
            call struc_eval_io( v_optic%unrng,rw,n,range1,range2,val)
          endif
        case ('window')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_optic%window,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_optic%window)/kind(v_optic%window)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_optic%window)/kind(v_optic%window)
            call struc_eval_io( v_optic%window,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uoptic_io: failed to find ',name
          stop
        end select
        end subroutine struc_uoptic_io
       
C uordn 11
       
C offe=-1 is only supported
        subroutine uordn(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_ordn):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'uordn offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uordn
       
        subroutine struc_uordn_io(v_ordn,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_ordn):: v_ordn
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call uordn(v_ordn,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_ordn%size,rw,1,range1,range2,val)
          endif
        case ('efre')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ordn%efre,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ordn%efre)/kind(v_ordn%efre)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%efre)/kind(v_ordn%efre)
            call struc_eval_io( v_ordn%efre,rw,n,range1,range2,val)
          endif
        case ('mode')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%mode,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%mode)/kind(v_ordn%mode)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%mode)/kind(v_ordn%mode)
            call struc_eval_io( v_ordn%mode,rw,n,range1,range2,val)
          endif
        case ('ncl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%ncl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%ncl)/kind(v_ordn%ncl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%ncl)/kind(v_ordn%ncl)
            call struc_eval_io( v_ordn%ncl,rw,n,range1,range2,val)
          endif
        case ('ndofh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%ndofh,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%ndofh)/kind(v_ordn%ndofh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%ndofh)/kind(v_ordn%ndofh)
            call struc_eval_io( v_ordn%ndofh,rw,n,range1,range2,val)
          endif
        case ('oclp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%oclp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%oclp)/kind(v_ordn%oclp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%oclp)/kind(v_ordn%oclp)
            call struc_eval_io( v_ordn%oclp,rw,n,range1,range2,val)
          endif
        case ('oclssl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%oclssl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%oclssl)/kind(v_ordn%oclssl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%oclssl)/kind(v_ordn%oclssl)
            call struc_eval_io( v_ordn%oclssl,rw,n,range1,range2,val)
          endif
        case ('oiaxg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%oiaxg,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%oiaxg)/kind(v_ordn%oiaxg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%oiaxg)/kind(v_ordn%oiaxg)
            call struc_eval_io( v_ordn%oiaxg,rw,n,range1,range2,val)
          endif
        case ('omapgv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%omapgv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%omapgv)/kind(v_ordn%omapgv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%omapgv)/kind(v_ordn%omapgv)
            call struc_eval_io( v_ordn%omapgv,rw,n,range1,range2,val)
          endif
        case ('ooffch')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_ordn%ooffch,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_ordn%ooffch)/kind(v_ordn%ooffch)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%ooffch)/kind(v_ordn%ooffch)
            call struc_eval_io( v_ordn%ooffch,rw,n,range1,range2,val)
          endif
        case ('rmaxg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_ordn%rmaxg,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_ordn%rmaxg)/kind(v_ordn%rmaxg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_ordn%rmaxg)/kind(v_ordn%rmaxg)
            call struc_eval_io( v_ordn%rmaxg,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uordn_io: failed to find ',name
          stop
        end select
        end subroutine struc_uordn_io
       
C upot 58
       
C offe=-1 is only supported
        subroutine upot(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_pot):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'upot offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine upot
       
        subroutine struc_upot_io(v_pot,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_pot):: v_pot
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call upot(v_pot,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_pot%size,rw,1,range1,range2,val)
          endif
        case ('bfield')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_pot%bfield,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_pot%bfield)/kind(v_pot%bfield)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%bfield)/kind(v_pot%bfield)
            call struc_eval_io( v_pot%bfield,rw,n,range1,range2,val)
          endif
        case ('nlma')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%nlma,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%nlma)/kind(v_pot%nlma)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%nlma)/kind(v_pot%nlma)
            call struc_eval_io( v_pot%nlma,rw,n,range1,range2,val)
          endif
        case ('nlml')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%nlml,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%nlml)/kind(v_pot%nlml)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%nlml)/kind(v_pot%nlml)
            call struc_eval_io( v_pot%nlml,rw,n,range1,range2,val)
          endif
        case ('nrhos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%nrhos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%nrhos)/kind(v_pot%nrhos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%nrhos)/kind(v_pot%nrhos)
            call struc_eval_io( v_pot%nrhos,rw,n,range1,range2,val)
          endif
        case ('oaamom')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oaamom,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oaamom)/kind(v_pot%oaamom)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oaamom)/kind(v_pot%oaamom)
            call struc_eval_io( v_pot%oaamom,rw,n,range1,range2,val)
          endif
        case ('obxc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%obxc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%obxc)/kind(v_pot%obxc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%obxc)/kind(v_pot%obxc)
            call struc_eval_io( v_pot%obxc,rw,n,range1,range2,val)
          endif
        case ('odddpf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%odddpf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%odddpf)/kind(v_pot%odddpf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%odddpf)/kind(v_pot%odddpf)
            call struc_eval_io( v_pot%odddpf,rw,n,range1,range2,val)
          endif
        case ('oddpf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oddpf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oddpf)/kind(v_pot%oddpf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oddpf)/kind(v_pot%oddpf)
            call struc_eval_io( v_pot%oddpf,rw,n,range1,range2,val)
          endif
        case ('oddpfr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oddpfr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oddpfr)/kind(v_pot%oddpfr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oddpfr)/kind(v_pot%oddpfr)
            call struc_eval_io( v_pot%oddpfr,rw,n,range1,range2,val)
          endif
        case ('odel')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%odel,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%odel)/kind(v_pot%odel)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%odel)/kind(v_pot%odel)
            call struc_eval_io( v_pot%odel,rw,n,range1,range2,val)
          endif
        case ('odpf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%odpf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%odpf)/kind(v_pot%odpf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%odpf)/kind(v_pot%odpf)
            call struc_eval_io( v_pot%odpf,rw,n,range1,range2,val)
          endif
        case ('odpfr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%odpfr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%odpfr)/kind(v_pot%odpfr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%odpfr)/kind(v_pot%odpfr)
            call struc_eval_io( v_pot%odpfr,rw,n,range1,range2,val)
          endif
        case ('ofes')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ofes,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ofes)/kind(v_pot%ofes)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ofes)/kind(v_pot%ofes)
            call struc_eval_io( v_pot%ofes,rw,n,range1,range2,val)
          endif
        case ('ogma')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ogma,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ogma)/kind(v_pot%ogma)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ogma)/kind(v_pot%ogma)
            call struc_eval_io( v_pot%ogma,rw,n,range1,range2,val)
          endif
        case ('ogmar')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ogmar,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ogmar)/kind(v_pot%ogmar)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ogmar)/kind(v_pot%ogmar)
            call struc_eval_io( v_pot%ogmar,rw,n,range1,range2,val)
          endif
        case ('ogrrme')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ogrrme,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ogrrme)/kind(v_pot%ogrrme)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ogrrme)/kind(v_pot%ogrrme)
            call struc_eval_io( v_pot%ogrrme,rw,n,range1,range2,val)
          endif
        case ('ohab')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ohab,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ohab)/kind(v_pot%ohab)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ohab)/kind(v_pot%ohab)
            call struc_eval_io( v_pot%ohab,rw,n,range1,range2,val)
          endif
        case ('oivso')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oivso,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oivso)/kind(v_pot%oivso)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oivso)/kind(v_pot%oivso)
            call struc_eval_io( v_pot%oivso,rw,n,range1,range2,val)
          endif
        case ('omad')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%omad,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%omad)/kind(v_pot%omad)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%omad)/kind(v_pot%omad)
            call struc_eval_io( v_pot%omad,rw,n,range1,range2,val)
          endif
        case ('oorhat')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oorhat,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oorhat)/kind(v_pot%oorhat)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oorhat)/kind(v_pot%oorhat)
            call struc_eval_io( v_pot%oorhat,rw,n,range1,range2,val)
          endif
        case ('opalp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opalp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opalp)/kind(v_pot%opalp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opalp)/kind(v_pot%opalp)
            call struc_eval_io( v_pot%opalp,rw,n,range1,range2,val)
          endif
        case ('opapg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opapg,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opapg)/kind(v_pot%opapg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opapg)/kind(v_pot%opapg)
            call struc_eval_io( v_pot%opapg,rw,n,range1,range2,val)
          endif
        case ('opdel')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opdel,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opdel)/kind(v_pot%opdel)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opdel)/kind(v_pot%opdel)
            call struc_eval_io( v_pot%opdel,rw,n,range1,range2,val)
          endif
        case ('opf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opf)/kind(v_pot%opf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opf)/kind(v_pot%opf)
            call struc_eval_io( v_pot%opf,rw,n,range1,range2,val)
          endif
        case ('opfnc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opfnc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opfnc)/kind(v_pot%opfnc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opfnc)/kind(v_pot%opfnc)
            call struc_eval_io( v_pot%opfnc,rw,n,range1,range2,val)
          endif
        case ('opfr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opfr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opfr)/kind(v_pot%opfr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opfr)/kind(v_pot%opfr)
            call struc_eval_io( v_pot%opfr,rw,n,range1,range2,val)
          endif
        case ('opmpol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opmpol,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opmpol)/kind(v_pot%opmpol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opmpol)/kind(v_pot%opmpol)
            call struc_eval_io( v_pot%opmpol,rw,n,range1,range2,val)
          endif
        case ('opnu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opnu,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opnu)/kind(v_pot%opnu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opnu)/kind(v_pot%opnu)
            call struc_eval_io( v_pot%opnu,rw,n,range1,range2,val)
          endif
        case ('opp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opp)/kind(v_pot%opp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opp)/kind(v_pot%opp)
            call struc_eval_io( v_pot%opp,rw,n,range1,range2,val)
          endif
        case ('oppi')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oppi,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oppi)/kind(v_pot%oppi)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oppi)/kind(v_pot%oppi)
            call struc_eval_io( v_pot%oppi,rw,n,range1,range2,val)
          endif
        case ('oppn')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oppn,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oppn)/kind(v_pot%oppn)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oppn)/kind(v_pot%oppn)
            call struc_eval_io( v_pot%oppn,rw,n,range1,range2,val)
          endif
        case ('opprel')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opprel,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opprel)/kind(v_pot%opprel)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opprel)/kind(v_pot%opprel)
            call struc_eval_io( v_pot%opprel,rw,n,range1,range2,val)
          endif
        case ('opti')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%opti,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%opti)/kind(v_pot%opti)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%opti)/kind(v_pot%opti)
            call struc_eval_io( v_pot%opti,rw,n,range1,range2,val)
          endif
        case ('oqc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oqc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oqc)/kind(v_pot%oqc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oqc)/kind(v_pot%oqc)
            call struc_eval_io( v_pot%oqc,rw,n,range1,range2,val)
          endif
        case ('oqmom')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oqmom,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oqmom)/kind(v_pot%oqmom)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oqmom)/kind(v_pot%oqmom)
            call struc_eval_io( v_pot%oqmom,rw,n,range1,range2,val)
          endif
        case ('oqnu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oqnu,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oqnu)/kind(v_pot%oqnu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oqnu)/kind(v_pot%oqnu)
            call struc_eval_io( v_pot%oqnu,rw,n,range1,range2,val)
          endif
        case ('oqpp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oqpp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oqpp)/kind(v_pot%oqpp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oqpp)/kind(v_pot%oqpp)
            call struc_eval_io( v_pot%oqpp,rw,n,range1,range2,val)
          endif
        case ('oqt')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oqt,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oqt)/kind(v_pot%oqt)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oqt)/kind(v_pot%oqt)
            call struc_eval_io( v_pot%oqt,rw,n,range1,range2,val)
          endif
        case ('orhos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%orhos,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%orhos)/kind(v_pot%orhos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%orhos)/kind(v_pot%orhos)
            call struc_eval_io( v_pot%orhos,rw,n,range1,range2,val)
          endif
        case ('orhrmx')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%orhrmx,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%orhrmx)/kind(v_pot%orhrmx)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%orhrmx)/kind(v_pot%orhrmx)
            call struc_eval_io( v_pot%orhrmx,rw,n,range1,range2,val)
          endif
        case ('osab')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%osab,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%osab)/kind(v_pot%osab)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%osab)/kind(v_pot%osab)
            call struc_eval_io( v_pot%osab,rw,n,range1,range2,val)
          endif
        case ('osgw')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%osgw,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%osgw)/kind(v_pot%osgw)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%osgw)/kind(v_pot%osgw)
            call struc_eval_io( v_pot%osgw,rw,n,range1,range2,val)
          endif
        case ('osmpot')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%osmpot,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%osmpot)/kind(v_pot%osmpot)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%osmpot)/kind(v_pot%osmpot)
            call struc_eval_io( v_pot%osmpot,rw,n,range1,range2,val)
          endif
        case ('osmrho')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%osmrho,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%osmrho)/kind(v_pot%osmrho)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%osmrho)/kind(v_pot%osmrho)
            call struc_eval_io( v_pot%osmrho,rw,n,range1,range2,val)
          endif
        case ('osop')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%osop,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%osop)/kind(v_pot%osop)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%osop)/kind(v_pot%osop)
            call struc_eval_io( v_pot%osop,rw,n,range1,range2,val)
          endif
        case ('osoptc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%osoptc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%osoptc)/kind(v_pot%osoptc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%osoptc)/kind(v_pot%osoptc)
            call struc_eval_io( v_pot%osoptc,rw,n,range1,range2,val)
          endif
        case ('osrout')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%osrout,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%osrout)/kind(v_pot%osrout)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%osrout)/kind(v_pot%osrout)
            call struc_eval_io( v_pot%osrout,rw,n,range1,range2,val)
          endif
        case ('otau')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%otau,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%otau)/kind(v_pot%otau)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%otau)/kind(v_pot%otau)
            call struc_eval_io( v_pot%otau,rw,n,range1,range2,val)
          endif
        case ('ovab')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ovab,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ovab)/kind(v_pot%ovab)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ovab)/kind(v_pot%ovab)
            call struc_eval_io( v_pot%ovab,rw,n,range1,range2,val)
          endif
        case ('ovdif')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ovdif,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ovdif)/kind(v_pot%ovdif)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ovdif)/kind(v_pot%ovdif)
            call struc_eval_io( v_pot%ovdif,rw,n,range1,range2,val)
          endif
        case ('oves')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%oves,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%oves)/kind(v_pot%oves)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%oves)/kind(v_pot%oves)
            call struc_eval_io( v_pot%oves,rw,n,range1,range2,val)
          endif
        case ('ovintr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ovintr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ovintr)/kind(v_pot%ovintr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ovintr)/kind(v_pot%ovintr)
            call struc_eval_io( v_pot%ovintr,rw,n,range1,range2,val)
          endif
        case ('ovrmax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ovrmax,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ovrmax)/kind(v_pot%ovrmax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ovrmax)/kind(v_pot%ovrmax)
            call struc_eval_io( v_pot%ovrmax,rw,n,range1,range2,val)
          endif
        case ('ovshf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_pot%ovshf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_pot%ovshf)/kind(v_pot%ovshf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%ovshf)/kind(v_pot%ovshf)
            call struc_eval_io( v_pot%ovshf,rw,n,range1,range2,val)
          endif
        case ('vconst')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_pot%vconst,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_pot%vconst)/kind(v_pot%vconst)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%vconst)/kind(v_pot%vconst)
            call struc_eval_io( v_pot%vconst,rw,n,range1,range2,val)
          endif
        case ('vmtz')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_pot%vmtz,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_pot%vmtz)/kind(v_pot%vmtz)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%vmtz)/kind(v_pot%vmtz)
            call struc_eval_io( v_pot%vmtz,rw,n,range1,range2,val)
          endif
        case ('vmtz0')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_pot%vmtz0,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_pot%vmtz0)/kind(v_pot%vmtz0)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_pot%vmtz0)/kind(v_pot%vmtz0)
            call struc_eval_io( v_pot%vmtz0,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_upot_io: failed to find ',name
          stop
        end select
        end subroutine struc_upot_io
       
C uspec 74
       
C offe=-1 is only supported
        subroutine uspec(struc,offe,nli,is,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_spec):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'uspec offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uspec
       
        subroutine struc_uspec_io(v_spec,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_spec):: v_spec
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call uspec(v_spec,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_spec%size,rw,1,range1,range2,val)
          endif
        case ('z')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%z,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%z)/kind(v_spec%z)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%z)/kind(v_spec%z)
            call struc_eval_io( v_spec%z,rw,n,range1,range2,val)
          endif
        case ('mass')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%mass,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%mass)/kind(v_spec%mass)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%mass)/kind(v_spec%mass)
            call struc_eval_io( v_spec%mass,rw,n,range1,range2,val)
          endif
        case ('rmt')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rmt,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rmt)/kind(v_spec%rmt)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rmt)/kind(v_spec%rmt)
            call struc_eval_io( v_spec%rmt,rw,n,range1,range2,val)
          endif
        case ('ntorb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%ntorb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%ntorb)/kind(v_spec%ntorb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%ntorb)/kind(v_spec%ntorb)
            call struc_eval_io( v_spec%ntorb,rw,n,range1,range2,val)
          endif
        case ('naug')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%naug,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%naug)/kind(v_spec%naug)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%naug)/kind(v_spec%naug)
            call struc_eval_io( v_spec%naug,rw,n,range1,range2,val)
          endif
        case ('norb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%norb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%norb)/kind(v_spec%norb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%norb)/kind(v_spec%norb)
            call struc_eval_io( v_spec%norb,rw,n,range1,range2,val)
          endif
        case ('rsmfa')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rsmfa,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rsmfa)/kind(v_spec%rsmfa)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rsmfa)/kind(v_spec%rsmfa)
            call struc_eval_io( v_spec%rsmfa,rw,n,range1,range2,val)
          endif
        case ('rsma')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rsma,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rsma)/kind(v_spec%rsma)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rsma)/kind(v_spec%rsma)
            call struc_eval_io( v_spec%rsma,rw,n,range1,range2,val)
          endif
        case ('rg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rg,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rg)/kind(v_spec%rg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rg)/kind(v_spec%rg)
            call struc_eval_io( v_spec%rg,rw,n,range1,range2,val)
          endif
        case ('lmxa')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%lmxa,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%lmxa)/kind(v_spec%lmxa)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%lmxa)/kind(v_spec%lmxa)
            call struc_eval_io( v_spec%lmxa,rw,n,range1,range2,val)
          endif
        case ('kmxh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%kmxh,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%kmxh)/kind(v_spec%kmxh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%kmxh)/kind(v_spec%kmxh)
            call struc_eval_io( v_spec%kmxh,rw,n,range1,range2,val)
          endif
        case ('lmxl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%lmxl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%lmxl)/kind(v_spec%lmxl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%lmxl)/kind(v_spec%lmxl)
            call struc_eval_io( v_spec%lmxl,rw,n,range1,range2,val)
          endif
        case ('kmxt')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%kmxt,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%kmxt)/kind(v_spec%kmxt)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%kmxt)/kind(v_spec%kmxt)
            call struc_eval_io( v_spec%kmxt,rw,n,range1,range2,val)
          endif
        case ('rsmv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rsmv,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rsmv)/kind(v_spec%rsmv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rsmv)/kind(v_spec%rsmv)
            call struc_eval_io( v_spec%rsmv,rw,n,range1,range2,val)
          endif
        case ('norp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%norp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%norp)/kind(v_spec%norp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%norp)/kind(v_spec%norp)
            call struc_eval_io( v_spec%norp,rw,n,range1,range2,val)
          endif
        case ('coreh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%coreh,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_spec%coreh)/kind(v_spec%coreh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%coreh)/kind(v_spec%coreh)
            call struc_eval_io( v_spec%coreh,rw,n,range1,range2,val)
          endif
        case ('coreq')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%coreq,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%coreq)/kind(v_spec%coreq)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%coreq)/kind(v_spec%coreq)
            call struc_eval_io( v_spec%coreq,rw,n,range1,range2,val)
          endif
        case ('a')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%a,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%a)/kind(v_spec%a)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%a)/kind(v_spec%a)
            call struc_eval_io( v_spec%a,rw,n,range1,range2,val)
          endif
        case ('nr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%nr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%nr)/kind(v_spec%nr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%nr)/kind(v_spec%nr)
            call struc_eval_io( v_spec%nr,rw,n,range1,range2,val)
          endif
        case ('eref')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%eref,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%eref)/kind(v_spec%eref)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%eref)/kind(v_spec%eref)
            call struc_eval_io( v_spec%eref,rw,n,range1,range2,val)
          endif
        case ('etf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%etf,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%etf)/kind(v_spec%etf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%etf)/kind(v_spec%etf)
            call struc_eval_io( v_spec%etf,rw,n,range1,range2,val)
          endif
        case ('beta')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%beta,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%beta)/kind(v_spec%beta)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%beta)/kind(v_spec%beta)
            call struc_eval_io( v_spec%beta,rw,n,range1,range2,val)
          endif
        case ('lfoca')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%lfoca,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%lfoca)/kind(v_spec%lfoca)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%lfoca)/kind(v_spec%lfoca)
            call struc_eval_io( v_spec%lfoca,rw,n,range1,range2,val)
          endif
        case ('ctail')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%ctail,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%ctail)/kind(v_spec%ctail)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%ctail)/kind(v_spec%ctail)
            call struc_eval_io( v_spec%ctail,rw,n,range1,range2,val)
          endif
        case ('etail')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%etail,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%etail)/kind(v_spec%etail)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%etail)/kind(v_spec%etail)
            call struc_eval_io( v_spec%etail,rw,n,range1,range2,val)
          endif
        case ('name')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%name,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_spec%name)/kind(v_spec%name)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%name)/kind(v_spec%name)
            call struc_eval_io( v_spec%name,rw,n,range1,range2,val)
          endif
        case ('orhoc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%orhoc,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%orhoc)/kind(v_spec%orhoc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%orhoc)/kind(v_spec%orhoc)
            call struc_eval_io( v_spec%orhoc,rw,n,range1,range2,val)
          endif
        case ('stc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%stc,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%stc)/kind(v_spec%stc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%stc)/kind(v_spec%stc)
            call struc_eval_io( v_spec%stc,rw,n,range1,range2,val)
          endif
        case ('lmxb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%lmxb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%lmxb)/kind(v_spec%lmxb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%lmxb)/kind(v_spec%lmxb)
            call struc_eval_io( v_spec%lmxb,rw,n,range1,range2,val)
          endif
        case ('lmxf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%lmxf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%lmxf)/kind(v_spec%lmxf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%lmxf)/kind(v_spec%lmxf)
            call struc_eval_io( v_spec%lmxf,rw,n,range1,range2,val)
          endif
        case ('stni')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%stni,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%stni)/kind(v_spec%stni)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%stni)/kind(v_spec%stni)
            call struc_eval_io( v_spec%stni,rw,n,range1,range2,val)
          endif
        case ('stnm')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%stnm,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%stnm)/kind(v_spec%stnm)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%stnm)/kind(v_spec%stnm)
            call struc_eval_io( v_spec%stnm,rw,n,range1,range2,val)
          endif
        case ('rham')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rham,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rham)/kind(v_spec%rham)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rham)/kind(v_spec%rham)
            call struc_eval_io( v_spec%rham,rw,n,range1,range2,val)
          endif
        case ('rfoca')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rfoca,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rfoca)/kind(v_spec%rfoca)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rfoca)/kind(v_spec%rfoca)
            call struc_eval_io( v_spec%rfoca,rw,n,range1,range2,val)
          endif
        case ('dv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%dv,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%dv)/kind(v_spec%dv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%dv)/kind(v_spec%dv)
            call struc_eval_io( v_spec%dv,rw,n,range1,range2,val)
          endif
        case ('mxcst')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%mxcst,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%mxcst)/kind(v_spec%mxcst)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%mxcst)/kind(v_spec%mxcst)
            call struc_eval_io( v_spec%mxcst,rw,n,range1,range2,val)
          endif
        case ('group')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%group,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%group)/kind(v_spec%group)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%group)/kind(v_spec%group)
            call struc_eval_io( v_spec%group,rw,n,range1,range2,val)
          endif
        case ('grp2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%grp2,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%grp2)/kind(v_spec%grp2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%grp2)/kind(v_spec%grp2)
            call struc_eval_io( v_spec%grp2,rw,n,range1,range2,val)
          endif
        case ('nxi')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%nxi,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%nxi)/kind(v_spec%nxi)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%nxi)/kind(v_spec%nxi)
            call struc_eval_io( v_spec%nxi,rw,n,range1,range2,val)
          endif
        case ('qc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%qc,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%qc)/kind(v_spec%qc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%qc)/kind(v_spec%qc)
            call struc_eval_io( v_spec%qc,rw,n,range1,range2,val)
          endif
        case ('lmxpb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%lmxpb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%lmxpb)/kind(v_spec%lmxpb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%lmxpb)/kind(v_spec%lmxpb)
            call struc_eval_io( v_spec%lmxpb,rw,n,range1,range2,val)
          endif
        case ('pb1')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%pb1,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_spec%pb1)/kind(v_spec%pb1)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%pb1)/kind(v_spec%pb1)
            call struc_eval_io( v_spec%pb1,rw,n,range1,range2,val)
          endif
        case ('pb2')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%pb2,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_spec%pb2)/kind(v_spec%pb2)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%pb2)/kind(v_spec%pb2)
            call struc_eval_io( v_spec%pb2,rw,n,range1,range2,val)
          endif
        case ('colxbs')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%colxbs,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%colxbs)/kind(v_spec%colxbs)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%colxbs)/kind(v_spec%colxbs)
            call struc_eval_io( v_spec%colxbs,rw,n,range1,range2,val)
          endif
        case ('lxi')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%lxi,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%lxi)/kind(v_spec%lxi)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%lxi)/kind(v_spec%lxi)
            call struc_eval_io( v_spec%lxi,rw,n,range1,range2,val)
          endif
        case ('radxbs')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%radxbs,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%radxbs)/kind(v_spec%radxbs)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%radxbs)/kind(v_spec%radxbs)
            call struc_eval_io( v_spec%radxbs,rw,n,range1,range2,val)
          endif
        case ('rcut')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rcut,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rcut)/kind(v_spec%rcut)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rcut)/kind(v_spec%rcut)
            call struc_eval_io( v_spec%rcut,rw,n,range1,range2,val)
          endif
        case ('rint')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rint,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rint)/kind(v_spec%rint)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rint)/kind(v_spec%rint)
            call struc_eval_io( v_spec%rint,rw,n,range1,range2,val)
          endif
        case ('eh3')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%eh3,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%eh3)/kind(v_spec%eh3)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%eh3)/kind(v_spec%eh3)
            call struc_eval_io( v_spec%eh3,rw,n,range1,range2,val)
          endif
        case ('rs3')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rs3,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rs3)/kind(v_spec%rs3)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rs3)/kind(v_spec%rs3)
            call struc_eval_io( v_spec%rs3,rw,n,range1,range2,val)
          endif
        case ('vmtz')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%vmtz,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%vmtz)/kind(v_spec%vmtz)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%vmtz)/kind(v_spec%vmtz)
            call struc_eval_io( v_spec%vmtz,rw,n,range1,range2,val)
          endif
        case ('kmxv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%kmxv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%kmxv)/kind(v_spec%kmxv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%kmxv)/kind(v_spec%kmxv)
            call struc_eval_io( v_spec%kmxv,rw,n,range1,range2,val)
          endif
        case ('rcfa')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%rcfa,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%rcfa)/kind(v_spec%rcfa)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%rcfa)/kind(v_spec%rcfa)
            call struc_eval_io( v_spec%rcfa,rw,n,range1,range2,val)
          endif
        case ('p')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%p,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%p)/kind(v_spec%p)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%p)/kind(v_spec%p)
            call struc_eval_io( v_spec%p,rw,n,range1,range2,val)
          endif
        case ('q')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%q,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%q)/kind(v_spec%q)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%q)/kind(v_spec%q)
            call struc_eval_io( v_spec%q,rw,n,range1,range2,val)
          endif
        case ('alpha')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%alpha,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%alpha)/kind(v_spec%alpha)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%alpha)/kind(v_spec%alpha)
            call struc_eval_io( v_spec%alpha,rw,n,range1,range2,val)
          endif
        case ('idmod')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%idmod,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%idmod)/kind(v_spec%idmod)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%idmod)/kind(v_spec%idmod)
            call struc_eval_io( v_spec%idmod,rw,n,range1,range2,val)
          endif
        case ('idxdn')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%idxdn,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%idxdn)/kind(v_spec%idxdn)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%idxdn)/kind(v_spec%idxdn)
            call struc_eval_io( v_spec%idxdn,rw,n,range1,range2,val)
          endif
        case ('hcr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%hcr,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%hcr)/kind(v_spec%hcr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%hcr)/kind(v_spec%hcr)
            call struc_eval_io( v_spec%hcr,rw,n,range1,range2,val)
          endif
        case ('exi')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%exi,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%exi)/kind(v_spec%exi)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%exi)/kind(v_spec%exi)
            call struc_eval_io( v_spec%exi,rw,n,range1,range2,val)
          endif
        case ('ngcut')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%ngcut,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%ngcut)/kind(v_spec%ngcut)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%ngcut)/kind(v_spec%ngcut)
            call struc_eval_io( v_spec%ngcut,rw,n,range1,range2,val)
          endif
        case ('qpol')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%qpol,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%qpol)/kind(v_spec%qpol)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%qpol)/kind(v_spec%qpol)
            call struc_eval_io( v_spec%qpol,rw,n,range1,range2,val)
          endif
        case ('chfa')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%chfa,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%chfa)/kind(v_spec%chfa)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%chfa)/kind(v_spec%chfa)
            call struc_eval_io( v_spec%chfa,rw,n,range1,range2,val)
          endif
        case ('orbp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%orbp,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%orbp)/kind(v_spec%orbp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%orbp)/kind(v_spec%orbp)
            call struc_eval_io( v_spec%orbp,rw,n,range1,range2,val)
          endif
        case ('enu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%enu,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%enu)/kind(v_spec%enu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%enu)/kind(v_spec%enu)
            call struc_eval_io( v_spec%enu,rw,n,range1,range2,val)
          endif
        case ('pz')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%pz,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%pz)/kind(v_spec%pz)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%pz)/kind(v_spec%pz)
            call struc_eval_io( v_spec%pz,rw,n,range1,range2,val)
          endif
        case ('idu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%idu,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%idu)/kind(v_spec%idu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%idu)/kind(v_spec%idu)
            call struc_eval_io( v_spec%idu,rw,n,range1,range2,val)
          endif
        case ('uh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%uh,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%uh)/kind(v_spec%uh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%uh)/kind(v_spec%uh)
            call struc_eval_io( v_spec%uh,rw,n,range1,range2,val)
          endif
        case ('jh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%jh,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%jh)/kind(v_spec%jh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%jh)/kind(v_spec%jh)
            call struc_eval_io( v_spec%jh,rw,n,range1,range2,val)
          endif
        case ('iq1')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%iq1,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%iq1)/kind(v_spec%iq1)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%iq1)/kind(v_spec%iq1)
            call struc_eval_io( v_spec%iq1,rw,n,range1,range2,val)
          endif
        case ('ivso')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_spec%ivso,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_spec%ivso)/kind(v_spec%ivso)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%ivso)/kind(v_spec%ivso)
            call struc_eval_io( v_spec%ivso,rw,n,range1,range2,val)
          endif
        case ('ehvl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%ehvl,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%ehvl)/kind(v_spec%ehvl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%ehvl)/kind(v_spec%ehvl)
            call struc_eval_io( v_spec%ehvl,rw,n,range1,range2,val)
          endif
        case ('vso')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_spec%vso,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_spec%vso)/kind(v_spec%vso)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_spec%vso)/kind(v_spec%vso)
            call struc_eval_io( v_spec%vso,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uspec_io: failed to find ',name
          stop
        end select
        end subroutine struc_uspec_io
       
C usite 27
       
C offe=-1 is only supported
        subroutine usite(struc,offe,nli,is,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_site):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'usite offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine usite
       
        subroutine struc_usite_io(v_site,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_site):: v_site
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call usite(v_site,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_site%size,rw,1,range1,range2,val)
          endif
        case ('bfield')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%bfield,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%bfield)/kind(v_site%bfield)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%bfield)/kind(v_site%bfield)
            call struc_eval_io( v_site%bfield,rw,n,range1,range2,val)
          endif
        case ('clabel')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%clabel,rw,1,range1,range2,v_ri%r)
            v_ri%t= 1
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 1
             v_ri%i= sizeof(v_site%clabel)/kind(v_site%clabel)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%clabel)/kind(v_site%clabel)
            call struc_eval_io( v_site%clabel,rw,n,range1,range2,val)
          endif
        case ('class')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%class,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%class)/kind(v_site%class)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%class)/kind(v_site%class)
            call struc_eval_io( v_site%class,rw,n,range1,range2,val)
          endif
        case ('cli')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%cli,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%cli)/kind(v_site%cli)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%cli)/kind(v_site%cli)
            call struc_eval_io( v_site%cli,rw,n,range1,range2,val)
          endif
        case ('delta')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%delta,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%delta)/kind(v_site%delta)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%delta)/kind(v_site%delta)
            call struc_eval_io( v_site%delta,rw,n,range1,range2,val)
          endif
        case ('dpole')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%dpole,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%dpole)/kind(v_site%dpole)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%dpole)/kind(v_site%dpole)
            call struc_eval_io( v_site%dpole,rw,n,range1,range2,val)
          endif
        case ('eula')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%eula,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%eula)/kind(v_site%eula)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%eula)/kind(v_site%eula)
            call struc_eval_io( v_site%eula,rw,n,range1,range2,val)
          endif
        case ('force')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%force,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%force)/kind(v_site%force)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%force)/kind(v_site%force)
            call struc_eval_io( v_site%force,rw,n,range1,range2,val)
          endif
        case ('mpole')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%mpole,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%mpole)/kind(v_site%mpole)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%mpole)/kind(v_site%mpole)
            call struc_eval_io( v_site%mpole,rw,n,range1,range2,val)
          endif
        case ('ndelta')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%ndelta,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%ndelta)/kind(v_site%ndelta)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%ndelta)/kind(v_site%ndelta)
            call struc_eval_io( v_site%ndelta,rw,n,range1,range2,val)
          endif
        case ('norb')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%norb,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%norb)/kind(v_site%norb)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%norb)/kind(v_site%norb)
            call struc_eval_io( v_site%norb,rw,n,range1,range2,val)
          endif
        case ('offh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%offh,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%offh)/kind(v_site%offh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%offh)/kind(v_site%offh)
            call struc_eval_io( v_site%offh,rw,n,range1,range2,val)
          endif
        case ('orho')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%orho,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%orho)/kind(v_site%orho)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%orho)/kind(v_site%orho)
            call struc_eval_io( v_site%orho,rw,n,range1,range2,val)
          endif
        case ('ov0')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%ov0,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%ov0)/kind(v_site%ov0)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%ov0)/kind(v_site%ov0)
            call struc_eval_io( v_site%ov0,rw,n,range1,range2,val)
          endif
        case ('ov1')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%ov1,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%ov1)/kind(v_site%ov1)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%ov1)/kind(v_site%ov1)
            call struc_eval_io( v_site%ov1,rw,n,range1,range2,val)
          endif
        case ('pl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%pl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%pl)/kind(v_site%pl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%pl)/kind(v_site%pl)
            call struc_eval_io( v_site%pl,rw,n,range1,range2,val)
          endif
        case ('plv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%plv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%plv)/kind(v_site%plv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%plv)/kind(v_site%plv)
            call struc_eval_io( v_site%plv,rw,n,range1,range2,val)
          endif
        case ('pnu')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%pnu,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%pnu)/kind(v_site%pnu)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%pnu)/kind(v_site%pnu)
            call struc_eval_io( v_site%pnu,rw,n,range1,range2,val)
          endif
        case ('pos')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%pos,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%pos)/kind(v_site%pos)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%pos)/kind(v_site%pos)
            call struc_eval_io( v_site%pos,rw,n,range1,range2,val)
          endif
        case ('pos0')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%pos0,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%pos0)/kind(v_site%pos0)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%pos0)/kind(v_site%pos0)
            call struc_eval_io( v_site%pos0,rw,n,range1,range2,val)
          endif
        case ('pz')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%pz,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%pz)/kind(v_site%pz)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%pz)/kind(v_site%pz)
            call struc_eval_io( v_site%pz,rw,n,range1,range2,val)
          endif
        case ('relax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%relax,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%relax)/kind(v_site%relax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%relax)/kind(v_site%relax)
            call struc_eval_io( v_site%relax,rw,n,range1,range2,val)
          endif
        case ('sid')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%sid,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%sid)/kind(v_site%sid)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%sid)/kind(v_site%sid)
            call struc_eval_io( v_site%sid,rw,n,range1,range2,val)
          endif
        case ('spec')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_site%spec,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_site%spec)/kind(v_site%spec)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%spec)/kind(v_site%spec)
            call struc_eval_io( v_site%spec,rw,n,range1,range2,val)
          endif
        case ('vel')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%vel,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%vel)/kind(v_site%vel)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%vel)/kind(v_site%vel)
            call struc_eval_io( v_site%vel,rw,n,range1,range2,val)
          endif
        case ('vshft')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_site%vshft,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_site%vshft)/kind(v_site%vshft)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_site%vshft)/kind(v_site%vshft)
            call struc_eval_io( v_site%vshft,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_usite_io: failed to find ',name
          stop
        end select
        end subroutine struc_usite_io
       
C ustr 30
       
C offe=-1 is only supported
        subroutine ustr(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_str):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'ustr offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ustr
       
        subroutine struc_ustr_io(v_str,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_str):: v_str
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call ustr(v_str,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_str%size,rw,1,range1,range2,val)
          endif
        case ('adec')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%adec,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%adec)/kind(v_str%adec)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%adec)/kind(v_str%adec)
            call struc_eval_io( v_str%adec,rw,n,range1,range2,val)
          endif
        case ('amode')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%amode,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%amode)/kind(v_str%amode)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%amode)/kind(v_str%amode)
            call struc_eval_io( v_str%amode,rw,n,range1,range2,val)
          endif
        case ('drwats')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%drwats,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%drwats)/kind(v_str%drwats)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%drwats)/kind(v_str%drwats)
            call struc_eval_io( v_str%drwats,rw,n,range1,range2,val)
          endif
        case ('iinv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%iinv,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%iinv)/kind(v_str%iinv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%iinv)/kind(v_str%iinv)
            call struc_eval_io( v_str%iinv,rw,n,range1,range2,val)
          endif
        case ('ivl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%ivl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%ivl)/kind(v_str%ivl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%ivl)/kind(v_str%ivl)
            call struc_eval_io( v_str%ivl,rw,n,range1,range2,val)
          endif
        case ('kaps')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%kaps,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%kaps)/kind(v_str%kaps)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%kaps)/kind(v_str%kaps)
            call struc_eval_io( v_str%kaps,rw,n,range1,range2,val)
          endif
        case ('lequiv')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%lequiv,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%lequiv)/kind(v_str%lequiv)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%lequiv)/kind(v_str%lequiv)
            call struc_eval_io( v_str%lequiv,rw,n,range1,range2,val)
          endif
        case ('lmaxw')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%lmaxw,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%lmaxw)/kind(v_str%lmaxw)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%lmaxw)/kind(v_str%lmaxw)
            call struc_eval_io( v_str%lmaxw,rw,n,range1,range2,val)
          endif
        case ('loka')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%loka,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%loka)/kind(v_str%loka)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%loka)/kind(v_str%loka)
            call struc_eval_io( v_str%loka,rw,n,range1,range2,val)
          endif
        case ('lshow')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%lshow,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%lshow)/kind(v_str%lshow)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%lshow)/kind(v_str%lshow)
            call struc_eval_io( v_str%lshow,rw,n,range1,range2,val)
          endif
        case ('mxnbr')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%mxnbr,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%mxnbr)/kind(v_str%mxnbr)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%mxnbr)/kind(v_str%mxnbr)
            call struc_eval_io( v_str%mxnbr,rw,n,range1,range2,val)
          endif
        case ('nalf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%nalf,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%nalf)/kind(v_str%nalf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%nalf)/kind(v_str%nalf)
            call struc_eval_io( v_str%nalf,rw,n,range1,range2,val)
          endif
        case ('nbisi')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%nbisi,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%nbisi)/kind(v_str%nbisi)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%nbisi)/kind(v_str%nbisi)
            call struc_eval_io( v_str%nbisi,rw,n,range1,range2,val)
          endif
        case ('ncupl')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%ncupl,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%ncupl)/kind(v_str%ncupl)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%ncupl)/kind(v_str%ncupl)
            call struc_eval_io( v_str%ncupl,rw,n,range1,range2,val)
          endif
        case ('ndust')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%ndust,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%ndust)/kind(v_str%ndust)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%ndust)/kind(v_str%ndust)
            call struc_eval_io( v_str%ndust,rw,n,range1,range2,val)
          endif
        case ('nitab')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%nitab,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%nitab)/kind(v_str%nitab)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%nitab)/kind(v_str%nitab)
            call struc_eval_io( v_str%nitab,rw,n,range1,range2,val)
          endif
        case ('nkaps')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%nkaps,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%nkaps)/kind(v_str%nkaps)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%nkaps)/kind(v_str%nkaps)
            call struc_eval_io( v_str%nkaps,rw,n,range1,range2,val)
          endif
        case ('nttab')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%nttab,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%nttab)/kind(v_str%nttab)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%nttab)/kind(v_str%nttab)
            call struc_eval_io( v_str%nttab,rw,n,range1,range2,val)
          endif
        case ('oadot')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%oadot,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%oadot)/kind(v_str%oadot)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%oadot)/kind(v_str%oadot)
            call struc_eval_io( v_str%oadot,rw,n,range1,range2,val)
          endif
        case ('oalp')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%oalp,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%oalp)/kind(v_str%oalp)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%oalp)/kind(v_str%oalp)
            call struc_eval_io( v_str%oalp,rw,n,range1,range2,val)
          endif
        case ('oiax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%oiax,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%oiax)/kind(v_str%oiax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%oiax)/kind(v_str%oiax)
            call struc_eval_io( v_str%oiax,rw,n,range1,range2,val)
          endif
        case ('ontab')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%ontab,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%ontab)/kind(v_str%ontab)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%ontab)/kind(v_str%ontab)
            call struc_eval_io( v_str%ontab,rw,n,range1,range2,val)
          endif
        case ('os')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%os,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%os)/kind(v_str%os)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%os)/kind(v_str%os)
            call struc_eval_io( v_str%os,rw,n,range1,range2,val)
          endif
        case ('osdot')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_str%osdot,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_str%osdot)/kind(v_str%osdot)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%osdot)/kind(v_str%osdot)
            call struc_eval_io( v_str%osdot,rw,n,range1,range2,val)
          endif
        case ('rfit')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%rfit,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%rfit)/kind(v_str%rfit)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%rfit)/kind(v_str%rfit)
            call struc_eval_io( v_str%rfit,rw,n,range1,range2,val)
          endif
        case ('rmax')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%rmax,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%rmax)/kind(v_str%rmax)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%rmax)/kind(v_str%rmax)
            call struc_eval_io( v_str%rmax,rw,n,range1,range2,val)
          endif
        case ('rmaxg')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%rmaxg,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%rmaxg)/kind(v_str%rmaxg)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%rmaxg)/kind(v_str%rmaxg)
            call struc_eval_io( v_str%rmaxg,rw,n,range1,range2,val)
          endif
        case ('skmsh')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%skmsh,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%skmsh)/kind(v_str%skmsh)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%skmsh)/kind(v_str%skmsh)
            call struc_eval_io( v_str%skmsh,rw,n,range1,range2,val)
          endif
        case ('wztcf')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_str%wztcf,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_str%wztcf)/kind(v_str%wztcf)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_str%wztcf)/kind(v_str%wztcf)
            call struc_eval_io( v_str%wztcf,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ustr_io: failed to find ',name
          stop
        end select
        end subroutine struc_ustr_io
       
C utb 7
       
C offe=-1 is only supported
        subroutine utb(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_tb):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
C  +2 is margin
C next call will overwrite nsize if rw='u'
          val=sizeof(struc)/8+2
          if (mod(val,2).eq.1) val=val+1
          call dcopy(val,dble(nli),0,struc,1)
          struc%size=val
        else
          write(*,*) 'utb offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine utb
       
        subroutine struc_utb_io(v_tb,name,arr,rw,range1,range2,val,v_ri)
C use val or v_ri
C arr='1': use v_ri and only 1 component if it is vector, used for dgets*, igets*
C arr='s': use v_ri and return type and size of varname
C arr='a': use val and pack/unpack
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_tb):: v_tb
        character(*),intent(in):: name
        character,intent(in):: rw,arr
        integer:: range1,range2
        integer:: val
        type(s_ri):: v_ri
        integer:: n,v_iv(10)
        select case (trim(adjustl(name)))
        case ('size')
          if (rw.eq.'p') then
          v_iv(1)=-1
          call utb(v_tb,v_iv,0,0,0,0)
          else
          call struc_eval_io( v_tb%size,rw,1,range1,range2,val)
          endif
        case ('alam')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_tb%alam,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_tb%alam)/kind(v_tb%alam)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_tb%alam)/kind(v_tb%alam)
            call struc_eval_io( v_tb%alam,rw,n,range1,range2,val)
          endif
        case ('alsc')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_tb%alsc,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_tb%alsc)/kind(v_tb%alsc)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_tb%alsc)/kind(v_tb%alsc)
            call struc_eval_io( v_tb%alsc,rw,n,range1,range2,val)
          endif
        case ('ebfit')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_tb%ebfit,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_tb%ebfit)/kind(v_tb%ebfit)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_tb%ebfit)/kind(v_tb%ebfit)
            call struc_eval_io( v_tb%ebfit,rw,n,range1,range2,val)
          endif
        case ('fmode')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_tb%fmode,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_tb%fmode)/kind(v_tb%fmode)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_tb%fmode)/kind(v_tb%fmode)
            call struc_eval_io( v_tb%fmode,rw,n,range1,range2,val)
          endif
        case ('nbfit')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_i8_realbody( v_tb%nbfit,rw,1,range1,range2,v_ri%i)
            v_ri%t=2
            v_ri%r= v_ri%i
           else if (arr.eq.'s') then
             v_ri%t= 2
             v_ri%i= sizeof(v_tb%nbfit)/kind(v_tb%nbfit)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_tb%nbfit)/kind(v_tb%nbfit)
            call struc_eval_io( v_tb%nbfit,rw,n,range1,range2,val)
          endif
        case ('rmfit')
          if (arr.eq.'1') then
            if ( range2.eq.-1 ) range2 = range1
            call struc_eval_io_r8_realbody( v_tb%rmfit,rw,1,range1,range2,v_ri%r)
            v_ri%t= 4
            v_ri%i= nint(v_ri%r)
           else if (arr.eq.'s') then
             v_ri%t= 4
             v_ri%i= sizeof(v_tb%rmfit)/kind(v_tb%rmfit)
             v_ri%r= v_ri%i
          else
            n= sizeof(v_tb%rmfit)/kind(v_tb%rmfit)
            call struc_eval_io( v_tb%rmfit,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_utb_io: failed to find ',name
          stop
        end select
        end subroutine struc_utb_io
       
        subroutine struc_packupack_val1(subname,name,struc,arr,rw,range1,range2,value,v_ri)
        use m_struc_def
        implicit none
        real(8):: struc
        integer:: range1,range2
        integer:: value
        character(*):: subname,name
        character:: rw,arr
        type(s_ri):: v_ri
        select case (trim(adjustl(subname)))
          case ('bz')
            call struc_ubz_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('array')
            call struc_uarray_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('ctrl')
            call struc_uctrl_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('gw')
            call struc_ugw_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('ham')
            call struc_uham_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('lat')
            call struc_ulat_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('mix')
            call struc_umix_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('move')
            call struc_umove_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('optic')
            call struc_uoptic_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('ordn')
            call struc_uordn_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('pot')
            call struc_upot_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('spec')
            call struc_uspec_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('site')
            call struc_usite_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('str')
            call struc_ustr_io(struc,name,arr,rw,range1,range2,value,v_ri)
          case ('tb')
            call struc_utb_io(struc,name,arr,rw,range1,range2,value,v_ri)
        case default
          write(*,*)'struc_packupack_val1: failed to find str=',subname
          stop
        end select
        end subroutine struc_packupack_val1
       
