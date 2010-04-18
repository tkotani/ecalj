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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_array%ozos,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uarray_io: failed to find ',name
          stop
        end select
        end subroutine struc_uarray_io
       
C uctrl 56
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_bz%dosw)
             v_ri%r= v_ri%i
          else
             n= size(v_bz%dosw)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_bz%lcond)
             v_ri%r= v_ri%i
          else
             n= size(v_bz%lcond)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_bz%lshft)
             v_ri%r= v_ri%i
          else
             n= size(v_bz%lshft)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_bz%nkabc)
             v_ri%r= v_ri%i
          else
             n= size(v_bz%nkabc)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_bz%semsh)
             v_ri%r= v_ri%i
          else
             n= size(v_bz%semsh)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_bz%zval,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ubz_io: failed to find ',name
          stop
        end select
        end subroutine struc_ubz_io
       
C uarray 32
       
C offe=-1 is only supported
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
             v_ri%i= size(v_ctrl%defm)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%defm)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ctrl%lpgf)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%lpgf)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ctrl%lstonr)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%lstonr)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ctrl%mdprm)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%mdprm)
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
             v_ri%i= size(v_ctrl%modep)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%modep)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ctrl%nesabc)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%nesabc)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ctrl%omax1)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%omax1)
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
             v_ri%i= size(v_ctrl%omax2)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%omax2)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ctrl%sdprm)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%sdprm)
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
             v_ri%i= size(v_ctrl%sdxsi)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%sdxsi)
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
             v_ri%i= size(v_ctrl%smalit)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%smalit)
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
             v_ri%i= size(v_ctrl%tol)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%tol)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ctrl%zbak)
             v_ri%r= v_ri%i
          else
             n= size(v_ctrl%zbak)
            call struc_eval_io( v_ctrl%zbak,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uctrl_io: failed to find ',name
          stop
        end select
        end subroutine struc_uctrl_io
       
C ugw 19
       
C offe=-1 is only supported
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
             v_ri%i= size(v_gw%delre)
             v_ri%r= v_ri%i
          else
             n= size(v_gw%delre)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_gw%nkabc)
             v_ri%r= v_ri%i
          else
             n= size(v_gw%nkabc)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_gw%rdsig,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ugw_io: failed to find ',name
          stop
        end select
        end subroutine struc_ugw_io
       
C uham 56
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ham%dabc)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%dabc)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ham%eterms)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%eterms)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ham%kmto)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%kmto)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ham%ldham)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%ldham)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ham%pmax)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%pmax)
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
             v_ri%i= size(v_ham%pmin)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%pmin)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ham%qpoff)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%qpoff)
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
             v_ri%i= size(v_ham%qss)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%qss)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ham%sigp)
             v_ri%r= v_ri%i
          else
             n= size(v_ham%sigp)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_ham%udiag,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uham_io: failed to find ',name
          stop
        end select
        end subroutine struc_uham_io
       
C ulat 43
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_lat%dist)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%dist)
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
             v_ri%i= size(v_lat%gam)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%gam)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_lat%nabc)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%nabc)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_lat%plat)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%plat)
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
             v_ri%i= size(v_lat%plat0)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%plat0)
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
             v_ri%i= size(v_lat%plat2)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%plat2)
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
             v_ri%i= size(v_lat%plate)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%plate)
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
             v_ri%i= size(v_lat%platl)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%platl)
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
             v_ri%i= size(v_lat%platr)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%platr)
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
             v_ri%i= size(v_lat%qlat)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%qlat)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_lat%slat)
             v_ri%r= v_ri%i
          else
             n= size(v_lat%slat)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_lat%vol,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ulat_io: failed to find ',name
          stop
        end select
        end subroutine struc_ulat_io
       
C umix 23
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_mix%r)
             v_ri%r= v_ri%i
          else
             n= size(v_mix%r)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_mix%tj)
             v_ri%r= v_ri%i
          else
             n= size(v_mix%tj)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_mix%w)
             v_ri%r= v_ri%i
          else
             n= size(v_mix%w)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_mix%wc,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_umix_io: failed to find ',name
          stop
        end select
        end subroutine struc_umix_io
       
C umove 11
       
C offe=-1 is only supported
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
             v_ri%i= size(v_move%ct)
             v_ri%r= v_ri%i
          else
             n= size(v_move%ct)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_move%modt)
             v_ri%r= v_ri%i
          else
             n= size(v_move%modt)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_move%prmint)
             v_ri%r= v_ri%i
          else
             n= size(v_move%prmint)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_move%tstot,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_umove_io: failed to find ',name
          stop
        end select
        end subroutine struc_umove_io
       
C uoptic 13
       
C offe=-1 is only supported
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
             v_ri%i= size(v_optic%axes)
             v_ri%r= v_ri%i
          else
             n= size(v_optic%axes)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_optic%ocrng)
             v_ri%r= v_ri%i
          else
             n= size(v_optic%ocrng)
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
             v_ri%i= size(v_optic%unrng)
             v_ri%r= v_ri%i
          else
             n= size(v_optic%unrng)
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
             v_ri%i= size(v_optic%window)
             v_ri%r= v_ri%i
          else
             n= size(v_optic%window)
            call struc_eval_io( v_optic%window,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uoptic_io: failed to find ',name
          stop
        end select
        end subroutine struc_uoptic_io
       
C uordn 11
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_ordn%rmaxg)
             v_ri%r= v_ri%i
          else
             n= size(v_ordn%rmaxg)
            call struc_eval_io( v_ordn%rmaxg,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uordn_io: failed to find ',name
          stop
        end select
        end subroutine struc_uordn_io
       
C upot 58
       
C offe=-1 is only supported
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
             v_ri%i= size(v_pot%bfield)
             v_ri%r= v_ri%i
          else
             n= size(v_pot%bfield)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_pot%vconst)
             v_ri%r= v_ri%i
          else
             n= size(v_pot%vconst)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_pot%vmtz0,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_upot_io: failed to find ',name
          stop
        end select
        end subroutine struc_upot_io
       
C uspec 74
       
C offe=-1 is only supported
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
             v_ri%i= size(v_site%bfield)
             v_ri%r= v_ri%i
          else
             n= size(v_site%bfield)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_site%delta)
             v_ri%r= v_ri%i
          else
             n= size(v_site%delta)
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
             v_ri%i= size(v_site%dpole)
             v_ri%r= v_ri%i
          else
             n= size(v_site%dpole)
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
             v_ri%i= size(v_site%eula)
             v_ri%r= v_ri%i
          else
             n= size(v_site%eula)
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
             v_ri%i= size(v_site%force)
             v_ri%r= v_ri%i
          else
             n= size(v_site%force)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_site%pnu)
             v_ri%r= v_ri%i
          else
             n= size(v_site%pnu)
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
             v_ri%i= size(v_site%pos)
             v_ri%r= v_ri%i
          else
             n= size(v_site%pos)
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
             v_ri%i= size(v_site%pos0)
             v_ri%r= v_ri%i
          else
             n= size(v_site%pos0)
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
             v_ri%i= size(v_site%pz)
             v_ri%r= v_ri%i
          else
             n= size(v_site%pz)
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
             v_ri%i= size(v_site%relax)
             v_ri%r= v_ri%i
          else
             n= size(v_site%relax)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_site%vel)
             v_ri%r= v_ri%i
          else
             n= size(v_site%vel)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_site%vshft,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_usite_io: failed to find ',name
          stop
        end select
        end subroutine struc_usite_io
       
C ustr 30
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_spec%coreq)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%coreq)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_spec%colxbs)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%colxbs)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_spec%rcfa)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%rcfa)
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
             v_ri%i= size(v_spec%p)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%p)
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
             v_ri%i= size(v_spec%q)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%q)
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
             v_ri%i= size(v_spec%alpha)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%alpha)
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
             v_ri%i= size(v_spec%idmod)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%idmod)
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
             v_ri%i= size(v_spec%idxdn)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%idxdn)
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
             v_ri%i= size(v_spec%hcr)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%hcr)
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
             v_ri%i= size(v_spec%exi)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%exi)
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
             v_ri%i= size(v_spec%ngcut)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%ngcut)
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
             v_ri%i= size(v_spec%qpol)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%qpol)
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
             v_ri%i= size(v_spec%chfa)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%chfa)
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
             v_ri%i= size(v_spec%orbp)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%orbp)
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
             v_ri%i= size(v_spec%enu)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%enu)
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
             v_ri%i= size(v_spec%pz)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%pz)
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
             v_ri%i= size(v_spec%idu)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%idu)
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
             v_ri%i= size(v_spec%uh)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%uh)
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
             v_ri%i= size(v_spec%jh)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%jh)
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
             v_ri%i= size(v_spec%iq1)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%iq1)
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
             v_ri%i= size(v_spec%ivso)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%ivso)
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
             v_ri%i= size(v_spec%ehvl)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%ehvl)
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
             v_ri%i= size(v_spec%vso)
             v_ri%r= v_ri%i
          else
             n= size(v_spec%vso)
            call struc_eval_io( v_spec%vso,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_uspec_io: failed to find ',name
          stop
        end select
        end subroutine struc_uspec_io
       
C usite 27
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_str%iinv)
             v_ri%r= v_ri%i
          else
             n= size(v_str%iinv)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_str%kaps)
             v_ri%r= v_ri%i
          else
             n= size(v_str%kaps)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_str%nbisi)
             v_ri%r= v_ri%i
          else
             n= size(v_str%nbisi)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_str%skmsh)
             v_ri%r= v_ri%i
          else
             n= size(v_str%skmsh)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
            call struc_eval_io( v_str%wztcf,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_ustr_io: failed to find ',name
          stop
        end select
        end subroutine struc_ustr_io
       
C utb 7
       
C offe=-1 is only supported
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_tb%ebfit)
             v_ri%r= v_ri%i
          else
             n= size(v_tb%ebfit)
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
             v_ri%i= 1
             v_ri%r= v_ri%i
          else
             n= 1
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
             v_ri%i= size(v_tb%nbfit)
             v_ri%r= v_ri%i
          else
             n= size(v_tb%nbfit)
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
             v_ri%i= size(v_tb%rmfit)
             v_ri%r= v_ri%i
          else
             n= size(v_tb%rmfit)
            call struc_eval_io( v_tb%rmfit,rw,n,range1,range2,val)
          endif
        case default
          write(*,*) 'struc_utb_io: failed to find ',name
          stop
        end select
        end subroutine struc_utb_io
       
        subroutine uarray_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_array):: struc
Cki        integer:: uarray_size
        struc%size=uarray_size()
        struc%ncl=0
        struc%nclasp=0
        struc%nofgl=0
        struc%nofgr=0
        struc%npadl=0
        struc%npadr=0
        struc%npl=0
        struc%oclabl=0
        struc%ocllst=0
        struc%oclp=0
        struc%oclssl=0
        struc%ogroup=0
        struc%ogrp2=0
        struc%ohave=0
        struc%oics=0
        struc%oifrlx=0
        struc%oipc=0
        struc%oipcp=0
        struc%oips=0
        struc%omxcst=0
        struc%onrc=0
        struc%onrcp=0
        struc%opgfsl=0
        struc%opgfvl=0
        struc%opgord=0
        struc%opgplp=0
        struc%opos=0
        struc%ormax=0
        struc%osordn=0
        struc%ostni=0
        struc%ozos=0
        end subroutine uarray_init
       
        subroutine uarray_show(struc)
        use m_struc_def
        implicit none
        type(s_array):: struc
        write(*,*) struc
        end subroutine uarray_show
       
        subroutine uarray(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_array):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call uarray_init(struc)
        else
          write(*,*) 'uarray offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uarray
       
        subroutine ubz_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_bz):: struc
Cki        integer:: ubz_size
        struc%size=ubz_size()
        struc%def=0.0d0
        struc%dosw(:)=0.0d0
        struc%ef=0.0d0
        struc%efmax=0.0d0
        struc%fsmom=0.0d0
        struc%lcond(:)=0.0d0
        struc%lio=0
        struc%lmet=0
        struc%lmull=0
        struc%lopt=0
        struc%lshft(:)=0
        struc%n=0
        struc%ndos=0
        struc%nevmx=0
        struc%nkabc(:)=0
        struc%nkp=0
        struc%ntet=0
        struc%odos=0
        struc%oidtet=0
        struc%oipq=0
        struc%opdos=0
        struc%oqp=0
        struc%ostar=0
        struc%owtkp=0
        struc%range=0.0d0
        struc%semsh(:)=0.0d0
        struc%stnr=0.0d0
        struc%w=0.0d0
        struc%zval=0.0d0
        end subroutine ubz_init
       
        subroutine ubz_show(struc)
        use m_struc_def
        implicit none
        type(s_bz):: struc
        write(*,*) struc
        end subroutine ubz_show
       
        subroutine ubz(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_bz):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call ubz_init(struc)
        else
          write(*,*) 'ubz offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ubz
       
        subroutine uctrl_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_ctrl):: struc
Cki        integer:: uctrl_size
        struc%size=uctrl_size()
        struc%defm(:)=0.0d0
        struc%elin=0.0d0
        struc%lasa=0
        struc%lbas=0
        struc%lcd=0
        struc%lcgf=0
        struc%ldos=0
        struc%lfp=0
        struc%lfrce=0
        struc%lgen3=0
        struc%lham=0
        struc%lmet=0
        struc%lncol=0
        struc%loptc=0
        struc%lordn=0
        struc%lpgf(:)=0
        struc%lqp=0
        struc%lrel=0
        struc%lrs=0
        struc%lscr=0
        struc%lstonr(:)=0
        struc%lstr=0
        struc%lsx=0
        struc%ltb=0
        struc%lves=0
        struc%lxcf=0
        struc%maxit=0
        struc%mdprm(:)=0.0d0
        struc%modep(:)=0
        struc%nbas=0
        struc%nbasp=0
        struc%nclass=0
        struc%nesabc(:)=0
        struc%nitmv=0
        struc%nl=0
        struc%nmap=0
        struc%npl=0
        struc%nsite=0
        struc%nspec=0
        struc%nspin=0
        struc%nvario=0
        struc%omax1(:)=0.0d0
        struc%omax2(:)=0.0d0
        struc%pfloat=0
        struc%quit=0
        struc%rmaxes=0.0d0
        struc%rmines=0.0d0
        struc%sclwsr=0.0d0
        struc%sdmod=0
        struc%sdprm(:)=0.0d0
        struc%sdxsi(:)=0.0d0
        struc%smalit(:)=0
        struc%tol(:)=0.0d0
        struc%wsrmax=0.0d0
        struc%zbak(:)=0.0d0
        end subroutine uctrl_init
       
        subroutine uctrl_show(struc)
        use m_struc_def
        implicit none
        type(s_ctrl):: struc
        write(*,*) struc
        end subroutine uctrl_show
       
        subroutine uctrl(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_ctrl):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call uctrl_init(struc)
        else
          write(*,*) 'uctrl offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uctrl
       
        subroutine ugw_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_gw):: struc
Cki        integer:: ugw_size
        struc%size=ugw_size()
        struc%delre(:)=0.0d0
        struc%deltaw=0.0d0
        struc%deltax=0.0d0
        struc%ecuts=0.0d0
        struc%gcutb=0.0d0
        struc%gcutx=0.0d0
        struc%gsmear=0.0d0
        struc%lgw=0
        struc%lmxpb=0
        struc%mksig=0
        struc%nband=0
        struc%nime=0
        struc%nkabc(:)=0
        struc%pb1=0.0d0
        struc%pb2=0.0d0
        struc%pbtol=0.0d0
        struc%qoffp=0.0d0
        struc%rdsig=0
        end subroutine ugw_init
       
        subroutine ugw_show(struc)
        use m_struc_def
        implicit none
        type(s_gw):: struc
        write(*,*) struc
        end subroutine ugw_show
       
        subroutine ugw(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_gw):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call ugw_init(struc)
        else
          write(*,*) 'ugw offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ugw
       
        subroutine uham_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_ham):: struc
Cki        integer:: uham_size
        struc%size=uham_size()
        struc%alfsi=0.0d0
        struc%amgm=0.0d0
        struc%bandw=0
        struc%dabc(:)=0.0d0
        struc%ehf=0.0d0
        struc%ehk=0.0d0
        struc%elind=0.0d0
        struc%eterms(:)=0.0d0
        struc%hord=0
        struc%kmto(:)=0.0d0
        struc%lasa=0
        struc%ldham(:)=0
        struc%lgen3=0
        struc%lham=0
        struc%lmaxu=0
        struc%lmxax=0
        struc%lncol=0
        struc%lsig=0
        struc%ltb=0
        struc%lxcf=0
        struc%nbf=0
        struc%ndham=0
        struc%ndhrs=0
        struc%ndofH=0
        struc%neula=0
        struc%nkaph=0
        struc%nlibu=0
        struc%nmto=0
        struc%npwmin=0
        struc%npwpad=0
        struc%nqsig=0
        struc%obdots=0
        struc%oeula=0
        struc%ohrs=0
        struc%oiaxs=0
        struc%oindxo=0
        struc%olmxa=0
        struc%omagf=0
        struc%ontabs=0
        struc%ooffH=0
        struc%oqsig=0
        struc%oveps=0.0d0
        struc%pmax(:)=0.0d0
        struc%pmin(:)=0.0d0
        struc%pwemax=0.0d0
        struc%pwemin=0.0d0
        struc%pwmode=0
        struc%qpoff(:)=0.0d0
        struc%qss(:)=0.0d0
        struc%rsrnge=0.0d0
        struc%rsstol=0.0d0
        struc%seref=0.0d0
        struc%sigp(:)=0.0d0
        struc%thrpv=0.0d0
        struc%udiag=0
        end subroutine uham_init
       
        subroutine uham_show(struc)
        use m_struc_def
        implicit none
        type(s_ham):: struc
        write(*,*) struc
        end subroutine uham_show
       
        subroutine uham(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_ham):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call uham_init(struc)
        else
          write(*,*) 'uham offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uham
       
        subroutine ulat_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_lat):: struc
Cki        integer:: ulat_size
        struc%size=ulat_size()
        struc%alat=0.0d0
        struc%as=0.0d0
        struc%avw=0.0d0
        struc%awald=0.0d0
        struc%dist(:)=0.0d0
        struc%gam(:)=0.0d0
        struc%gmax=0.0d0
        struc%ldist=0
        struc%nabc(:)=0
        struc%ng=0
        struc%nkd=0
        struc%nkdmx=0
        struc%nkq=0
        struc%nkqmx=0
        struc%npgrp=0
        struc%nsgrp=0
        struc%oag=0
        struc%obgv=0
        struc%ocg=0
        struc%ocy=0
        struc%odlv=0
        struc%ogv=0
        struc%oidxcg=0
        struc%oips0=0
        struc%oistab=0
        struc%ojcg=0
        struc%okv=0
        struc%opos=0
        struc%oqlv=0
        struc%osymgr=0
        struc%plat(:)=0.0d0
        struc%plat0(:)=0.0d0
        struc%plat2(:)=0.0d0
        struc%plate(:)=0.0d0
        struc%platl(:)=0.0d0
        struc%platr(:)=0.0d0
        struc%qlat(:)=0.0d0
        struc%rpad=0.0d0
        struc%slat(:)=0.0d0
        struc%tol=0.0d0
        struc%tolft=0.0d0
        struc%vol=0.0d0
        end subroutine ulat_init
       
        subroutine ulat_show(struc)
        use m_struc_def
        implicit none
        type(s_lat):: struc
        write(*,*) struc
        end subroutine ulat_show
       
        subroutine ulat(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_lat):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call ulat_init(struc)
        else
          write(*,*) 'ulat offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ulat
       
        subroutine umix_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_mix):: struc
Cki        integer:: umix_size
        struc%size=umix_size()
        struc%b=0.0d0
        struc%bl=0.0d0
        struc%bv=0.0d0
        struc%elind=0.0d0
        struc%fn=0.0d0
        struc%kill=0
        struc%lxpot=0
        struc%mmix=0
        struc%mode=0
        struc%model=0
        struc%n=0
        struc%nitu=0
        struc%nmix=0
        struc%nsave=0
        struc%r(:)=0.0d0
        struc%rms1=0.0d0
        struc%rms2=0.0d0
        struc%tj(:)=0.0d0
        struc%tolu=0.0d0
        struc%umix=0.0d0
        struc%w(:)=0.0d0
        struc%wc=0.0d0
        end subroutine umix_init
       
        subroutine umix_show(struc)
        use m_struc_def
        implicit none
        type(s_mix):: struc
        write(*,*) struc
        end subroutine umix_show
       
        subroutine umix(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_mix):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call umix_init(struc)
        else
          write(*,*) 'umix offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine umix
       
        subroutine umove_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_move):: struc
Cki        integer:: umove_size
        struc%size=umove_size()
        struc%ct(:)=0.0d0
        struc%gyro=0.0d0
        struc%kt=0.0d0
        struc%modt(:)=0
        struc%nmodt=0
        struc%prmint(:)=0.0d0
        struc%tnow=0.0d0
        struc%ts=0.0d0
        struc%tsequ=0.0d0
        struc%tstot=0.0d0
        end subroutine umove_init
       
        subroutine umove_show(struc)
        use m_struc_def
        implicit none
        type(s_move):: struc
        write(*,*) struc
        end subroutine umove_show
       
        subroutine umove(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_move):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call umove_init(struc)
        else
          write(*,*) 'umove offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine umove
       
        subroutine uoptic_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_optic):: struc
Cki        integer:: uoptic_size
        struc%size=uoptic_size()
        struc%axes(:)=0
        struc%cll=0
        struc%cln=0
        struc%cls=0
        struc%clsite=0
        struc%esciss=0.0d0
        struc%mode=0
        struc%nchi2=0
        struc%ne=0
        struc%ocrng(:)=0
        struc%unrng(:)=0
        struc%window(:)=0.0d0
        end subroutine uoptic_init
       
        subroutine uoptic_show(struc)
        use m_struc_def
        implicit none
        type(s_optic):: struc
        write(*,*) struc
        end subroutine uoptic_show
       
        subroutine uoptic(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_optic):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call uoptic_init(struc)
        else
          write(*,*) 'uoptic offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uoptic
       
        subroutine uordn_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_ordn):: struc
Cki        integer:: uordn_size
        struc%size=uordn_size()
        struc%efre=0.0d0
        struc%mode=0
        struc%ncl=0
        struc%ndofh=0
        struc%oclp=0
        struc%oclssl=0
        struc%oiaxg=0
        struc%omapgv=0
        struc%ooffch=0
        struc%rmaxg(:)=0.0d0
        end subroutine uordn_init
       
        subroutine uordn_show(struc)
        use m_struc_def
        implicit none
        type(s_ordn):: struc
        write(*,*) struc
        end subroutine uordn_show
       
        subroutine uordn(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_ordn):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call uordn_init(struc)
        else
          write(*,*) 'uordn offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uordn
       
        subroutine upot_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_pot):: struc
Cki        integer:: upot_size
        struc%size=upot_size()
        struc%bfield(:)=0.0d0
        struc%nlma=0
        struc%nlml=0
        struc%nrhos=0
        struc%oaamom=0
        struc%obxc=0
        struc%odddpf=0
        struc%oddpf=0
        struc%oddpfr=0
        struc%odel=0
        struc%odpf=0
        struc%odpfr=0
        struc%ofes=0
        struc%ogma=0
        struc%ogmar=0
        struc%ogrrme=0
        struc%ohab=0
        struc%oivso=0
        struc%omad=0
        struc%oorhat=0
        struc%opalp=0
        struc%opapg=0
        struc%opdel=0
        struc%opf=0
        struc%opfnc=0
        struc%opfr=0
        struc%opmpol=0
        struc%opnu=0
        struc%opp=0
        struc%oppi=0
        struc%oppn=0
        struc%opprel=0
        struc%opti=0
        struc%oqc=0
        struc%oqmom=0
        struc%oqnu=0
        struc%oqpp=0
        struc%oqt=0
        struc%orhos=0
        struc%orhrmx=0
        struc%osab=0
        struc%osgw=0
        struc%osmpot=0
        struc%osmrho=0
        struc%osop=0
        struc%osoptc=0
        struc%osrout=0
        struc%otau=0
        struc%ovab=0
        struc%ovdif=0
        struc%oves=0
        struc%ovintr=0
        struc%ovrmax=0
        struc%ovshf=0
        struc%vconst(:)=0.0d0
        struc%vmtz=0.0d0
        struc%vmtz0=0.0d0
        end subroutine upot_init
       
        subroutine upot_show(struc)
        use m_struc_def
        implicit none
        type(s_pot):: struc
        write(*,*) struc
        end subroutine upot_show
       
        subroutine upot(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_pot):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call upot_init(struc)
        else
          write(*,*) 'upot offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine upot
       
        subroutine usite_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_site):: struc
Cki        integer:: usite_size
        struc%size=usite_size()
        struc%bfield(:)=0.0d0
        struc%clabel=0.0d0
        struc%class=0
        struc%cli=0
        struc%delta(:)=0.0d0
        struc%dpole(:)=0.0d0
        struc%eula(:)=0.0d0
        struc%force(:)=0.0d0
        struc%mpole=0.0d0
        struc%ndelta=0
        struc%norb=0
        struc%offh=0
        struc%orho=0
        struc%ov0=0
        struc%ov1=0
        struc%pl=0
        struc%plv=0
        struc%pnu(:)=0.0d0
        struc%pos(:)=0.0d0
        struc%pos0(:)=0.0d0
        struc%pz(:)=0.0d0
        struc%relax(:)=0
        struc%sid=0
        struc%spec=0
        struc%vel(:)=0.0d0
        struc%vshft=0.0d0
        end subroutine usite_init
       
        subroutine usite_show(struc)
        use m_struc_def
        implicit none
        type(s_site):: struc
        write(*,*) struc
        end subroutine usite_show
       
        subroutine usite(struc,offe,nli,is,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_site):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call usite_init(struc)
        else
          write(*,*) 'usite offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine usite
       
        subroutine uspec_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_spec):: struc
Cki        integer:: uspec_size
        struc%size=uspec_size()
        struc%z=0.0d0
        struc%mass=0.0d0
        struc%rmt=0.0d0
        struc%ntorb=0
        struc%naug=0
        struc%norb=0
        struc%rsmfa=0.0d0
        struc%rsma=0.0d0
        struc%rg=0.0d0
        struc%lmxa=0
        struc%kmxh=0
        struc%lmxl=0
        struc%kmxt=0
        struc%rsmv=0.0d0
        struc%norp=0
        struc%coreh=0.0d0
        struc%coreq(:)=0.0d0
        struc%a=0.0d0
        struc%nr=0
        struc%eref=0.0d0
        struc%etf=0.0d0
        struc%beta=0.0d0
        struc%lfoca=0
        struc%ctail=0.0d0
        struc%etail=0.0d0
        struc%name=0.0d0
        struc%orhoc=0
        struc%stc=0.0d0
        struc%lmxb=0
        struc%lmxf=0
        struc%stni=0.0d0
        struc%stnm=0.0d0
        struc%rham=0.0d0
        struc%rfoca=0.0d0
        struc%dv=0.0d0
        struc%mxcst=0
        struc%group=0
        struc%grp2=0
        struc%nxi=0
        struc%qc=0.0d0
        struc%lmxpb=0
        struc%pb1=0.0d0
        struc%pb2=0.0d0
        struc%colxbs(:)=0.0d0
        struc%lxi=0
        struc%radxbs=0.0d0
        struc%rcut=0.0d0
        struc%rint=0.0d0
        struc%eh3=0.0d0
        struc%rs3=0.0d0
        struc%vmtz=0.0d0
        struc%kmxv=0
        struc%rcfa(:)=0.0d0
        struc%p(:)=0.0d0
        struc%q(:)=0.0d0
        struc%alpha(:)=0.0d0
        struc%idmod(:)=0
        struc%idxdn(:)=0
        struc%hcr(:)=0.0d0
        struc%exi(:)=0.0d0
        struc%ngcut(:)=0
        struc%qpol(:)=0.0d0
        struc%chfa(:)=0.0d0
        struc%orbp(:)=0.0d0
        struc%enu(:)=0.0d0
        struc%pz(:)=0.0d0
        struc%idu(:)=0
        struc%uh(:)=0.0d0
        struc%jh(:)=0.0d0
        struc%iq1(:)=0
        struc%ivso(:)=0
        struc%ehvl(:)=0.0d0
        struc%vso(:)=0.0d0
        end subroutine uspec_init
       
        subroutine uspec_show(struc)
        use m_struc_def
        implicit none
        type(s_spec):: struc
        write(*,*) struc
        end subroutine uspec_show
       
        subroutine uspec(struc,offe,nli,is,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_spec):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call uspec_init(struc)
        else
          write(*,*) 'uspec offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine uspec
       
        subroutine ustr_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_str):: struc
Cki        integer:: ustr_size
        struc%size=ustr_size()
        struc%adec=0.0d0
        struc%amode=0
        struc%drwats=0.0d0
        struc%iinv(:)=0.0d0
        struc%ivl=0
        struc%kaps(:)=0.0d0
        struc%lequiv=0
        struc%lmaxw=0
        struc%loka=0
        struc%lshow=0
        struc%mxnbr=0
        struc%nalf=0
        struc%nbisi(:)=0
        struc%ncupl=0
        struc%ndust=0
        struc%nitab=0
        struc%nkaps=0
        struc%nttab=0
        struc%oadot=0
        struc%oalp=0
        struc%oiax=0
        struc%ontab=0
        struc%os=0
        struc%osdot=0
        struc%rfit=0.0d0
        struc%rmax=0.0d0
        struc%rmaxg=0.0d0
        struc%skmsh(:)=0.0d0
        struc%wztcf=0.0d0
        end subroutine ustr_init
       
        subroutine ustr_show(struc)
        use m_struc_def
        implicit none
        type(s_str):: struc
        write(*,*) struc
        end subroutine ustr_show
       
        subroutine ustr(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_str):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call ustr_init(struc)
        else
          write(*,*) 'ustr offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine ustr
       
        subroutine utb_init(struc)
        use m_struc_def
        use m_struc_func
        implicit none
        type(s_tb):: struc
Cki        integer:: utb_size
        struc%size=utb_size()
        struc%alam=0.0d0
        struc%alsc=0.0d0
        struc%ebfit(:)=0.0d0
        struc%fmode=0
        struc%nbfit(:)=0
        struc%rmfit(:)=0.0d0
        end subroutine utb_init
       
        subroutine utb_show(struc)
        use m_struc_def
        implicit none
        type(s_tb):: struc
        write(*,*) struc
        end subroutine utb_show
       
        subroutine utb(struc,offe,nli,off,cast,nelt)
        use m_struc_def
        implicit none
        type(s_tb):: struc
        integer:: offe(1),nli,is,off,cast,nelt
        integer::val
        if (offe(1).eq.-1) then
          call utb_init(struc)
        else
          write(*,*) 'utb offe=',offe(1),' not suppoted'
          stop
        endif
        end subroutine utb
       
