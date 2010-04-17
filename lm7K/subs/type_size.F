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

        subroutine uarray_init(struc)
        use m_struc_def
         use m_struc_func
        implicit none
        type(s_array):: struc
        integer:: uarray_size
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
        integer:: ubz_size
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
        integer:: uctrl_size
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
        integer:: ugw_size
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
        integer:: uham_size
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
        integer:: ulat_size
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
        integer:: umix_size
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
        integer:: umove_size
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
        integer:: uoptic_size
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
        integer:: uordn_size
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
        integer:: upot_size
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
        integer:: usite_size
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
        integer:: uspec_size
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
        integer:: ustr_size
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
        integer:: utb_size
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
       
