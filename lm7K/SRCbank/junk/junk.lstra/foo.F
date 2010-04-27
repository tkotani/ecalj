        module m_struc_def
       
        type s_ri
         integer:: t
         integer(8):: i
         real(8):: r
        end type
       
        type s_bz
        real(8):: size
        real(8) :: def 
        real(8) :: dosw(2) 
        real(8) :: ef 
        real(8) :: efmax 
        real(8) :: fsmom 
        real(8) :: lcond(4) 
        integer(8) :: lio 
        integer(8) :: lmet 
        integer(8) :: lmull 
        integer(8) :: lopt 
        integer(8) :: lshft(3) 
        integer(8) :: n 
        integer(8) :: ndos 
        integer(8) :: nevmx 
        integer(8) :: nkabc(3) 
        integer(8) :: nkp 
        integer(8) :: ntet 
        integer(8) :: odos 
        integer(8) :: oidtet 
        integer(8) :: oipq 
        integer(8) :: opdos 
        integer(8) :: oqp 
        integer(8) :: ostar 
        integer(8) :: owtkp 
        real(8) :: range 
        real(8) :: semsh(10) 
        real(8) :: stnr 
        real(8) :: w 
        real(8) :: zval 
        end type s_bz
       
        type s_array
        real(8):: size
        integer(8) :: ncl 
        integer(8) :: nclasp 
        integer(8) :: nofgl 
        integer(8) :: nofgr 
        integer(8) :: npadl 
        integer(8) :: npadr 
        integer(8) :: npl 
        integer(8) :: oclabl 
        integer(8) :: ocllst 
        integer(8) :: oclp 
        integer(8) :: oclssl 
        integer(8) :: ogroup 
        integer(8) :: ogrp2 
        integer(8) :: ohave 
        integer(8) :: oics 
        integer(8) :: oifrlx 
        integer(8) :: oipc 
        integer(8) :: oipcp 
        integer(8) :: oips 
        integer(8) :: omxcst 
        integer(8) :: onrc 
        integer(8) :: onrcp 
        integer(8) :: opgfsl 
        integer(8) :: opgfvl 
        integer(8) :: opgord 
        integer(8) :: opgplp 
        integer(8) :: opos 
        integer(8) :: ormax 
        integer(8) :: osordn 
        integer(8) :: ostni 
        integer(8) :: ozos 
        end type s_array
       
        type s_ctrl
        real(8):: size
        real(8) :: defm(6) 
        real(8) :: elin 
        integer(8) :: lasa 
        integer(8) :: lbas 
        integer(8) :: lcd 
        integer(8) :: lcgf 
        integer(8) :: ldos 
        integer(8) :: lfp 
        integer(8) :: lfrce 
        integer(8) :: lgen3 
        integer(8) :: lham 
        integer(8) :: lmet 
        integer(8) :: lncol 
        integer(8) :: loptc 
        integer(8) :: lordn 
        integer(8) :: lpgf(2) 
        integer(8) :: lqp 
        integer(8) :: lrel 
        integer(8) :: lrs 
        integer(8) :: lscr 
        integer(8) :: lstonr(3) 
        integer(8) :: lstr 
        integer(8) :: lsx 
        integer(8) :: ltb 
        integer(8) :: lves 
        integer(8) :: lxcf 
        integer(8) :: maxit 
        real(8) :: mdprm(6) 
        integer(8) :: modep(3) 
        integer(8) :: nbas 
        integer(8) :: nbasp 
        integer(8) :: nclass 
        integer(8) :: nesabc(3) 
        integer(8) :: nitmv 
        integer(8) :: nl 
        integer(8) :: nmap 
        integer(8) :: npl 
        integer(8) :: nsite 
        integer(8) :: nspec 
        integer(8) :: nspin 
        integer(8) :: nvario 
        real(8) :: omax1(3) 
        real(8) :: omax2(3) 
        integer(8) :: pfloat 
        integer(8) :: quit 
        real(8) :: rmaxes 
        real(8) :: rmines 
        real(8) :: sclwsr 
        integer(8) :: sdmod 
        real(8) :: sdprm(5) 
        real(8) :: sdxsi(4) 
        integer(8) :: smalit(2) 
        real(8) :: tol(3) 
        real(8) :: wsrmax 
        real(8) :: zbak(2) 
        end type s_ctrl
       
        type s_gw
        real(8):: size
        real(8) :: delre(2) 
        real(8) :: deltaw 
        real(8) :: deltax 
        real(8) :: ecuts 
        real(8) :: gcutb 
        real(8) :: gcutx 
        real(8) :: gsmear 
        integer(8) :: lgw 
        integer(8) :: lmxpb 
        integer(8) :: mksig 
        integer(8) :: nband 
        integer(8) :: nime 
        integer(8) :: nkabc(3) 
        real(8) :: pb1 !string
        real(8) :: pb2 !string
        real(8) :: pbtol 
        real(8) :: qoffp 
        integer(8) :: rdsig 
        end type s_gw
       
        type s_ham
        real(8):: size
        real(8) :: alfsi 
        real(8) :: amgm 
        integer(8) :: bandw 
        real(8) :: dabc(3) 
        real(8) :: ehf 
        real(8) :: ehk 
        real(8) :: elind 
        real(8) :: eterms(20) 
        integer(8) :: hord 
        real(8) :: kmto(6) 
        integer(8) :: lasa 
        integer(8) :: ldham(16) 
        integer(8) :: lgen3 
        integer(8) :: lham 
        integer(8) :: lmaxu 
        integer(8) :: lmxax 
        integer(8) :: lncol 
        integer(8) :: lsig 
        integer(8) :: ltb 
        integer(8) :: lxcf 
        integer(8) :: nbf 
        integer(8) :: ndham 
        integer(8) :: ndhrs 
        integer(8) :: ndofH 
        integer(8) :: neula 
        integer(8) :: nkaph 
        integer(8) :: nlibu 
        integer(8) :: nmto 
        integer(8) :: npwmin 
        integer(8) :: npwpad 
        integer(8) :: nqsig 
        integer(8) :: obdots 
        integer(8) :: oeula 
        integer(8) :: ohrs 
        integer(8) :: oiaxs 
        integer(8) :: oindxo 
        integer(8) :: olmxa 
        integer(8) :: omagf 
        integer(8) :: ontabs 
        integer(8) :: ooffH 
        integer(8) :: oqsig 
        real(8) :: oveps 
        real(8) :: pmax(10) 
        real(8) :: pmin(10) 
        real(8) :: pwemax 
        real(8) :: pwemin 
        integer(8) :: pwmode 
        real(8) :: qpoff(3) 
        real(8) :: qss(4) 
        real(8) :: rsrnge 
        real(8) :: rsstol 
        real(8) :: seref 
        real(8) :: sigp(10) 
        real(8) :: thrpv 
        integer(8) :: udiag 
        end type s_ham
       
        type s_lat
        real(8):: size
        real(8) :: alat 
        real(8) :: as 
        real(8) :: avw 
        real(8) :: awald 
        real(8) :: dist(9) 
        real(8) :: gam(4) 
        real(8) :: gmax 
        integer(8) :: ldist 
        integer(8) :: nabc(3) 
        integer(8) :: ng 
        integer(8) :: nkd 
        integer(8) :: nkdmx 
        integer(8) :: nkq 
        integer(8) :: nkqmx 
        integer(8) :: npgrp 
        integer(8) :: nsgrp 
        integer(8) :: oag 
        integer(8) :: obgv 
        integer(8) :: ocg 
        integer(8) :: ocy 
        integer(8) :: odlv 
        integer(8) :: ogv 
        integer(8) :: oidxcg 
        integer(8) :: oips0 
        integer(8) :: oistab 
        integer(8) :: ojcg 
        integer(8) :: okv 
        integer(8) :: opos 
        integer(8) :: oqlv 
        integer(8) :: osymgr 
        real(8) :: plat(9) 
        real(8) :: plat0(9) 
        real(8) :: plat2(9) 
        real(8) :: plate(9) 
        real(8) :: platl(9) 
        real(8) :: platr(9) 
        real(8) :: qlat(9) 
        real(8) :: rpad 
        real(8) :: slat(9) 
        real(8) :: tol 
        real(8) :: tolft 
        real(8) :: vol 
        end type s_lat
       
        type s_mix
        real(8):: size
        real(8) :: b 
        real(8) :: bl 
        real(8) :: bv 
        real(8) :: elind 
        real(8) :: fn !string
        integer(8) :: kill 
        integer(8) :: lxpot 
        integer(8) :: mmix 
        integer(8) :: mode 
        integer(8) :: model 
        integer(8) :: n 
        integer(8) :: nitu 
        integer(8) :: nmix 
        integer(8) :: nsave 
        real(8) :: r(3) !string
        real(8) :: rms1 
        real(8) :: rms2 
        real(8) :: tj(10) 
        real(8) :: tolu 
        real(8) :: umix 
        real(8) :: w(3) 
        real(8) :: wc 
        end type s_mix
       
        type s_move
        real(8):: size
        real(8) :: ct(3) 
        real(8) :: gyro 
        real(8) :: kt 
        integer(8) :: modt(3) 
        integer(8) :: nmodt 
        real(8) :: prmint(20) 
        real(8) :: tnow 
        real(8) :: ts 
        real(8) :: tsequ 
        real(8) :: tstot 
        end type s_move
       
        type s_optic
        real(8):: size
        integer(8) :: axes(18) 
        integer(8) :: cll 
        integer(8) :: cln 
        integer(8) :: cls 
        integer(8) :: clsite 
        real(8) :: esciss 
        integer(8) :: mode 
        integer(8) :: nchi2 
        integer(8) :: ne 
        integer(8) :: ocrng(2) 
        integer(8) :: unrng(2) 
        real(8) :: window(2) 
        end type s_optic
       
        type s_ordn
        real(8):: size
        real(8) :: efre 
        integer(8) :: mode 
        integer(8) :: ncl 
        integer(8) :: ndofh 
        integer(8) :: oclp 
        integer(8) :: oclssl 
        integer(8) :: oiaxg 
        integer(8) :: omapgv 
        integer(8) :: ooffch 
        real(8) :: rmaxg(3) 
        end type s_ordn
       
        type s_pot
        real(8):: size
        real(8) :: bfield(4) 
        integer(8) :: nlma 
        integer(8) :: nlml 
        integer(8) :: nrhos 
        integer(8) :: oaamom 
        integer(8) :: obxc 
        integer(8) :: odddpf 
        integer(8) :: oddpf 
        integer(8) :: oddpfr 
        integer(8) :: odel 
        integer(8) :: odpf 
        integer(8) :: odpfr 
        integer(8) :: ofes 
        integer(8) :: ogma 
        integer(8) :: ogmar 
        integer(8) :: ogrrme 
        integer(8) :: ohab 
        integer(8) :: oivso 
        integer(8) :: omad 
        integer(8) :: oorhat 
        integer(8) :: opalp 
        integer(8) :: opapg 
        integer(8) :: opdel 
        integer(8) :: opf 
        integer(8) :: opfnc 
        integer(8) :: opfr 
        integer(8) :: opmpol 
        integer(8) :: opnu 
        integer(8) :: opp 
        integer(8) :: oppi 
        integer(8) :: oppn 
        integer(8) :: opprel 
        integer(8) :: opti 
        integer(8) :: oqc 
        integer(8) :: oqmom 
        integer(8) :: oqnu 
        integer(8) :: oqpp 
        integer(8) :: oqt 
        integer(8) :: orhos 
        integer(8) :: orhrmx 
        integer(8) :: osab 
        integer(8) :: osgw 
        integer(8) :: osmpot 
        integer(8) :: osmrho 
        integer(8) :: osop 
        integer(8) :: osoptc 
        integer(8) :: osrout 
        integer(8) :: otau 
        integer(8) :: ovab 
        integer(8) :: ovdif 
        integer(8) :: oves 
        integer(8) :: ovintr 
        integer(8) :: ovrmax 
        integer(8) :: ovshf 
        real(8) :: vconst(3) 
        real(8) :: vmtz 
        real(8) :: vmtz0 
        end type s_pot
       
        end module m_struc_def
