Cgetarg...info...           structure ['supot', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['supot', 'lat', 'slat']
Cgetarg...info...           structure ['supot', 'pot', 'spot']
Cgetarg...info...           use_to_add ['supot', 'm_struc_def']
#define FP 1
      subroutine supot(mode,sctrl,slat,spot)
       
       use m_struc_def  !Cgetarg

C- Initialization for potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make Madelung matrix for monopoles (ASA)
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: nbas nbasp nspin lcd,4
Ci   slat  :struct containing information about the lattice
Ci     Elts read: alat vol opos awald nkd nkq odlv oqlv plat nabc gmax
Ci                nsgrp osymgr oag
Ci     Stored:    ng ogv okv gmax nabc oips0 obgv
Ci   spot  :struct containing information about the potential
Ci     Stored:    omad osmrho osmpot oorhat
Co Outputs
Cr Remarks
Cr   Permanent arrays are created to handle the potential, depending
Cr   on the mode
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode
Cgetarg       double precision sctrl(1),slat(1),spot(1)
       type(s_ctrl)::sctrl
       type(s_lat)::slat
       type(s_pot)::spot

C ... Local parameters
      integer nbas,nbasp,nbaspp,nsp,omad,opos,nkd,nkq,odlv,oqlv,igets,
     .  ngabc(3),n1,n2,n3,k1,k2,k3,ngmx,ogv,okv,ng,osmrho,osmpot,oorhat,
     .  oips0,obgv,ngrp,og,oag,iprint
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(9),gmax,xx
C ... Heap
      integer w(1)
      common /w/ w
      call info0(50,1,0,' supot : potential setup ... ')

Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nbasp , integer(8) :: sctrl%nbasp 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg       call upack('ctrl nbas nbasp nspin',sctrl,nbas,nbasp,nsp,0,0)
       
       nbas=sctrl%nbas
       nbasp=sctrl%nbasp
       nsp=sctrl%nspin

Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg...info...           integer :: opos , integer(8) :: slat%opos 1
Cgetarg       call upack('lat alat vol opos',slat,alat,vol,opos,0,0)
       
       alat=slat%alat
       vol=slat%vol
       opos=slat%opos

Cgetarg...info...           real(8) :: awald , real(8) :: slat%awald 1
Cgetarg...info...           integer :: nkd , integer(8) :: slat%nkd 1
Cgetarg...info...           integer :: nkq , integer(8) :: slat%nkq 1
Cgetarg...info...           integer :: odlv , integer(8) :: slat%odlv 1
Cgetarg...info...           integer :: oqlv , integer(8) :: slat%oqlv 1
Cgetarg       call upack('lat awald nkd nkq odlv oqlv',slat,awald,nkd,nkq,odlv,
Cgetarg      .  oqlv)
       
       awald=slat%awald
       nkd=slat%nkd
       nkq=slat%nkq
       odlv=slat%odlv
       oqlv=slat%oqlv


      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
        call defdr(omad,nbaspp**2)
        call madmat(nbaspp,w(opos),awald,alat,vol,w(odlv),nkd,
     .    w(oqlv),nkq,w(omad))
Cgetarg...info...           integer :: omad , integer(8) :: spot%omad 1
Cgetarg         call pack1('pot omad',spot,omad)
         
         spot%omad=omad 

      endif

#if FP
C --- Setup for FT charge density, potential representation ---
Cgetarg       if (igets('ctrl lcd,4',sctrl) .ne. 0) then
       if ( iand(4,int(sctrl%lcd)) .ne. 0 ) then 

Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(9) , real(8) :: slat%plat 9
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           real(8) :: gmax , real(8) :: slat%gmax 1
Cgetarg         call upack('lat alat plat nabc gmax',slat,alat,plat,ngabc,gmax,
Cgetarg      .    0)
         
         alat=slat%alat
         i_copy_size=size(slat%plat) 
         call dcopy(i_copy_size,slat%plat,1,plat,1) 
         i_copy_size=size(slat%nabc) 
         call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
         gmax=slat%gmax


C   ... Generate energy cutoff gmax or n1..n3
c        call pshpr(iprint()-0)
        call mshsiz(alat,plat,0,gmax,ngabc,ngmx)
c        call fftz30(n1,n2,n3,k1,k2,k3)
        k1 = n1
        k2 = n2
        k3 = n3

C   ... Make list of lattice vectors within cutoff
        call defrr(ogv, ngmx*3)
        call defi (okv, ngmx*3)
        call gvlist(alat,plat,w,n1,n2,n3,gmax,8,ngmx,ng,w(okv),
     .    w(ogv),xx,xx)
        if (ng .ne. ngmx) call rx('supot: bug in gvlist')

Cgetarg...info...           integer :: ng , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           integer :: okv , integer(8) :: slat%okv 1
Cgetarg...info...           real(8) :: gmax , real(8) :: slat%gmax 1
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg         call pack5('lat ng ogv okv gmax nabc',slat,ng,ogv,okv,gmax,
Cgetarg      .    ngabc)
         
         slat%ng=ng 
         slat%ogv=ogv 
         slat%okv=okv 
         slat%gmax=gmax 
         i_copy_size=size(slat%nabc) 
         call ii8copy(i_copy_size,ngabc,1,slat%nabc,1) 

C       call prmx('gv',w(ogv),ngmx,ng,3)

        call defcc(osmrho, -k1*k2*k3*nsp)
        call defcc(osmpot, -k1*k2*k3*nsp)
        call defi(oorhat,3*nbaspp)
Cgetarg...info...           integer :: osmrho , integer(8) :: spot%osmrho 1
Cgetarg...info...           integer :: osmpot , integer(8) :: spot%osmpot 1
Cgetarg...info...           integer :: oorhat , integer(8) :: spot%oorhat 1
Cgetarg         call pack5('pot osmrho osmpot oorhat',spot,osmrho,osmpot,oorhat,
Cgetarg      .    0,0)
         
         spot%osmrho=osmrho 
         spot%osmpot=osmpot 
         spot%oorhat=oorhat 


C   ... Setup for symmetrization
        call defi (oips0,   ng)
        call defcc(obgv,    ng)
Cgetarg...info...           integer :: ngrp , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: og , integer(8) :: slat%osymgr 1
Cgetarg...info...           integer :: oag , integer(8) :: slat%oag 1
Cgetarg         call upack('lat nsgrp osymgr oag',slat,ngrp,og,oag,0,0)
         
         ngrp=slat%nsgrp
         og=slat%osymgr
         oag=slat%oag

        call sgvsym(ngrp,w(og),w(oag),ng,w(ogv),w(oips0),w(obgv))
Cgetarg...info...           integer :: oips0 , integer(8) :: slat%oips0 1
Cgetarg...info...           integer :: obgv , integer(8) :: slat%obgv 1
Cgetarg         call pack2('lat oips0 obgv',slat,oips0,obgv)
         
         slat%oips0=oips0 
         slat%obgv=obgv 


      endif
#endif
Cgetarg       end
       end subroutine supot 


