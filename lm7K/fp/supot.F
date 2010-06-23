Cchp1 cmdvar= real(8),rv_p_,oqlv
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
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_pot)::spot

C ... Local parameters
Cchp1       integer nbas,nbasp,nbaspp,nsp,omad,opos,nkd,nkq,odlv,oqlv,igets,
Cchp1      .ngabc(3),n1,n2,n3,k1,k2,k3,ngmx,ogv,okv,ng,osmrho,osmpot,oorhat,
Cchp1      .oips0,obgv,ngrp,og,oag,iprint
       integer:: nbas , nbasp , nbaspp , nsp , omad , opos , nkd , nkq 
     . , odlv , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx 
     . , ogv , okv , ng , osmrho , osmpot , oorhat , oips0 , obgv , 
     . ngrp , og , oag , iprint 
       real(8),pointer :: rv_p_oqlv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(9),gmax,xx
C ... Heap
      integer w(1)
      common /w/ w
      call info0(50,1,0,' supot : potential setup ... ')


      nbas=sctrl%nbas
      nbasp=sctrl%nbasp
      nsp=sctrl%nspin


      alat=slat%alat
      vol=slat%vol
      opos=slat%opos


      awald=slat%awald
      nkd=slat%nkd
      nkq=slat%nkq
      odlv=slat%odlv
Cchp1       oqlv=slat%oqlv
       rv_p_oqlv => slat%rv_p_oqlv 



      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
        call defdr(omad,nbaspp**2)
Cchp1         call madmat(nbaspp,w(opos),awald,alat,vol,w(odlv),nkd,
Cchp1      .  w(oqlv),nkq,w(omad))
         call madmat ( nbaspp , w ( opos ) , awald , alat , vol , w ( 
     .   odlv ) , nkd , rv_p_oqlv , nkq , w ( omad ) ) 


        spot%omad=omad

      endif

#if FP
C --- Setup for FT charge density, potential representation ---
      if ( iand(4,int(sctrl%lcd)) .ne. 0 ) then


        alat=slat%alat
        i_copy_size=size(slat%plat)
        call dcopy(i_copy_size,slat%plat,1,plat,1)
        i_copy_size=size(slat%nabc)
        call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
        gmax=slat%gmax


C   ... Generate energy cutoff gmax or n1..n3
c        call pshpr(iprint()-0)
        call mshsiz(alat,plat,0,gmax,ngabc,ngmx)
        call fftz30(n1,n2,n3,k1,k2,k3)

C   ... Make list of lattice vectors within cutoff
        call defrr(ogv, ngmx*3)
        call defi (okv, ngmx*3)
        call gvlist(alat,plat,w,n1,n2,n3,gmax,8,ngmx,ng,w(okv),
     .  w(ogv),xx,xx)
        if (ng .ne. ngmx) call rx('supot: bug in gvlist')


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

        spot%osmrho=osmrho
        spot%osmpot=osmpot
        spot%oorhat=oorhat


C   ... Setup for symmetrization
        call defi (oips0,   ng)
        call defcc(obgv,    ng)

        ngrp=slat%nsgrp
        og=slat%osymgr
        oag=slat%oag

        call sgvsym(ngrp,w(og),w(oag),ng,w(ogv),w(oips0),w(obgv))

        slat%oips0=oips0
        slat%obgv=obgv


      endif
#endif
      end subroutine supot


