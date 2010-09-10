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
                   integer:: nbas , nbasp , nbaspp , nsp , nkd , nkq , igets , ngabc(3) 
     .             , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ng , ngrp , iprint 
                   real(8),pointer :: rv_p_opos(:)

                  real(8),pointer :: rv_p_og(:)

                 integer,pointer :: iv_p_okv(:)

                real(8),pointer :: rv_p_ogv(:)

               real(8),pointer :: rv_p_odlv(:)

              complex(8),pointer :: zv_p_obgv(:)

             real(8),pointer :: rv_p_oag(:)

            type(s_rv1),pointer :: sv_p_oorhat(:)

           real(8),pointer :: rv_p_omad(:)

          complex(8),pointer :: zv_p_osmrho(:)

         complex(8),pointer :: zv_p_osmpot(:)

        integer,pointer :: iv_p_oips0(:)

       real(8),pointer :: rv_p_oqlv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(9),gmax,xx
C ... Heap
       integer ::iwdummy 

      call info0(50,1,0,' supot : potential setup ... ')


      nbas=sctrl%nbas
      nbasp=sctrl%nbasp
      nsp=sctrl%nspin


      alat=slat%alat
      vol=slat%vol
       rv_p_opos => slat%rv_p_opos 



      awald=slat%awald
      nkd=slat%nkd
      nkq=slat%nkq
       rv_p_odlv => slat%rv_p_odlv 

       rv_p_oqlv => slat%rv_p_oqlv 



      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
         allocate(rv_p_omad(nbaspp**2))
         rv_p_omad(:)=0.0d0

            call madmat ( nbaspp , rv_p_opos , awald , alat , vol , rv_p_odlv 
     .      , nkd , rv_p_oqlv , nkq , rv_p_omad ) 





         spot%rv_p_omad => rv_p_omad 


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
         allocate(rv_p_ogv(abs(ngmx*3)))
c         if (ngmx*3<0) rv_p_ogv(:)=0.0d0

         allocate(iv_p_okv(abs(ngmx*3)))
c         if (ngmx*3<0) iv_p_okv(:)=0

           call gvlist ( alat , plat , iwdummy , n1 , n2 , n3 , gmax , 8 
     .     , ngmx , ng , iv_p_okv , rv_p_ogv , xx , xx ) 



        if (ng .ne. ngmx) call rx('supot: bug in gvlist')


        slat%ng=ng
         slat%rv_p_ogv => rv_p_ogv 

         slat%iv_p_okv => iv_p_okv 

        slat%gmax=gmax
        i_copy_size=size(slat%nabc)
        call ii8copy(i_copy_size,ngabc,1,slat%nabc,1)

C       call prmx('gv',w(ogv),ngmx,ng,3)

         allocate(zv_p_osmrho(k1*k2*k3*nsp))
c         if (-k1*k2*k3*nsp<0) zv_p_osmrho(:)=0.0d0
         zv_p_osmrho(:)=0d0

         allocate(zv_p_osmpot(k1*k2*k3*nsp))
         zv_p_osmpot(:)=0.0d0

         allocate(sv_p_oorhat(3*nbaspp))
cki         sv_p_oorhat(:)=0
         do i=1,3*nbaspp
            nullify(sv_p_oorhat(i)%v)
         enddo


         spot%zv_p_osmrho => zv_p_osmrho 

         spot%zv_p_osmpot => zv_p_osmpot 

         spot%sv_p_oorhat => sv_p_oorhat 



C   ... Setup for symmetrization
         allocate(iv_p_oips0(ng))
         iv_p_oips0(:)=0.0d0

         allocate(zv_p_obgv(ng))
         zv_p_obgv(:)=0.0d0


        ngrp=slat%nsgrp
         rv_p_og => slat%rv_p_osymgr 

         rv_p_oag => slat%rv_p_oag 


             call sgvsym ( ngrp , rv_p_og , rv_p_oag , ng , rv_p_ogv , iv_p_oips0 
     .       , zv_p_obgv ) 






         slat%iv_p_oips0 => iv_p_oips0 

         slat%zv_p_obgv => zv_p_obgv 



      endif
#endif
      end subroutine supot


