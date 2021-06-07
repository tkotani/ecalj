      module m_lmfgw
      use m_struc_def,only: s_rv1
      integer,protected:: jobgw=-999
      type(s_rv1),allocatable,protected  :: sv_p_oppix(:,:)
      type(s_rv1),allocatable,protected  :: sv_p_otaux(:,:)
      type(s_rv1),allocatable,protected  :: sv_p_osigx(:,:)
      complex(8),allocatable,protected   :: spotx(:,:,:,:)
      
      contains
      subroutine m_lmfgw_init()
      use m_density,only: zv_a_osmrho=>osmrho,sv_a_oorhat=>orhoat
      use m_lmfinit,only:nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nab,nppn,sspec=>v_sspec,ssite=>v_ssite,
     & nkap0,alat=>lat_alat
      use m_supot,only: k1,k2,k3
      use m_MPItk,only: master_mpi,mlog
      use m_suham,only: nchan=>pot_nlma, nvl=>pot_nlml,ham_ldham
      use m_lattic,only: qlat=>lat_qlat, plat=>lat_plat
      use m_suham,only: iv_a_oindxo,ham_ldham, ndham=>ham_ndham
      use m_mkpot,only: mkpot
      integer:: i,lfrzw,i1,i_copy_size,ib,ifildima,ificlass,ipqn,is,isp,konf,l,ldim,lmaxa,ifilnlaindx,
     & lmxax,ifinlaindx,ndima,npqn,mxint
      integer:: ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb,nat,ipb(nbas)
      real(8) ::qbz,vconst,vrmt(nbas),qval,qsc
      complex(8),allocatable:: vorb(:,:,:,:)
      real(8) ,allocatable :: fes1_rv(:), gpot0_rv(:), hab_rv(:),ppnl_rv(:,:,:,:), qmom_rv(:)
      real(8) ,allocatable :: sab_rv(:), vab_rv(:),vval_rv(:),sav_rv(:),bas(:,:)
      integer,allocatable:: lmxa(:),konft(:,:,:),ipc(:),ipcx(:)
      real(8) ::xx(5),pnu(n0,2),pnz(n0,2)
      character lsym(0:n0-1)*1, lorb(3)*1, dig(9)*1, strn4*4
      integer:: idxdn(n0,nkap0),fopnx
      data lsym /'s','p','d','f','g','5','6','7','8','9'/
      data lorb /'p','d','l'/
      data dig /'1','2','3','4','5','6','7','8','9'/
      
c      qbg = ctrl_zbak(1) !homogenious background charge
      lfrzw = 0
      if(ham_frzwf) lfrzw = 1  !freeze all augmentation wave
      if(master_mpi) then
         write(stdo,*)
         write(stdo,*) '=== lmfgw-MPIK: Choose one of following jobs: ==='
         write(stdo,*) '   0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima'
         write(stdo,*) '   1 : GW setup mode'
         write(stdo,*) ' job 0 or 1?'
         read (5,*) jobgw
      endif
      call mpibc1(jobgw,1,2, mlog,'bndfp','jobgw')
      if(.not.( jobgw==0.or.jobgw==1)) call rx(' lmf-MPIK is only for jobgw=0,1 now')
      if(jobgw/=0) then
         write(stdo,"(a)")' Making one-particle potential without XC part ...'
         allocate( sv_p_osigx(3,nbas), sv_p_otaux(3,nbas), sv_p_oppix(3,nbas))
         call dfaugm(sv_p_osigx, sv_p_otaux, sv_p_oppix) !for sig,tau,ppi without XC(LDA)
         allocate( spotx(k1,k2,k3,nsp)) !smooth potential without XC
         allocate( vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
         allocate(  qmom_rv(nvl))
         allocate( gpot0_rv(nvl))
         allocate( vval_rv(nchan))
         allocate(  hab_rv(nab*n0*nsp*nbas))
         allocate(  vab_rv(nab*n0*nsp*nbas))
         allocate(  sab_rv(nab*n0*nsp*nbas))
         allocate( ppnl_rv(nppn,n0,nsp,nbas))
         allocate( fes1_rv(3*nbas))
         spotx=0d0
         i = 1 + 10*lfrzw + 100 ! Adding 100 means excluding XC(LDA) part. nolxc=T
!     We obtain  sv_p_osigx , sv_p_otaux , sv_p_oppix, smpotx  without XC by the 'call mkpot'.
!     Other output are overwritten by the 2nd call of mkpot 
         call mkpot (0,         !nbas , ssite , sspec ,   0 , lcplxp , k1 , k2 , k3 ,
     i        zv_a_osmrho , sv_a_oorhat , qbg  
     o        , spotx, qmom_rv , vconst , vrmt , sv_p_osigx , sv_p_otaux , sv_p_oppix 
     o        , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o        , vval_rv , fes1_rv , i , vorb ) 
         deallocate(vorb,qmom_rv,gpot0_rv,vval_rv,hab_rv,vab_rv,sab_rv,ppnl_rv,fes1_rv)
      endif
!!-------------------------------      
      if(jobgw/=0) return
      if(.not.master_mpi) return
      
!!=====================================================================
!! only for jobgw=0      
!!  Count number of atoms : exclude floating orbitals
c      alat=lat_alat
c      plat =lat_plat
      nat = 0
      do  i = 1, nbas
        is = ssite(i)%spec
        lmaxa = int(sspec(is)%lmxa)
        if (lmaxa .gt. -1) then
          nat = nat + 1
        endif
        ipb(i) = nat
      enddo
      ldim=ham_ldham(1)
!! ndima
      ndima = 0
      lmxax = -1
      do  ib = 1, nbas
        is=ssite(ib)%spec
        lmaxa=sspec(is)%lmxa
        call dcopy(size(sspec(is)%pz),sspec(is)%pz,1,pnz,1)
        lmxax = max(lmxax,lmaxa)
        if (lmaxa .gt. -1) then
          do  l = 0, lmaxa
            npqn = 2
            if (pnz(l+1,1) .ne. 0) npqn = 3
            ndima = ndima + npqn*(2*l+1)
          enddo
        endif
      enddo
!! === Make files SYMOPS, LATTC, CLASS, NLAindx, ldima ===
        call info(30,1,1,' Creating files SYMOPS, LATTC, CLASS, NLAindx, ldima',0,0)
        allocate(lmxa(nbas),bas(3,nbas)) 
!! ldima
        open(newunit= ifildima,file='ldima') 
        do  ib = 1, nbas
          lmxa(ib) = sspec(int(ssite(ib)%spec))%lmxa  !int(sspec(ips(ib))%lmxa)
          if (lmxa(ib) .gt. -1) then
            call orbl ( ib , 0 , ldim , iv_a_oindxo , norb , ltab , ktab , xx , offl , i1 )
            write(ifildima,"(3i10)") i1
          endif
       enddo
       close(ifildima)
!! LATTC
        lmxax = mxint(nbas,lmxa)
        allocate(konft(0:lmxax,nbas,nsp))
        do ib = 1, nbas
           is = ssite(ib)%spec
           call dcopy(size(ssite(ib)%pnu), ssite(ib)%pnu,1,pnu,1)
           call dcopy(size(ssite(ib)%pz),  ssite(ib)%pz, 1,pnz,1)
           do  isp = 1, nsp
              do  l  = 0, lmxa(ib)
                 konft(l,ib,isp) = pnu(l+1,isp)
                 if( mod(pnz(l+1,isp),10d0) .lt. pnu(l+1,isp) .and. pnz(l+1,isp) .gt. 0) then
                    konft(l,ib,isp) = mod(pnz(l+1,isp),10d0)
                 endif    
              enddo
           enddo
        enddo
        call wlattc(alat,plat,nbas,nat,ipb,lmxax,lmxa,nsp,konft)
!! NLAindx
        ifinlaindx = fopnx('NLAindx',2,2,-1)
        rewind ifinlaindx
        write(ifinlaindx,'(''----NLAindx start---------------''/I6)') ndima
        ndima = 0
!!    This loop order is backwardly compatible with prior versions
        do  ipqn = 1, 3
          do  ib = 1, nbas
            is =  ssite(ib)%spec
            lmaxa=sspec(is)%lmxa
            i_copy_size=size(ssite(ib)%pnu)
            call dcopy  (i_copy_size ,ssite(ib)%pnu,1,pnu,1)
            i_copy_size=size(ssite(ib)%pnu)
            call dcopy  (i_copy_size,ssite(ib)%pz, 1,pnz,1)
            call icopy(sspec(is)%idxdn,sspec(is)%idxdn,1,idxdn,1)
            if (lmaxa .gt. -1) then
              do  l = 0, lmaxa
                npqn = 2
                if (pnz(l+1,1) .ne. 0) npqn = 3
                if (ipqn .le. npqn) then
                  konf = pnu(l+1,1)
                  if (ipqn .eq. 3) konf = mod(pnz(l+1,1),10d0)
                  strn4 = dig(konf)//lsym(l)//'_'//lorb(ipqn)
                  if (idxdn(l+1,1) .eq. 1 .or. idxdn(l+1,2) .eq. 1)
     .            call chcase(0,1,strn4(2:2))
                  write(ifinlaindx,'(i6,i3,i4,i6,4x,a)')
     .            ipqn, l, ipb(ib), ndima, strn4
                  ndima = ndima + (2*l+1)
                endif
              enddo
            endif
          enddo
        enddo
        call fclr(' ',ifinlaindx)
!! CLASS
        ificlass = fopnx('CLASS',2,2,-1)
        rewind ificlass
        allocate(ipc(nbas),ipcx(nbas))
        do ib=1,nbas
           ipc(ib) = ssite(ib)%class !call spackv_array_copy_i8_i('u',ssite(ib)%class,1,ib,ipc)
        enddo
        call pvsug1(nbas,lmxa,ipc,ipcx)
        do  i = 1, nbas
           if(lmxa(i) .gt. -1) then
              write(ificlass,'(2I4)') ipb(i), ipcx(i)
           endif
        enddo
      deallocate(ipc,ipcx)
      call fclr(' ',ificlass)
      call fexit(0,1,' OK! '//'lmfgw mode=0 generated LATTC CLASS NLAindx ldima',0)
!!================================================      
      end subroutine
      end module
