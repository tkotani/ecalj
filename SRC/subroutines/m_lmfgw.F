      module m_lmfgw
      public:: m_lmfgw_init
      private
      contains
      subroutine m_lmfgw_init() !writing files for GW calculaiton
      use m_lmfinit,only:nbas,stdo,qbg=>zbak,ham_frzwf,lmaxu,nsp,nlibu,n0,nab,nppn,sspec=>v_sspec,ssite=>v_ssite,
     & nkap0,alat=>lat_alat,ham_ldham,iprmb
      use m_supot,only: k1,k2,k3
      use m_MPItk,only: master_mpi,mlog
      use m_lmfinit,only: nchan=>pot_nlma, nvl=>pot_nlml
      use m_lattic,only: qlat=>lat_qlat, plat=>lat_plat
      use m_suham,only: ndham=>ham_ndham
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
      integer:: fopnx
      data lsym /'s','p','d','f','g','5','6','7','8','9'/
      data lorb /'p','d','l'/
      data dig /'1','2','3','4','5','6','7','8','9'/
c     integer:: idxdn(n0,nkap0)
c      qbg = ctrl_zbak(1) !homogenious background charge
!!=====================================================================
c!! jobgw=1      
c      if(jobgw/=0) then
c        call m_mkpot_novxc() !Get osigx,otaux oppix spotx (without XC(LDA))
c        return
c      endif
!!=====================================================================
!! jobgw=0      
!!  Count number of atoms : exclude floating orbitals
      if(.not.master_mpi) return
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
            call orbl ( ib , 0 , ldim , iprmb , norb , ltab , ktab , xx , offl , i1 )
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
        open(newunit=ifinlaindx,file='NLAindx')
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
c            call icopy(sspec(is)%idxdn,sspec(is)%idxdn,1,idxdn,1)
            if (lmaxa .gt. -1) then
              do  l = 0, lmaxa
                npqn = 2
                if (pnz(l+1,1) .ne. 0) npqn = 3
                if (ipqn .le. npqn) then
                  konf = pnu(l+1,1)
                  if (ipqn .eq. 3) konf = mod(pnz(l+1,1),10d0)
                  strn4 = dig(konf)//lsym(l)//'_'//lorb(ipqn)
c                  if (idxdn(l+1,1) .eq. 1 .or. idxdn(l+1,2) .eq. 1)
c     .            call chcase(0,1,strn4(2:2))
                  write(ifinlaindx,'(i6,i3,i4,i6,4x,a)')
     .            ipqn, l, ipb(ib), ndima, strn4
                  ndima = ndima + (2*l+1)
                endif
              enddo
            endif
          enddo
        enddo
        close(ifinlaindx)
!! CLASS
c        ificlass = fopnx('CLASS',2,2,-1)
        open(newunit=ificlass,file="CLASS")
c        rewind ificlass
        allocate(ipc(nbas),ipcx(nbas))
        do ib=1,nbas
           ipc(ib) = ssite(ib)%class 
        enddo
        call pvsug1(nbas,lmxa,ipc,ipcx)
        do  i = 1, nbas
           if(lmxa(i) .gt. -1) then
              write(ificlass,'(2I4)') ipb(i), ipcx(i)
           endif
        enddo
        deallocate(ipc,ipcx)
        close(ificlass)
c      call fclr(' ',ificlass)
c      call fexit(0,1,' OK! '//'lmfgw mode=0 generated LATTC CLASS NLAindx ldima',0)
!!================================================      
      end subroutine
      
C      subroutine fmain
CC- debugs pvsug1
C      implicit none
C
C      integer nbas
C      parameter (nbas=8)
C      integer ipc(nbas),lmxa(nbas),ipcx(nbas)
C      data ipc/3,1,4,3,2,6,5,1/
C      data lmxa/-1,1,4,-1,2,6,-1,1/
C
C      call pvsug1(nbas,lmxa,ipc,ipcx)
C
C      end
      subroutine pvsug1(nbas,lmxa,ipc,ipcx)

C- Remakes class table, expunging classes with floating orbitals.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   lmxa  :augmentation l-cutoff
Ci   ipc   :class table: site ib belongs to class ipc(ib)
Co Inputs/Outputs
Co   ipcx  :expunged class table: classes with lmxa=-1 are expunged
Co         :and the remaining classes are sequentially renumbered
Co         :preserving the order of the remaining classes.
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Mar 07
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,lmxa(nbas),ipc(nbas),ipcx(nbas)
C ... Local parameters
      integer i,ip,ic,ipskip,nelim
      integer prm(nbas)

      call ivshel(1,nbas,ipc,prm,.true.)
C      do  i = 1, nbas
C        prm(i) = prm(i)+1
C      enddo

C     nelim = number of classes that have been eliminated so far
      nelim = 0
C     ipskip is the number of the last class that was skipped.
C     Multiple occurrences of a skipped class must still only
C     reduce the net number of classes by one.
C     We must avoid incrementing nelim when multiple sites
C     correspond to a skipped class.
      ipskip = 0
C     Loop over sites in order of increasing class index, ip
      do  i = 1, nbas
        ip = prm(i)+1
        ic = ipc(ip)
C       Test whether this class should be purged
        if (lmxa(ip) .lt. 0) then
          if (ipskip .ne. ic) nelim = nelim+1
          ipskip = ic
          ipcx(ip) = -1
        else
          ipcx(ip) = ic - nelim
        endif
      enddo
      end subroutine pvsug1

      subroutine wlattc(alat,plat,nbas,nat,ipb,lmxax,lmxa,nsp,konf)
C- Write LATTC file
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nbas  :size of basis, including floating orbitals sites
Ci   nat   :size of true basis (exclude floating orbitals sites)
Ci   ipb   :index to true basis (excluding floating orbitals)
Ci         :given site index including those orbitals
Ci   lmxax :global maximum of lmxa
Ci   lmxa  :augmentation l-cutoff
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   konf  :principal quantum numbers
Co Outputs
Co   file LATTC is written to disk
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 Jul 05 handle sites with lmxa=-1 -> no augmentation
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nat,nsp,lmxax,lmxa(nbas),konf(0:lmxax,nbas,nsp)
      integer ipb(nbas)
      real(8) :: plat(3,3),alat
      integer ifi,ib,isp,fopnx
      open(newunit=ifi,file='LATTC')
      write(ifi,"(e24.16)") alat
      write(ifi,"(3e24.16)") plat(1:3,1)
      write(ifi,"(3e24.16)") plat(1:3,2)
      write(ifi,"(3e24.16)") plat(1:3,3)
      write(ifi,*) ' -1d10 ! This is dummy. True QpGcut_psi is in GWIN0'
      write(ifi,*) ' ------------------------------------------- '
      write(ifi,"(2i4,' ! nbas lmxax (max l for argumentaion)')")
     .nat,lmxax
      write(ifi,*) ' ------------------------------------------- '
      do  isp = 1, nsp
        write(ifi,"(' -- ibas lmxa konf(s) konf(p) konf(d)... '
     .  , ' isp=',2i2)")isp
        do  ib = 1, nbas
          if (lmxa(ib) .gt. -1) then
            write(ifi,"('   ',99i4)")
     .      ipb(ib),lmxa(ib),konf(0:lmxa(ib),ib,isp)
          endif
        enddo
      enddo
      close(ifi)
      end subroutine wlattc
      
      end module m_lmfgw
