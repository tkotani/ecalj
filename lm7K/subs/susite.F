Cgetarg...info...           structure ['susite', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['susite', 'site', 'v_ssite']
Cgetarg...info...           structure ['susite', 'spec', 'sspec']
Cgetarg...info...           structure ['susite', 'ham', 'sham']
Cgetarg...info...           structure ['susite', 'lat', 'slat']
Cgetarg...info...           structure ['susite', 'array', 'sarray']
Cgetarg...info...           structure ['susite', 'pot', 'spot']
Cgetarg...info...           structure ['pvsub2', 'site', 'ssite']
      module m_susite_func
      contains
      subroutine susite(sctrl,sham,spot,sarray,slat,sspec,v_ssite)
C- Sets up permanent arrays related to basis
C takao have removed things not related to lmf.
      use m_struc_def
      implicit none
      type(s_site),pointer:: v_ssite(:)
      type(s_ctrl)::sctrl
      type(s_ham):: sham
      type(s_pot):: spot
      type(s_array)::sarray
      type(s_lat):: slat
      type(s_spec),pointer:: sspec(:)
      
C ... Local parameters
c      integer procid,master,mpipid
c      logical lgors,pass1,lorder,mlog
      integer :: nsite,nbas,nl,ib, i,i_spacks,i_copy_size,i_spackv
      integer :: oips,opos

cigets,,nbasp,ncl,npl,i,ib, !npadl,npadr,
c     .  ix(3),neul,nbf,opgplp,opgord,opos,oips,ovshft,opgfsl,
c     .  opgfvl,owk,oeula,omagf,oaamom,oclssl,ocllst,oclord,oclp,nclp,
c     .  igetss,iprint
c      parameter (nclp=9)
c      double precicsion plat(9),platl(9),platr(9),pgfn(3),xsi(5)
c      double precision dval,dglob,seref,dgetss
      double precision seref
      character*8 spid

c      real(8),allocatable:: woeula(:),womagf(:)
C ... Heap
      integer w(1)
      common /w/ w

Cgetarg...info...           integer :: nsite , integer(8) :: sctrl%nsite 1
Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nl , integer(8) :: sctrl%nl 1
Cgetarg       call upack('ctrl nsite nbas nl',sctrl,nsite,nbas,nl,0,0)
       nsite=sctrl%nsite
       nbas=sctrl%nbas
       nl=sctrl%nl

C     call rxx(nsite.ne.nbas,'susite not set up for nsite.ne.nbas')
Cki      ssize = nint(dval(v_ssite,1))
c      ssize=v_ssite(1)%size
c      pass1 = .false.
c      nbasp = nbas
c      npl   = nbas
c      master = 0
c      procid = mpipid(1)
c      mlog = .false.

       print *,'susite=',nsite,nbas,nl
C ... Pack species labels; make sumeref
      seref = 0
      do  ib = 1, nsite
Cgetarg         i = igetss('site spec',ib,v_ssite)
         print *,'vvv000 i=',v_ssite(ib)%spec
         i = int(v_ssite(ib)%spec) 
         print *,'vvv111 i=',i

C       Skip over sites with no species (e.g. multipoles)
        if (i .gt. 0) then
Cgetarg           call spacks(0,'spec name',sspec,spid,i,i)
           do i_spacks=i,i 
           call spacks_copy('u',sspec(i_spacks)%name,i,i,spid,i_spacks)
           enddo

Cgetarg           call spacks(1,'site clabel',v_ssite,spid,ib,ib)
           do i_spacks=ib,ib 
           call spacks_copy('p',v_ssite(i_spacks)%clabel,ib,ib,spid,i_spacks)
           enddo

Cgetarg           seref = seref + dgetss('spec eref',i,sspec)
           seref = seref + (sspec(i)%eref) 

        endif
      enddo
Cgetarg...info...           real(8) :: seref , real(8) :: sham%seref 1
Cgetarg       call pack1('ham seref',sham,seref)
      sham%seref=seref 

C ... Start again if ssite was reallocated
c   10 pass1 = .not. pass1
Cgetarg...info...           real(8) :: plat(9) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: platl(9) , real(8) :: slat%platl 9
Cgetarg...info...           real(8) :: platr(9) , real(8) :: slat%platr 9
Cgetarg       call upack('lat plat platl platr',slat,plat,platl,platr,0,0)
       
c       i_copy_size=size(slat%plat) 
c       call dcopy(i_copy_size,slat%plat,1,plat,1) 
c       i_copy_size=size(slat%platl) 
c       call dcopy(i_copy_size,slat%platl,1,platl,1) 
c       i_copy_size=size(slat%platr) 
c       call dcopy(i_copy_size,slat%platr,1,platr,1) 

      print *,'vvv21'
      call defrr(opos,3*nsite)
Cgetarg       call spackv(10,'site pos',v_ssite,1,nsite,w(opos))
      i_copy_size=size(v_ssite(1)%pos)
      do i_spackv=1,nsite 
      call spackv_array_copy_r8_r8('u',v_ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
      enddo

      print *,'vvv22'
      call defi(oips,nsite)
Cgetarg       call spackv(10,'site spec',v_ssite,1,nsite,w(oips))
      i_copy_size=1; 
      do i_spackv=1,nsite 
      call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,w(oips))
      enddo


C$$$C --- Order-N: group all sites in a crystal into clusters ---
C$$$Cgetarg       if (lgors('ctrl lordn,-1',sctrl)) then
C$$$       if ( iand(-1,int(sctrl%lordn)) .ne.0 ) then 

C$$$        call defi(oclssl, nsite)
C$$$        call defi(ocllst, nsite+1)
C$$$        call defi(oclp,  -(nsite+1)*nclp)
C$$$        call defi(oclord, nsite)
C$$$        call clsset(11,slat,nsite,w(opos),ncl,w(oclssl),w(ocllst),
C$$$     .    w(oclord),w(oclp))

C$$$C   ... Reorder the site structure according to clord
C$$$        call defdr(owk, ssize*nsite)
C$$$        call dvprm(ssize,nsite,v_ssite,w(owk),w(oclord),.true.)
C$$$        call rlse(owk)
C$$$        call dvset(v_ssite,1,1,dble(ssize))
C$$$C   ... Unpack arrays, now permuted; shorten pos
C$$$        ix(1) = 1
C$$$        ix(2) = 1
C$$$        ix(3) = 1
C$$$Cgetarg         call spackv(10,'site pos',v_ssite,1,nsite,w(opos))
C$$$         i_copy_size=size(v_ssite(1)%pos)
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_r8_r8('u',v_ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
C$$$         enddo

C$$$        call shorps(nbas,plat,ix,w(opos),w(opos))
C$$$Cgetarg         call spackv(11,'site pos',v_ssite,1,nsite,w(opos))
C$$$         i_copy_size=size(v_ssite(1)%pos)
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_r8_r8('p',v_ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
C$$$         enddo

C$$$Cgetarg         call spackv(10,'site spec',v_ssite,1,nsite,w(oips))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,w(oips))
C$$$         enddo

C$$$C       call shstru('site',v_ssite,1,nsite)
C$$$        call redfi(oclp,   nclp*(ncl+1))
C$$$Cgetarg...info...           integer :: ncl , integer(8) :: sarray%ncl 1
C$$$Cgetarg...info...           integer :: ocllst , integer(8) :: sarray%ocllst 1
C$$$Cgetarg...info...           integer :: oclp , integer(8) :: sarray%oclp 1
C$$$Cgetarg...info...           integer :: oclssl , integer(8) :: sarray%oclssl 1
C$$$Cgetarg         call pack5('array ncl ocllst oclp oclssl',sarray,ncl,ocllst,
C$$$Cgetarg      .    oclp,oclssl,0)
         
C$$$         sarray%ncl=ncl 
C$$$         sarray%ocllst=ocllst 
C$$$         sarray%oclp=oclp 
C$$$         sarray%oclssl=oclssl 

C$$$      endif

C --- Check and order principal layers, and sites by PL ---
Cgetarg       lpgf = igets('ctrl lpgf',sctrl)
c       lpgf = int(sctrl%lpgf(1)) 

C$$$      if (lpgf .ne. 0) then
C$$$        call defi(opgfsl,-nsite)
C$$$Cgetarg         call spackv(10,'site pl',v_ssite,1,nsite,w(opgfsl))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%pl,i_copy_size,i_spackv+1-1,w(opgfsl))
C$$$         enddo
C$$$        call defi(opgfvl,-nsite)
C$$$Cgetarg         call spackv(10,'site plv',v_ssite,1,nsite,w(opgfvl))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%plv,i_copy_size,i_spackv+1-1,w(opgfvl))
C$$$         enddo
C$$$        call defdr(ovshft, -8-2*nbasp)
C$$$        call defi(opgplp,  6*(npl+2))
C$$$C        print *, 'debuging susite'
C$$$C        call ivset(w(opgplp),1,6*(npl+2),-99)
C$$$C   ... no accumulated sites in -1 PL
C$$$        w(opgplp) = 0
C$$$        call defi(opgord,  2*nbas)
C$$$        if (pass1) call pshpr(min(iprint(),1))
C$$$C       if (pass1) call pshpr(0)
C$$$C       call prmx('pos before pgfset',w(opos),3,3,nbas)
C$$$        lorder = .not. pass1
C$$$        call pgfset(sspec,nbas,w(opos),plat,lorder,.true.,w(oips),
C$$$     .    w(ovshft),w(opgfsl),w(opgfvl),w(opgord),pgfn,npl,npadl,npadr,
C$$$     .    w(opgplp))
C$$$        nbasp = nbas + (npadl + npadr)
C$$$        ix(1) = dglob('nbasp',dble(nbasp),1)
C$$$        ix(2) = dglob('npl',dble(npl),1)
C$$$        nsite = nbas + 2*(npadl + npadr)
C$$$        if (pass1) then
C$$$          call poppr
C$$$Cki          call redfrr(v_ssite,nsite*ssize)
C$$$           allocate(v_ssite2(nsite))
C$$$           nsite_old=size(v_ssite)
C$$$           v_ssite2(1:nsite_old)=v_ssite(1:nsite_old)
C$$$           deallocate(v_ssite)
C$$$           allocate(v_ssite(nsite))
C$$$           do i=1,nsite; v_ssite(i)%size=usite_size(); enddo 
C$$$           v_ssite(1:nsite_old)=v_ssite2(1:nsite_old)
C$$$           deallocate(v_ssite2)
C$$$C         Zero out enlarged last part of it
C$$$          call defrr(owk,-nsite*ssize)
C$$$          call dcopy(nbas*ssize,v_ssite,1,w(owk),1)
C$$$          call dcopy(nsite*ssize,w(owk),1,v_ssite,1)
C$$$          call rlse(owk)
C$$$          goto 10
C$$$        endif
C$$$        call defdr(owk, ssize*nbas)
C$$$C   ... Reset pos,pfgsl, since they may be shifted by pgfset
C$$$Cgetarg         call spackv(11,'site pos',v_ssite,1,nsite,w(opos))
C$$$         i_copy_size=size(v_ssite(1)%pos)
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_r8_r8('p',v_ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
C$$$         enddo

C$$$Cgetarg         call spackv(11,'site pl',v_ssite,1,nsite,w(opgfsl))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('p',v_ssite(i_spackv)%pl,i_copy_size,i_spackv+1-1,w(opgfsl))
C$$$         enddo

C$$$C   ... Reorder the site structure according to pgord
C$$$        do  12  i = 1, nbas
C$$$   12   w(opgord+i-1) = w(opgord+i-1)-1
C$$$C       call awrit2('%n:1i',' ',100,6,nbas,w(opgord))
C$$$C       call prmx('ssite before perm',v_ssite,ssize,ssize,nbas)
C$$$        call dvperm(ssize,nbas,v_ssite,w(owk),w(opgord),.true.)
C$$$        call dvset(v_ssite,1,1,dble(ssize))
C$$$C       call prmx('ssite after perm',v_ssite,ssize,ssize,nbas)
C$$$C   ... Doubly pad ssite
C$$$        call pgfpp(nbas,npadl,npadr,ssize,1,.false.,.false.,
C$$$     .    v_ssite,v_ssite)
C$$$C   ... Unpack arrays, now permuted and padded
C$$$Cgetarg         call spackv(10,'site pl',v_ssite,1,nsite,w(opgfsl))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%pl,i_copy_size,i_spackv+1-1,w(opgfsl))
C$$$         enddo

C$$$Cgetarg         call spackv(10,'site plv',v_ssite,1,nsite,w(opgfvl))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%plv,i_copy_size,i_spackv+1-1,w(opgfvl))
C$$$         enddo

C$$$Cgetarg         call spackv(10,'site pos',v_ssite,1,nsite,w(opos))
C$$$         i_copy_size=size(v_ssite(1)%pos)
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_r8_r8('u',v_ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
C$$$         enddo

C$$$Cgetarg         call spackv(10,'site spec',v_ssite,1,nsite,w(oips))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,w(oips))
C$$$         enddo


C$$$C   ... Shift doubly padded bas; repack
C$$$        if (nbasp .gt. nbas) then
C$$$        call pgbasp(nbas,npadl,npadr,w(opos),plat,platl,platr,w(owk))
C$$$        call dscal(9,2d0,platl,1)
C$$$        call pgbasp(nbasp,npadl,npadr,w(owk),plat,platl,platr,w(opos))
C$$$        call dscal(9,.5d0,platl,1)
C$$$Cgetarg         call spackv(11,'site pos',v_ssite,1,nsite,w(opos))
C$$$         i_copy_size=size(v_ssite(1)%pos)
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_r8_r8('p',v_ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
C$$$         enddo

C$$$C       call prmx('pos, double pad',w(opos),3,3,nsite)
C$$$C       call prmx('ssite, double pad',v_ssite,ssize,ssize,nsite)
C$$$C   ... Shift pgfsl in padding layers, and repack
C$$$        call pvsub1(w(opgfsl),nbas,npadl,npadr)
C$$$Cgetarg         call spackv(11,'site pl',v_ssite,1,nsite,w(opgfsl))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('p',v_ssite(i_spackv)%pl,i_copy_size,i_spackv+1-1,w(opgfsl))
C$$$         enddo

C$$$C        call awrit2('%n:1i',' ',100,6,nsite,w(opgfsl))
C$$$C        call awrit2('%n:1i',' ',100,6,nbas,w(opgfvl))
C$$$C   ... Pick up padded ips,pgfsl,pgfvl
C$$$Cgetarg         call spackv(10,'site spec',v_ssite,1,nsite,w(oips))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,w(oips))
C$$$         enddo

C$$$Cgetarg         call spackv(10,'site plv',v_ssite,1,nsite,w(opgfvl))
C$$$         i_copy_size=1; 
C$$$         do i_spackv=1,nsite 
C$$$         call spackv_array_copy_i8_i('u',v_ssite(i_spackv)%plv,i_copy_size,i_spackv+1-1,w(opgfvl))
C$$$         enddo

C$$$        call redfi(opgplp,   6*(npl+2))
C$$$Cgetarg...info...           integer :: opgfsl , integer(8) :: sarray%opgfsl 1
C$$$Cgetarg...info...           integer :: opgfvl , integer(8) :: sarray%opgfvl 1
C$$$Cgetarg...info...           integer :: opgplp , integer(8) :: sarray%opgplp 1
C$$$Cgetarg         call pack5('array opgfsl opgfvl opgplp',sarray,opgfsl,opgfvl,
C$$$Cgetarg      .    opgplp,0,0)
         
C$$$         sarray%opgfsl=opgfsl 
C$$$         sarray%opgfvl=opgfvl 
C$$$         sarray%opgplp=opgplp 

C$$$Cgetarg...info...           integer :: npl , integer(8) :: sctrl%npl 1
C$$$Cgetarg         call pack1('ctrl npl',sctrl,npl)
         
C$$$         sctrl%npl=npl 

C$$$Cgetarg...info...           integer :: ovshft , integer(8) :: spot%ovshf 1
C$$$Cgetarg         call pack1('pot ovshf',spot,ovshft)
         
C$$$         spot%ovshf=ovshft 

C$$$C       call shstru('site',v_ssite,1,nsite)
C$$$        endif
C$$$        call rlse(opgord)
C$$$Cgetarg...info...           integer :: nsite , integer(8) :: sctrl%nsite 1
C$$$Cgetarg...info...           integer :: nbasp , integer(8) :: sctrl%nbasp 1
C$$$Cgetarg         call pack5('ctrl nsite nbasp',sctrl,nsite,nbasp,0,0,0)
         
C$$$         sctrl%nsite=nsite 
C$$$         sctrl%nbasp=nbasp 

C$$$      else
C   ... No partitioning of hamiltonian: make global pgplp
c        npadl = 0
c        npadr = 0
c        call defi(opgplp, -12)
c        w(opgplp+0) = -1
c        w(opgplp+6) =  nbas
c        call defdr(ovshft, -8-2*nbas)
Cgetarg...info...           integer :: opgplp , integer(8) :: sarray%opgplp 1
Cgetarg         call pack5('array opgfsl opgfvl opgplp',sarray,1,1,opgplp,0,0)
         
c         sarray%opgfsl=1 
c         sarray%opgfvl=1 
c         sarray%opgplp=opgplp 

Cgetarg         call pack1('ctrl npl',sctrl,1)
         
         sctrl%npl=1 

Cgetarg...info...           integer :: ovshft , integer(8) :: spot%ovshf 1
Cgetarg         call pack1('pot ovshf',spot,ovshft)
         
c         spot%ovshf=ovshft 

C$$$      endif

C --- Euler angles and external magnetic fields ---
Cgetarg       if (lgors('ctrl lncol,11',sctrl)) then
C$$$      if ( iand(11,int(sctrl%lncol)) .ne.0 ) then 

C$$$        call defdr(oeula, -nbasp*nl**2*3+3)
C$$$        call dpzero(xsi,5)
C$$$        if (procid .eq. master) then
C$$$          allocate(woeula(nbasp*nl**2*3+3))
C$$$          call dcopy(nbasp*nl**2*3+3, w(oeula),1,woeula,1)
C$$$CKi          call pvsub2(0,v_ssite,nbas,nbasp,nl,w(oeula),neul,xsi)
C$$$          call pvsub2(0,v_ssite,nbas,nbasp,nl,woeula,neul,xsi)
C$$$          call dcopy(nbasp*nl**2*3+3, woeula,1,w(oeula),1)
C$$$          deallocate(woeula)
C$$$        endif
C$$$        call mpibc1(neul,1,2,mlog,'susite','neul')
C$$$        call mpibc1(xsi,5,4,mlog,'susite','xsi')
C$$$        call redfrr(oeula, nbasp*neul*3+3)
C$$$        call mpibc1(w(oeula),nbas*neul*3,4,mlog,'susite','eula')
C$$$Cgetarg...info...           real(8) :: xsi(5) , real(8) :: sctrl%sdxsi 4
C$$$Cgetarg         call pack1('ctrl sdxsi',sctrl,xsi)
         
C$$$         i_copy_size=size(sctrl%sdxsi) 
C$$$         call dcopy(i_copy_size,xsi,1,sctrl%sdxsi,1) 

C$$$Cgetarg...info...           integer :: neul , integer(8) :: sham%neula 1
C$$$Cgetarg...info...           integer :: oeula , integer(8) :: sham%oeula 1
C$$$Cgetarg         call pack2('ham neula oeula',sham,neul,oeula)
         
C$$$         sham%neula=neul 
C$$$         sham%oeula=oeula 


C$$$        nbf = 1
C$$$        omagf = 1
C$$$Cgetarg         if (lgors('ctrl lncol,8',sctrl)) then
C$$$       if ( iand(8,int(sctrl%lncol)) .ne.0 ) then 

C$$$          call defdr(omagf, -nbasp*nl**2*3)
C$$$          if (procid .eq. master) then
C$$$            allocate(womagf(nbasp*nl**2*3))
C$$$            call dcopy(nbasp*nl**2*3, w(omagf),1,womagf,1)
C$$$CKi            call pvsub2(1,v_ssite,nbas,nbasp,nl,w(omagf),nbf,xsi)
C$$$            call pvsub2(1,v_ssite,nbas,nbasp,nl,womagf,nbf,xsi)
C$$$            call dcopy(nbasp*nl**2*3, womagf,1,w(omagf),1)
C$$$            deallocate(womagf)
C$$$          endif
C$$$          call mpibc1(nbf,1,2,mlog,'susite','nbf')
C$$$          call mpibc1(w(omagf),nbas*neul*3,4,mlog,'susite','magf')
C$$$          call redfrr(omagf, nbasp*nbf*3)
C$$$C         call prmx('b field',w(omagf),nbasp*nbf,nbasp*nbf,3)
C$$$        endif
C$$$Cgetarg...info...           integer :: nbf , integer(8) :: sham%nbf 1
C$$$Cgetarg...info...           integer :: omagf , integer(8) :: sham%omagf 1
C$$$Cgetarg         call pack2('ham nbf omagf',sham,nbf,omagf)
         
C$$$         sham%nbf=nbf 
C$$$         sham%omagf=omagf 


C$$$        call defdr(oaamom, -nbasp)
C$$$Cgetarg...info...           integer :: oaamom , integer(8) :: spot%oaamom 1
C$$$Cgetarg         call pack1('pot oaamom',spot,oaamom)
         
C$$$         spot%oaamom=oaamom 

C$$$      else
c        call defdr(oeula,1)
Cgetarg...info...           integer :: oeula , integer(8) :: sham%oeula 1
Cgetarg         call pack2('ham neula oeula',sham,0,oeula)
         
c       sham%neula=0 
c       sham%oeula=oeula 

c      endif

Cgetarg...info...           integer :: opos , integer(8) :: slat%opos 1
Cgetarg       call pack5('lat opos',slat,opos,0,0,0,0)
       
       slat%opos=opos 

Cgetarg...info...           integer :: opos , integer(8) :: sarray%opos 1
Cgetarg...info...           integer :: oips , integer(8) :: sarray%oips 1
Cgetarg...info...           integer :: npadl , integer(8) :: sarray%npadl 1
Cgetarg...info...           integer :: npadr , integer(8) :: sarray%npadr 1
Cgetarg       call pack5('array opos oips npadl npadr',sarray,opos,oips,npadl,
Cgetarg      .  npadr,0)
       
       sarray%opos=opos 
       sarray%oips=oips 
c       sarray%npadl=npadl 
c       sarray%npadr=npadr 

C     call pack5('ctrl nsite nbasp',sctrl,nsite,nbasp,0,0,0)
C     call shstru('site',v_ssite,1,nsite)

      end subroutine susite
C$$$      subroutine pvsub1(pgfsl,nbas,npadl,npadr)

C$$$C- Sets PL indices for padding layers
C$$$C     implicit none
C$$$      integer nbas,npadl,npadr,pgfsl(nbas)
C$$$      integer nbasp,i

C$$$      nbasp = nbas + (npadl + npadr)
C$$$      do  10  i = 1, npadl
C$$$      pgfsl(nbas+i) = pgfsl(nbas+i) - 1
C$$$   10 pgfsl(nbasp+i) = pgfsl(nbasp+i) - 2
C$$$      do  20  i = 1, npadr
C$$$      pgfsl(nbas+npadl+i) = pgfsl(nbas+npadl+i) + 1
C$$$   20 pgfsl(nbasp+npadl+i) = pgfsl(nbasp+npadl+i) + 2

C$$$      end subroutine pvsub1
C$$$      subroutine pvsub2(mode,ssite,nbas,nbasp,nl,eula,neul,xsi)

C$$$C- Unpacks Euler angles or magnetic fields; reads from disk
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   mode  :0 for euler angles, 1 for bfield
C$$$Ci   ssite :struct for site-specific information; see routine usite
C$$$Ci     Elts read:
C$$$Ci     Stored:    eula
C$$$Ci     Passed to: spackv
C$$$Ci   nbas  :size of basis
C$$$Ci   nbasp :size of padded basis (layer programs)
C$$$Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
C$$$Ci   nl    :(global maximum l) + 1
C$$$Co Outputs
C$$$Co   eula  :(mode 0) Euler angles for noncollinear spins
C$$$Co         :(mode 1) Magnetic field read
C$$$Co   neul  :1 if Euler angles (Bfield) are l-independent, nl otherwise
C$$$Co   xsi   :global deamon parameters for spin dynamics (not read mode 1)
C$$$Cl Local variables
C$$$Cl         :
C$$$Cr Remarks
C$$$Cr
C$$$Cu Updates
C$$$Cu   29 Jan 03 Added mode to also read b-field
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$C ... Passed parameters
C$$$      use m_struc_def
C$$$      integer mode,nbas,nbasp,nl,neul
C$$$CKi      double precision ssite(1),eula(nbasp,*)
C$$$      double precision eula(nbasp,*)
C$$$      type(s_site)::ssite(:)
C$$$C ... Local parameters
C$$$      integer oeula,ifi,fopn,fxst
C$$$      double precision xsi(5)
C$$$      character name*6
C$$$C ... Heap
C$$$      integer w(1)
C$$$      common /w/ w

C$$$      call isanrg(mode,0,1,'pvsub2','mode',.true.)
C$$$      if (mode .eq. 0) name = 'eula  '
C$$$      if (mode .eq. 1) name = 'bfield'

C$$$      neul = 1
C$$$      call defr(oeula,3*nbasp)
C$$$Cki      call spackv(10,'site '//name,ssite,1,nbasp,w(oeula))
C$$$      if (mode.eq.0) then
C$$$Cgetarg       call spackv(10,'site eula',ssite,1,nbasp,w(oeula))
C$$$       i_copy_size=size(ssite(1)%eula)
C$$$       do i_spackv=1,nbasp 
C$$$       call spackv_array_copy_r8_r8('u',ssite(i_spackv)%eula,i_copy_size,i_spackv+1-1,w(oeula))
C$$$       enddo

C$$$      else if (mode.eq.1) then
C$$$Cgetarg       call spackv(10,'site bfield',ssite,1,nbasp,w(oeula))
C$$$       i_copy_size=size(ssite(1)%bfield)
C$$$       do i_spackv=1,nbasp 
C$$$       call spackv_array_copy_r8_r8('u',ssite(i_spackv)%bfield,i_copy_size,i_spackv+1-1,w(oeula))
C$$$       enddo

C$$$      endif
C$$$C     call prmx(name,w(oeula),3,3,nbasp)
C$$$      call dmcpy(w(oeula),1,3,eula,nbasp,1,nbasp,3)
C$$$C     call prmx(name//' transposed',eula,nbasp,nbasp,3)

C$$$C ... Read angles from disk, if available
C$$$      if (fxst(name) .eq. 1) then
C$$$        ifi = fopn(name)
C$$$        rewind ifi
C$$$        if (mode .eq. 0) call ioeula(nbasp,nl,eula,neul,xsi,ifi)
C$$$        if (mode .eq. 1) call iomagf(nbasp,nl,eula,neul,ifi)
C$$$        call fclose(ifi)
C$$$      endif
C$$$C     call prmx(name//' yet again',eula,nbasp*neul,nbasp*neul,3)

C$$$C     At present, no noncollinearity in padding layers
C$$$      if (nbasp .gt. nbas)
C$$$     .  call dmcpy(0d0,0,0,eula(1+nbas,1),nbasp,1,nbasp-nbas,neul*3)

C$$$C     Repack if neula is 1
C$$$      if (neul .eq. 1) then
C$$$        call dmcpy(eula,nbasp,1,w(oeula),1,3,nbasp,3)
C$$$CKi        call spackv(11,'site '//name,ssite,1,nbasp,w(oeula))
C$$$        if (mode.eq.0) then
C$$$Cgetarg         call spackv(11,'site eula',ssite,1,nbasp,w(oeula))
C$$$         i_copy_size=size(ssite(1)%eula)
C$$$         do i_spackv=1,nbasp 
C$$$         call spackv_array_copy_r8_r8('p',ssite(i_spackv)%eula,i_copy_size,i_spackv+1-1,w(oeula))
C$$$         enddo

C$$$        else
C$$$Cgetarg         call spackv(11,'site bfield',ssite,1,nbasp,w(oeula))
C$$$         i_copy_size=size(ssite(1)%bfield)
C$$$         do i_spackv=1,nbasp 
C$$$         call spackv_array_copy_r8_r8('p',ssite(i_spackv)%bfield,i_copy_size,i_spackv+1-1,w(oeula))
C$$$         enddo

C$$$        endif
C$$$      endif

C$$$      call rlse(oeula)

C$$$      end subroutine pvsub2

      end module m_susite_func

