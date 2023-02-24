module m_sxcf_count !job scheduler for self-energy calculation. icount mechanism
  use m_readeigen,only: Readeval
  use m_itq,only: ntq,nbandmx
  use m_genallcf_v3,only: nlmto,nspin,nctot,niw,ecore
  use m_read_bzdata,only: qibz,qbz,wk=>wbz,nqibz,nqbz,wklm,lxklm,nq0i, wqt=>wt,q0i, irk
  use m_readfreq_r,only: freq_r, nw_i,nw,freqx,wx=>wwx,nblochpmx,mrecl,expa_,npm,nprecx
  use m_readhbe,only: nband,mrecg
  use m_hamindex,only: ngrp
  use m_mpi,only: MPI__sxcf_rankdivider
  use m_ftox
  implicit none
  public sxcf_scz_count
  !=== Job scheduler ==============
  integer,public:: ncount 
  integer,allocatable,public:: ispc(:),kxc(:),irotc(:),ipc(:),krc(:),nstateMax(:),nstti(:),nstte(:)
  integer,allocatable,public:: nwxic(:), nwxc(:), nt_maxc(:)
  !=========================================================
contains
  subroutine sxcf_scz_count(ef,esmr,exchange,ixc,nspinmx) 
    intent(in)              ef,esmr,exchange,ixc,nspinmx
    logical :: exchange
    integer :: isp,nspinmx,jobsw 
!    integer :: nbandmx(nqibz,nspinmx)
    real(8) :: ef,esmr
    real(8):: ebmx
    complex(8),pointer::zsec(:,:)
    complex(8),pointer::ww(:,:)
    integer,allocatable :: ifrcw(:),ifrcwi(:)
    integer :: ip, it, itp, i, ix, kx, irot, kr
    integer :: nt0p, nt0m,nstate , nbmax, ntqxx 
    integer :: nt,ixs,iw,ivc,ifvcoud,ngb0
    integer :: ifwd,nrot,nwp,ierr 
    integer :: iqini,iqend
    integer :: invr,ia,nn,ntp0,no,itpp,nrec,itini,itend,nbmxe
    integer :: iwp,nwxi,nwx,iir, igb1,igb2,ix0,iii
    integer :: invrot,nocc,nlmtobnd,nt0,verbose,ififr, istate,  nt_max ,noccx
    real(8) :: ekc(nctot+nband),ekq(nband), det, q(3) !,ua_
    real(8) :: wtt,wfac,we!,esmrx
    real(8) :: qvv(3),eq(nband),omega(ntq),quu(3),freqw,ratio
    real(8) :: qibz_k(3),qbz_kr(3),vc,omega0,omg
    complex(8),allocatable,target:: zwz(:,:,:),zw(:,:)
    real(8), parameter :: wfaccut=1d-8,tolq=1d-8
    complex(8), parameter :: img=(0d0,1d0)
    character(10) :: i2char
    real(8)::polinta, wfacx, wfacx2, weavx2, wexx,ua2_(niw),freqw1,q_r(3),qk(3)
    logical,parameter :: debug=.false.,timemix=.true.
    logical ::   oncew, onceww !, eibz4sig  
    real(8),allocatable:: we_(:,:),wfac_(:,:)
    complex(8),allocatable:: w3p(:),wtff(:)
    logical:: tote=.false.!, hermitianW
    real(8),allocatable:: vcoud_(:),wfft(:)
    logical:: iprx,cmdopt0
    integer:: ixx,ixc,icount,ndivmx
    real(8),parameter:: pi=4d0*datan(1d0), fpi=4d0*pi, tpi=8d0*datan(1d0),ddw=10d0
    integer:: kxold,nccc,icount0
    complex(8),allocatable:: zmelc(:,:,:)
    integer,allocatable::ndiv(:),nstatei(:,:),nstatee(:,:),irkip(:,:,:,:)
    integer:: job
    !!----------------------------------------------------------------
    if(npm==2) call rx('sxcf_fal2_sc: npm=2 need to be examined')
    if(ixc==3.and.nctot==0) return
    rankdivider: block !  We divide irkip_all into irkip for nodes. irkip is dependent on rank.
      !  Total number of none zero irkip for all ranks is the number of nonzero irkip_all
      integer:: irkip_all(nspinmx,nqibz,ngrp,nqibz),iqq,is
      do is = 1,nspinmx
         do iqq=1,nqibz
            irkip_all(is,:,:,iqq)=irk
         enddo
      enddo
      allocate(    irkip(nspinmx,nqibz,ngrp,nqibz)) ! nrkip is weight correspoinding to irkip for a node.
      call MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nqibz,  irkip)
    endblock rankdivider
    PreIcountBlock: Block!Get nstateMax(ncount),ndiv(icount),nstatei(j,icount),nstatee(j,icount)
      integer:: ndivide,nstateavl,nnn,nloadav,nrem,idiv,j
      integer,allocatable:: nload(:)
      ncount=count(irkip/=0)
      allocate(nstateMax(ncount))
      iqini = 1
      iqend = nqibz             
      icount=0
      kxloop:do kx = iqini,iqend !quick loop
         isploop: do isp = 1,nspinmx !empty run to get index for icount ordering
            if(sum(irkip(isp,kx,:,:))==0) cycle ! next kx
            irotloop: do irot = 1,ngrp !over rotations irot ===
               if(sum(irkip(isp,kx,irot,:))==0) cycle ! next ip
               iqloop: do 1150 ip = 1,nqibz         
                  kr = irkip(isp,kx,irot,ip) ! index for rotated kr in the FBZ
                  if(kr==0) cycle
                  icount=icount+1
                  q(1:3)= qibz(1:3,ip)
                  qbz_kr= qbz (:,kr)     !rotated qbz vector. 
                  qk =  q - qbz_kr        
                  ekq = readeval(qk, isp) 
                  ekc(nctot+1:nctot+nband) = ekq (1:nband)
                  nt0p = count(ekq<ef+ddw*esmr) +nctot 
                  if(exchange) then
                     nstateMax(icount) = nt0p
                  else   
                     ebmx=1d10 !this is needed probably for filling 1d99 for ekc(i) above boundary.
                     nbmxe = count(ekc<ebmx)-nctot !nocc (ekc,ebmx,nstatetot)-nctot!
                     nbmax  = min(nband,nbmxe) 
                     nstateMax(icount) = nctot + nbmax ! = nstate for the case of correlation
                  endif
1150           enddo iqloop
            enddo irotloop
         enddo isploop
      enddo kxloop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      ! ndivide = 2 
      nstateavl = 16  ! middle stats are batched by nstateavl.
      !nstateavl=max(sum(nstatemax)/(ncount*ndivide),1)
      if(ixc==3) nstateavl= maxval(nstatemax)
      ! size of average load of middle states (in G)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(ndiv(ncount))
      ndiv = (nstatemax-1)/nstateavl + 1  !number of division for middle states.
      ndivmx = maxval(ndiv)
      allocate(nstatei(ndivmx,ncount),nstatee(ndivmx,ncount),nload(ndivmx))
      do icount=1,ncount
         nnn = nstatemax(icount)    !total number of middle states for given icount
         ndiv(icount) = (nnn-1)/nstateavl + 1  !number of division for icount
         nloadav = nnn/ndiv(icount) !number of average load of middle states
         nrem = nnn - ndiv(icount)*nloadav !remnant count
         nload(1:nrem) = nloadav+1
         nload(nrem+1:ndiv(icount)) = nloadav !write(6,ftox)'nload=',nload(1:ndiv(icount))
         nstatei(:,icount)= [(sum(nload(1:idiv-1))+1,idiv=1,ndiv(icount))]!init index for(idiv,icount)
         nstatee(:,icount)= [(sum(nload(1:idiv)),    idiv=1,ndiv(icount))]!end  index
      enddo
      !      do icount=1,ncount
      !      do j=1,ndiv(icount)
      !        write(6,ftox)'nnnx icou ndiv=',icount,j,nstatei(j,icount),nstatee(j,icount),nstatemax(icount)
      !      enddo
      !      enddo
      deallocate(nload,nstatemax)
    EndBlock PreIcountBlock
    write(6,*)'nnn init ncount=',ncount
    ncount = ncount*ndivmx
    ! icount mechanism for sum in MAINicountloop 3030
    IcountBlock: Block !quick loop to gather index sets for main loop
      integer:: idiv
      !      ncount=count(irkip/=0)
      allocate(ispc(ncount),kxc(ncount),irotc(ncount),ipc(ncount),krc(ncount))
      allocate(nwxic(ncount), nwxc(ncount), nt_maxc(ncount),nstateMax(ncount))
      allocate(nstti(ncount),nstte(ncount))
      iqini = 1
      iqend = nqibz             !no sum for offset-Gamma points.
      icount=0
      icount0=0
      do 130 kx = iqini,iqend !this is empty run to get index for icount ordering
         do 120 isp = 1,nspinmx 
            if(sum(irkip(isp,kx,:,:))==0) cycle ! next kx
            do 140 irot = 1,ngrp !over rotations irot ===
               if(sum(irkip(isp,kx,irot,:))==0) cycle ! next ip
               do 150 ip = 1,nqibz         
                  kr = irkip(isp,kx,irot,ip) ! index for rotated kr in the FBZ
                  if(kr==0) cycle
                  icount0=icount0+1
                  do idiv=1,ndiv(icount0) !icount loop have further division of middle states by ndiv
                     icount=icount+1
                     nstti(icount)=nstatei(idiv,icount0) ! [nstti,nstte] specify range of middle states.
                     nstte(icount)=nstatee(idiv,icount0) !
                     ispc(icount)=isp!icount specify isp,kx,irot,iq. (kx,irot) gives kr in the all FZ.
                     kxc(icount)=kx
                     irotc(icount)=irot
                     ipc(icount)=ip
                     krc(icount)=kr
                     qibz_k = qibz(:,kx)
                     q(1:3)= qibz(1:3,ip)
                     eq = readeval(q,isp)
                     omega(1:ntq) = eq(1:ntq)
                     qbz_kr= qbz (:,kr)     !rotated qbz vector. 
                     qk =  q - qbz_kr        
                     ekq = readeval(qk, isp) 
                     ekc(nctot+1:nctot+nband) = ekq (1:nband)
                     nt0  = count(ekc<ef) 
                     nt0p = count(ekq<ef+ddw*esmr) +nctot 
                     nt0m = count(ekq<ef-ddw*esmr) +nctot
                     ntqxx = nbandmx(ip,isp) ! ntqxx is number of bands for <i|sigma|j>.
                     !write(6,*) icount, ispc(icount),kxc(icount),' irot ',irot,ip,kr
                     if(exchange) then
                        nstateMax(icount) = nt0p
                     else   
                        ebmx=1d10 !this is because filling 1d99 for ekc(i) above boundary.
                        nbmxe = count(ekc<ebmx)-nctot !nocc (ekc,ebmx,nstatetot)-nctot!
                        nbmax  = min(nband,nbmxe) 
                        nstateMax(icount) = nctot + nbmax ! = nstate for the case of correlation
                        call get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstateMax(icount),freq_r,&
                             nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,nwxi,nwx,nt_max)
                        !Get index nwxi nwx nt_max.
                        ! get_nwx is not written clearly, but works and not time-consuming.
                        nwxic(icount)=nwxi 
                        nwxc(icount)=nwx
                        nt_maxc(icount)=nt_max
                     endif
                  enddo
150            enddo
140         enddo
120      enddo
130   enddo
      if(icount0/=count(irkip/=0)) call rx('sxcf: icount/=count(irkip/=0)')
      ncount=icount
    EndBlock IcountBlock
    write(6,*)'sxcf_scz_count: nnn dev  ncount=',ncount
  end subroutine sxcf_scz_count
  subroutine get_nwx(omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,&
       nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug,&
       nwxi,nwx,nt_max)
    implicit none
    intent(in)::     omega,ntq,ntqxx,nt0p,nt0m,nstate,freq_r,&
       nw_i,nw,esmr,ef,ekc,wfaccut,nctot,nband,debug
    intent(out)::     &
       nwxi,nwx,nt_max
    !> Determine indexes of a range for calculation. !! It is better to clean this up...
    integer:: nctot,nw_i,nw,nstate,nt0p,nt0m,ntq, nband,ntqxx
    real(8):: esmr,ef,ekc(nctot+nband),wfaccut,freq_r(nw_i:nw)
    real(8):: wfac,wfacx2,we,weavx2,esmrx,wexx
    real(8),pointer::omg
    real(8),target:: omega(ntq)
    integer:: nt_max,nwxi,nwx,itp,it,itini,itend,iwp,ixs=-9999,ixsmin,ixsmx,verbose
    logical::debug
    !!     maximum ixs reqired.
    ixsmx =0
    ixsmin=0
    do 301 itp = 1,ntqxx
       omg => omega(itp) 
       if (omg < ef) then
          itini= 1
          itend= nt0p
       else
          itini= nt0m+1
          itend= nstate
       endif
       do 311 it=itini,itend
          esmrx = esmr
          if(it<=nctot) esmrx = 0d0
          wfac = wfacx2(omg,ef, ekc(it),esmrx)
          if(wfac<wfaccut) cycle !Gaussian case
          we = .5d0*(weavx2(omg,ef,ekc(it),esmr)-omg)
          if(it<=nctot) then
             if(wfac>wfaccut) call rx( "sxcf: it<=nctot.and.wfac/=0")
          endif
          do iwp = 1,nw
             ixs=iwp
             if(freq_r(iwp)>abs(we)) exit
          enddo
          if(ixs>ixsmx  .and. omg>=ef ) ixsmx  = ixs
          if(ixs>ixsmin .and. omg< ef ) ixsmin = ixs
          wexx  = we
          if(ixs+1 > nw) then
             write (*,*)'nw_i ixsmin wexx',nw_i,ixsmin,wexx,' omg ekc(it) ef ', omg,ekc(it),ef
             call rx( ' sxcf 222: |w-e| out of range')
          endif
311    enddo !continue
301 enddo !continue               
    if(nw_i==0) then          !time reversal
       nwxi = 0
       nwx  = max(ixsmx+1,ixsmin+1)
    else                      !no time revarsal working?
       nwxi = -ixsmin-1
       nwx  =  ixsmx+1
    endif
    if(nwx > nw .or. nwxi < nw_i ) call rx( ' get_nwx : |w-e| > max(w)')
    nt_max=nt0p               !initial nt_max
    do 401 itp = 1,ntqxx
       omg => omega(itp)
       if (omg > ef) then
          do  it = nt0m+1,nstate ! nt0m corresponds to efm
             wfac = wfacx2 (ef,omg, ekc(it),esmr)
             if(wfac>wfaccut) then
                if (it > nt_max) nt_max=it ! nt_max is  unocc. state
             endif               ! that ekc(it>nt_max)-omega > 0
          enddo                 ! so it > nt_max does not contribute to omega pole integral
       endif
401 enddo !continue               
  end subroutine get_nwx
end module m_sxcf_count

