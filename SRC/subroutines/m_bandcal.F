      module m_bandcal
      use m_struc_def
      use m_suham,only: ndhamx=>ham_ndhamx,nspx=>ham_nspx 
      use m_qplist,only: nkp
      use m_lmfinit, only: nsp,nlibu,lmaxu,nbas,nl
      use m_subzi, only: numq
      use m_mkqp,only: ntet=> bz_ntet
      use m_qplist,only: qplist
      use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_set
      use m_lmfinit,only: lrsig=>ham_lsig, lso,ham_scaledsigma,lekkl,
     &     lmet=>bz_lmet,stdo,nbas,epsovl=>ham_oveps,nspc,plbnd,lfrce=>ctrl_lfrce,
     &      pwmode=>ham_pwmode,pwemax=>ham_pwemax,stdl
      use m_MPItk,only: mlog, master_mpi, procid,strprocid, numprocs=>nsize, mlog_MPIiq
      use m_subzi, only: m_subzi_init, nevmx,lswtk,rv_a_owtkb,rv_a_oswtk
      use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate,k1,k2,k3,
     o    zv_a_osmpot, qmom_rv, vconst, vrmt, sv_p_osig, sv_p_otau, sv_p_oppi 
     o  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o     , vval_rv , fes1_rv
      use m_seneinterp,only: senex,getsenex,dsene
      use m_procar,only: m_procar_init,m_procar_closeprocar
      use m_clsmode,only: m_clsmode_set1
      
!! outputs of bandcal ---------------------------
      real(8),protected::  sumqv(3,2),sumev(3,2)
      integer,allocatable,protected:: ndimhx_(:),nev_(:),nevls(:,:)
      real(8),allocatable,protected:: evlall(:,:,:),frcband(:,:), orbtm_rv(:,:,:)
      complex(8),allocatable,protected::  smrho_out(:),dmatu(:,:,:,:)
      type(s_rv1),allocatable,protected:: sv_p_oeqkkl(:,:), sv_p_oqkkl(:,:)
!! ------------------------------------------------     
      logical,private:: debug,sigmamode,call_m_bandcal_2nd,procaron,writeham
      logical,private:: dmatuinit=.true.
      contains
! Set up Hamiltonian, diagonalization and accumulation to srout and so on.
      subroutine m_bandcal_init(iqini,iqend,lrout,ef0,ifih,lwtkb) !,kpproc)
c     &     evlall,dos_rv,dmatu,orbtm_rv,smrho_out,sumqv,sumev,sv_p_oqkkl,sv_p_oeqkkl, frc, ndimhx_,nev_)
      implicit none
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:)!,ovlmtoi(:,:),ovliovl(:,:) 
      integer:: iq,nmx,ispinit,isp,jsp,nev,ifih,lwtkb,iqini,iqend,lrout,ifig
      real(8):: qp(3),ef0,def=0d0,xv(3)
      real(8),allocatable:: evl(:,:)
      complex(8),allocatable :: t_zv(:,:)
      integer::  pidorigin,iprint,i,ibas
      logical:: ltet,cmdopt0,dmatuinit=.true.
      if(master_mpi) write(stdo,*)' mmmmm m_bandcal_init'
      
      sigmamode = mod(lrsig,10)/=0
      writeham= cmdopt0('--writeham')
      PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
      debug    = cmdopt0('--debugbndfp')
c      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') !fullmesh stop just after do 2010 iq loop.
!! pdos mode = PROCARon .and. fullmesh (--mkprocar and --fullmesh. See job_pdos).      
      ltet = ntet>0
      if(plbnd==0.and.lso/=0.and.lwtkb==0)call rx('metal weights required to get orb.moment')
      
      if(lso/=0) allocate(orbtm_rv(nl,nsp,nbas)) !for spin-orbit coupling
      if(lso/=0) orbtm_rv=0d0
      
      allocate( evlall(ndhamx,nspx,nkp))
      if(lfrce>0) then
         allocate( frcband(3,1:nbas*numq)) !force for band
         frcband  = 0d0
      endif   
      allocate( ndimhx_(nkp),nev_(nkp),nevls(nkp,nspx))
      ndimhx_=0
      nev_   =0                
      nevls  =0
      
      if(nlibu>0.and.dmatuinit) then
        allocate( dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu))
        dmatuinit=.false.
      endif  
      
      if(lrout/=0) then
         allocate( sv_p_oeqkkl(3,nbas), sv_p_oqkkl(3,nbas))
         call dfqkkl( sv_p_oqkkl ) !zero clear 
         if(lekkl==1) call dfqkkl( sv_p_oeqkkl )!zero clear
         allocate( smrho_out(k1*k2*k3*numq*nsp) )
         smrho_out = 0d0
      endif

      call_m_bandcal_2nd =.false. 
      if(plbnd==0) call_m_bandcal_2nd= (lmet>=0 .and. lrout>0 )
      if(call_m_bandcal_2nd) open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
      allocate( evl(ndhamx,nspx))
      
!! These are accumulation varivables 
      if(nlibu>0)  dmatu=0d0    !density matrix for U initialization
      sumev = 0d0
      sumqv = 0d0
      if (lswtk==1)  rv_a_oswtk=0d0
      do 2010 iq = iqini, iqend !This is a big iq loop
         qp = qplist(:,iq)
         if(debug) print *,' do 2010 iq procid=',iq,procid,iq,iqini,iqend
         if(iq==iqini) call mlog_MPIiq(iq,iqini,iqend)
         call m_Igv2x_set(qp)   ! Get napw and so on for given qp
         nmx=min(nevmx,ndimhx)!nmx: max number of eigenfunctions we will obtain. Smaller is faster.
         if(lso==1) then
            allocate(hamm(ndimh,ndimh,4),ovlm(ndimh,ndimh,4)) !spin offdiagonal included
            ispinit = 2         !if lso, nsp=2
         else   
            if(lso/=1) allocate(hamm(ndimh,ndimh,1),ovlm(ndimh,ndimh,1)) !only for one spin
            ispinit=1
         endif   
         do 2005 isp = ispinit,nsp
            jsp = isp
            if(lso==1) jsp = 1
!! hambl calls augmbl.  augmbl calculates C_kl(=array b in the program)
!! finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s.
!! Note that F2=Hankel head at a site + Hankel tail contributions from the other site.
            
!! == Set up Hamiltonian by hambl. ==============
!!    Hamiltonian: hamm(1:ndimh,1:ndimh,3) means off-diagonal section when SO=1.
!!    Overlap matrix: ovlm 
!!    senex:  Sigma-Vxc
!! ==========================================
!! Generate sene(Sigma-Vxc) for given sfz.
!! Determine interpolated self-energy sene at qp from sfz. 
!! sigmat = Sigma-Vxc is generated in a basis of ndimsig (usually MTOs only) 
!!     ... Bloch transform sigm(RS)-sigm(k). :RS means realspace
!! Main input  => ham_iv_a_oiaxs,ham_rv_a_ohrs
!! Main output sene. See m_seneinterp
            hamm=0d0
            ovlm=0d0
            if(lso==1) then !L.S case
              call hambl(lso,1,qp,   zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,
     &              hamm(1,1,   1),ovlm(1,1, 1),hamm(1,1, 3)) !isp=1
              call hambl(lso,2,qp,   zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,
     &              hamm(1,1,   2),ovlm(1,1, 2),hamm(1,1, 3)) !isp=2 note readme in hambl
              call sopert2 ( hamm , hamm ) !SO case 
              call sopert2 ( ovlm , ovlm )
              if(sigmamode) then
                 call getsenex(lrsig,qp, 1, ndimh,ovlm(1,1, 1))
                 hamm(:,:,1) = hamm(:,:,1) + ham_scaledsigma*senex !senex_up= Vxc(QSGW)-Vxc(LDA)
                 call getsenex(lrsig,qp, 2, ndimh,ovlm(1,1, 2))
                 hamm(:,:,2) = hamm(:,:,2) + ham_scaledsigma*senex !senex_dn= Vxc(QSGW)-Vxc(LDA)
              endif
            else ! No SO or Lz.Sz case
              call hambl(lso,isp,qp, zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,
     &           hamm(1,1, 1),ovlm(1,1, 1),hamm(1,1, 3)) ! hamm(1,1,3) is dummy
              if(sigmamode) then !!Add  Vxc(QSGW)-Vxc
                 call getsenex(lrsig,qp,isp,ndimh,ovlm(1,1, 1))
                 hamm(:,:, 1) = hamm(:,:, 1) + ham_scaledsigma*senex !senex= Vxc(QSGW)-Vxc(LDA)
                 call dsene()
              endif
            endif 
            
            if(iprint()>=30) 
     &      write(stdo,'(" bndfp: kpt ",i5," of ",i5, " k=",3f8.4," ndimh = nmto+napw = ",3i5,f13.5)')
     &           iq,nkp,qp,ndimh,ndimh-napw,napw
            if(writeham) then 
                write(stdo,"(a,3f9.5)") "Hamiltonian: Writing hamm and ovlm for qp= ",qp
                write(ifih) qp,ndimhx,lso,epsovl,jsp
                if(lso==1) then !L.S case 
                  write(ifih) ovlm(1:ndimhx,1:ndimhx,1:4) 
                  write(ifih) hamm(1:ndimhx,1:ndimhx,1:4)
                else  !nspc=1 spin diagonal case
                  write(ifih) ovlm(1:ndimhx,1:ndimhx,1)
                  write(ifih) hamm(1:ndimhx,1:ndimhx,1)
                endif
            endif
             
             allocate(t_zv(ndimhx,nmx))
!! == Diagonalize Hamiltonian ==
!! ndimhx: dimension of Hamitonian
!! hamm:Hamiltonian, ovlm: overlap matrix
!! t_zv:eigenfunciton. evl: eigenvalue.
!! nmx: input, number of requested eigenvalues(functions). If nmx=0, no eigenfunctions but all eigenvalues.
!! nev: out number of obtained eigenfvalues(funcitons)
             call zhev_tk4(ndimhx, hamm, ovlm, nmx, nev, evl(1, jsp ), t_zv, epsovl)
             if(writeham.and.master_mpi) call prtev(t_zv, ndimhx , evl(1, jsp ) , nmx , nev )
             
             if(call_m_bandcal_2nd) then
                write(ifig) nev,nmx,ndimhx
                write(ifig) evl(1:nev,jsp)
                write(ifig) t_zv(1:ndimhx,1:nmx)
             endif
           
!! nev: number of eigenvalues 
            evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
            nevls(iq,jsp) = nev  !nov2014 isp and jsp is confusing...
            evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)
            if(master_mpi.and. epsovl>=1.000001d-14.and.plbnd/=0) then
              write(stdo,"(' : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
            endif
            
            if(plbnd==0.and.lwtkb/=-1) then !lwtkb=-1,0,1
               if(lso/=0)  call mkorbm(jsp, nev, iq,qp, t_zv,  orbtm_rv) !Orbital magnetic moment 
               if(nlibu>0.and.nev>0.and.lwtkb==0) call rx('metal weights required for LDA+U calculation')
               if(nlibu>0.and.nev>0) call mkdmtu(jsp, iq, qp, nev, t_zv,  dmatu)
               if(cmdopt0('--cls'))  call m_clsmode_set1(nmx,jsp,iq,qp,nev,t_zv) !all inputs
            endif
            
            if(lrout/=0.and.lwtkb>=0) then
c               call readindensitymodesetting() !dummy
               !Accumulate output density and sampling DOS. To make DOS when lmet=4
               call addrbl (jsp, qp, lrout
     .              , iq , lfrce,  zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi
     &              , t_zv,evl,nev,  ef0,def  !emin,emax,ndos, dos_rv
     o              , smrho_out, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frcband)
            endif     
            if(PROCARon) call m_procar_init(iq,isp,ef0,evl,ndimh,jsp,qp,nev,t_zv,ndimhx,nmx)
            if(allocated(t_zv)) deallocate(t_zv)
 2005    continue !== end loop over isp (main loop in parallel mode)==
         if(allocated(hamm)) deallocate(hamm,ovlm)
         ndimhx_(iq)= ndimhx
         nev_(iq)   = nev
 2010 continue                  !end of iq loop
!! ... Average forces so net force on system is zero (APW case)
      if (pwemax>0.and.mod(pwmode,10)>0.and.lfrce/=0) then
         do i = 1, 3
           xv(i)=sum(frcband(i,1:nbas))/nbas
        enddo
        do  ibas= 1, nbas
            frcband(:,ibas) = frcband(:,ibas) - xv(:)
         enddo
      endif
      if(PROCARon) call m_procar_closeprocar()
      if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid
      if(call_m_bandcal_2nd) close(ifig)
      deallocate(evl)
      end subroutine
!!
      subroutine m_bandcal_2nd(iqini,iqend,lrout,ef0) !,emin,emax)!,ndos)
      implicit none
      integer:: iq,nmx,ispinit,isp,jsp,nev,iqini,iqend,lrout,ifig,i,ibas
      real(8):: qp(3),ef0,def=0d0,xv(3)
      real(8),allocatable:: evl(:,:)
      complex(8),allocatable :: t_zv(:,:)
      logical:: cmdopt0
      if(master_mpi) write(stdo,*)' mmmmm m_bandcal_2nd'
!accumulation      
      call dfqkkl( sv_p_oqkkl ) !zero clear 
      if(lekkl==1) call dfqkkl( sv_p_oeqkkl ) !zero clear
      if (lfrce>0)  frcband  = 0d0
      if (lswtk==1) rv_a_oswtk=0d0
      if(lso/=0) orbtm_rv=0d0
      if(allocated(smrho_out)) deallocate(smrho_out)
      allocate( smrho_out(k1*k2*k3*numq*nsp) )
      smrho_out = 0d0
      sumev = 0d0 
      sumqv = 0d0
      allocate(evl(ndhamx,nspx))
      open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
      do 12010 iq = iqini, iqend !This is a big iq loop
         qp = qplist(:,iq)
         call m_Igv2x_set(qp)   ! Get napw and so on for given qp
         if(lso==1) ispinit = 2
         if(lso/=1) ispinit=1
         do 12005 isp = ispinit,nsp
            jsp = isp
            if(lso==1) jsp = 1
            read(ifig) nev,nmx  !ndimhx <---supplied by m_Igv2x_set
            if (allocated(t_zv)) deallocate(t_zv)
            allocate(t_zv(ndimhx,nmx))
            read(ifig) evl(1:nev,jsp)
            read(ifig) t_zv(1:ndimhx,1:nmx)
            evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
            if( lso/=0)            call mkorbm(jsp, nev, iq,qp, t_zv,  orbtm_rv)
            if( nlibu>0.and.nev>0) call mkdmtu(jsp, iq,qp, nev, t_zv,  dmatu)
            if( cmdopt0('--cls'))  call m_clsmode_set1(nmx,jsp,iq,qp,nev,t_zv) !all inputs
            if( lrout/=0 ) then
c               call readindensitymodesetting() !dummy
               call addrbl (jsp, qp, lrout
     .              , iq , lfrce,  zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi
     &              , t_zv,evl,nev,  ef0,def !emin,emax,ndos, dos_rv
     o              , smrho_out, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frcband)
            endif     
            if(allocated(t_zv)) deallocate(t_zv)
12005    continue 
12010 continue    
!! ... Average forces so net force on system is zero (APW case)
      if (pwemax>0.and.mod(pwmode,10)>0.and.lfrce/=0) then
         do i = 1, 3
           xv(i)=sum(frcband(i,1:nbas))/nbas
        enddo
        do  ibas= 1, nbas
            frcband(:,ibas) = frcband(:,ibas) - xv(:)
         enddo
      endif
      close(ifig)
      deallocate(evl)
      end subroutine

      subroutine m_bandcal_clean()
      if (allocated(orbtm_rv)) deallocate(orbtm_rv)
c      if (allocated(dos_rv))   deallocate(dos_rv)
      if (allocated(smrho_out)) deallocate(smrho_out)
      if (allocated(frcband))  deallocate(frcband)
      if (allocated(ndimhx_))  deallocate(ndimhx_,nev_,nevls)
      if(allocated(sv_p_oqkkl)) deallocate( sv_p_oqkkl)
      if(allocated(sv_p_oeqkkl))deallocate( sv_p_oeqkkl)
      deallocate(evlall)
      end

!     !  Allreduce density-related quantities
      subroutine m_bandcal_allreduce(lrout,lwtkb)
      integer:: lrout,nnn,ib,i,lwtkb
      if(debug) print *,'goto m_bandcal_allreduce'
      if (lrout .ne. 0) then
         nnn=size(sumqv);    call mpibc2_real(sumqv,nnn,'bndfp_sumqv') 
         nnn=size(sumev);    call mpibc2_real(sumev,nnn,'bndfp_sumev')
         nnn=size(smrho_out); call mpibc2_complex(smrho_out,nnn,'bndfp_smrho')
         if (lswtk==1) then
            nnn=size(rv_a_oswtk);  call mpibc2_complex(rv_a_oswtk,'bndfp_swtk') 
         endif 
         do  ib = 1, nbas
         do  i = 1, 3
            if(allocated(sv_p_oqkkl(i,ib)%v)) then
               nnn = size(sv_p_oqkkl(i,ib)%v)
               if(nnn>0) call mpibc2_real(sv_p_oqkkl(i,ib)%v,nnn,'bndfp_qkkl')
            endif   
            if(lekkl==1.and.allocated(sv_p_oeqkkl(i,ib)%v)) then
               nnn = size(sv_p_oeqkkl(i,ib)%v)
               if(nnn>0) call mpibc2_real(sv_p_oeqkkl(i,ib)%v,nnn,'bndfp_eqkkl') 
            endif
         enddo
         enddo
c         if( ndos>0 ) nnn=size(dos_rv)
c         if( ndos>0 ) call mpibc2_real(dos_rv,nnn,'bndfp_dos')
         if(lfrce/=0) nnn=size(frcband)
         if(lfrce/=0) call mpibc2_real(frcband,nnn,'bndfp_frcband')        
         if(nlibu>0)  nnn=size(dmatu)
         if(nlibu>0)  call mpibc2_complex(dmatu,nnn,'bndfp_dmatu') 
         if(lso/=0.and.lwtkb/=-1) nnn=size(orbtm_rv)
         if(lso/=0.and.lwtkb/=-1) call mpibc2_real(orbtm_rv,nnn,'bndfp_orbtm') 
      endif
      end subroutine
      
!!------------------------------
c$$$      subroutine m_bandcal_dosw(lwtkb,lrout,  dosw,evtop,ecbot) !goto99)
c$$$      ! lwtkb,lrout,def are used only for
c$$$      use m_lmfinit,only:ctrl_zbak,bz_w
c$$$      intent(in)   ::           lwtkb,lrout
c$$$      intent(inout)::                             dosw !input is just initial guess
c$$$      intent(out)  ::                                  evtop,ecbot
c$$$c      intent(out)::                                  goto99 
c$$$!!  ===   Repeat loop for printout and goto 99 ===
c$$$!!   jsp=isp in the collinear case; jsp=1 in the noncollinear
c$$$!!     Thus jsp should be used in place of isp
c$$$!     !     isp serves as a flag for the noncollinear case
c$$$!     ! block10
c$$$      logical:: ltet !goto99,
c$$$      real(8):: ef00,ef0,dosw(2),qp(3),qbg,evl(ndhamx,nspx),dum,ebot,ecbot,evtop
c$$$      integer:: iq,ipr,lwtkb,isp,jsp,nev_iq,i,lrout,iprint
c$$$      character(10):: i2char
c$$$c      goto99=.false.
c$$$      ltet = ntet>0
c$$$      qbg = ctrl_zbak(1) !homogenious background charge
c$$$      ipr=iprint()
c$$$      evtop = -99999
c$$$      ecbot = -evtop
c$$$      ebot = 1000d0
c$$$      do  iq = 1, nkp
c$$$         qp=qplist(:,iq)
c$$$         do isp = 1, nspx
c$$$            jsp = isp
c$$$            evl(1:ndhamx,jsp) = evlall(1:ndhamx,jsp,iq)
c$$$            nev_iq    = nev_(iq) 
c$$$            if(plbnd==0.and.ipr>=10 .and. iq==1) write (stdl,"('fp evl',8f8.4)") (evl(i,jsp),i=1,nev_iq)
c$$$            ebot = dmin1(ebot,evl(1,jsp))
c$$$            i = max(1,nint(qval-qbg)/(3-nspc))
c$$$            evtop = max(evtop,evl(i,jsp))
c$$$            ecbot = min(ecbot,evl(i+1,jsp))
c$$$            if (lmet==0 .and. iq==1 .and. jsp==1) ef0 = evtop
c$$$            if(debug) print *,'eeeee44444444444 plbnd=',plbnd
c$$$            if(plbnd==0.and.lwtkb/=-1) then
c$$$                  if (iq .eq. 1 .and. jsp .eq. nsp ) then !
c$$$                     ef00 = ef0
c$$$                     call fixef0(qval-qbg,jsp,1,nev_iq,ndhamx,evl,dosw,ef0)
c$$$!!        :on output, dosw is revised if ebot<dosw(1) or dosw(2)<ef0
c$$$                     if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
c$$$     .                    lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
c$$$                          if (master_mpi) write(stdo,"(a)") 
c$$$     .                       ' ... Fermi level reset in second spin channel ... restart band pass'
c$$$                          call rx0('tk think ecalj is going though not maintained branch).')
c$$$                         !goto99=.true.
c$$$                        exit    !this was the case of make co test at ecalj/TestInstall/ (did I remove this?)
c$$$                     endif
c$$$                  endif
c$$$            endif
c$$$!!     Check for cases when nevmx is too small : i=2 => fatal error
c$$$            if(plbnd==0.and.lwtkb/=-1) then
c$$$                  i = 0
c$$$                  if (nevmx.ge.0 .and. lmet .ne. 0) then
c$$$                     dum = evl(max(nev_iq,1),jsp)
c$$$                     if (.not. ltet .and. ef0+5*bz_w .gt. dum) i=2
c$$$c                     if (lmet.eq.4 .and. ef0+def+5*bz_w .gt.dum)i=2
c$$$                  endif
c$$$                  if(i .eq. 2) then
c$$$                     write(stdo,"(a,f13.5,a, f13.5)")
c$$$     &                 'evl(nev='//trim(i2char(nev_iq))//')=',
c$$$     &                 evl(max(nev_iq,1),jsp),' but ef0=',ef0
c$$$                     call rx('bndfp:... restart with larger nevmx: bndfp')
c$$$                  endif
c$$$                  
c$$$            endif
c$$$         enddo                  ! end second loop over isp
c$$$      enddo                     !end second loop over iq
c$$$      end subroutine
      
!     ! eterms are alwasy broadcasted. comment out if conditions jan2012


      subroutine m_bandcal_symsmrho()
      call symsmrho(smrho_out)
      end
      
      end module
