      module m_bandcal
      
      contains
! Set up Hamiltonian, diagonalization and accumulation to srout and so on.
      subroutine m_bandcal_init( iqini, iqend,sigmamode,lrsig,lfrce,lrout,
     &     ef0,def,emin,emax,ndos,nevls,ifih,plbnd,lwtkb,kpproc,gotosecondpass, 
     &     evlall,dos_rv,dmatu,orbtm_rv,srout_zv,sumqv,sumev,sv_p_oqkkl,sv_p_oeqkkl, frc, ndimhx_,nev_)
      use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_set
      use m_MPItk,only: mlog, master_mpi, procid,strprocid, numprocs=>nsize, mlog_MPIiq
      use m_subzi, only: m_subzi_init, nevmx,numq,lswtk,rv_a_owtkb,rv_a_oswtk,ldos
      use m_lmfinit, only: lso,nsp,ham_scaledsigma,nlibu,lmaxu,
     & lmet=>bz_lmet,stdo,nbas,epsovl=>ham_oveps,nspc
      use m_mkpot,only: m_mkpot_init,m_mkpot_deallocate,
     o    zv_a_osmpot, qmom_rv, vconst, vrmt, sv_p_osig, sv_p_otau, sv_p_oppi 
     o  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     o     , vval_rv , fes1_rv
      use m_struc_def
      use m_seneinterp,only: senex,getsenex,dsene
      use m_qplist,only: nkp, qplist
      use m_suham,only: iv_a_oindxo, nchan=>pot_nlma, nvl=>pot_nlml, ndham=>ham_ndham,ham_ldham
     & ,ndhamx=>ham_ndhamx,nspx=>ham_nspx 
      use m_procar,only: m_procar_init,m_procar_closeprocar
      use m_clsmode,only: m_clsmode_set1
      implicit none
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),ovlmtoi(:,:),ovliovl(:,:) 
      logical:: debug,sigmamode,fullmesh,gotosecondpass,procaron,writeham,cmdopt0
      integer:: iq,nmx,ispinit,isp,lrsig,jsp,nev,ifih,plbnd,lwtkb,lfrce,ndos,kpproc(*),iqini,iqend,lrout,ifig
      real(8):: qp(3),ef0,def,emin,emax,sumqv(*),sumev(*)
      real(8),allocatable:: evl(:,:)
      complex(8),allocatable :: t_zv(:,:)
      complex(8)::  dmatu(*),srout_zv(*)
      real(8)::  dos_rv(*),frc(*),orbtm_rv(*)
      type(s_rv1) :: sv_p_oeqkkl(*)
      type(s_rv1) :: sv_p_oqkkl(*)
      real(8):: evlall(ndhamx,nspx,nkp)
      integer::  ndimhx_(nkp),nev_(nkp),nevls(nkp,nspx),pidorigin,iprint

      allocate(evl(ndhamx,nspx))
      writeham= cmdopt0('--writeham')
      PROCARon = cmdopt0('--mkprocar') !write PROCAR(vasp format).
      debug    = cmdopt0('--debugbndfp')
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') !fullmesh stop just after do 2010 iq loop.
      if(gotosecondpass) open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
      do 2010 iq = iqini, iqend  !This is a big iq loop
         qp = qplist(:,iq)
         if(debug) print *,' do 2010 iq procid=',iq,procid,iq,iqini,iqend
         call mlog_MPIiq(iq,kpproc)
         call m_Igv2x_set(qp)   ! Get napw and so on for given qp
         nmx=min(nevmx,ndimhx)!nmx: max number of eigenfunctions we will obtain. Smaller is faster.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cC... range check qqq is in [ikc1/nk1, (ikc1+1)/nk1],... qpara is [0,1],[0,1],[0,1]
c      print *,' --- goto qparam--- '
c      call qparam(qp,qsmesh,nk1,nk2,nk3,plat, ikc,qpara,igshift)
c      ! eight corners -->  qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
c      call interparam(qpara, wgc)
c      cycle
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
!!note kino. hambl calls augmbl.  augmbl calculates C_kl(=array b in the program)
!! finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s.
!! Note that F2=Hankel head at a site + Hankel tail contributions from the other site.
            
!! == Set up Hamiltonian by hambl. ==============
!!    Hamiltonian: hamm(1:ndimh,1:ndimh,3) means off-diagonal section when SO=1.
!!    Overlap matrix: ovlm 
!!    senex:  Sigma-Vxc
!! ==========================================
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
            else !no SO or Lz.Sz case
              call hambl(lso,isp,qp, zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi,
     &           hamm(1,1,1),ovlm(1,1,1),hamm(1,1,3)) ! hamm(1,1,3) is dummy
!! Generate sene(Sigma-Vxc) for given sfz.
!! Determine interpolated self-energy sene at qp from sfz. 
!! sigmat = Sigma-Vxc is generated in a basis of ndimsig (usually MTOs only) 
!!     ... Bloch transform sigm(RS)-sigm(k). :RS means realspace
!! Main input  => ham_iv_a_oiaxs,ham_rv_a_ohrs
!! Main output sene. See m_seneinterp
              if(sigmamode) then !!Add  Vxc(QSGW)-Vxc
                 call getsenex(lrsig,qp,isp,ndimh,ovlm(1,1,1))
                 hamm(:,:,1) = hamm(:,:,1) + ham_scaledsigma*senex !senex= Vxc(QSGW)-Vxc(LDA)
                 call dsene()
              endif
            endif 
            
            if(iprint()>=30) 
     &           write(stdo,'(" bndfp: kpt ",i5," of ",i5, " k=",3f8.4," ndimh = nmto+napw = ",3i5,f13.5)')
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
!! nmx: input, number of requested eigenvalues(funcitons). If nmx=0, no eigenfunctions but all eigenvalues.
!! nev: out number of obtained eigenfvalues(funcitons)
             call zhev_tk4(ndimhx, hamm, ovlm, nmx, nev, evl(1, jsp ), t_zv, epsovl)
             if(writeham.and.master_mpi) call prtev(t_zv, ndimhx , evl(1, jsp ) , nmx , nev )
             
             if(gotosecondpass) then
                write(ifig) nev,nmx,ndimhx
                write(ifig) evl(1:nev,jsp)
                write(ifig) t_zv(1:ndimhx,1:nmx)
             endif
           
!! nev: number of eigenvalues 
            evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
            nevls(iq,jsp) = nev  !nov2014 isp and jsp is confusing...
            if( master_mpi.and. epsovl>=1.000001d-14.and.plbnd/=0) then
              write(stdo,"(' : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
            endif
            evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)
            
!!   ... Save data for this qp to moments file
            if(plbnd==0.and.lwtkb/=-1) then
!!   --- Orbital magnetic moment (requires lso) ---
               if (lso/=0) then
                  if (lwtkb .eq. 0) call rx('metal weights required for orb. moment calculation')
                  call mkorbm(jsp, nev, iq,qp, t_zv,  orbtm_rv)
               endif
              if (nlibu>0.and.nev>0) then !density matrix dmatu for LDA+U 
                  if (lwtkb .eq. 0) call rx('metal weights required for LDA+U calculation')
                  call mkdmtu(jsp, iq, qp, nev, t_zv,  dmatu)
               endif
               if(cmdopt0('--cls')) then !!! Core-level spectroscopy
                  call m_clsmode_set1(nmx,jsp,iq,qp,nev,t_zv) !all inputs
               elseif(lrout.ne.0 .or. lmet.eq.4) then
                  call readindensitymodesetting() !dummy
                 if(lwtkb>=0) then !Accumulate output density and sampling DOS. To make DOS when lmet=4
                   call addrbl (jsp, qp, lrout
     .             , iq , lfrce,  zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi
     &             , t_zv,evl,nev,  ef0,def,emin,emax,ndos
     o             , dos_rv, srout_zv, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frc)
                 endif     
               endif
            endif
            if(PROCARon) call m_procar_init(fullmesh,iq,isp,nspx,nspc,ef0,evl,ndimh,jsp,qp,nev,t_zv,ndimhx,nmx)
            if(allocated(t_zv)) deallocate(t_zv)
 2005    continue !== end loop over isp (main loop in parallel mode)==
         if(allocated(hamm)) deallocate(hamm,ovlm)
         ndimhx_(iq)= ndimhx
         nev_(iq)   = nev
 2010 continue                  !end of iq loop
      if(PROCARon) call m_procar_closeprocar()
      if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid
      if(gotosecondpass) close(ifig)
      if(writeham) close(ifih)
      end subroutine
      end module
