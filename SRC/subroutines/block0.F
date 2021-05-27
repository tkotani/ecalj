      ix99=1
      if(lrout/=0) then
         call dfqkkl( sv_p_oqkkl ) !zero clear 
         if(lekkl==1) call dfqkkl( sv_p_oeqkkl )!zero clear
         if(allocated(fh_rv)) deallocate(srout_zv,fh_rv,fes2_rv)
         allocate( srout_zv(k1*k2*k3*numq*nsp) )
         allocate( fh_rv(3*nbas)  )
         allocate( fes2_rv(3*nbas))
         srout_zv = 0d0
         fh_rv= 0d0  
      endif
      ebot = 1000d0
      sumev = 0d0
      sumqv = 0d0
      if (ldos==1 .or. lmet>0) dos_rv=0d0
      if (lfrce>0) frc  = 0d0
      if (lswtk==1) rv_a_oswtk=0d0
!! Mullikan mode was here. removed for simplicity
!! For insulator, valence band top and conduction band bottom
      evtop = -99999
      ecbot = -evtop
      if(lso/=0) orbtm_rv=0d0
      if(gotosecondpass) open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
      
!     ! Start k and isp loop 2010 and 2005
      do 12010 iq = iqini, iqend !This is a big iq loop
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
c            allocate(hamm(ndimh,ndimh,4),ovlm(ndimh,ndimh,4)) !spin offdiagonal included
            ispinit = 2         !if lso, nsp=2
         else   
c            if(lso/=1) allocate(hamm(ndimh,ndimh,1),ovlm(ndimh,ndimh,1)) !only for one spin
            ispinit=1
         endif   
         do 12005 isp = ispinit,nsp
            jsp = isp
            if(lso==1) jsp = 1
c            if(gotosecondpass.and.ix99==2) then   ! second pass
              if (allocated(t_zv)) deallocate(t_zv)
              allocate(t_zv(ndimhx,nmx))
              read(ifig) nev,nmx
              read(ifig) evl(1:nev,jsp)
              read(ifig) t_zv(1:ndimhx,1:nmx)
c            endif
           
!! nev: number of eigenvalues 
            evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
            nevls(iq,jsp) = nev  !nov2014 isp and jsp is confusing...
            if( master_mpi.and. epsovl>=1.000001d-14.and.plbnd/=0) then
              write(stdo,"(' : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
            endif
            evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)
            
!!   ... Save data for this qp to moments file
            if(plbnd==0) then
!!   --- Orbital magnetic moment (requires lso) ---
               if (lso/=0.and.lwtkb/=-1) then
                  if (lwtkb .eq. 0) call rx('metal weights required for orb. moment calculation')
                  allocate(auso_zv(nlmax*ndhamx*3*nsp*nbas))
                  auso_zv=0.0d0
                  call makusq(0 , nbas,0, nev, jsp,1,qp,t_zv, auso_zv )
                  call mkorbm(jsp, nev, iq,                   auso_zv, orbtm_rv )
                  deallocate(auso_zv)
               endif
c$$$!! Mulliken analysis was here
!! Make new density matrix dmatu for LDA+U (implementation of LDA+U is a little complicated).
               if (lwtkb .ne. -1) then
                  if (nlibu>0.and.nev>0) then
                     if (lwtkb .eq. 0) call rx('metal weights required for LDA+U calculation')
                     allocate(ausu_rv(2*nlmax*ndhamx*3*nsp*nbas))
                     ausu_rv=0.0d0
                     call makusq(0, nbas,0 , nev,  jsp, 1, qp, t_zv, ausu_rv )
                     call mkdmtu(jsp, iq, nev, ausu_rv , dmatu)
                     if (allocated(ausu_rv)) deallocate(ausu_rv)
                  endif
               endif
!! Core-level spectroscopy
               if(cmdopt0('--cls').and.lwtkb.ne.-1) then
                  call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
                  call m_clsmode_set1(nmx,jsp,iq,qp,nev,t_zv) !all inputs
               elseif (lwtkb.ne.-1 .and. (lrout.ne.0 .or. lmet.eq.4)) then
                  call readindensitymodesetting() !dummy
!! Accumulate output density and sampling DOS
!! Even if no output rho, still call addrbl to make DOS when lmet=4
                 if(lwtkb>=0) then
                   call addrbl (jsp, qp, lrout
     .             , iq , lfrce,  zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi
     &             , t_zv,evl,nev,  ef0,def,emin,emax,ndos
     o             , dos_rv, srout_zv, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frc)
                 endif     
               endif
            endif
            if(PROCARon) call m_procar_init(fullmesh,iq,isp,nspx,nspc,ef0,evl,ndimh,jsp,qp,nev,t_zv,ndimhx,nmx)
            if(allocated(t_zv)) deallocate(t_zv)
12005    continue !== end loop over isp (main loop in parallel mode)==
c         if(allocated(hamm)) deallocate(hamm,ovlm)
         ndimhx_(iq)= ndimhx
         nev_(iq)   = nev
12010 continue                  !end of iq loop
c      if(PROCARon) m_procar_closeprocar()
      if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid
      if(gotosecondpass) close(ifig)
