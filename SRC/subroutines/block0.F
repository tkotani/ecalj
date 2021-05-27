!     ! ix99=2 second path
      print *,'second path gggggg999 goto ix99 lwtkb=',lwtkb
      if(lrout/=0) then
         call dfqkkl( sv_p_oqkkl ) !zero clear 
         if(lekkl==1) call dfqkkl( sv_p_oeqkkl )!zero clear
         srout_zv = 0d0
         fh_rv= 0d0
      endif
      sumev = 0d0
      sumqv = 0d0
      if (ldos==1 .or. lmet>0) dos_rv=0d0
      if (lfrce>0)  frc  = 0d0
      if (lswtk==1) rv_a_oswtk=0d0
      if(lso/=0)    orbtm_rv=0d0
      if(gotosecondpass) open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
      
!! Start k and isp loop 2010 and 2005
      do 12010 iq = iqini, iqend !This is a big iq loop
         qp = qplist(:,iq)
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
            read(ifig) nev,nmx
            if (allocated(t_zv)) deallocate(t_zv)
            allocate(t_zv(ndimhx,nmx))
            read(ifig) evl(1:nev,jsp)
            read(ifig) t_zv(1:ndimhx,1:nmx)
c            endif
            
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
c            if(PROCARon) call m_procar_init(fullmesh,iq,isp,nspx,nspc,ef0,evl,ndimh,jsp,qp,nev,t_zv,ndimhx,nmx)
            if(allocated(t_zv)) deallocate(t_zv)
12005    continue !== end loop over isp (main loop in parallel mode)==
c         if(allocated(hamm)) deallocate(hamm,ovlm)
12010 continue                  !end of iq loop
      if(gotosecondpass) close(ifig)
