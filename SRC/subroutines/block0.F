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
      evtop = -99999 !! For insulator, valence band top and conduction band bottom
      ecbot = -evtop
      if(lso/=0) orbtm_rv=0d0
      open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
!! Start k and isp loop 2010 and 2005
      do 12010 iq = iqini, iqend !This is a big iq loop
         qp = qplist(:,iq)
         call m_Igv2x_set(qp)   ! Get napw and so on for given qp
         if(lso==1) ispinit = 2
         if(lso/=1) ispinit=1
         do 12005 isp = ispinit,nsp
            jsp = isp
            if(lso==1) jsp = 1
            if (allocated(t_zv)) deallocate(t_zv)
              read(ifig) nev,nmx
              allocate(t_zv(ndimhx,nmx))
              read(ifig) evl(1:nev,jsp)
              read(ifig) t_zv(1:ndimhx,1:nmx)
            evl(nev+1:ndhamx,jsp)=1d99 !to skip these data
            if(plbnd==0.and.lwtkb.ne.-1) then
               if (lso/=0) then
                  if (lwtkb .eq. 0) call rx('metal weights required for orb. moment calculation')
                  call mkorbm(jsp, nev, iq,qp, t_zv,  orbtm_rv)
               endif
               if (nlibu>0.and.nev>0) then
                     if (lwtkb .eq. 0) call rx('metal weights required for LDA+U calculation')
                     call mkdmtu(jsp, iq,qp, nev, t_zv,  dmatu)
               endif
               if(cmdopt0('--cls')) then
                  call m_clsmode_set1(nmx,jsp,iq,qp,nev,t_zv) !all inputs
               elseif (lrout.ne.0 .or. lmet.eq.4) then
                  call readindensitymodesetting() !dummy
                  if(lwtkb>=0) then
                   call addrbl (jsp, qp, lrout
     .             , iq , lfrce,  zv_a_osmpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi
     &             , t_zv,evl,nev,  ef0,def,emin,emax,ndos
     o             , dos_rv, srout_zv, sumqv, sumev, sv_p_oqkkl,sv_p_oeqkkl, frc)
                  endif     
               endif
            endif
            if(allocated(t_zv)) deallocate(t_zv)
12005    continue !== end loop over isp (main loop in parallel mode)==
         ndimhx_(iq)= ndimhx
         nev_(iq)   = nev
12010 continue                  !end of iq loop
      close(ifig)
