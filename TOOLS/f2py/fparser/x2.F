      subroutine bndfp(nbas,nsp,nlibu,lmaxu,lldau,
     .  ssite,sspec,slat,sctrl,sham,spot,sbz,sstrn,
     .  ndham,leks,lrout,lfrce,lpnu,dmxp,iter,maxit,evl,frc,
     .  dmatu,vorb)
      use m_hamindex
        if (ldim .ne. ndham) then !june2009 takao
           write(6,"(a,2i5)")" ldim(dim of lmto) = ",ldim,
     &     " ndham(max of ham dimension)= ",ndham
          write(6,*)"go to rdsigm for PMT mode..." 
        endif
      end

