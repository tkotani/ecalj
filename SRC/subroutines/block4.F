!!  Allreduce density-related quantities
      if(debug) print *,'goto all reduce'
      if (lrout .ne. 0) then
         nnn=size(sumqv);    call mpibc2_real(sumqv,nnn,'bndfp_sumqv') 
         nnn=size(sumev);    call mpibc2_real(sumev,nnn,'bndfp_sumev')
         nnn=size(srout_zv); call mpibc2_complex(srout_zv,nnn,'bndfp_smrho')
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
         if( ndos>0 ) nnn=size(dos_rv)
         if( ndos>0 ) call mpibc2_real(dos_rv,nnn,'bndfp_dos')
         if(lfrce/=0) nnn=size(frc)
         if(lfrce/=0) call mpibc2_real(frc,nnn,'bndfp_frc')        
         if(nlibu>0)  nnn=size(dmatu)
         if(nlibu>0)  call mpibc2_complex(dmatu,nnn,'bndfp_dmatu') 
         if(lso/=0.and.lwtkb/=-1) nnn=size(orbtm_rv)
         if(lso/=0.and.lwtkb/=-1) call mpibc2_real(orbtm_rv,nnn,'bndfp_orbtm') 
      endif
!! eterms are alwasy broadcasted. comment out if conditions jan2012
      eterms=ham_eterms ! call dcopy(size(ham_eterms),ham_eterms,1,eterms,1)
      call mpibc2_real(eterms(19),1,'bndfp_rhosig')
      ham_eterms=eterms ! call dcopy(size(ham_eterms),eterms,1,ham_eterms,1)
      entime = MPI_WTIME()
      call info2(20,0,0,' MPI broadcast took %;1d sec', (entime-sttime),0)



      
