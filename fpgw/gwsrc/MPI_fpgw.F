      subroutine MPI_sxcf_rankdivider(nrank,nspinmx,nqibz,ngrp,nq,irkip_all, irank,irkip)
      implicit none
      integer:: nspinmx,nqibz,ngrp,nq, nrank,irank
      integer:: irkip(nspinmx,nqibz,ngrp,nq)
      integer:: irkip_all(nspinmx,nqibz,ngrp,nq)
      if(nrank==1.and.irank==0) then
         irkip=irkip_all
      else
         stop 'not implemented yet'
      endif   
      end
