      subroutine phispinsym_ssite_set()
      use m_lmfinit,only: nbas,nsp,sspec=>v_sspec,n0,ssite=>v_ssite
      use m_MPItk,only:master_mpi
      use m_struc_def
      implicit none
      integer:: ib,i_copy_size,is,lmxa,l
c      type(s_site):: ssite(*) ssite(:) did not work (2021jun)
      real(8):: pnu(n0,2),pnz(n0,2),pmean
      if(master_mpi)write(6,*)' --phispinsym use spin-averaged potential for phi and phidot' 
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      ib=1
c      i_copy_size=size(ssite(ib)%pnu)
c         print *,'111 ibas pnu=',ib,i_copy_size
c         print *,'222 ibas pnu=',ssite(ib)%pnu(1:10)
c         print *,'333 ibas pnu=',ssite(ib)%pnu(11:20)
c         call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
c         stop 'xxxxxxxxxxxxxxxx'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
      do ib = 1,nbas
         i_copy_size=size(ssite(ib)%pnu)
         call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
         i_copy_size=size(ssite(ib)%pz)
         call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
         is   = ssite(ib)%spec
         lmxa = sspec(is)%lmxa
         do l=0,lmxa
            pmean = sum(pnu(l+1,1:nsp))/nsp
            if(master_mpi.and.nsp==2) write(6,"('  ibas l=',2i3,' pnu=',2f10.5,' -->',f10.5)")
     &           ib,l,pnu(l+1,1:nsp),pmean
            pnu(l+1,1:nsp) = pmean
            if(pnz(l+1,1)/=0 ) then
               pmean = sum(pnz(l+1,1:nsp))/nsp
               if(master_mpi.and.nsp==2)write(6,"('  ibas l=',2i3,' pnz=',2f10.5,' -->',f10.5)")
     &              ib,l,pnz(l+1,1:nsp),pmean
               pnz(l+1,1:nsp) = pmean
            endif
         enddo
         i_copy_size=size(ssite(ib)%pnu)
         call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
         i_copy_size=size(ssite(ib)%pz)
         call dcopy(i_copy_size,pnz,1,ssite(ib)%pz,1)
      enddo
      end
