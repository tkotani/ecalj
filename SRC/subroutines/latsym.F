      subroutine latsym(prgnam)
      use m_lmfinit,only: v_ssite,nbas,sstrnsymg,ctrl_noinv,
     &     symgaf,iv_a_oips,slabl_,mxspec,nsite,procid,master
! v_ssite%pos modified

      use m_lattic,only: rv_a_opos,m_lattic_init
      use m_mksym,only: mksymaf,mksym
!! Set Symmetry. Driver for calling mksymaf and mksym
!-------------      
!!note: Antiferro symmetric self-energy mode. Feb2021 recovered
!! To obtain self-consistency, it may be useful to keep AF condition during iteration.
!! We need to set SYMGRPAF (still developing...)
! sample for nio
!     SYMGRPAF i:(1,1,1) ! translation + inversion makes antiferro symmetry.
!     SYMGRP r3d   !this keeps spin axis. (probably larger symmetry allowed for SO=0)
!STRUC   ALAT={a} PLAT= 0.5 0.5 1.0  0.5 1.0 0.5  1.0 0.5 0.5
!        NBAS= 4  NSPEC=3
!SITE    ATOM=Niup POS=  .0   .0   .0   AF=1    <--- AF symmetric pair
!        ATOM=Nidn POS= 1.0  1.0  1.0   AF=-1   <---
!        ATOM=O POS=  .5   .5   .5
!        ATOM=O POS= 1.5  1.5  1.5
!------------------------------------------------------------
      implicit none
      character,intent(in)::  prgnam*(*)
      integer:: ibas,lc,j,iprint
      integer,parameter::recln=511
      logical ::cmdopt,ipr10
      character strn*(recln),strn2*(recln),outs(recln)
!! lattice set up      
      call m_lattic_init() 
      do ibas=1,nbas
        v_ssite(ibas)%pos = rv_a_opos(:,ibas)
      enddo
!! mksym mksymaf generate symmetry operations; split species into classes  ---
!! mksym may correct atomic position slightly (not tested well).
      strn = 'find'
      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
      if (cmdopt('--nosym',7,0,outs).or.cmdopt('--pdos',6,0,outs)) strn = ' '
!! when lmfgw, 1st digit of lc is zero--> no inversion added in mksym.
      lc = 10
      if (.not. prgnam == 'LMFGWD') lc = 12 ! add inversion for (not lmfgw)
      if (prgnam=='LMFA'.or.prgnam=='LMCHK') then
         lc=20
         if(ctrl_noinv /=1 ) lc = lc+1  !+1 means add inversion 
      endif
      if(len_trim(symgaf)>0) lc=12 ! inversion allowed for AF case. Good?
      if(.not.(prgnam=='LMFA'.or.prgnam=='LMCHK')) ipr10= iprint()>10 !this is only for master
      if(len_trim(symgaf)>0) then
        if(ipr10) then
           write(6,*)
           write(6,"(a)")       ' AF: ======================================== '
           write(6,"(a)")       ' AF: Antiferro mode: SPGGRAF='//trim(symgaf)
           write(6,"(a)")       ' AF:  (neglct waring in GENSYM) '
           do j=1,nbas
              write(6,"(a,2i3)") ' AF:  ibas,AF=',j,v_ssite(j)%iantiferro
           enddo
        endif
        strn2=trim(strn)//' '//trim(symgaf)
        call mksymaf(v_ssite,iv_a_oips,nbas,procid==master,strn2,lc,slabl_,mxspec,nsite)
        if(ipr10) write(6,"(a)") ' AF: ===== end of AF section================= '
        if(ipr10) write(6,"(a)")       
      endif
      if(procid==master) call pshpr(60)
      call mksym(lc,slabl_,strn,v_ssite,iv_a_oips) 
      if(procid==master) call poppr
      end subroutine latsym
