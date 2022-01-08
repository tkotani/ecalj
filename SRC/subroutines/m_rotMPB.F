      module m_rotMPB2 !WARN. rotMPB looks work OK. Used in eibz4x0()=T in switch.F  However,
!!---------------------------------------------------
!! == Mixed product basis rotator. ==
!! a curious strange segmention fault occurs in a Si444 case when we test a code qout= matmul(symops(:,:,igx),qin). 
!! I don't know why. 
!! This has developed for eibz mode jun2012, not tested completely. 
!! In future, this routine can be for some application.
!!takao2012jun started
      use m_pbindex,only: norbt, ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
     &     max_ibas_tbl,max_l_tbl,max_k_tbl,max_offset_tbl
      use m_hamindex, only: qlat,plat,invgx,
     &     miat,tiat,shtvg,symops,nbas,ngrp
      use m_readqgcou,only: imxc,ngvecc,qtt_,nqnum,ngc,ngveccrev
      implicit none

      public:: RotMPB2
      private
      integer:: lxxa
      real(8),allocatable:: dlmm(:,:,:,:)

      contains
      subroutine RotMPB2(nbloch,ngbb,qin,igx,igxt,ginv
     o ,zrotm) !zcousqr=Rotate_igx(zcousq) igxt=-1 means timereversal case.
      intent(in)::       nbloch,ngbb,qin,igx,igxt,ginv
      intent(out)::
     o  zrotm 
!! --- zrotm(J,J') = <Mbar^k_J| \hat{A}^k_i Mbar^k_J'>. ---
!! See Eq.(51) around in PRB81 125102(2010). 
!! Eq. (51) can be written as
!!   P_IJ= \sum_i T_alpha_i [ zrotm_i^dagger (I,I') P'_I'J' zrom_i(J'J) ],
!! where P'_I'J' is ths sum not yet symmetrized. 

!!
!! Eexactrly speaking, we insert conversion matrix  between Enu basis and M_I basis.
!! It is zcousq(or zzr) in x0kf_v4h.F
!!
!! input qin = q
!! \hat{A}^k_i  is specified by symops(:,:,igx),and igxt (-1 for time-reversal).
!! Note that k= \hat{A}^k_i(k) (S_A^k)
      real(8):: qin(3),ginv(3,3),platt(3,3),qout(3),qu(3),sss
      integer:: igx,igxt,ngbb,iqin,iqout,ngcx,nbloch,nl,i,j
      complex(8):: zrotm(ngbb,ngbb)
      integer,save:: init=1
      logical ::debug=.false.
      real(8):: tolq=1d-8
      if(debug) write(6,*)' rotMPB2:'
      if(init==1) then
        lxxa= 2*max_l_tbl
        nl= lxxa +1 
        allocate(dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa,ngrp))
        call rotdlmm(symops,ngrp,nl,dlmm)
        init=0
      endif   
      zrotm=0d0
      call rotmto2(qin,nbloch,ngbb,
     i  norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl,
     i  symops(:,:,igx),shtvg(:,igx),dlmm(:,:,:,igx),lxxa,miat(:,igx),tiat(:,:,igx),igxt,nbas, 
     o  zrotm(1:nbloch,1:nbloch))
      if(debug) write(6,*)' end of rotmto2 '
      call iqindx2(qin,ginv,qtt_,nqnum, iqin,qu) !iqindx2 is slow. If required, speed it up in a manner as iqindx2_.
      if(debug) write(6,*)' end of iqindx2rotmto2 '
Cstop2rx 2013.08.09 kino      if(sum(abs(qin-qu))>1d-6) stop 'rotMPB:qin is not included in QGcou'
      if(sum(abs(qin-qu))>tolq) call rx( 'rotMPB2:qin is not included in QGcou')
      sss=1d0
      if(igxt==-1) sss=-1d0
      if(debug) write(6,*)' iqindx2 111111111111 sss=',sss,igx,qin
      call iqindx2(sss*matmul(symops(:,:,igx),qin),ginv,qtt_,nqnum, iqout,qout)
      if(debug) write(6,*)' iqindx2 111111111111222222222222'
!!
      ngcx=ngc(iqin)
      if(ngcx/=ngc(iqout).or.ngcx/=ngbb-nbloch) then
         write(6,*)' ngbb qin  =',ngbb,qin
         write(6,*)' iqin iqout=',iqin,iqout
         write(6,*)' qout=',qout
         write(6,*)' ngc(iqin),ngc(iqout)=',ngc(iqin),ngc(iqout)
         write(6,*)' ngcx nbloch ngbb=',ngcx,nbloch,ngbb
         write(6,*)' ngc(iqout)=',ngc(iqout)
Cstop2rx 2013.08.09 kino         if(ngcx/=ngbb-nbloch) stop 'rotMPB2:ngc(iqin)/=ngbb-nbloch'
         if(ngcx/=ngbb-nbloch) call rx( 'rotMPB2:ngc(iqin)/=ngbb-nbloch')
Cstop2rx 2013.08.09 kino         stop 'rotMPB2:ngc(iqin)/=ngc(iqout)'
         call rx( 'rotMPB2:ngc(iqin)/=ngc(iqout)')
      endif   
      platt=transpose(plat)
c      write(6,*)' iqin iqout igx=',iqin,iqout,igx
c      write(6,*)' goto  rotipw ',sum(zcousq(nbloch+1:nbloch+ngcx,1:ngbb))
c      write(6,*)'qin  =',qin,iqin,igx
c      write(6,*)'qout =',qout,iqout
c      write(6,*)'ngcx ngbb=',ngcx,ngbb
c      write(6,*)'platt=',platt
c      write(6,*)' qlat=',qlat
      call rotipw2(qin,qout,ngcx,ngbb,
     &  platt,qlat,symops(1,1,igx),ngvecc(1,1,iqin),ngveccrev(:,:,:,iqout),shtvg(:,igx),igxt,imxc,
     o  zrotm(nbloch+1:nbloch+ngcx,nbloch+1:nbloch+ngcx) )
c      if(igxt<0) then !time-reversal
c        zcousqr = dconjg(zcousqr)
c      endif  
      end subroutine rotMPB2
      end module m_rotMPB2 
