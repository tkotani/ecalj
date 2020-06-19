!! ------------------------------------------------------
!> Return eigenvalus and eigenfunctions for given q and isp.
!! qtt(1:3, nqtt)  :q-vector in full BZ (no symmetry) in QGpsi, QGcou
!! qtti(1:3,nqi) :eivenvalues, eigenvectors are calculated only for them.
!!                  See lmfgw (q-vector with irr flag in QGpsi).
!!
!!* we have to call init_foobar to call readeval, readcphi, readgeig.
!!----------------
      module m_readeigen
      use m_iqindx_qtt,only: Iqindx2_, Init_iqindx_qtt
      use m_hamindex,only:   Readhamindex, ngpmx, nqtt, nqi, nqnum, qtt,iqimap,
     &     iqmap,igmap,shtvg,qlat,symops,plat,invgx,
     &     miat,tiat,dlmm,shtvg,symops,lxxa,nbas
      implicit none
      public:: Readeval, Readgeig, Readcphi, Init_readeigen, Init_readeigen2, Lowesteval,
     &     Init_readeigen_mlw_noeval, Init_readeigen_phi_noeval , Onoff_write_pkm4crpa,
     &     Readcphiw,Readgeigw

      private
      integer:: ifgeigW,ifcphiW,norbtx,imx
      complex(8),allocatable:: geigW(:,:,:,:),cphiW(:,:,:,:)
      integer(4)::nwf_ev
      complex(8),allocatable:: evud_w(:,:,:,:)
      real(8),allocatable:: eval_w(:,:,:)      
      logical:: evalwan=.true.
      real(8),allocatable,private:: evud(:,:,:), ginv(:,:)
      logical,private:: init=.true.,init2=.true.,keepeig
      integer(4),allocatable,private:: ngp(:)
      integer(4),private:: nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   nnnn, nsp,ldim2,mrecg,ifcphi,ifgeig
      complex(8),allocatable,private:: geig(:,:,:,:),cphi(:,:,:,:)
      real(8),private:: leval
      integer(4),private :: nwf,ifcphi_mlw,ifgeig_mlw, nqixx
      complex(8),allocatable,private:: geig_mlw(:,:,:,:),
     &           cphi_mlw(:,:,:,:)
      logical,private:: debug=.false.
      integer,allocatable,private:: ngvecp(:,:,:), ngvecprev(:,:,:,:)
      integer,allocatable,private:: l_tbl(:),k_tbl(:),ibas_tbl(:)
     & ,offset_tbl(:),offset_rev_tbl(:,:,:)
      logical,private:: Wpkm4crpa=.false.

      contains
      subroutine onoff_write_pkm4crpa(lll)
      logical:: lll
      Wpkm4crpa=lll
      end subroutine

      subroutine readgeigW(q,ngp_in,isp, qu,geigen)
      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: geigen(ngp_in,nwf)
      if(init2) call rx( 'readgeig_mlw: modele is not initialized yet')
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
      if(ngp_in < ngp(iq)) then
        write(6,*)'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        call rx( 'readgeig_mlw: ngpmx<ngp(iq)')
      endif
      if(keepeig) then
        geigen(1:ngp(iq),1:nwf) = geigW(1:ngp(iq),1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifgeigW, rec=ikpisp) geigen(1:ngpmx,1:nwf)
      endif
      end subroutine

      subroutine readcphiW(q,ldim2_in,isp,  qu,cphif)
      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: cphif(ldim2_in,nwf)
      if(init2) call rx( 'readcphi_mlw: modele is not initialized yet')
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.
      if(keepeig) then
        cphif(1:ldim2,1:nwf) = cphiW(1:ldim2,1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifcphi_mlw, rec=ikpisp) cphif(1:ldim2,1:nwf)
      endif
      end subroutine


!> Return ev(1:nband) for given q(1:3) and isp 
      subroutine readeval(q,isp, ev)
      integer(4), intent(in) :: isp
      real(8), intent(in) :: q(3)
      real(8), intent(out) :: ev(nband)
      integer(4):: iq,iqindx,i
      real(8):: qu(3)
      if(init) call rx( 'readeigen: modele is not initialized yet')
c      iq = iqindx(q, ginv,qtt,nqtt)
       call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
      if(debug) then
        write(6,*)'iq iqimap(iq)=',iq,iqimap(iq)
        write(6,"('iq iqimap(iq) q=',2i8,3f13.5)")iq,iqimap(iq),q
      endif
      ev(1:nband) = evud(1:nband,iqimap(iq),isp) !iqimap is given in suham.F/gen_hamindex
      if(debug) then
        write(6,"(9f9.4)")ev(1:9)
      endif
      end subroutine

      subroutine findemaxmin_r(emax,emin)
      real(8):: emax,emin,eee
      integer:: ib,iq,isp
      emax=-1d9
      do ib=1,nband
         do iq=1,nqixx
            do isp=1,nsp
               eee= evud(ib,iq,isp)
               if(eee>Emax.and.eee<1d9) Emax=eee
            enddo
         enddo
      enddo
      emin = minval(evud)
      end subroutine

!> Return ev(1:nband) for given q(1:3) and isp 
      subroutine readgeig(q,ngp_in,isp, qu,geigen)
      real(8),intent(in) :: q(3)
      integer,intent(in) :: isp,ngp_in
      real(8),intent(out) :: qu(3)
      complex(8), intent(out) :: geigen(ngp_in,nband)

      integer:: iq,iqindx,ikpisp,napw,iqq,nnn(3),ig,igg,ig2,iqi,igxt
      real(8)   :: ddd(3),platt(3,3),qpg(3),qpgr(3),qtarget(3),qout(3),qin(3)
      complex(8):: geigenr(ngp_in,nband),img=(0d0,1d0),img2pi

      img2pi=2d0*4d0*datan(1d0)*img
      platt=transpose(plat) !this is inverse of qlat
      if(init2) call rx( 'readgeig: modele is not initialized yet')
!! nqtt=nqnum. (aug2012. It was nqtt/2 before aug2012 when we use shorbz). See lm7k/fp/suham.F
c      call iqindx2(q,ginv,qtt,nqtt/2, iq, qu) !qu is used q. q-qu is a G vector.
      call iqindx2_(q, iq, qu) !qu is used q. q-qu is a G vector.
      if(debug) write(6,*)' readgeig:xxx iq=',iq
      iqq=iqmap(iq)
      iqi=iqimap(iq)
      igg=igmap(iq)
      qtarget=qtt(:,iq) ! iqq is mapped to qtarget=qu=qtt(:,iq)
      if(debug) then
        write(6,*)' ddddd ginv=',ginv
        write(6,*)' ddddd readgeig: iq iqq igg=',iq,iqq,igg
        write(6,*)' ddddd q      =',q
        write(6,*)' ddddd qu      =',qu
        write(6,*)' ddddd qtarget=',qtarget
        write(6,*)' ddddd q (iqq)=',qtt(:,iqq)
        write(6,"(a,3i5,3f10.4,2i5)")' ngp(iq) ngp(iqq)=',iq,iqq,igg,q,ngp(iq),ngp(iqq)
      endif
!!  qtt(iqq) is rotated to qtt(iq) by sympos(  ,igg).
      if(ngp(iq)==0) return
      if(ngp(iq)/=ngp(iqq)) then
        write(6,*)' ddddd ginv=',ginv
        write(6,*)' ddddd readgeig: iq iqq igg=',iq,iqq,igg
        write(6,*)' ddddd q      =',q
        write(6,*)' ddddd qu      =',qu
        write(6,*)' ddddd qtarget=',qtarget
        write(6,*)' ddddd q (iqq)=',qtt(:,iqq)
        write(6,"(a,3i5,3f10.4,2i5)")' ngp(iq) ngp(iqq)=',iq,iqq,igg,q,ngp(iq),ngp(iqq)
        call rx( 'readgeig:x2 ngp(iq)/=ngp(iqq)')
      endif  
      if(debug) write(6,*)' ngp=',ngp(iq)
      if(keepeig) then
        geigenr(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iqi,isp)
      else
        ikpisp= isp + nsp*(iqi-1)
        read(ifgeig, rec=ikpisp) geigenr(1:ngpmx,1:nband)
      endif
      if(ngp_in < ngp(iq)) then
        write(6,*)'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
        call rx( 'readgeig: ngpmx<ngp(iq)')
      endif
!! new rotation of IPW. june2012
!!   qinput: qtt(:,iqq)  ---> qtarget: qtt(:,iq) ( G-vector difference from symops*qtt(:,iqq) )
      igxt=1 !not timereversal
      call rotipw(qtt(:,iqq),qtt(:,iq),geigenr(1:ngp(iqq),1:nband),ngp(iqq),nband,
     &  platt,qlat,symops(:,:,igg),ngvecp(:,:,iqq),ngvecprev(:,:,:,iq),shtvg(:,igg),igxt,imx,
     &  geigen(1:ngp(iq),1:nband))
      return
!! Followings are not used now...  june2012
!! We can use the followings instead of rotipw (both should give the same result).
!! ----------
!! Rotation  See mto part of lm7K/subs/rotwv.F
!! NOTE that qrot = matmul(symops(:,:,igg),qtt(:,iqq))
!! rotation wave
c      ddd = shtvg(:,invgx(igg))
cccccccccccccccccccccccccccccccccc
c      write(6,*)' xddddd1= ',-matmul(symops(:,:,igg),ddd)
c      write(6,*)' xddddd2= ',shtvg(:,igg)
c      write(6,*)' xddddd   '
ccccccccccccccccccccccccccccccccccc
c      if(debug) then
c        write(6,*)' ddddd=',igg,ddd
c        write(6,"('ddddd iq =',2i5,3f13.5)") iq,  ngp(iq), qtt(:,iqq)
c        write(6,"('ddddd iqq=',2i5,3f13.5,2x,3f13.5)") iqq, ngp(iq), q,qu
c        write(6,*)' ddddd '
c      endif
      do ig = 1,ngp(iqq)
        qpg = qtt(:,iqq) + matmul( qlat(:,:),ngvecp(:,ig,iqq)) !q+G
        qpgr = matmul(symops(:,:,igg),qpg)             !rotated q+G
        call conv2int( matmul(platt,qpgr-qtarget),nnn,3 )
        ig2 = ngvecprev(nnn(1),nnn(2),nnn(3),iq)
        geigen(ig2,:)= geigenr(ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
      enddo
      end subroutine readgeig

!!
      subroutine readcphi(q,ldim2_in,isp,  qu,cphif)
!!-- return mto part of eigenfunction for given q(1:3) and isp
      real(8), intent(in) :: q(3)
      integer, intent(in)  :: ldim2_in
      integer, intent(in)  :: isp
      real(8), intent(out)  :: qu(3)
      complex(8), intent(out)  :: cphif(ldim2_in,nband)

      integer:: iq,iqindx,ikpisp,iqq,iorb,ibaso,ibas,k,l,
     & ini1,ini2,iend1,iend2,igg,ig,iqi,i,igxt
      real(8)   :: qrot(3) ,qout(3)
      complex(8):: phase,cphifr(ldim2_in,nband),phaseatom !takao 1->*->nband
!      complex(8),parameter:: img=(0d0,1d0),img2pi=2d0*4d0*datan(1d0)*img ! MIZUHO-IR
      complex(8),parameter:: img=(0d0,1d0) ! MIZUHO-IR

      complex(8):: img2pi ! MIZUHO-IR
      img2pi = 2d0*4d0*datan(1d0)*img ! MIZUHO-IR
Cstop2rx 2013.08.09 kino      if(init2) stop 'readcphi: modele is not initialized yet'
      if(init2) call rx( 'readcphi: modele is not initialized yet')
!! nqtt=nqnum. (aug2012. It was nqtt/2 before aug2012 when we use shorbz). See lm7k/fp/suham.F
c      call iqindx2(q,ginv,qtt,nqtt/2, iq, qu) !qu is used q. q-qu is a G vector.
      call iqindx2_(q, iq, qu) !qu is used q. q-qu= G vectors. qu=qtt(:,iq)
!!  For given iq, iqq = iqmap(iq)
!! qtt(:,iq)= matmul(sympos(  ,igg),qtt(:,iqq))  and  qtt(:,iqq) = qtti(:,iqi)
      igg=igmap(iq)
      iqq=iqmap(iq)  ! mapped from qtt(:,iqq) to qtt(:,iq);
      ! qtt(:,iq)=matmul(sym(igg),qtt(:,iqq))+some Gvector(see iqindx2 above)
      iqi=iqimap(iq) ! iqi is index for irr.=1 (cphi calculated. See lqg4gw and lm7K/subs/suham.F:L200).
      ! qtt(:,iqq) = qtti(:,iqi) is satisfied.
      ! qtti(:,iqi) are only when eigenfunctions are calculated).
      if(keepeig) then
        cphifr(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iqi,isp)
      else
        ikpisp= isp + nsp*(iqi-1)
        read(ifcphi, rec=ikpisp) cphifr(1:ldim2,1:nband)
      endif
      if(debug) write(6,"('readcphi:: xxx sum of cphifr=',3i4,4d23.16)")ldim2_in,ldim2,norbtx,
     &                sum(cphifr(1:ldim2,1:nband)),sum(abs(cphifr(1:ldim2,1:nband)))
!!
      igxt=1 !not timereversal
      call rotmto(qtt(:,iqq),cphifr,ldim2,nband,
     i  norbtx,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,
     i  maxval(ibas_tbl),maxval(l_tbl),maxval(k_tbl),
     i  ginv,symops(1,1,igg),shtvg(:,igg),dlmm(:,:,:,igg),lxxa,miat(:,igg),tiat(:,:,igg),igxt,nbas, 
     o  cphif)
      return

c$$$!! Followings are not used now...  june2012
c$$$!! We can use the followings instead of rotmto (both should give the same result).
c$$$!! ----------
c$$$!! Rotation  See mto part of lm7K/subs/rotwv.F
c$$$!! NOTE that qrot = matmul(symops(:,:,igg),qtt(:,iqq))
c$$$!! (then MTO itself is rotated. a little diffent from here).
c$$$!! qrot= symops(igg) qtt
c$$$      qrot=matmul(symops(:,:,igg),qtt(:,iqq)) !qrot-qtt(:,iq) is some G vector.
c$$$      phaseq   = exp( -img2pi*sum( qrot*shtvg(:,igg)) )
c$$$cddd
c$$$c      phaseq   = exp( img2pi*sum( qtt(:,iqq)*shtvg(:,invgx(igg)) ) )
c$$$cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c$$$cccccccccccccccccccccccccccccccccc
c$$$c      write(6,*)' xddddd1= ', sum( qtt(:,iqq)*shtvg(:,invgx(igg)) )
c$$$c      write(6,*)' xddddd2= ',-sum( qrot*shtvg(:,igg))
c$$$c      write(6,*)' xddddd diff=   ',sum( qtt(:,iqq)*shtvg(:,invgx(igg)) )+sum( qrot*shtvg(:,igg))
c$$$ccccccccccccccccccccccccccccccccccc
c$$$      ibaso=-999
c$$$      do iorb=1,norbt !orbital-blocks are specified by ibas, l, and k.
c$$$        !Index of Hamiltonian is devided into these blocks.
c$$$        ibas = ibas_tbl(iorb)
c$$$        if(ibas/=ibaso) then
c$$$          phase = exp( -img2pi*sum(qrot*(tiat(:,ibas,igg)-shtvg(:,igg))) )
c$$$        endif
c$$$        ibaso=ibas
c$$$        l   = l_tbl(iorb)
c$$$        k   = k_tbl(iorb)
c$$$        ini1 = offset_tbl(iorb)+1
c$$$        iend1 = ini1+2*l
c$$$        ini2 = offset_rev_tbl(miat(ibas,igg),l,k)+1
c$$$        iend2 = ini2+2*l
c$$$        cphif(ini2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),cphifr(ini1:iend1,:))
c$$$     &                       *phase*phaseq
c$$$        if(debug) then
c$$$          write(6,*)ibaso,ibas,phase
c$$$          write(6,*)ini1,iend1,ini2,iend2
c$$$          write(6,*)'sumd  dlmm =',igg,l,sum(abs(dlmm(-l:l,-l:l,l,igg)))
c$$$          write(6,*)'sumcr cphifr=',l,k,igg,sum(abs(cphifr(ini1:iend1,:)))
c$$$          write(6,*)'sumcf cphif=',l,k,igg,sum(abs(cphif(ini2:iend2,:)))
c$$$        endif
c$$$      enddo
c$$$ccccccccccccccccccccccccccccccccc
c$$$c       write(6,"('readcphi:: xxxx sum of cphifr=',2i4,4d23.16)")ldim2_in,ldim2,
c$$$c     &                sum(cphif(1:ldim2,1:nband)),sum(abs(cphif(1:ldim2,1:nband)))
c$$$c       stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxx takao test xxxxxxxxxxxxxxxxx'
c$$$cccccccccccccccccccccccccccccccccc
c$$$ccccccccccccccccccccccc
c$$$c      write(6,"('rrcc0:',i3,3f8.4)")iq,qu
c$$$c      do i=1,1
c$$$c      write(6,"('rrccb:',5(2x,2f8.4))")
c$$$c     &  cphif(14,i)/cphif(1,i) !*abs(cphif(1,i))
c$$$c      enddo
c$$$c      cphif(5:13, :)=0d0
c$$$c      cphif(14,:)=0d0
c$$$c      cphif(18:26,:)=0d0
c$$$
c$$$c      cphif(1:8,:)=0d0
c$$$c      cphif(19:26,:)=0d0
c$$$c      cphif(1,:)=0d0
c$$$c      cphif(19,:)=0d0
c$$$cccccccccccccccccccccccccccccccccccc
      end subroutine readcphi


c$$$c------------------------------
c$$$      subroutine readgeig(q,ngp_in,isp, qu,geigen)
c$$$! return ev(1:nband) for given q(1:3) and isp
c$$$      implicit none
c$$$      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
c$$$      real(8)   :: q(3),qu(3)
c$$$      complex(8):: geigen(ngp_in,nband)
c$$$      if(init2) then
c$$$        stop 'readgeig: modele is not initialized yet'
c$$$      endif
c$$$      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.
c$$$      if(ngp_in < ngp(iq)) then
c$$$        write(6,*)'readgeig: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
c$$$        stop 'readgeig: ngpmx<ngp(iq)'
c$$$      endif
c$$$! for ngp=0 Aug2005
c$$$      if(ngp(iq)==0) return
c$$$      if(keepeig) then
c$$$        geigen(1:ngp(iq),1:nband) = geig(1:ngp(iq),1:nband,iq,isp)
c$$$      else
c$$$        ikpisp= isp + nsp*(iq-1)
c$$$        read(ifgeig, rec=ikpisp) geigen(1:ngpmx,1:nband)
c$$$      endif
c$$$      end subroutine readeig
c$$$c------------------------------
c$$$      subroutine readcphi(q,ldim2_in,isp,  qu,cphif)
c$$$! return ev(1:nband) for given q(1:3) and isp
c$$$      implicit none
c$$$      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
c$$$      real(8)   :: q(3),qu(3)
c$$$      complex(8):: cphif(ldim2_in,nband) !takao 1->*->nband
c$$$      if(init2) then
c$$$        stop 'readcphi: modele is not initialized yet'
c$$$      endif
c$$$      call iqindx2(q,ginv,qtt,nqtt, iq, qu) !qu is used q.  q-qu= G vectors.
c$$$      if(keepeig) then
c$$$        cphif(1:ldim2,1:nband) = cphi(1:ldim2,1:nband,iq,isp)
c$$$      else
c$$$        ikpisp= isp + nsp*(iq-1)
c$$$        read(ifcphi, rec=ikpisp) cphif(1:ldim2,1:nband)
c$$$      endif
c$$$      end subroutine readcphi

c------------------------------
      subroutine readgeig_mlw(q,ngp_in,isp, qu,geigen)
! return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ngp_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: geigen(ngp_in,nwf)
      if(init2) then
Cstop2rx 2013.08.09 kino        stop 'readgeig_mlw: modele is not initialized yet'
        call rx( 'readgeig_mlw: modele is not initialized yet')
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      write(6,*)' readgeig: q=',q
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.

      if(ngp_in < ngp(iq)) then
        write(6,*)'readgeig_mlw: ngpmx<ngp(iq)',iq,ngpmx,ngp(iq),q
Cstop2rx 2013.08.09 kino        stop 'readgeig_mlw: ngpmx<ngp(iq)'
        call rx( 'readgeig_mlw: ngpmx<ngp(iq)')
      endif
      if(keepeig) then
        geigen(1:ngp(iq),1:nwf) = geig_mlw(1:ngp(iq),1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifgeig_mlw, rec=ikpisp) geigen(1:ngpmx,1:nwf)
      endif
      end subroutine

c------------------------------
      subroutine readcphi_mlw(q,ldim2_in,isp,  qu,cphif)
!!-- return ev(1:nband) for given q(1:3) and isp
      implicit none
      integer(4):: isp,iq,iqindx,ldim2_in,ikpisp
      real(8)   :: q(3),qu(3)
      complex(8):: cphif(ldim2_in,nwf)
      if(init2) then
Cstop2rx 2013.08.09 kino        stop 'readcphi_mlw: modele is not initialized yet'
        call rx( 'readcphi_mlw: modele is not initialized yet')
      endif
c      iq = iqindx(q, ginv,qtt,nqtt)
c      write(6,*)' readgcphi_mlw: q=',q
      call iqindx2_(q, iq, qu) !qu is used q.  q-qu= G vectors.

      if(keepeig) then
        cphif(1:ldim2,1:nwf) = cphi(1:ldim2,1:nwf,iq,isp)
      else
        ikpisp= isp + nsp*(iq-1)
        read(ifcphi_mlw, rec=ikpisp) cphif(1:ldim2,1:nwf)
      endif
      end subroutine

ccccccccccccccccccccccccccc
      subroutine init_readeigen(ginv_in,nsp_in,nband_in,mrece_in)
!!-- initialization. Save QpGpsi EVU EVD to arrays.--
      integer:: iq,is,ifiqg,nnnn,ikp,nsp_in, 
     & nband_in,mrece_in,iopen,iclose,isx,mrecb_in,ik,ib,verbose,ifile_handle
      integer:: ifev,nband_ev, nqi_, nsp_ev ,ngpmx_ ,nqnum_
      real(8):: ginv_in(3,3)
      real(8):: QpGcut_psi
      real(8),allocatable:: qtt_(:,:),qtti_(:,:)

      write(6,*) 'init_readeigen:'
      allocate(ginv(3,3))
      ginv  = ginv_in
      nsp   = nsp_in
      if(nsp<0.or.nsp>2) call rx( 'init_reaeigen:nsp wrong')
      nband = nband_in
      mrece = mrece_in
c      write(6,*) 'init_readeigen:1'

!!--- store requied index
      call readhamindex()
c      write(*,*)ngrp,nbas,kxx,lxx
      write(*,*)'nqi=',nqi!,nqnum
c      write(*,*)'imx=',imx,ngpmx
c      do iq=1,nqtt
c        write(*,*) qtt(:,iq)
c      enddo

      call init_iqindx_qtt(ginv)

!!--- Readin QGpsi --- 
      ifiqg  = ifile_handle()
      open(ifiqg ,file='QGpsi',form='unformatted')
      read(ifiqg) nqnum_ , ngpmx_, QpGcut_psi, nnnn,nqi_ ,imx 
      write(6,*)'read(ifiqg)', nqnum , ngpmx_, QpGcut_psi, nnnn,nqi
      if(nqi  /=  nqi_) call rx( 'init_readeigen:nqi/=nqi_ 11111')
      if(nqnum/=  nqnum_) call rx( 'init_readeigen:nqnum/=nqnum_ 11111')
      if(ngpmx_/=ngpmx) call rx('ngpmx error: 1111111 readeigen')
      allocate( qtt_(3,nqnum),ngp(nqnum) )
      allocate( ngvecp(3,ngpmx,nqnum))
      allocate( ngvecprev(-imx:imx,-imx:imx,-imx:imx,nqnum) )
      do ikp = 1,nqnum
        read (ifiqg) qtt_(:,ikp), ngp(ikp)
        read (ifiqg) ngvecp(1:3, 1:ngp(ikp),ikp),ngvecprev(-imx:imx,-imx:imx,-imx:imx,ikp)
      enddo
      close(ifiqg)
      deallocate(qtt_)

!!--- Read eigenvalue ---
      ifev= iopen('EValue',0,-1,0)    !feb2012takao
      read(ifev) nband_ev, nqi_, nsp_ev
      write(6,*)'read EValue: nband_ev,nqi,nsp_ev', nband_ev, nqi_, nsp_ev
      if(nband_ev/=nband) call rx( 'init_readeigen:nband_ev/=nband')
      if(nsp_ev  /=  nsp) call rx( 'init_readeigen:nsp_ev/=nsp')
      if(nqi  /=  nqi_) call rx( 'init_readeigen:nqi/=nqi_')
      nqixx=nqi
      allocate(evud(nband,nqi_,nsp),qtti_(3,nqi_))
      read(ifev) qtti_(1:3,1:nqi_) !qtti
      read(ifev) evud(1:nband, 1:nqi, 1:nsp )
      close(ifev)
      if(debug) then
        do is= 1,nsp
          do ik= 1,nqi
            do ib= 1,nband
              write(6,"('ib ik e=',3i5,f13.5,2x,3f9.4)") ib,ik,is,evud(ib,ik,is), qtti_(1:3,ik)
            enddo
          enddo
        enddo
        if(debug) write(6,*)'init_readeigen:end'
      endif

c comment out mar2012takao
c      ifev(1) = iopen('EVU', 0,0,mrece)
c      if(nsp==2) ifev(2) = iopen('EVD', 0,0,mrece)
C$$$      do is= 1,nsp
C$$$        call rwdd (ifev(is),   !read eigenvalues
C$$$     d              nband,nqi,
C$$$     o              evud(1,1,is))
C$$$c        write(6,*) 'evud is=',is,minval(evud(:,:,is))
C$$$        if(debug) then
C$$$          do ik=1,nqtt
C$$$            do ib=1,nband
C$$$              write(6,"('ib ik e=',2i5,f13.5)") ib,ik,evud(ib,ik,is)
C$$$            enddo
C$$$          enddo
C$$$        endif
C$$$      enddo
C$$$c      stop
C$$$      isx = iclose('EVU')
C$$$      if(nsp==2) isx = iclose('EVD')
      leval= minval(evud)
      init=.false.
      end subroutine init_readeigen

      real(8) function lowesteval()
      lowesteval=leval
      end function
c-------------------------------------------------------------
      subroutine init_readeigen2(mrecb_in,ldim2_in,mrecg_in)
c      use m_hamindex,only: nqi
c this should be called after init_readgeigen
      implicit none
      integer(4):: iq,is,ifiqg,nnnn,ikp,ldim2_in, 
     & iopen,iclose,isx,mrecb_in,ikpisp,mrecg_in,verbose,ifoc,
     & i1,i2,i3,i4,i5,iorb,iorbold
      logical :: keepeigen

      mrecb = mrecb_in
      ldim2 = ldim2_in
      mrecg = mrecg_in
      keepeig = keepeigen()

      init2=.false.

!! === readin @MNLA_CPHI for rotation of MTO part of eigenfunction cphi ===
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      read(ifoc,*)
      norbtx=0
c      write(*,*) ' goto readin @MNLA_CPHI: norbtx=',norbtx
      do
        read(ifoc,*,end=106) i1,i2,i3,i4,i5,iorb
c        write(6,*) i1,i2,i3,i4,i5,iorb
        if(iorb>norbtx) norbtx=iorb
      enddo
 106  continue
      write(6,*) ' end of readin @MNLA_CPHI: norbtx=',norbtx
      rewind ifoc
      read(ifoc,*)
      allocate(l_tbl(norbtx),k_tbl(norbtx),ibas_tbl(norbtx),offset_tbl(norbtx))
      iorbold=0

      do
        read(ifoc,*,end=107)i1,i2,i3,i4,i5,iorb
        if(iorbold/=iorb) then
          k_tbl(iorb)=i2
          l_tbl(iorb)=i3
          ibas_tbl(iorb)=i4
          offset_tbl(iorb)=i5-1
          iorbold=iorb
        endif
      enddo
 107  continue

      allocate(offset_rev_tbl(maxval(ibas_tbl),0:maxval(l_tbl),maxval(k_tbl)))
      offset_rev_tbl=-99999
      do iorb=1,norbtx
        offset_rev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))= offset_tbl(iorb)
      enddo
      write(6,*)'init_readeigen2: table for rotation of cphi'
c      do iorb=1,norbtx
c        write(*,*)ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb),
c     &   offset_tbl(iorb),
c     &  offset_rev_tbl(ibas_tbl(iorb),l_tbl(iorb),k_tbl(iorb))
c      enddo
      ifoc=iclose('@MNLA_CPHI')
c      stop 'init2 xxxxxxxxxxxxx'
c --- Readin CPHI and GEIG
      if(Keepeig     ) write(6,*)' KeepEigen=T; readin geig and cphi into m_readeigen'
      if(.not.Keepeig) write(6,*)' KeepEigen=F; not keep geig and cphi in m_readeigen'
      ifgeig = iopen('GEIG',0,0,mrecg)
      ifcphi = iopen('CPHI',0,0,mrecb) ! Coeeficients of eigenfunction
c      write(6,*)' ifgeig ifcphi=',ifgeig,ifcphi
c      write(6,*)' ldim2,nband,nqi,nsp',ldim2,ngpmx,nband,nqi,nsp
      if(.not.keepeig) return
c      ifgeig = iopen('GEIG',0,-1,0)
c      ifcphi = iopen('CPHI',0,-1,0) ! Coeeficients of eigenfunction
      allocate(geig(ngpmx,nband,nqi,nsp))
      allocate(cphi(ldim2,nband,nqi,nsp))
      write(6,*)' size geig=',ngpmx,nband,nqi,nsp,ldim2,size(geig),size(cphi)
      do ikp= 1,nqi
        do is= 1,nsp
c        read(ifgeig) geig(1:ngpmx,1:nband,ikp,is) !geig(1:ngp(ikp),1:nband,ikp,1:is)
c        read(ifcphi) cphi(1:ldim2,1:nband,ikp,is)
          ikpisp= is + nsp*(ikp-1)
          if(ngpmx/=0) read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is) !add ngpmx/=0 Aug2005
ccccccccccccccccccccccc
c          if(verbose()>40)
c     &       write(6,'("Readin geig sumcheck ",i5,3f9.4,"  ",i5,2d13.6)')
c     &      ikp,qtti(:,ikp),ngp(ikp),sum(geig(1:ngp(ikp),1:nband,ikp,is))
cccccccccccccccccccccc
c        read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nband,ikp,is)
          read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nband,ikp,is)
        enddo
      enddo
      write(6,*)'eeeee end of init_readeigen2 eeeee for keepeig=T'
      isx= iclose('CPHI')
      isx= iclose('GEIG')
      end subroutine




c-------------------------------------------------------------
      subroutine init_readeigen_mlw(nwf_o)
c allocate cphi_mlw and geig_mlw
c for hx0fp0_mlw (ixc=112 mode)
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,
     & iqbz,iqbz2,nwf2,ib,iwf,iwf2,iko_ix,iko_fx,in,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      write(6,*)' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
Cstop2rx 2013.08.09 kino          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) call rx( "init_readeigen_mlw: nqbz2 error")
        if (iko_ix2.ne.iko_ix)  call rx( "init_readeigen_mlw: iko_ix2 error")
        if (iko_fx2.ne.iko_fx)  call rx( "init_readeigen_mlw: iko_fx2 error")
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
            if (iq0i2 .ne. iq0i) call rx( 'init_readeigen_mlw: iq0i error')
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c allocate geig_mlw and cphi_mlw
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqtt
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
              enddo
            enddo
          else
c
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
          endif
        enddo
      enddo
      deallocate(dnk,uum)


      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig_mlw(ngpmx,nwf,nqtt,nsp))
        allocate(cphi_mlw(ldim2,nwf,nqtt,nsp))
        geig_mlw = 0d0
        cphi_mlw = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig_mlw(:,iwf,ikp,is) = geig_mlw(:,iwf,ikp,is) + 
     &           geig(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi_mlw(:,iwf,ikp,is) = cphi_mlw(:,iwf,ikp,is) + 
     &           cphi(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

          enddo
        enddo


      else
        ifcphi_mlw = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_mlw = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_mlw, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_mlw, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
      endif
      deallocate(cbwf)
c replace dimensions
      nwf_o = nwf
c ---
      end subroutine

c$$$c-------------------------------------------------------------
c$$$      subroutine init_readeigen_mlw_noeval1D(
c$$$     &                 nwf_o,nband_o,mrecb_o,mrecg_o)
c$$$      use m_hamindex
c$$$
c$$$c replace cphi and geig
c$$$c for hwmat
c$$$c this should be called after init_readgeigen2
c$$$      implicit none
c$$$      integer(4):: iq,is,ifiqg,ikp,
c$$$     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
c$$$     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
c$$$     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
c$$$     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2,
c$$$     & iset,nset,iset2,nset2,ini
c$$$      real(8):: q(3),rnorm,cnorm
c$$$      real(8),allocatable :: eval(:,:,:)
c$$$      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
c$$$     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
c$$$     &                          geig3(:,:),cphi3(:,:),
c$$$     &                          geig4(:,:),cphi4(:,:),
c$$$     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)
c$$$
c$$$
c$$$      logical :: keepeigen
c$$$      keepeig = keepeigen()
c$$$
c$$$      write(6,*)' init_readeigen_mlw_noeval1D'
c$$$
c$$$
c$$$c --- Readin MLWU/D and UU1dU/D
c$$$      do is = 1,nsp
c$$$
c$$$c fileopen
c$$$        if (is.eq.1) then
c$$$          ifmlw  = iopen('MLWU',0,0,0)
c$$$          ifuu   = iopen('UU1dU',0,0,0)
c$$$        else
c$$$          ifmlw  = iopen('MLWD',0,0,0)
c$$$          ifuu   = iopen('UU1dD',0,0,0)
c$$$        endif
c$$$
c$$$c read UU1dU/D
c$$$        read(ifuu)nqbz,nq0i,nwf,nset
c$$$
c$$$        iko_ix = 1
c$$$        iko_fx = nwf
c$$$Cstop2rx 2013.08.09 kino        if (nqtt.ne.nset*nqbz) stop
c$$$        if (nqtt.ne.nset*nqbz) call rx( ''//
c$$$     >      "init_readeigen_nowval1D: nqtt error")
c$$$
c$$$c nqbz mesh-points
c$$$
c$$$        if (is.eq.1) allocate(cbwf(nwf,nwf,nqtt,nsp))
c$$$
c$$$        do iset = 1,nset
c$$$          ini = (iset-1)*nqbz
c$$$          read(ifmlw)nqbz2,nwf2,iset2,nset2
c$$$Cstop2rx 2013.08.09 kino          if (nqbz.ne.nqbz2) stop
c$$$          if (nqbz.ne.nqbz2) call rx( ''//
c$$$     >         "init_readeigen_nowval1D: nqbz error")
c$$$Cstop2rx 2013.08.09 kino          if (nwf.ne.nwf2) stop
c$$$          if (nwf.ne.nwf2) call rx( ''//
c$$$     >         "init_readeigen_nowval1D: nwf error")
c$$$Cstop2rx 2013.08.09 kino          if (iset.ne.iset2) stop
c$$$          if (iset.ne.iset2) call rx( ''//
c$$$     >         "init_readeigen_nowval1D: iset error")
c$$$Cstop2rx 2013.08.09 kino          if (nset.ne.nset2) stop
c$$$          if (nset.ne.nset2) call rx( ''//
c$$$     >         "init_readeigen_nowval1D: nset error")
c$$$          do iqbz = 1,nqbz
c$$$            read(ifmlw)iqbz2,q(1:3)
c$$$            if (iqbz2 .ne. iqbz)
c$$$Cstop2rx 2013.08.09 kino     &           stop 'init_readeigen_mlw_noeval1D: iqbz error'
c$$$     &           call rx( 'init_readeigen_mlw_noeval1D: iqbz error')
c$$$            read(ifmlw)cbwf(1:nwf,1:nwf,ini+iqbz,is)
c$$$          enddo
c$$$        enddo
c$$$
c$$$c fileclose
c$$$        if (is.eq.1) then
c$$$          ifmlw  = iclose('MLWU')
c$$$          ifuu   = iclose('UU1dU')
c$$$        else
c$$$          ifmlw  = iclose('MLWD')
c$$$          ifuu   = iclose('UU1dD')
c$$$        endif
c$$$
c$$$c end of is-loop
c$$$      enddo
c$$$
c$$$
c$$$c replace evud
c$$$      deallocate(evud)
c$$$      allocate(evud(nwf,nqtt,nsp))
c$$$      evud = 0d0
c$$$
c$$$
c$$$      mrecb_o = mrecb * nwf / nband
c$$$      mrecg_o = mrecg * nwf / nband
c$$$
c$$$
c$$$      if(keepeig) then
c$$$
c$$$        allocate(geig2(ngpmx,nband,nqtt,nsp))
c$$$        allocate(cphi2(ldim2,nband,nqtt,nsp))
c$$$        geig2 = geig
c$$$        cphi2 = cphi
c$$$        deallocate(geig,cphi)
c$$$        allocate(geig(ngpmx,nwf,nqtt,nsp))
c$$$        allocate(cphi(ldim2,nwf,nqtt,nsp))
c$$$        geig = 0d0
c$$$        cphi = 0d0
c$$$        do ikp= 1,nqtt
c$$$          do is= 1,nsp
c$$$            do iwf= 1,nwf
c$$$              do ib= iko_ix,iko_fx
c$$$                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
c$$$     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
c$$$                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
c$$$     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
c$$$              enddo
c$$$            enddo
c$$$
c$$$          enddo
c$$$        enddo
c$$$
c$$$
c$$$        deallocate(geig2,cphi2)
c$$$
c$$$      else
c$$$        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
c$$$        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)
c$$$
c$$$        allocate(geig3(ngpmx,nwf))
c$$$        allocate(cphi3(ldim2,nwf))
c$$$        allocate(geig4(ngpmx,nband))
c$$$        allocate(cphi4(ldim2,nband))
c$$$        do ikp= 1,nqtt
c$$$          do is= 1,nsp
c$$$            ikpisp= is + nsp*(ikp-1)
c$$$            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
c$$$            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)
c$$$
c$$$            geig3 = 0d0
c$$$            cphi3 = 0d0
c$$$            do iwf= 1,nwf
c$$$              do ib= iko_ix,iko_fx
c$$$                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
c$$$                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
c$$$              enddo
c$$$            enddo
c$$$
c$$$            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
c$$$            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
c$$$          enddo
c$$$        enddo
c$$$
c$$$        deallocate(geig3,geig4,cphi3,cphi4)
c$$$        isx= iclose('CPHI')
c$$$        isx= iclose('GEIG')
c$$$        isx= iclose('CPHI.mlw')
c$$$        isx= iclose('GEIG.mlw')
c$$$        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
c$$$        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
c$$$      endif
c$$$
c$$$
c$$$      deallocate(cbwf)
c$$$
c$$$
c$$$c replace dimensions
c$$$      nband = nwf
c$$$      nwf_o = nwf
c$$$      nband_o = nband
c$$$      mrecb = mrecb_o
c$$$      mrecg = mrecg_o
c$$$
c$$$c ---
c$$$
c$$$      end subroutine
c$$$


c-------------------------------------------------------------
      subroutine init_readeigen_mlw_noeval(nwf_o,nband_o,mrecb_o,mrecg_o)
c      use m_hamindex,only: nqi,nqnum,qtti

c replace cphi and geig
c for hwmat
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm,qu(3),tolq=1d-8
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:),cphi2(:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:)
      logical :: keepeigen
      integer:: ikpx,ifi,ifile_handle
      character*(8):: fname

      keepeig = keepeigen()
      write(6,*)' init_readeigen_mlw_noeval'

c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
Cstop2rx 2013.08.09 kino          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1)  allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) call rx( "init_readeigen_mlw: nqbz2 error")
        if (iko_ix2.ne.iko_ix) call rx( "init_readeigen_mlw: iko_ix2 error")
        if (iko_fx2.ne.iko_fx) call rx( "init_readeigen_mlw: iko_fx2 error")
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
Cstop2rx 2013.08.09 kino            if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
            if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
Cstop2rx 2013.08.09 kino            if (iq0i2 .ne. iq0i) stop 'init_readeigen_mlw: iq0i error'
            if (iq0i2 .ne. iq0i) call rx( 'init_readeigen_mlw: iq0i error')
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo


c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqi,nsp)) !nqtt
      evud = 0d0

c replace geig and cphi
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp)) 
      cbwf = 0d0
!!
      if(Wpkm4crpa) then
         fname='pkm4crpa'
         ifi=ifile_handle()
         open(ifi,file=fname,form='formatted',status='unknown')
         write(ifi,"('== p_km^alpha in PRB83,121101 ! weight in l-subspace ==')")
         write(ifi,"('( = c^sigma_km in book of 45th IFFK by Ersoy)')")
         write(ifi,"(8i8)") nqtt,nwf,nsp,iko_ix,iko_fx
         write(ifi,"('       |pkm|**2          ib      iq     is       q(1:3)')")
      endif

      do ikp = 1,nqtt
        iqbz = mod(ikp,nqbz) 
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = dnk(ib,iwf,iqbz,is)
              enddo
            enddo
          else
c
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                cbwf(ib,iwf,ikp,is) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
          endif
!! --- write pkm4crpa
          if(Wpkm4crpa) then
             do ib = iko_ix,iko_fx
                write(ifi,"(f19.15, 3i8, 3f13.6 )")
     &            sum(abs(cbwf(ib,1:nwf,ikp,is))**2),ib, ikp, is, qtt(1:3,ikp)
             enddo
          endif

ctm norm check
c         do iwf  = 1,nwf
c         do iwf2 = 1,nwf
c           rnorm = 0d0
c           cnorm = 0d0
c           do ib = iko_ix,iko_fx
c              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c           enddo
c           do ib = 1,nwf
c              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c           enddo
c           if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c           write(7700,"(4i5,2f12.6)")is,ikp,iwf,iwf2,rnorm,cnorm
c         enddo
c         enddo
c         write(7300,"(5i5)")is,ikp,iko_ix,iko_fx,nwf
c         write(7300,*)cbwf(:,:,ikp,is)
        enddo
      enddo
      if(Wpkm4crpa) then
         close(ifi)
      endif

      deallocate(dnk,uum)
      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband

      if(keepeig) then
        write(6,*)' xxx nband=',nband
        allocate(geig2(ngpmx,nband))  !nqtt -->nqi
        allocate(cphi2(ldim2,nband))
c$$$        allocate(geig2(ngpmx,nband,nqtt,nsp))  !nqtt -->nqi
c$$$        allocate(cphi2(ldim2,nband,nqtt,nsp))
c$$$        do ikp= 1,nqtt !nqtt --> nqi
c$$$        do is= 1,nsp
c$$$           call readgeig(qtt(:,ikp),ngpmx,is, qu,geig2(1,1,ikp,is))
c$$$           if(sum(abs(qtt(:,ikp)-qu))>1d-8) call rx('init_readeigen_mlw_noeval 1111')
c$$$           call readcphi(qtt(:,ikp),ldim2,is, qu,cphi2(1,1,ikp,is))
c$$$           if(sum(abs(qtt(:,ikp)-qu))>1d-8) call rx('init_readeigen_mlw_noeval 2222')
c$$$           write(6,"('xxx ikp=',i5,3f9.5,2d13.5)")ikp,qtt(:,ikp)
c$$$     &     ,sum(abs(geig2(:,:,ikp,is))),sum(abs(cphi2(:,:,ikp,is)))
c$$$        enddo
c$$$        enddo
c$$$        write(6,*)' xxxx size geig=',size(geig),size(cphi)  
c$$$        write(6,*)' xxxx size geig2=',size(geig2),size(cphi2)
c$$$        write(6,*)' size geig2=',ngpmx,nband,nqi,nsp,ldim2
ccccccccccccccccccccccccccccccccccccccccccc
c        deallocate(geig,cphi)
ccccccccccccccccccccccccccccccccccccccccccc
        allocate(geigW(ngpmx,nwf,nqtt,nsp))
        allocate(cphiW(ldim2,nwf,nqtt,nsp))
        geigW = 0d0
        cphiW = 0d0
        do ikp= 1,nqtt ! nqi
        do is= 1,nsp
           if(debug) write(6,"(' ikp=',i5,3f10.5)") ikp,qtt(:,ikp)
           call readgeig(qtt(:,ikp),ngpmx,is, qu,geig2)
           if(debug)print *,'qqqqqq1',qu
           if(debug)print *,'qqqqqq2',qtt(:,ikp)
           if(sum(abs(qtt(:,ikp)-qu))>tolq) call rx('init_readeigen_mlw_noeval 1111')
           call readcphi(qtt(:,ikp),ldim2,is, qu,cphi2)
           if(sum(abs(qtt(:,ikp)-qu))>tolq) call rx('init_readeigen_mlw_noeval 2222')
c           call readgeig(qtti(:,ikp),ngpmx,is, qu,geig2)
c           if(sum(abs(qtti(:,ikp)-qu))>1d-8) call rx('init_readeigen_mlw_noeval 1111')
c           call readcphi(qtti(:,ikp),ldim2,is, qu,cphi2)
c           if(sum(abs(qtti(:,ikp)-qu))>1d-8) call rx('init_readeigen_mlw_noeval 2222')
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
cc                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
cc     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikpx,is)
cc                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
cc     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikpx,is)
c                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
c                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                geigW(:,iwf,ikp,is) = geigW(:,iwf,ikp,is) + geig2(:,ib)*cbwf(ib,iwf,ikp,is)
                cphiW(:,iwf,ikp,is) = cphiW(:,iwf,ikp,is) + cphi2(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo
check write
c            do iwf  = 1,nwf
c            do iwf2 = 1,nwf
c               rnorm = 0d0
c               cnorm = 0d0
c               do ib = 1,ldim2
c                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*
c     &                                   cphi(ib,iwf2,ikp,is))
c               enddo
c               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
c            enddo
c            enddo
c            write(7500,*)ikp,ldim2,nwf
c            write(7500,*)cphi(:,:,ikp,is)
          enddo
        enddo
        deallocate(geig2,cphi2,geig,cphi)
      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)
        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)
            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo
            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo
        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeigW = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphiW = iopen('CPHI.mlw',0,0,mrecb_o)
      endif
      deallocate(cbwf)
c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o
      end subroutine init_readeigen_mlw_noeval


c-------------------------------------------------------------
      subroutine init_readeigen_phi_noeval(nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi and geig
c for hwmat_phi
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbz2,nqbze,nqbze2,
     & iqbz,iqbz2,nwf2,nsp2,nlmto2,ngpmx2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,ifdim
      real(8):: q(3),rnorm,cnorm
      complex(8),allocatable :: geig2(:,:),cphi2(:,:)

      logical :: keepeigen
      keepeig = keepeigen()

      write(6,*)' init_readeigen_phi_noeval'

      ifdim = iopen('PHIG.d',1,0,0)
      read(ifdim,*) nsp2,nqbz2,nwf2,nlmto2,ngpmx2
      isx= iclose('PHIG.d')
      if(nsp2.ne.nsp) then
        write(6,*)'nsp,nsp2',nsp,nsp2
Cstop2rx 2013.08.09 kino        stop 'init_readeigen_phi: ns'
        call rx( 'init_readeigen_phi: ns')
      endif
      if(nlmto2.ne.ldim2) then
        write(6,*)'nlmto,nlmto2',ldim2,nlmto2
Cstop2rx 2013.08.09 kino        stop 'init_readeigen_phi: nlmto'
        call rx( 'init_readeigen_phi: nlmto')
      endif
      if(ngpmx2.ne.ngpmx) then
        write(6,*)'ngpmx,ngpmx2',ngpmx,ngpmx2
Cstop2rx 2013.08.09 kino        stop 'init_readeigen_phi: ngpmx'
        call rx( 'init_readeigen_phi: ngpmx')
      endif
      nwf = nwf2
      nqbz=nqbz2

c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = 0d0

c replace geig and cphi

      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband

      isx= iclose('CPHI')
      isx= iclose('GEIG')
      ifgeig = iopen('GEIGg',0,0,mrecg_o)
      ifcphi = iopen('CPHIg',0,0,mrecb_o)

      if(keepeig) then
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        do ikp= 1,nqtt
          iqbz=ikp
c         iqbz=mod(ikp,nqbz)
c         if (iqbz.eq.0) iqbz=nqbz
          do is= 1,nsp
            ikpisp = is + nsp*(iqbz-1)
            read(ifgeig, rec=ikpisp) geig(1:ngpmx,1:nwf,ikp,is)
            read(ifcphi, rec=ikpisp) cphi(1:ldim2,1:nwf,ikp,is)
          enddo
        enddo
      else
c         allocate(geig2(ngpmx,nwf),cphi2(ldim2,nwf))
c         do ikp= 1,nqtt
c           if (ikp.le.nqbz) cycle
c           iqbz=mod(ikp,nqbz)
c           if (iqbz.eq.0) iqbz=nqbz
c         do is= 1,nsp
c           ikpisp = is + nsp*(iqbz-1)
c           read(ifgeig, rec=ikpisp) geig2(1:ngpmx,1:nwf)
c           read(ifcphi, rec=ikpisp) cphi2(1:ldim2,1:nwf)
c           ikpisp = is + nsp*(ikp-1)
c           write(ifgeig, rec=ikpisp) geig2(1:ngpmx,1:nwf)
c           write(ifcphi, rec=ikpisp) cphi2(1:ldim2,1:nwf)
c         enddo
c         enddo
c         deallocate(geig2,cphi2)
      endif
c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o
c ---
      end subroutine



c-------------------------------------------------------------
      subroutine init_readeigen_mlw_diag(nwf_o,nband_o,mrecb_o,mrecg_o)
c replace cphi, geig and evud
c for hx0fp0_mlw, ixc=11 mode
c this should be called after init_readgeigen2
      implicit none
      integer(4):: iq,is,ifiqg,ikp,
     & iopen,iclose,isx,mrecb_o,ikpisp,mrecg_o,
     & nwf_o,nband_o,ifmlw,ifmlwe,nqbz,nqbze,nqbze2,iqbz,iqbz2,nwf2,
     & ib,iwf,iwf2,iko_ix,iko_fx,in,ifcphi_o,ifgeig_o,
     & ifuu,nqbz2,nq0i,iko_ix2,iko_fx2,iq0i,iq0i2,j1,j2
      real(8):: q(3),rnorm,cnorm
      real(8),allocatable :: eval(:,:,:)
      complex(8),allocatable :: dnk(:,:,:,:),evec(:,:,:,:),
     &                          geig2(:,:,:,:),cphi2(:,:,:,:),
     &                          geig3(:,:),cphi3(:,:),
     &                          geig4(:,:),cphi4(:,:),
     &                          cbwf(:,:,:,:),uum(:,:,:,:,:),
     &                          dnk2(:,:)


      logical :: keepeigen
      keepeig = keepeigen()

      write(6,*)' init_readeigen_mlw'


c --- Readin MLWU/D, MLWEU/D, and UUq0U/D
      do is = 1,nsp

c fileopen
        if (is.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
          ifmlwe = iopen('MLWEU',0,0,0)
          ifuu   = iopen('UUq0U',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
          ifmlwe = iopen('MLWED',0,0,0)
          ifuu   = iopen('UUq0D',0,0,0)
        endif

c nqbz mesh-points
        read(ifmlw)nqbz,nwf,iko_ix,iko_fx
        if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp), dnk2(iko_ix:iko_fx,nwf))
        do iqbz = 1,nqbz
          read(ifmlw)iqbz2,q(1:3)
Cstop2rx 2013.08.09 kino          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
        enddo

        read(ifuu)
        read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
        if (is.eq.1) allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nsp))
        if (nqbz2.ne.nqbz) call rx( "init_readeigen_mlw: nqbz2 error")
        if (iko_ix2.ne.iko_ix) call rx( "init_readeigen_mlw: iko_ix2 error")
        if (iko_fx2.ne.iko_fx) call rx( "init_readeigen_mlw: iko_fx2 error")
        do iqbz = 1,nqbz
          do iq0i =1,nq0i
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
            if (iq0i2 .ne. iq0i) call rx( 'init_readeigen_mlw: iq0i error')
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
          enddo
        enddo

c nqbze diag-points
        read(ifmlw)nqbze,nwf2
Cstop2rx 2013.08.09 kino        if (nwf2 .ne. nwf) stop 'init_readeigen_mlw: nwf error'
        if (nwf2 .ne. nwf) call rx( 'init_readeigen_mlw: nwf error')
        if (nqbze.ne.nqtt .and. nqbze.ne.nqtt-1)
Cstop2rx 2013.08.09 kino     &      stop 'init_readeigen_mlw: nqbze error'
     &      call rx( 'init_readeigen_mlw: nqbze error')
        read(ifmlwe)nqbze2,nwf2
Cstop2rx 2013.08.09 kino        if (nwf2 .ne. nwf) stop 'init_readeigen_mlw: nwf error'
        if (nwf2 .ne. nwf) call rx( 'init_readeigen_mlw: nwf error')
Cstop2rx 2013.08.09 kino        if (nqbze2.ne.nqbze) stop 'init_readeigen_mlw: nqbze2 error'
        if (nqbze2.ne.nqbze) call rx( 'init_readeigen_mlw: nqbze2 error')
        if (is.eq.1)allocate(evec(nwf,nwf,nqtt,nsp),eval(nwf,nqtt,nsp))
        evec = 0d0
        eval = 0d0
        do iqbz = 1,nqbze
          read(ifmlw)iqbz2,q(1:3)
Cstop2rx 2013.08.09 kino          if (iqbz2 .ne. iqbz) stop 'init_readeigen_mlw: iqbz error'
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen_mlw: iqbz error')
          read(ifmlw)evec(1:nwf,1:nwf,iqbz,is)

          read(ifmlwe)iqbz2,q(1:3)
Cstop2rx 2013.08.09 kino          if (iqbz2 .ne. iqbz) stop 'init_readeigen: iqbz error'
          if (iqbz2 .ne. iqbz) call rx( 'init_readeigen: iqbz error')
          read(ifmlwe)eval(1:nwf,iqbz,is)
        enddo

c fileclose
        if (is.eq.1) then
          ifmlw  = iclose('MLWU')
          ifmlwe = iclose('MLWEU')
          ifuu   = iclose('UUq0U')
        else
          ifmlw  = iclose('MLWD')
          ifmlwe = iclose('MLWED')
          ifuu   = iclose('UUq0D')
        endif

c end of is-loop
      enddo


c replace evud
      deallocate(evud)
      allocate(evud(nwf,nqtt,nsp))
      evud = eval
      deallocate(eval)

check write
c      do iq = 1,nqbze
c      write(6800,*)'***',iq
c      write(6800,*)evec(:,:,iq,1)
c      do iwf = 1,nwf
c         write(6800,*)evud(iwf,iq,1)
c      enddo
c      enddo

c replace geig and cphi
      allocate(cbwf(iko_ix:iko_fx,nwf,nqtt,nsp))
      cbwf = 0d0
      do ikp = 1,nqbze
        iqbz = mod(ikp,nqbz)
        if (iqbz .eq. 0) iqbz = nqbz
        iq0i = (ikp - iqbz)/nqbz
        do is= 1,nsp
c
          if (iq0i.eq.0) then
            do ib = iko_ix,iko_fx
              do in = 1,nwf
                do iwf= 1,nwf
                  cbwf(ib,iwf,ikp,is) = cbwf(ib,iwf,ikp,is) + 
     &         dnk(ib,in,iqbz,is)*evec(in,iwf,ikp,is)
                enddo
              enddo
            enddo
c
          else
c   <psi(k+q0,n) | psi(k+q0,m)^B>
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^B>
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^B>
c
c psi^B : bloch fn. corresponding to maxloc Wannier fn.
            dnk2 = 0
            do ib = iko_ix,iko_fx
              do iwf= 1,nwf
                dnk2(ib,iwf) = 
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
              enddo
            enddo
c
            do ib = iko_ix,iko_fx
              do in = 1,nwf
                do iwf= 1,nwf
                  cbwf(ib,iwf,ikp,is) = cbwf(ib,iwf,ikp,is) + 
     &         dnk2(ib,in)*evec(in,iwf,ikp,is)
                enddo
              enddo
            enddo
          endif


ctm norm check
c         do iwf  = 1,nwf
c         do iwf2 = 1,nwf
c           rnorm = 0d0
c           cnorm = 0d0
c           do ib = iko_ix,iko_fx
c              rnorm = rnorm + dreal(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(cbwf(ib,iwf,ikp,is))
c     &                                     *cbwf(ib,iwf2,ikp,is))
c              rnorm = rnorm + dreal(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c              cnorm = cnorm + dimag(dconjg(dnk(ib,iwf,iqbz,is))
c     &                                    *dnk(ib,iwf2,iqbz,is))
c           enddo
c           do ib = 1,nwf
c              rnorm = rnorm + dreal(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c              cnorm = cnorm + dimag(dconjg(evec(ib,iwf,ikp,is))
c     &                                    *evec(ib,iwf2,ikp,is))
c           enddo
c           if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c           write(7700,"(4i5,2f12.6)")is,ikp,iwf,iwf2,rnorm,cnorm
c         enddo
c         enddo
c         write(7300,"(5i5)")is,ikp,iko_ix,iko_fx,nwf
c         write(7300,*)cbwf(:,:,ikp,is)
        enddo
      enddo
      deallocate(dnk,evec,uum,dnk2)

      mrecb_o = mrecb * nwf / nband
      mrecg_o = mrecg * nwf / nband


      if(keepeig) then

        allocate(geig2(ngpmx,nband,nqtt,nsp))
        allocate(cphi2(ldim2,nband,nqtt,nsp))
        geig2 = geig
        cphi2 = cphi
        deallocate(geig,cphi)
        allocate(geig(ngpmx,nwf,nqtt,nsp))
        allocate(cphi(ldim2,nwf,nqtt,nsp))
        geig = 0d0
        cphi = 0d0
        do ikp= 1,nqtt
          do is= 1,nsp
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig(:,iwf,ikp,is) = geig(:,iwf,ikp,is) + 
     &           geig2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
                cphi(:,iwf,ikp,is) = cphi(:,iwf,ikp,is) + 
     &           cphi2(:,ib,ikp,is)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

check write
c            do iwf  = 1,nwf
c            do iwf2 = 1,nwf
c               rnorm = 0d0
c               cnorm = 0d0
c               do ib = 1,ldim2
c                  rnorm = rnorm + dreal(dconjg(cphi(ib,iwf,ikp,is))*
c     &                                   cphi(ib,iwf2,ikp,is))
c               enddo
c               if (iwf.eq.iwf2) rnorm = rnorm - 1d0
c               write(7600,"(4i5,f12.6)")is,ikp,iwf,iwf2,rnorm
c            enddo
c            enddo
c            write(7500,*)ikp,ldim2,nwf
c            write(7500,*)cphi(:,:,ikp,is)
          enddo
        enddo


        deallocate(geig2,cphi2)

      else
        ifcphi_o = iopen('CPHI.mlw', 0,-1,mrecb_o)
        ifgeig_o = iopen('GEIG.mlw', 0,-1,mrecg_o)

        allocate(geig3(ngpmx,nwf))
        allocate(cphi3(ldim2,nwf))
        allocate(geig4(ngpmx,nband))
        allocate(cphi4(ldim2,nband))
        do ikp= 1,nqtt
          do is= 1,nsp
            ikpisp= is + nsp*(ikp-1)
            read(ifgeig, rec=ikpisp) geig4(1:ngpmx,1:nband)
            read(ifcphi, rec=ikpisp) cphi4(1:ldim2,1:nband)

            geig3 = 0d0
            cphi3 = 0d0
            do iwf= 1,nwf
              do ib= iko_ix,iko_fx
                geig3(:,iwf) = geig3(:,iwf) +  geig4(:,ib)*cbwf(ib,iwf,ikp,is)
                cphi3(:,iwf) = cphi3(:,iwf) +  cphi4(:,ib)*cbwf(ib,iwf,ikp,is)
              enddo
            enddo

            write(ifgeig_o, rec=ikpisp) geig3(1:ngpmx,1:nwf)
            write(ifcphi_o, rec=ikpisp) cphi3(1:ldim2,1:nwf)
          enddo
        enddo

        deallocate(geig3,geig4,cphi3,cphi4)
        isx= iclose('CPHI')
        isx= iclose('GEIG')
        isx= iclose('CPHI.mlw')
        isx= iclose('GEIG.mlw')
        ifgeig = iopen('GEIG.mlw',0,0,mrecg_o)
        ifcphi = iopen('CPHI.mlw',0,0,mrecb_o)
      endif

      deallocate(cbwf)

c replace dimensions
      nband = nwf
      nwf_o = nwf
      nband_o = nband
      mrecb = mrecb_o
      mrecg = mrecg_o
c ---
      end subroutine
      end module m_readeigen

      
c-----------------------------------------------------------------------
      module m_iqindx_wan
      implicit none
!! To accelarate inverse mapping, q to 'integer index',
!! we prepare integer index of q, its resolution is given by epsd.
!!      use m_hamindex, only: qtt,nqtt
!      use m_iqindx_qtt, only: init_iqindx_qtt
      public:: Iqindx2_wan
      
      private
      integer,allocatable,private:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
      real(8),private:: epsd=1d-7 !key parameter to map to integer index
      integer,private:: nkey(3)
      real(8),private:: ginv(3,3)
      real(8),private,allocatable::qtt(:,:)
      integer,private::nqtt
      contains
!------------------------------------------------------
      subroutine iqindx2_wan(q,  iqindx,qu)
      intent(in) ::          q
      intent(out) ::             iqindx,qu
!! ginv is the inverse of plat (premitive translation vector).
!! Use kk1,kk2,kk3,nkey(1:3),iqkkk to get iqindx.
      real(8):: q(3),qu(3)
      integer:: iqindx
      integer:: i_out, iq,iqx ,kkk3(3),ik1,ik2,ik3
      real(8):: qx(3),qzz(3)
      logical::debug=.true., init=.true.
      if (init) then
         call wan_getqbz()      !!! get nqtt, qtt(3,nqtt),ginv
         call init_iqindx_wan(ginv)
         init=.false.
      endif
!      print *,"iqindx2_wan:: nqtt=",nqtt
ccccccccccccccccccccc
      debug=.false.
      if(abs(q(1)+0.1d0)+abs(q(2)+0.1d0)<1d-3) then
         debug=.true.
      endif
ccccccccccccccccccccc
      if(debug) write(6,"(' iqindx2_: q=',3f20.15)") q
 !     print *,"ginv,qzz",ginv,qzz
      call rangedq(matmul(ginv,q), qzz)
      if(debug) write(6,"(' iqindx2_: q=',3f20.15)") qzz
!! we generate qzz integer index for qzz
      kkk3 = (qzz+0.5d0*epsd)/epsd 
      if(debug) write(6,*)'kkk3=',kkk3
      if(debug) write(6,*)'nkey=',nkey
      if(debug) write(6,*)'kk1=',kk1
      if(debug) write(6,*)'kk2=',kk2
      if(debug) write(6,*)'kk3=',kk3
      call tabkk(kkk3(1), kk1,nkey(1), ik1)
      call tabkk(kkk3(2), kk2,nkey(2), ik2)
      call tabkk(kkk3(3), kk3,nkey(3), ik3)
      if(debug) write(6,"(' 222222a q=',3i8,3f18.12)") kkk3,qzz
      if(debug) write(6,"(' 222222a ik1,ik2,ik3,q=',3i8,3f18.12)") ik1,ik2,ik3,q
      iqindx = iqkkk(ik1,ik2,ik3)
      if(debug) then
         do iqx=1,nqtt
         write(6,"(i5,3f13.5)")iqx,qtt(:,iqx)
         enddo
      endif
      qu =qtt(:,iqindx)
      if(debug) write(6,*) iqindx,qu
      end subroutine iqindx2_wan
!!-------------------
      subroutine wan_getqbz()
      implicit none
      integer:: ifwqb, ifile_handle
      ifwqb=ifile_handle()
      open(ifwqb,file="wanqbz",form='unformatted')
      read(ifwqb) ginv
      read(ifwqb) nqtt
      allocate(qtt(3,nqtt))
      read(ifwqb) qtt
      close(ifwqb)
      end subroutine wan_getqbz
c--------------------------------
      subroutine init_iqindx_wan(ginv_)
      intent(in):: ginv_
!! For magnon+wannier
!! === mapping of qtt ===
!! nkey, kk1,kk2,kk3, iqkkk are to get iqindx.
!!  q --> call rangedq(matmul(ginv,q), qx) ---> n= (qx+0.5*epsd)/epsd
!!       --->  ik1,ik2,ik3= tabkk(kkk,iqk,nkey) ---> iqkkk(ik1,ik2,ik3)
      real(8):: ginv_(3,3),qzz(3)
      real(8),allocatable:: qxx(:,:)
      integer:: isig,i,ix,kkk,kkk3(3),ik1,ik2,ik3,iq,ik
      integer,allocatable:: ieord(:)
      logical::debug=.false.
      ginv=ginv_
      allocate(ieord(nqtt))
      if(debug) write(6,"(a,2i5,20f9.4)")' iiiiii nqtt=',nqtt,size(qtt),ginv_(1:3,1:3),ginv(1:3,1:3)
      allocate(key(3,0:nqtt),qxx(3,nqtt))
      key=-99999
      do iq=1,nqtt
         call rangedq(matmul(ginv,qtt(:,iq)), qxx(:,iq))
c         write(6,"(a,i5,3f13.5,2x,3f13.5)") ' qqqttxx =',iq,qtt(:,iq),qxx(:,iq)
      enddo
c      write(6,*)'sssqqq=',sum(abs(qtt(1,1:nqtt))),      sum(abs(qtt(2,1:nqtt))),   sum(abs(qtt(3,1:nqtt)))
c
!! get key and nkey for each ix.
      key(:,0)=0 !dummy
      do ix =1,3
         call sortea(qxx(ix,:),ieord,nqtt,isig)
         ik=0
         do i=1,nqtt
            kkk=(qxx(ix,ieord(i))+0.5d0*epsd)/epsd  !kkk is digitized by 1/epsd
            if(i==1.or.key(ix,ik)<kkk) then
               ik=ik+1
               key(ix,ik) = kkk
            elseif (key(ix,ik)>kkk) then
               write(6,*)ix, ik,i, key(ix,ik), qxx(ix,ieord(i))
               call rx( 'iqindx: bug not sorted well')
            endif
         enddo
         nkey(ix)=ik
      enddo
      deallocate(ieord)
!!  key is reallocated. inverse mattping, iqkkk
      allocate( kk1(nkey(1)),kk2(nkey(2)),kk3(nkey(3)) )
      kk1(:) = key(1,1:nkey(1))
      kk2(:) = key(2,1:nkey(2))
      kk3(:) = key(3,1:nkey(3))
      deallocate(key)
      allocate( iqkkk(nkey(1),nkey(2),nkey(3)) )
      iqkkk=-99999
c      print *,'eeee kk1=',kk1
c      print *,'eeee kk2=',kk2
c      print *,'eeee kk3=',kk3
c      write(6,*)' initqqq nqtt=',nqtt
      do i=1,nqtt
         kkk3= (qxx(:,i)+0.5d0*epsd)/epsd !kkk is digitized by 1/epsd
         call tabkk(kkk3(1), kk1,nkey(1), ik1)
         call tabkk(kkk3(2), kk2,nkey(2), ik2)
         call tabkk(kkk3(3), kk3,nkey(3), ik3)
         iqkkk(ik1,ik2,ik3)=i
c         print *,'111111 kkk3=',kkk3,qxx(:,i)
c         write(6,"(' 111111111a ik1,ik2,ik3,iqxx=',3i8,3f18.12)") kkk3,qxx(:,i)
c         write(6,"(' 111111111a ik1,ik2,ik3,iqkkk=',4i8,3f18.12)") ik1,ik2,ik3,i,qtt(:,i)
      enddo
      deallocate(qxx)
      end subroutine init_iqindx_wan      
      end module m_iqindx_wan

c-----------------------------------------
c     taken from lm7k/subs/rotwv.F
      subroutine conv2int(a,ii,n)
      implicit none
      real(8), intent(in) :: a(n)
      integer, intent(in) :: n
      integer, intent(out) :: ii(n)
      integer::i
      real(8):: eps=1d-6,aa !eps=1d-8 until aug2012
      do i=1,n
        ii(i)=nint(a(i))
        if(abs(ii(i)-a(i))>eps)then
           write(6,*)' ii a eps=',ii(i),a(i),eps
           call rx( 'conv2int: a is not close to integer')
        endif   
      enddo
      end subroutine conv2int

      subroutine rx_(string) !error exit routine
      character*(*) string
      write(6,892) string
      write(6,890)
  890 format(' ---- Error exit')
  892 format(/' ---- ',a)
c ... write to log file also
      write(71,710) string
  710 format('++ ',a)
      call exit(-1)
      end
