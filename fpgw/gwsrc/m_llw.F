      module m_llw
      use m_readq0p,only: Readq0p,
     &     wqt,q0i,nq0i ,nq0iadd,ixyz
      use m_rdpp,only: nblochpmx
      use m_genallcf_v3,only: 
     &     nclass,natom,nspin,nl,nn,ngrp,
     &     nlmto,nlnmx, nctot,!niw, !nw_input=>nw,
     &     alat, deltaw,symgrp,clabl,iclass, !diw,dw,delta,
     &     invg, il, in, im, nlnm, 
     &     plat, pos, ecore, symgg , tpioa
      use m_freq,only: 
     &     frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,niw !output of getfreq
      use m_qbze,only: Setqbze,
     &     nqbze,nqibze,qbze,qibze
      use m_read_bzdata,only: Read_bzdata, !<--- 'call read_bzdata' sets up following data.
     &   ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,qbas,ginv,
     &   dq_,qbz,wbz,qibz,wibz,
     &     ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
      use m_readVcoud,only: 
     &     vcousq,zcousq,ngb
      use m_readq0p,only: Readq0p,
     &     wqt,q0i,nq0i ,nq0iadd,ixyz
      use m_rdpp,only: nbloch,mrecl
      implicit none
!-----------------------------
      public:: WVRllwR,WVIllwI,  MPI__sendllw,MPI__sendllw2
      complex(8),allocatable,protected,public:: llw(:,:), llwI(:,:)
      complex(8),allocatable,protected,public:: wmuk(:,:)
      logical,protected,public:: w4pmode
      integer,protected,public:: ngbq0
!-----------------------------      
      private
      complex(8),allocatable:: zw0(:,:),zw(:,:)
      complex(8),allocatable:: epstinv(:,:),epstilde(:,:)
      real(8),parameter:: pi=4d0*datan(1d0),fourpi = 4d0*pi
      contains
      
      subroutine WVRllwR(q,iq,zxq,nmbas1,nmbas2)
      use m_readqg,only: Readqg0
      intent(in)::       q,iq,    nmbas1,nmbas2 !zxq can be twiced when nspin=2
      integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllw
      integer:: nmbas1,nmbas2,ngc0,ifw4p,ifrcw
      real(8):: frr,q(3),vcou1,quu(3),eee
      logical::  localfieldcorrectionllw
      logical,save:: init=.true.
      complex(8):: zxq(nmbas1,nmbas2,nw_i:nw)
      character(5):: charnum5
!     !
      write(6,*)'WVRllwR init zxq dddxxx',nmbas1,nmbas2,nw_i,nq0i,nblochpmx !,zxq(1,1,iw) !, sum(abs(zxq))
      if(init) then !initialization related to w4pmode, zw, tpioa...
         allocate( llw(nw_i:nw,nq0i) )
         init=.false.
         llw= 1d99
         if(sum(ixyz)/=0) w4pmode=.true.
         if(w4pmode) then
            allocate( wmuk(2:nblochpmx,3))
            wmuk=1d99
         endif
!! we have to call WVRllwR beveor WVIllwI
c         tpioa=2d0*pi/alat
         allocate( zw(nblochpmx,nblochpmx) )
      endif   
!! ngb is q-dependent. released at the end of WVIllwi
      call readqg0('QGcou', (/0d0,0d0,0d0/),  quu,ngc0)
      ngbq0 = nbloch+ngc0
      allocate( zw0(ngb,ngb) ) 
      allocate( epstilde(ngb,ngb), epstinv(ngb,ngb))
      if (nspin == 1) zxq = 2d0*zxq !if paramagnetic, multiply x0 by 2
      nwmax = nw
      nwmin  = nw_i
      write(6, *)" === trace check for W-V === nwmin nwmax=",nwmin,nwmax
      write(6,*) 'qqqqqxx=',iq,q
      if(iq<=nqibz) then        !for mmmw
         open(newunit=ifrcw, file='WVR.'//charnum5(iq),form='unformatted',access='direct',recl=mrecl)
         do 1015 iw  = nwmin,nwmax
            frr= dsign(freq_r(abs(iw)),dble(iw))
c            imode = 1
            if(iq==1) then
               ix=1
               zw0(:,1)=0d0
               zw0(1,:)=0d0
            else
               ix=0
            endif
!     !  Eqs.(37),(38) in PRB81 125102 (Friedlich)
            do igb1=ix+1,ngb
               do igb2=ix+1,ngb
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
               enddo
            enddo
            epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
            call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
!     ! w4p writing eps
            if(iw==0.and.w4pmode) then 
!static epstinv is saved. For q=0 epstilde (mu=1 skipped). For q/=0 full matrix inversion.
!(ix=1 is set for q=0) 
               open(newunit=ifw4p,file='W4PHONON.'//charnum5(iq),form='unformatted')
               write(ifw4p) iq,q,ngb,ix !ix=0, or ix=1 for q=0 (iq=1)
               write(ifw4p) epstinv(ix+1:ngb,ix+1:ngb) 
               close(ifw4p)
            endif
            do igb1=1+ix,ngb
               do igb2=1+ix,ngb
                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
               enddo
            enddo
            zw(1:ngb,1:ngb) = zw0
            write(ifrcw, rec= iw-nw_i+1 ) zw !  WP = vsc-v
            call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
 1015    continue
         close(ifrcw)
      else  ! llw, Wing elements of W. See PRB81 125102
         iq0 = iq - nqibz 
         vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
c         print *,'dddvcou1=',iq,nqibz,iq0,q
c         print *,'dddvcou1=',vcou1,tpioa,fourpi,sum(q**2)
         do 1115 iw  = nwmin,nwmax
            frr= dsign(freq_r(abs(iw)),dble(iw))
c            write(6,*)'ddddxxx111a',iw,sum(abs(zxq(:,:,iw))),zxq(1,1,iw),'sss',sum(abs(vcousq(1:ngb)))
c            imode = 1
!! Full inversion to calculalte eps with LFC.
            if(localfieldcorrectionllw()) then
               ix=0
               eee=0d0
               do igb1=ix+1,ngb
               do igb2=ix+1,ngb
                  if(igb1==1.and.igb2==1) then
                     epstilde(igb1,igb2)= 1d0 - vcou1*zxq(1,1,iw)
                     cycle
                  endif
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) then
                     epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                  endif   
               enddo
               enddo
               epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
               call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
               if(iq0<=nq0i) llw(iw,iq0)= 1d0/epstinv(1,1)
!     ! Wing elements calculation july2016
!     ! We need check nqb is the same as that of q=0
               if(ixyz(iq0)/=0.and.iw==0) then
                  if(ngb/=ngbq0) then
                     write(6,*)q,iq0,ngb,ngbq0
                     call rx('hx0p0_sc: ngb/=ngbq0')
                  endif
                  wmuk(2:ngb,ixyz(iq0))=epstinv(1,2:ngb)/epstinv(1,1)
! this is dot(q(:)*w_mu(:,igb)). See PRB125102(2016) eq.(36)
               endif
            else
               if(iq0<=nq0i) llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw) 
            endif 
            if(iq0<=nq0i) write(6,"('epsWVR: iq iw_R omg(iw) eps(wFC) eps(woLFC) ',
     &           2i5,x,10(d13.6,2x,d13.6,x,d13.6,2x,d13.6,x,d13.6))")
     &           iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
 1115    continue               !iw
         open(newunit=ifllw,file='LLWR.'//charnum5(iq),form='unformatted')
         write(ifllw) llw
         close(ifllw)
      endif   
c      deallocate( zw0, epstilde, epstinv)
      end subroutine WVRllwR
      
!--------------------------------------------------------
      subroutine WVIllwi(q,iq,zxqi,nmbas1,nmbas2)
      intent(in)::       q,iq,     nmbas1,nmbas2 !zxqi can be twiced when nspin=2
      integer:: nmbas1,nmbas2
      integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllwi,ifrcwi
      real(8):: frr,q(3),vcou1
      logical::  localfieldcorrectionllw
      logical,save:: init=.true.
      complex(8):: zxqi(nmbas1,nmbas2,niw)
      character(5):: charnum5
      if(init) then
         allocate( llwI(niw,nq0i) )
         init=.false.
         llwI= 1d99
      endif   
c      allocate( zw0(ngb,ngb) ) 
c      allocate( epstilde(ngb,ngb), epstinv(ngb,ngb))
      write(6,*)'WVRllwI: init'
      if (nspin == 1) zxqi = 2d0*zxqi ! if paramagnetic, multiply x0 by 2
      if( iq<=nqibz ) then
        open(newunit=ifrcwi,file='WVI.'//charnum5(iq),form='unformatted',access='direct',recl=mrecl)
        do 1016 iw  = 1,niw
!!  Eqs.(37),(38) in PRB81 125102
           if(iq==1) then
              ix=1
              zw0(:,1)=0d0
              zw0(1,:)=0d0
           else
              ix=0
           endif
           do igb1=ix+1,ngb
              do igb2=ix+1,ngb
                 epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                 if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
              enddo
           enddo
           epstinv=epstilde
           call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
           do igb1=ix+1,ngb
              do igb2=ix+1,ngb
                 zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                 if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
              enddo
           enddo
           zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
           write(ifrcwi, rec= iw)  zw !  WP = vsc-v
           call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
 1016   continue
        close(ifrcwi)
      else
!! Full inversion to calculalte eps with LFC.
        iq0 = iq - nqibz 
        vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
        do 1116 iw  = 1,niw
           if(localfieldcorrectionllw()) then
              ix=0
              do igb1=ix+1,ngb
                 do igb2=ix+1,ngb
                    if(igb1==1.and.igb2==1) then
                       epstilde(igb1,igb2)= 1d0 - vcou1*zxqi(1,1,iw)
                       cycle
                    endif
                    epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                    if(igb1==igb2) then
                       epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                    endif   
                 enddo
              enddo
              epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
              call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
              if(iq0<=nq0i) llwI(iw,iq0)= 1d0/epstinv(1,1)
           else
              if(iq0<=nq0i) llwI(iw,iq0)=  1d0 -vcou1*zxqi(1,1,iw) 
           endif  
           if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)")
     &          iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
 1116   continue
        open(newunit=ifllwi,file='LLWI.'//charnum5(iq),form='unformatted')
        write(ifllwi) llwi
        close(ifllwi)
      endif
      deallocate(epstinv,epstilde,zw0)
      end subroutine WVILLWI

c----------------------------------------------------------
      subroutine MPI__sendllw2(iqxend) !for hx0fp0
      use m_mpi,only: MPI__hx0fp0_rankdivider2,MPI__task,MPI__Initialize,MPI__Finalize,MPI__root,
     &     MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size,
     &     MPI__ranktab,MPI__consoleout,MPI__barrier
      intent(in)::             iqxend
      integer:: iq0,dest,src,iq,iqxend
!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
      if(MPI__size/=1) then
        do iq=nqibz+1,iqxend
          iq0 = iq - nqibz
          if(MPI__ranktab(iq)/=0) then !jan2012
            if(MPI__ranktab(iq) == MPI__rank) then
              dest=0
              call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
              call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
            elseif(MPI__root) then
              src=MPI__ranktab(iq)
              call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
              call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
            endif
          endif
        enddo  
      endif
      end subroutine
c----------------------------------------------------------
      subroutine MPI__sendllw(iqxend) !for hx0fp0_sc
      use m_mpi,only: MPI__hx0fp0_rankdivider2Q,MPI__hx0fp0_rankdivider2S,
     &     MPI__Qtask,MPI__InitializeQSPBM,MPI__Finalize,MPI__root,
     &     MPI__Broadcast,MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__rank,MPI__size,
     &     MPI__Qranktab,MPI__consoleout,MPI__Ss,MPI__Se, MPI__allreducesumS,
     &     MPI__barrier, MPI__rankQ,MPI__rootQ,MPI__rootS
!     ! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
      intent(in)::            iqxend
      integer:: iq0,dest,src,iq,iqxend
      if(MPI__size/=1) then
        do iq=nqibz+1,iqxend
          iq0 = iq - nqibz
          if(MPI__Qranktab(iq)/=0) then !jan2012
            if(MPI__Qranktab(iq) == MPI__rankQ) then
              dest=0
              if(iq0<=nq0i) then
                 call MPI__DbleCOMPLEXsendQ(llw(nw_i,iq0),(nw-nw_i+1),dest)
                 call MPI__DbleCOMPLEXsendQ(llwI(1,iq0),niw,dest)
              endif
              if(ixyz(iq0)/=0) then
                call MPI__DbleCOMPLEXsendQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,dest)
              endif
            elseif(MPI__rootQ) then
              src=MPI__Qranktab(iq)
              if(iq0<=nq0i) then
                 call MPI__DbleCOMPLEXrecvQ(llw(nw_i,iq0),(nw-nw_i+1),src)
                 call MPI__DbleCOMPLEXrecvQ(llwI(1,iq0),niw,src)
              endif
              if(ixyz(iq0)/=0) then
                call MPI__DbleCOMPLEXrecvQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,src)
              endif
            endif
          endif
        enddo  
      endif
      end subroutine
      end module m_llw

          
C===================================================================
      subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
      implicit none
      integer:: nblochpmx,nbloch,ngb,iw,i,iq
      complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
      real(8):: freqq
      character*(*)::tagname
      trwv=0d0
      do i = 1,nbloch
        trwv = trwv + zw(i,i)
      enddo
      trwv2 = 0d0
      do i = 1,ngb
         trwv2 = trwv2 + zw(i,i)
      enddo  !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
      write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
      end


