! Get effective W(q=0,omega) for GW.
module m_w0w0i
  use m_llw,only: llw,llwI,wmuk,ngbq0,w4pmode
  use m_read_bzdata,only: q0i,nq0i, ixyz
  use m_genallcf_v3,only:  tpioa
  !! See the Friedlich's paper.
  implicit none

  public:: W0w0i
  complex(8),allocatable,protected,public :: w0(:),w0i(:),llmat(:,:)

  private
contains
  !----------------------------------------------------------
  subroutine finalizew4p()
    !      use m_w0w0i,only: llmat
    integer:: i,igb,ifw4p
    real(8):: qv(3,3)
    complex(8),allocatable:: wmu(:,:)
    !! Finalize w4phonon
    !! wmuk(ix)= matmul(wmu,qv) ==> wmu= matmul(wmuk,qvinv)
    do i=1,3
       qv(:,i)= tpioa*q0i(:,ixyz(i))
       qv(:,i)= qv(:,i)/sqrt(sum(qv(:,i)**2))
    enddo
    call matinv(3,qv)
    allocate( wmu(2:ngbq0,3) )
    do igb=2,ngbq0
       wmu(igb,:) =matmul(wmuk(igb,:),qv)
    enddo
    open(newunit=ifw4p,file='W4PHONON.HeadWing',form='unformatted')
    write(ifw4p) llmat(1:3,1:3),ngbq0 !for q~0
    write(ifw4p) wmu(2:ngbq0,1:3) !for q~0
    close(ifw4p)
    deallocate(wmu)
  end subroutine finalizew4p
  !-------------------------------------------------------------
  subroutine w0w0i(nw_i,nw,nq0i,niw,q0i) !llw,llwI,
    !! Get w0 and w0i (diagonal element at Gamma point) for given llw and llwi
    !! Outputs w0,w0i,llmat. See use m_w0w0i at the begining of this routine.
    use m_read_bzdata,only:  wbz,lxklm,dmlx,epinvq0i,epinv,wklm
    intent(in)::     nw_i,nw,nq0i,niw,q0i !llw,llwI,
    real(8)   :: q0i(1:3,1:nq0i)
    integer:: nw_i,nw,nq0i,nq0ix,niw,ifidmlx ,ifile_handle,i,ifw0w0i,ixc,nlxklm
    logical:: readw0w0itest
    complex(8):: llmat_dummy(3,3)
    write(6,*) 'w0w0i:'
    write(6,*)' ==== newaniso mode W(0) divergent part ==== '
    !! == W(0) divergent part ==  starting from llw(iw,iq0),llwI(iw,iq0)
    !! === <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
    allocate(w0(nw_i:nw),w0i(niw),llmat(3,3))
    write(6,*)' goto getw0 nq0i epinvq0i=',nq0i,epinvq0i,wklm
    !! wbz(1) is the weight for q=0 = 1/(n1*n2*n3)
    !! llmat is added. July2016. llw is calculated at iw=0, when nw_i<=0
    call getw0(llw, nw_i,nw,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0,  llmat)
    call getw0(llwI,1,niw  ,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0i, llmat_dummy)
    !$$$!! test mode
    !$$$      if(ixc/=1011) then
    !$$$c      ifw0w0i = ifile_handle()
    !$$$      open(newunit=ifw0w0i,file='W0W0I',form='unformatted')
    !$$$      write(ifw0w0i) nw_i,nw,niw,nq0i
    !$$$      write(ifw0w0i) llw(nw_i:nw,1:nq0i)
    !$$$      write(ifw0w0i) llwI(1:niw,1:nq0i)
    !$$$      write(ifw0w0i) w0(nw_i:nw)
    !$$$      write(ifw0w0i) w0i(1:niw)
    !$$$      close(ifw0w0i)
    !$$$      endif
    do i=nw_i,nw
       write(6,"('w0 =',i4,2f13.4)")i,w0(i)
    enddo
    do i=1,niw
       write(6,"('w0i=',i4,2f13.4)")i,w0i(i)
    enddo
    !! modivy files WVR and WVI
    call ModifyWV0()
    !!
    if(w4pmode) call FinalizeW4p() !W for phonon mode finalized.
  end subroutine w0w0i
  !-------------------------------------------------
  subroutine modifyWV0()
    use m_qbze,only: &
         nqbze,nqibze,qbze,qibze
    use m_rdpp,only: nblochpmx,mrecl
    use m_freq,only: niw ,nw,nw_i
    integer:: ifrcwx,iq,ircw,iw,nini,nend,mreclx
    real(8)::q(3)
    complex(8),allocatable:: zw(:,:)
    character(10):: i2char
    mreclx=mrecl
    !! Read WVR and WVI at Gamma point, and give correct W(0) (averaged in the Gamma cell, where
    !! Gamma cell) is the micro cell of BZ including Gamma point).
    !! === w0,w0i are stored to zw for q=0 ===
    !! === w_ks*wk are stored to zw for iq >nqibz ===
    ! We assume iq=1 is for rank=0
    allocate( zw(nblochpmx,nblochpmx) )
    iq = 1             !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
    q = qibze(:,iq)
    do ircw=1,2
       if (ircw==1) then
          nini=nw_i
          nend=nw
          open(newunit=ifrcwx,  file='WVR.'//i2char(iq), form='unformatted', &
               status='old',access='direct',recl=mreclx)
       elseif(ircw==2) then;  nini=1;      nend=niw;
          open(newunit=ifrcwx,  file='WVI.'//i2char(iq), form='unformatted', &
               status='old',access='direct',recl=mreclx)
       endif
       do iw=nini,nend
          read(ifrcwx, rec= iw-nini+1 ) zw !(1:ngb,1:ngb)
          if( iq==1 ) then
             if(ircw==1) zw(1,1) = w0(iw)
             if(ircw==2) zw(1,1) = w0i(iw)
          endif
          write(ifrcwx,rec=iw-nini+1) zw !(1:ngb,1:ngb)
       enddo
       close(ifrcwx)
    enddo
  end subroutine modifyWV0

end module m_w0w0i
