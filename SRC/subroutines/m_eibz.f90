module m_eibz
  use m_read_bzdata,only: Read_bzdata, nqbz,nqibz,n1,n2,n3,ginv, &
       dq_,qbz,wbz,qibz,wibz, &
       ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
  use m_genallcf_v3,only: Genallcf_v3, &
       nclass,natom,nspin,nl,nn, &
       nlmto,nlnmx, nctot, &
       alat, deltaw,clabl,iclass, &
       il, in, im, nlnm, &
       plat, pos, ecore
  use m_hamindex,only: symgg =>symops
  use m_qbze,only: Setqbze, &
       nqbze,nqibze,qbze,qibze
  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_rdpp,only: Rdpp, nxx,lx,nx,mdimx,nbloch,cgr,ppbrd
  use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
  use m_readqgcou,only: Readqgcou
  use m_hamindex,only: ngrp
  implicit none

  integer,allocatable,protected:: nwgt(:,:)
  integer,allocatable,protected:: neibz(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)

contains

  subroutine Seteibz(iqxini,iqxend,iprintx)
    logical:: eibzmode,eibz4x0,tiii,iprintx
    integer:: iqxini,iqxend,iqxendx,timereversal
    eibzmode = eibz4x0()
    !! === Use of symmetry. EIBZ procedure PRB81,125102 ===
    !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
    if(eibzmode) then
       !! commentout block inversion Use iqxendx=iqxend because of full inversion
       iqxendx=iqxend
       allocate( nwgt(nqbz,iqxini:iqxendx),igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx), &
            eibzsym(ngrp,-1:1,iqxini:iqxendx))
       iprintx=.false.
       write(6,*)
       write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal()
       call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz, &
            timereversal(),ginv,iprintx, &
            nwgt,igx,igxt,eibzsym,tiii)
       write(6,"('Used timeRevesal for EIBZ = ',l1)") tiii
       call cputid(0)
       call PBindex(natom,lx,2*(nl-1),nx) !this returns required index stored in arrays in m_pbindex.
       ! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
       call readqgcou() ! no input. Read QGcou and store date into variables.
    else
       iqxendx=iqxend
       allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx) &
            ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
    endif
  end subroutine Seteibz
end module m_eibz

!! ------------------------------------------------------------------
module m_eibzhs
  use m_read_bzdata,only: qbz,nqibz,ginv,irk,nqbz
  use m_hamindex,only: ngrp, symgg=>symops
  use m_readhbe,only: nband
  use m_rdpp,only:  nbloch
  implicit none
  !------------------------------
  integer,allocatable,protected:: nwgt(:,:)
  integer,allocatable,protected:: irkip_all(:,:,:,:),nrkip_all(:,:,:,:),eibzsym(:,:,:)
  logical:: tiii
  !------------------------------
contains
  subroutine Seteibzhs(nspinmx,nq,q,iprintx)
    intent(in)::         nspinmx,nq,q,iprintx
    integer ::   nspinmx,nq,iqxini,iqxend,iqq,is,kr,kx,igrp
    integer,allocatable:: neibz(:),nwgt(:,:),igx(:,:,:),igxt(:,:,:)
    logical:: eibz4sig,tiiiout,iprintx
    real(8):: q(3,nq)
    allocate(irkip_all(nspinmx,nqibz,ngrp,nq)) !this is global
    allocate(nrkip_all(nspinmx,nqibz,ngrp,nq)) !this is global
    !! eibz4sig() is EIBZ symmetrization or not...
    if(eibz4sig()) then
       allocate(nwgt(nqbz,1:nq),igx(ngrp*2,nqbz,nq))
       allocate(igxt(ngrp*2,nqbz,nq), eibzsym(ngrp,-1:1,nq))
       iqxini=1
       iqxend=nq
       !         write(6,"('TimeRevesal switch = ',l1)") timereversal()
       !         call eibzgen(nq,symgg,ngrp,q(:,iqxini:iqxend),
       !     &        iqxini,iqxend,qbz,nqbz,timereversal(),ginv,iprintx,
       !     o        nwgt,igx,igxt,eibzsym,tiii)
       !! Check timereversal is required for symmetrization operation or not. If tiii=timereversal=F is enforced,
       !! the symmetrization procedure in x0kf_v4h becomes a little time-consuming.
       tiii=.false.            !Enforce no time reversal. time reversal not yet...
       write(6,*)'NOTE:TimeReversal not yet implemented in hsfp0.sc.m.F'
       write(6,"('=== goto eibzgen === used timereversal=',l1)")tiii
       call eibzgen(nq,symgg,ngrp,q(:,iqxini:iqxend), &
            iqxini,iqxend,qbz,nqbz,tiii,ginv,iprintx, &
            nwgt,igx,igxt,eibzsym,tiiiout)
       nrkip_all=0
       irkip_all=0
       is=1                    ! not spin dependent
       do iqq=1,nq
          do kx=1,nqibz
             do igrp=1,ngrp
                kr = irk(kx,igrp) !ip_all(is,kx,igrp,iqq) !kr is index for qbz (for example, nonzero # of kr is 64 for 4x4x4)
                if(kr==0) cycle
                if(nwgt(kr,iqq)/=0) then
                   irkip_all(is,kx,igrp,iqq)= irk(kx,igrp)
                   nrkip_all(is,kx,igrp,iqq)= nwgt(kr,iqq)
                endif
             enddo
          enddo
       enddo
       !          do iqq=1,nq
       !             write(6,"('iq=',i4,' # of EIBZ: Used(TimeR 1 or -1)=',i3,'=',i3,'+',i3)")
       !      &           iqq,sum(eibzsym(:,:,iqq)),sum(eibzsym(:,1,iqq)),sum(eibzsym(:,-1,iqq))
       !             write(6,"('eibz: iqq sum(nrkip_all)=nqbz  ',i3,3f11.5,3i8)")
       !      &           iqq,q(:,iqq),sum(nrkip_all(is,:,:,iqq)),nqbz
       !             do kx=1,nqibz
       !                do igrp=1,ngrp
       !                   kr = irkip_all(is,kx,igrp,iqq) !kr is index for qbz
       !                   if(kr/=0) write(6,"('      ',i8,3f11.5,i8,2x,25(i4,i2))")
       !      &                 kr,qbz(:,kr),nrkip_all(is,kx,igrp,iqq)
       !      &                 ,(igx(i,kr,iqq),igxt(i,kr,iqq),i=1,nwgt(kr,iqq))
       !                enddo
       !             enddo
       !!     !   Probably partial group symmetrization is enough. But it may not reduce computational time so much.
       !         enddo
       if(nspinmx==2) then
          irkip_all(2,:,:,:)=irkip_all(1,:,:,:)
          nrkip_all(2,:,:,:)=nrkip_all(1,:,:,:)
       endif
    else                      ! not eibz4sig
       do is = 1,nspinmx
          do iqq=1,nq
             irkip_all(is,:,:,iqq)=irk
          enddo
       enddo
    endif
  end subroutine Seteibzhs
end module m_eibzhs