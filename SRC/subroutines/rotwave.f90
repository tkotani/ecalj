!! space-group rotation core of eigenfunctions, MTpartAPW expansion and MPB.
module m_rotwave
  public:: Rotmto,Rotmto2,Rotipw,Rotipw2,Rotwvigg
contains
  !!------------------------------------------------------
  subroutine rotwvigg(igg,q,qtarget,ndimh,napw_in,nband,evec,evecout,ierr)
    use m_hamindex,only: symops,invgx,miat,tiat,shtvg,qlat,plat,dlmm,ngrp,norbmto, &
         ibastab,ltab,ktab,offl,offlrev,getikt,igv2,igv2rev,napwk,nbas,pwmode
    use m_ftox
    implicit none
    intent(in)::        igg,q,qtarget,ndimh,napw_in,nband,evec
    intent(out)::                                              evecout,ierr
    !! ==  wave function rotator by space group operation.
    !! OUTPUT evecout, ierr
    !! NOTE:
    !! rotation of coefficients on PMT basis.
    !!  phi(r) = \sum_i evec(i,iband) |F_i> ==> Rotated[phi](r)=\sum_i evecout(i,iband) |F_i>  by sym(:,:,ig).
    !!  Rotated[phi](r)= phi[sym^-1(r)], where   sym(r)=r'= symops*r + shftvg.
    integer::ig,ndimh,napw_in,nband,ibaso,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr &
         ,igg,ikt2,ikt,l,ibas,ig2,k
    real(8)::q(3),gout(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3)
    complex(8):: evec(ndimh,nband),evecout(ndimh,nband),phase(nbas)
    real(8),parameter:: tolq=1d-4
    complex(8),parameter:: img=(0d0,1d0), img2pi=2*4d0*datan(1d0)*img
    platt = transpose(plat) !this is inverse of qlat
    ierr=1
    !! check q is really rotated to qtarget by symops(:,:,igg)
    call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,igg),q)) ), qx)
    if(sum(abs(qx))>tolq) then
       write(6,"(a,3f7.3,2x,3f7.3)")'  rotwvigg: qtarget is not a star of q',q,qtarget
       call rx( 'rotwvigg: qtarget is not symops(:,:,ig)*q')
    endif
    evecout = 0d0
    nlmto = ndimh-napw_in
    !! mto part
    if(nlmto/=0) then
       phase = [(exp(-img2pi*sum(qtarget*tiat(:,ibas,igg))),ibas=1,nbas)]
       do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          ibas = ibastab(iorb)
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+1
          iend1 = offl(iorb)+2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+1
          iend2 = offlrev(miat(ibas,igg),l,k)+2*l+1
          evecout(init2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),evec(init1:iend1,:))*phase(ibas)
       enddo
    endif
    !! apw part
    if(napw_in/=0) then
       ikt  = getikt(q)       !index for q
       ikt2 = getikt(qtarget) !index for qtarget
       if(napw_in /= napwk(ikt) ) then
          call rx_('rotwv: napw_in /= napw(ikt)')
       endif
       do ig = 1,napw_in
          if(pwmode>10) then
             qpg  = q + matmul( qlat(:,:),igv2(:,ig,ikt))  !q+G
             qpgr = matmul(symops(:,:,igg),qpg)            !rotated q+G
             nnn= nint(matmul(platt,qpgr-qtarget)) !integer representation of G= qpgr - qtarget
          else   
             block
               real(8):: gg(3),ggr(3) 
               gg  = matmul(qlat(:,:),igv2(:,ig,ikt))  !q+G
               ggr = matmul(symops(:,:,igg),gg)            !rotated G
               nnn = nint(matmul(platt,ggr)) !integer representation of G= qpgr - qtarget
             endblock
          endif
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2) !get index of G
          if(ig2>=999999) then
             block
               integer:: i1
             do i1=1,napwk(ikt)
                write(6,ftox)'yyy0 igv2', ftof(q,3),     ikt,i1, ' ',igv2(:,i1,ikt)
             enddo
             do i1=1,napwk(ikt2)
                write(6,ftox)'yyy1 igv2', ftof(qtarget,3),ikt2,i1,' ',igv2(:,i1,ikt2)
             enddo
             endblock
             write(6,ftox)'rotwvigg: q=',ftof( q,3),'qtarget=', ftof(qtarget,3)
             write(6,ftox)'rotwvigg  qr=',ftof(matmul(symops(:,:,igg),q),3)
             write(6,ftox)'rotwvigg: qpg=',ftof(qpg,3),'qpgr=', ftof(qpgr,3)
             write(6,ftox)'rorwvigg: igv2rev ikt2=',nnn(1),nnn(2),nnn(3)
             call rx('rotwvigg can not find index of mapped G vector ig2')
          endif
          evecout(nlmto+ig2,:)= evec(nlmto+ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
       enddo
    endif
    ierr=0
  end subroutine rotwvigg
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine rotmto(qin,nbloch,nband, &
       norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, &
       sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas,  cphiin, &
       cphiout)
    implicit none
    intent(in)::      qin,nbloch,nband, &
         norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, &
         sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas,  cphiin
    intent(out):: &
         cphiout
    !!== Rotation of a function of ProductBasis or MT part of wave function.
    !! (then MTO itself is rotated. a little diffent from here).
    !! qin  = qtti(:,iqi)
    !! cphiin cphi(1:ldim2,1:nband,iqi,isp) eigenfunction @ qtti(:,iqi)
    !! dlmm   dlmm(-l:l,-l:l,l,igg)
    !! sym   symops(:,:,igg)
    !! tiat   tiat(:,ibas,igg)
    !! shtvg  shtvg(:,igg)
    !! miat   miat(ibas,igg)
    !!  cphin (iqq) ---> cphiout (iq)
    ! input variables
    integer:: nbloch,nband,norbt,nbas, ibas_tbl(norbt),l_tbl(norbt),k_tbl(norbt),offset_tbl(norbt)
    integer:: miat(nbas),lxxa,ibas,ibaso,l,k,ini1,iend1,ini2,iend2,iorb
    integer:: max_ibas_tbl,max_l_tbl,max_k_tbl,igxt,offset_rev_tbl(max_ibas_tbl, 0:max_l_tbl, max_k_tbl)
    real(8) :: qrot(3), tiat(3,nbas),sym(3,3), qin(3),shtvg(3), dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa)
    complex(8):: cphiin(nbloch,1:nband)
    complex(8):: cphiout(nbloch,1:nband),phase(nbas)
    complex(8),parameter:: img=(0d0,1d0), img2pi = 2d0*4d0*datan(1d0)*img
    qrot = matmul(sym,qin)
    if(igxt==-1) qrot=-qrot !july2012takao
    phase = [(exp(-img2pi*sum(qrot*tiat(:,ibas))),ibas=1,nbas)]
    do iorb=1,norbt !orbital-blocks
       ibas = ibas_tbl(iorb)
       l = l_tbl(iorb)
       k = k_tbl(iorb)
       ini1 = offset_tbl(iorb)+1
       iend1 = ini1+2*l
       ini2 = offset_rev_tbl(miat(ibas),l,k)+1
       iend2 = ini2+2*l
       cphiout(ini2:iend2,:)= matmul(dlmm(-l:l,-l:l,l),cphiin(ini1:iend1,:))*phase(ibas)
    enddo
  end subroutine rotmto
  !! ------------------------------------------------------
  subroutine rotmto2(qin,nbloch,nband, &
       norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, &
       sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas, &
       zrotm)
    implicit none
    intent(in)::       qin,nbloch,nband, &
         norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, &
         sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas
    intent(inout):: &
         zrotm
    !!== Similar with rotmto but just for rotation matrix ==
    !! output: zrotm is the rotation matrix
    integer:: ibas,ibaso,nbloch,nband,norbt,iorb,nbas
    integer:: ibas_tbl(norbt),l_tbl(norbt),k_tbl(norbt),offset_tbl(norbt)
    integer:: l,k,ini1,iend1,ini2,iend2,miat(nbas),lxxa
    real(8):: tiat(3,nbas),sym(3,3), qin(3),qrot(3),shtvg(3) !shtvg
    real(8):: dlmm(-lxxa:lxxa,-lxxa:lxxa,0:lxxa)
    complex(8),parameter:: img=(0d0,1d0),img2pi = 2d0*4d0*datan(1d0)*img
    integer:: max_ibas_tbl,max_l_tbl,max_k_tbl,igxt
    integer::  offset_rev_tbl(max_ibas_tbl, 0:max_l_tbl, max_k_tbl),m1,m2,nrotm,irotm1,irotm2
    complex(8):: phase(nbas),zrotm(nbloch,nbloch)
    qrot = matmul(sym,qin)
    if(igxt==-1) qrot=-qrot !july2012takao
    phase = [(exp(-img2pi*sum(qrot*tiat(:,ibas))),ibas=1,nbas)]
    zrotm=0d0
    do iorb=1,norbt !orbital-blocks
       ibas = ibas_tbl(iorb)
       l    = l_tbl(iorb)
       k    = k_tbl(iorb)
       ini1  = offset_tbl(iorb)+1
       iend1 = ini1+2*l
       ini2  = offset_rev_tbl(miat(ibas),l,k)+1
       iend2 = ini2+2*l
       zrotm(ini2:iend2,ini1:iend1) = dlmm(-l:l,-l:l,l)*phase(ibas)
    enddo
  end subroutine rotmto2
  !! --------------------------------
  !> Rotation of Plane wave part. by sym
  !! Mapped from qtt(:,iqq) to qtt(:,iq)
  !!   qtt(:,iq)= matmul(sym(igg),qtt(:,iqq))+some G vector
  !!  geigenin (iqq) ---> geigenout (iq)
  subroutine rotipw(qin,qtarget,ngp,nband,platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx, &
       geigenin, &
       geigenout)
    implicit none
    intent(in)::      qin,qtarget,geigenin,ngp,nband,platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx
    integer :: ig,ig2,nnn(3)
    integer :: ngp,imx,igxt,nband, ngvecp(3,ngp),ngvecprev(-imx:imx,-imx:imx,-imx:imx)
    real(8) :: sym(3,3),qlat(3,3),platt(3,3),shtvg(3),qin(3),qtarget(3)
    real(8):: qpg(3),qpgr(3)
    complex(8) :: geigenin(ngp,nband), geigenout(ngp,nband)
    complex(8),parameter:: img=(0d0,1d0),img2pi = 2d0*4d0*datan(1d0)*img
    do ig = 1,ngp
       qpg = qin + matmul( qlat(:,:),ngvecp(:,ig)) ! q+G
       qpgr = matmul(sym,qpg)             !rotated q+G
       if(igxt==-1) qpgr=-qpgr ! xxxxxxxx need to check!
       nnn = nint( matmul(platt,qpgr-qtarget)) ! integer-representation of G=qpgr-qtarget
       ig2 = ngvecprev(nnn(1),nnn(2),nnn(3))   ! index for G
       geigenout(ig2,:)= geigenin(ig,:) * exp( -img2pi*sum(qpgr*shtvg) )
    enddo
  end subroutine rotipw
  !! --------------------------------
  subroutine rotipw2(qin,qtarget,ngp,nband,platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx, &
       zrotm)
    implicit none
    intent(in)::       qin,qtarget,ngp,nband,platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx
    !! similar with rotipw
    !! output: zrotm
    real(8):: sym(3,3),qlat(3,3),platt(3,3),shtvg(3)
    integer:: ngp,imx,nband, &
         ngvecp(3,ngp), ngvecprev(-imx:imx,-imx:imx,-imx:imx) 
    real(8):: qin(3),qpg(3),qpgr(3),qtarget(3)
    integer:: ig,ig2,nnn(3),igxt
    complex(8),parameter:: img=(0d0,1d0),img2pi = 2d0*4d0*datan(1d0)*img
    complex(8):: zrotm(ngp,ngp)
    do ig = 1,ngp
       qpg = qin + matmul( qlat(:,:),ngvecp(:,ig)) !q+G
       qpgr = matmul(sym,qpg)  !rotated q+G
       if(igxt==-1) qpgr=-qpgr
       nnn= nint( matmul(platt,qpgr-qtarget)) ! G= qpgr-qtarget
       ig2 = ngvecprev(nnn(1),nnn(2),nnn(3))  ! indeg of G
       zrotm(ig2,ig)= exp( -img2pi*sum(qpgr*shtvg) )
    enddo
  end subroutine rotipw2
end module m_rotwave
