module m_rotwave ! PMT wave funciton rotation ! space-group rotation of eigenfunctions, MTpartAPW expansion and MPB.
   use m_nvfortran,only:findloc
   public:: Rotmto,Rotmto2,Rotipw,Rotipw2,Rotevec
contains
   subroutine rotmatMTO(igg,q,qtarget,ndimh, rotmat) ! Rotation matrix for MTO
      !usage:
      !   qpr = matmul(symops(:,:,igrp),qp)
      !   call rotmatMTO(igrp,qp,qpr,ndimMTO, rotmat)
      !   hamm(1:ndimMTO,1:ndimMTO)=matmul(rotmat,matmul(hamm0(1:ndimMTO,1:ndimMTO),dconjg(transpose(rotmat))))
      !   ovlm(1:ndimMTO,1:ndimMTO)=matmul(rotmat,matmul(ovlm0(1:ndimMTO,1:ndimMTO),dconjg(transpose(rotmat))))
      !   use m_qplist,only: igv2qp,igv2revqp,napwkqp,qplist,nkp
      use m_mksym,only:   symops,miat,tiat,shtvg,dlmm,ngrp
      use m_lmfinit,only: norbmto,ibastab,ltab,ktab,offl,offlrev,nbas
      use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
      use m_lgunit,only: stdo
      use m_ftox
      implicit none
      intent(in)::     igg,q,qtarget,ndimh
      intent(out)::                          rotmat
      !   phi(r) = \sum_i evec(i,iband) |F_i> ==> Rotated[phi](r)=\sum_i evecout(i,iband) |F_i>  by sym(:,:,ig).
      !   Rotated[phi](r)= phi[sym^-1(r)], where   sym(r)=r'= symops*r + shftvg.
      !This comment is Checked at 2023-9-13
      integer::i,ig,ndimh,nband,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr,igg,ikt2,ikt,l,ibas,ig2,k,i1
      real(8)::q(3),gout(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3),qpgr2(3)
      complex(8):: phase(nbas),rotmat(ndimh,ndimh)
      real(8),parameter:: tolq=1d-4
      complex(8),parameter:: img=(0d0,1d0), img2pi=2*4d0*datan(1d0)*img
      character(256)::aaa
      logical::errig2=.false.
      platt = transpose(plat) ! inverse of qlat
      call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,igg),q)) ), qx) ! Check equivalence of q and qtarget
      if(sum(abs(qx))>tolq) then
         write(aaa,"(a,3f7.3,2x,3f7.3)")' 111 qtarget is not a star of q',q,qtarget
         call rx( '111 rotwvigg: qtarget is not symops(:,:,ig)*q'//trim(aaa))
      endif
      rotmat= 0d0
      nlmto = ndimh
      phase = [(exp(-img2pi*sum(qtarget*tiat(:,ibas,igg))), ibas=1,nbas)]
      OrbitalBlock: do iorb=1,norbmto
         ibas = ibastab(iorb)
         l   = ltab(iorb)
         k   = ktab(iorb)
         init1 = offl(iorb)+1
         iend1 = offl(iorb)+2*l+1
         init2 = offlrev(miat(ibas,igg),l,k)+1
         iend2 = offlrev(miat(ibas,igg),l,k)+2*l+1
         rotmat(init2:iend2,init1:iend1)=dlmm(-l:l,-l:l,l,igg)*phase(ibas)
      enddo OrbitalBlock     ! evecout(init2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),evec(init1:iend1,:))*phase(ibas)
   end subroutine rotmatMTO

   subroutine rotmatPMT_ProbablyOKbutNotCheckedYet(igg,q,qtarget,ndimh,napw, rotmat) ! Rotation matrix for PMT
      use m_qplist,only: igv2qp,igv2revqp,napwkqp,qplist,nkp
      use m_mksym,only:   symops,miat,tiat,shtvg,dlmm,ngrp
      use m_lmfinit,only: norbmto,ibastab,ltab,ktab,offl,offlrev,nbas
      use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
      use m_lgunit,only: stdo
      use m_ftox
      implicit none
      !   phi(r) = \sum_i evec(i,iband) |F_i> ==> Rotated[phi](r)=\sum_i evecout(i,iband) |F_i>  by sym(:,:,ig).
      !   Rotated[phi](r)= phi[sym^-1(r)], where   sym(r)=r'= symops*r + shftvg.
      !This comment is Checked at 2023-9-13
      integer::i,ig,ndimh,napw,nband,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr,igg,ikt2,ikt,l,ibas,ig2,k,i1
      real(8)::q(3),gout(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3),qpgr2(3)
      complex(8):: rotmat(ndimh,ndimh)
      real(8),parameter:: tolq=1d-4
      complex(8),parameter:: img=(0d0,1d0), img2pi=2*4d0*datan(1d0)*img
      character(256)::aaa
      logical::errig2=.false.
      nlmto = ndimh-napw
      rotmat=0d0
      call rotmatMTO(igg,q,qtarget,nlmto, rotmat(1:nlmto,1:nlmto)) ! Rotation matrix for MTO
      platt = transpose(plat) ! inverse of qlat
      call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,igg),q)) ), qx) ! Check equivalence of q and qtarget
      APWpart: if(napw/=0) then
         ikt  = findloc([(sum(abs(q      -qplist(:,i)))<1d-8,i=1,nkp)],value=.true.,dim=1)  !=index for q
         ikt2 = findloc([(sum(abs(qtarget-qplist(:,i)))<1d-8,i=1,nkp)],value=.true.,dim=1)
         !write(stdo,ftox)'rotevec: ikt q=',ikt,ftof(q,3),' ikt2 q=',ikt2,ftof(qtarget,3),'q-qtarget=',ftof(matmul(platt,q-qtarget),3)
         if(napw /= napwkqp(ikt) ) call rxii('rotevec: napw_in /= napw(ikt)',napw,napwkqp(ikt))
         igloop: do ig = 1,napw
            qpg  = q + matmul(qlat(:,:),igv2qp(:,ig,ikt))! q+G
            qpgr = matmul(symops(:,:,igg),qpg)           ! rotated. g(q+G) = q'+G' = qpgr
            nnn  = nint(matmul(platt,qpgr-qtarget))      ! integer representation for G'= qpgr - q'
            ig2 = igv2revqp(nnn(1),nnn(2),nnn(3),ikt2)   ! get index of G'
            !if(ig2>=999999) call rx('rotwave: q+G rotation error. (We have to set PWmode=11 for symgrpAF)')
            !evecout(nlmto+ig2,:)= evec(nlmto+ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
            rotmat(nlmto+ig2, nlmto+ig)=  exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
         enddo igloop
      endif APWpart
    endsubroutine rotmatPMT_ProbablyOKbutNotCheckedYet

   subroutine rotevec(igg,q,qtarget,ndimh,napw_in,nband,evec, evecout) ! Rotation of coefficients evec in PMT basis for q in qplist.
      use m_qplist,only: igv2qp,igv2revqp,napwkqp,qplist,nkp
      use m_mksym,only:   symops,miat,tiat,shtvg,dlmm,ngrp
      use m_lmfinit,only: norbmto,ibastab,ltab,ktab,offl,offlrev,nbas
      use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
      use m_lgunit,only: stdo
      use m_ftox
      implicit none
      intent(in)::     igg,q,qtarget,ndimh,napw_in,nband,evec
      intent(out)::                                            evecout
      !   phi(r) = \sum_i evec(i,iband) |F_i> ==> Rotated[phi](r)=\sum_i evecout(i,iband) |F_i>  by sym(:,:,ig).
      !   Rotated[phi](r)= phi[sym^-1(r)], where   sym(r)=r'= symops*r + shftvg.
      !This comment is Checked at 2023-9-13
      integer::i,ig,ndimh,napw_in,nband,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr,igg,ikt2,ikt,l,ibas,ig2,k,i1
      real(8)::q(3),gout(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3),qpgr2(3)
      complex(8):: evec(ndimh,nband),evecout(ndimh,nband),phase(nbas)
      real(8),parameter:: tolq=1d-4
      complex(8),parameter:: img=(0d0,1d0), img2pi=2*4d0*datan(1d0)*img
      character(256)::aaa
      logical::errig2=.false.
      platt = transpose(plat) ! inverse of qlat
      call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,igg),q)) ), qx) ! Check equivalence of q and qtarget
      if(sum(abs(qx))>tolq) then
         write(aaa,"(a,3f7.3,2x,3f7.3)")' qtarget is not a star of q',q,qtarget
         call rx( 'rotwvigg: qtarget is not symops(:,:,ig)*q'//trim(aaa))
      endif
      evecout = 0d0
      nlmto = ndimh-napw_in
      MTOpart: if(nlmto/=0) then
         phase = [(exp(-img2pi*sum(qtarget*tiat(:,ibas,igg))), ibas=1,nbas)]
         OrbitalBlock: do iorb=1,norbmto
            ibas = ibastab(iorb)
            l   = ltab(iorb)
            k   = ktab(iorb)
            init1 = offl(iorb)+1
            iend1 = offl(iorb)+2*l+1
            init2 = offlrev(miat(ibas,igg),l,k)+1
            iend2 = offlrev(miat(ibas,igg),l,k)+2*l+1
            !write(stdo,ftox)'iorb',iorb,'data from',init1,iend1,'to', init2,iend2,'ibas l k',ibas,l,k,&
            !     'igg=',igg,' miat',miat(ibas,igg),'tiat=',ftof(tiat(:,ibas,igg),3)
            evecout(init2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),evec(init1:iend1,:))*phase(ibas)
!          rotmat(init2:iend2,init1:iend1)=dlmm(-l:l,-l:-l,l,igg)*phase(ibas)
         enddo OrbitalBlock
      endif MTOpart
      APWpart: if(napw_in/=0) then
         ikt  = findloc([(sum(abs(q      -qplist(:,i)))<1d-8,i=1,nkp)],value=.true.,dim=1)  !=index for q
         ikt2 = findloc([(sum(abs(qtarget-qplist(:,i)))<1d-8,i=1,nkp)],value=.true.,dim=1)
         !write(stdo,ftox)'rotevec: ikt q=',ikt,ftof(q,3),' ikt2 q=',ikt2,ftof(qtarget,3),'q-qtarget=',ftof(matmul(platt,q-qtarget),3)
         if(napw_in /= napwkqp(ikt) ) call rxii('rotevec: napw_in /= napw(ikt)',napw_in,napwkqp(ikt))
         igloop: do ig = 1,napw_in
            qpg  = q + matmul(qlat(:,:),igv2qp(:,ig,ikt))! q+G
            qpgr = matmul(symops(:,:,igg),qpg)           ! rotated. g(q+G) = q'+G' = qpgr
            nnn  = nint(matmul(platt,qpgr-qtarget))      ! integer representation for G'= qpgr - q'
            ig2 = igv2revqp(nnn(1),nnn(2),nnn(3),ikt2)   ! get index of G'
            if(ig2>=999999) call rx('rotwave: q+G rotation error. (We have to set PWmode=11 for symgrpAF)')
            evecout(nlmto+ig2,:)= evec(nlmto+ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
         enddo igloop
      endif APWpart
   endsubroutine rotevec
   subroutine rotmto(qin,nbloch,nband,norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, &!Rotation of a function of ProductBasis or MT part of wave function. !    ! (then MTO itself is rotated. a little diffent from here).
      sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas,  cphiin, &
      cphiout)
      implicit none
      intent(in)::    qin,nbloch,nband,norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, &
         sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas,  cphiin
      intent(out):: &
         cphiout
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
   endsubroutine rotmto
   subroutine rotmto2(qin,nbloch,nband,norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, & !== Similar with rotmto but just for rotation matrix ==
      sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas, &
      zrotm)
      implicit none
      intent(in)::       qin,nbloch,nband, &
         norbt,ibas_tbl,l_tbl,k_tbl,offset_tbl,offset_rev_tbl,max_ibas_tbl,max_l_tbl,max_k_tbl, &
         sym,shtvg,dlmm,lxxa,miat,tiat,igxt,nbas
      intent(inout):: &
         zrotm
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
   endsubroutine rotmto2
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
         !       ig2 = findloc(ngvecp,val=nnn)
         geigenout(ig2,:)= geigenin(ig,:) * exp( -img2pi*sum(qpgr*shtvg) )
      enddo
   endsubroutine rotipw
   subroutine rotipw2(qin,qtarget,ngp,nband,platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx, zrotm)
      implicit none
      intent(in)::     qin,qtarget,ngp,nband,platt,qlat,sym,ngvecp,ngvecprev,shtvg,igxt,imx
      intent(out)::                                                                          zrotm
      real(8):: sym(3,3),qlat(3,3),platt(3,3),shtvg(3),qin(3),qpg(3),qpgr(3),qtarget(3)
      integer:: ngp,imx,nband, ngvecp(3,ngp), ngvecprev(-imx:imx,-imx:imx,-imx:imx),ig,ig2,nnn(3),igxt
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
   endsubroutine rotipw2
endmodule m_rotwave
