# matri element m_zmel

```
!! ----------------------------------------
cold  ntqxx--->nqmax
cold  nbmax -->nmmax
!!note: For usual correlation mode, I think nctot=0
!!note: For self-energy mode;   we calculate <iq1|\Sigma |iq2> , where iq1 and iq2 are in nqmax.
!!       nstate = nctot+nmmax
!!       allocate(zmelt(MPB,  intermediate phi nstate,  external state phi ntqxx))
!!       zmelt= < MPB     phi   | phi   > 
!!               <rkvec q-rkvec  |  q    >
 !                      cphim    | cphiq 
!                       ispm     | ispq
!            nctot+  nmini:nmmax | ncc + nqini:ntqxx
!                    middle state| end state
!
!!--- For dielectric funciton, we use irot=1 kvec=rkvec=q. We calulate \chi(q).
!!              q      rkvec     | q + rkvec  
!                    nkmin:nkmax | nkqmin:nkqmax
!                   (we fix nkmin=1)
!           or
!              nt0=nkmax-nkmin+1 | ntp0=nkqmax-nkqmin+1
!                      1:nt0     | 1:ntp0 
!                         occ    | unocc     
!                      (cphi_k   | cphi_kq !in x0kf)
!                    middle state| end state
!
!! NOTE: dimension
!!   nmtot = nctot+ nmmax-mnini+1
!!   nqtot = ncc  + ntqxx-nqini+1
!!   <q 1:ngb,      q-rkvec, 1:nmtot | rkvec, 1:nqtot>
!!   <end state,       middle state  |  MPB          >
!!   rkvec =mutmul(symops(:,:,irot),kvec)
!! ----------------------------------------
```


# Space group rotation of ecalj

* basic routines for spherical harmonics rotation.
/home/takao/ecalj/SRC/subroutines/rdpp.f90:      call rotcg(nl-1,symope,ngrp,cgr)
/home/takao/ecalj/SRC/subroutines/rotcg.f90:     call rotdlmm(symops, ng, 2*lmxax+1,dlmm)
rotcg --> rotdlmm      (rotation matrix of Ylm)
m_zmel->rdpp --> rotcg (rotated CG coefficients)

* eigenfunction rotation on the PMT basis
ecalj/SRC/main/hsfp0_sc.f90:   call rotwvigg(igrp, q(:,iqxx), q(:,iqxx), nhdim, ...)

* m_zmel matrix elements
  mt part  readcphif       -> readeigen ->rotmto
  ipw part drvmelp melpl2  -> readgeigf ->rotipw
  mptauf_zmel, ppbafp_v2->cgr --->ppbir
    * eigenfunction (MT+IPW expansion) rotation
  readeigen->rotipw,rotmto
 /home/takao/ecalj/SRC/subroutines/readeigen.f90: call rotipw(qtt(:,iqq),
 /home/takao/ecalj/SRC/subroutines/readeigen.f90: call rotmto(qtt(:,iqq),cphifr,ldim2,nband,
    * gvector rotation  rotgvec
 (the historical `ppbafp.fal.F` was merged into `ecalj/SRC/subroutines/m_zmel.f90`; rotgvec is invoked from there now.)

(Wannier-based matrix elements are implemented in `SRC/subroutines/hwmatK_MPI.f90`
— in principle only `readeigenW` and `cphieigW` differ from the standard path.)

* MPB rotation. symmetrization of x0
x0kf->rotMPB2->rotmto2,rotipw2
/home/takao/ecalj/SRC/subroutines/x0kf_v4h.f90:      call rotMPB2(nbloch,ngb,q,ig,itimer,ginv,zrotm)
/home/takao/ecalj/SRC/subroutines/m_rotMPB.f90:      call rotmto2(qin,nbloch,ngbb,
/home/takao/ecalj/SRC/subroutines/m_rotMPB.f90:      call rotipw2(qin,qout,ngcx,ngbb,
/home/takao/ecalj/SRC/subroutines/m_rotMPB.f90:      call rotdlmm(symops,ngrp,nl,dlmm)

* Self-energy Sigma rotation. Read sigm.* which is expanded by MTO basis only (neglecting IPW components)
/home/takao/ecalj/SRC/subroutines/rdsigm2.f90:   call rotsig(qin,q1,ndimh,napw_in,ldim,hq,gfbz(i1,i2,i3,:,:),ierr,iaf)






---------------------------
Old codes

/home/takao/ecalj/SRC/subroutines/m_hamindex.f90:    call rotdlmm(symops,ngrp, nl, dlmm)
/home/takao/ecalj/SRC/subroutines/m_hamindex.f90:    call rotdlmm(symops,ngrp, nl, dlmm)
/home/takao/ecalj/SRC/subroutines/m_q0p.f90:     call rotcg(lmxax,(/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/),1,cg)
/home/takao/ecalj/SRC/subroutines/mptauof.f90:   call rotcg(lmxax,(/1d0,0d0,0d0,0d0,1d0,0d0,0d
/home/takao/ecalj/SRC/wanniergw/hpsig.f90:         call rotcg(nl-1,symope,ngrpx,cg)
/home/takao/ecalj/SRC/wanniergw/hpsig_MPI.f90:      call rotcg(nl-1,symope,ngrpx,cg)
/home/takao/ecalj/SRC/wanniergw/huumat.f90:        call rotcg(nl-1,symope,ngrpx,cg)
/home/takao/ecalj/SRC/wanniergw/huumat_MPI.f90:     call rotcg(nl-1,symope,ngrpx,cg)
/home/takao/ecalj/SRC/wanniergw/hwmatK.f90:c        call rotgvec(symgg(:,:,irot), nqibz, 
/home/takao/ecalj/SRC/wanniergw/hwmatK_MPI.f90:c    call rotgvec(symgg(:,:,irot), nqibz, 


/home/takao/ecalj/SRC/subroutines/a2rotm.f90:      call rotma(phi,theta,angle,rotj)
/home/takao/ecalj/SRC/subroutines/chkdmu.f90:      call rotycs(-1,vorb,nbas,nsp,lmaxu,sspec,ssite,lldau)
/home/takao/ecalj/SRC/subroutines/sudmtu.f90:      call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,sspec,ssite,lldau)
/home/takao/ecalj/SRC/subroutines/mkplat.f90:        call rotmat(-1,.false.,nrot(m),mat,vecg)
/home/takao/ecalj/SRC/subroutines/rothrm.f90:        call rothrm(2,ndimh,iprmb,rotm,1,nbas,ndimh,uz,uz)
/home/takao/ecalj/SRC/subroutines/rothrm.f90:        call rothph(02,qpr,delT,ndimh,iprmb,1,nbas,ndimh,sq1)
/home/takao/ecalj/SRC/subroutines/rsmsym.f90:          call rotpnt(v,rv,g(1,ig))
/home/takao/ecalj/SRC/subroutines/symdmu.f90:        call rotspu(0,1,1,eula,1,u(1,1,ig))
/home/takao/ecalj/SRC/subroutines/symiax.f90:          call rotpnt(v,rv,g(1,ig))
/home/takao/ecalj/SRC/wanniergw/hmaxloc.f90:      call rot_hmnk(umnk,eunk,
```