subroutine hambl(isp,qin, smpot,vconst,sv_p_osig,sv_p_otau,sv_p_oppi, h, s)
  use m_lmfinit,only: nbas , sspec=>v_sspec,nsp
  use m_igv2x,only: napw, igvapwin=>igv2x, ndimh
  use m_supot,only: k1,k2,k3
  use m_struc_def,only: s_rv1,s_cv1
  use m_lattic,only:plat=>lat_plat,qlat=>lat_qlat
  use m_hsibl,only:hsibl
  use m_shortn3_qlat,only: shortn3_qlat
  !- Make the LDA/GGA Hamiltonian and overlap matrix for one k-point. No SOC contribution
  !i   isp   :spin index
  !i   qin    :Bloch vector (k-point).
  !i   k1,k2,k3 dimensions of smpot
  !i   smpot :smooth potential on uniform mesh (mkpot.f)
  !i   vconst:additional constant potential
  !i   osig,otau,oppi  augmentation matrices
  !i   ndimh :dimension of hamiltonian and ovrlap h,s
  !o Outputs
  !o   h     :Hamiltonian matrix
  !o   s     :overlap matrix
  !u See update log via git log after 2009
  !u   04 Jul 08 (T. Kotani) New PW addition to basis
  !u   03 Feb 05 (A. Chantis) calculate hso
  !u    1 Sep 04 Adapted to handle complex ppi; S.O. put into ppi
  !u   25 Aug 04 modifications for extended local orbitals
  !u   15 Jul 04 (Chantis) Add Lz.Sz spin-orbit coupling
  !u   10 Jan 03 Remove addition from hambl.  See hambls.f
  !u   10 Jan 03 put sigma back into Bloch transform only
  !u   14 Aug 02 Added overlap-only option and option for orthog sigm
  !u   20 Jul 02 Can add Bloch transform of sigma matrix to ham
  !u   18 May 00 Adapted from nfp mk_hamiltonian.f
  !     ----------------------------------------------------------------------
  implicit none
  integer:: mode,igvapw(3,napw),isp
  integer:: mode1,ipr,inn(3),ig,i
  type(s_cv1) :: sv_p_oppi(3,nbas)
  type(s_rv1) :: sv_p_otau(3,nbas)
  type(s_rv1) :: sv_p_osig(3,nbas)
  real(8):: q(3), vconst,vavg, qin(3),qq(3) 
  complex(8):: smpot(k1,k2,k3,nsp), h(ndimh,ndimh),s(ndimh,ndimh)
  call tcn('hambl')
  !xxx  qin=qinin+ [1d-6,2d-6,3d-6] !trick to avoid degeneracy oscillation.
  
  ! Hamiltonian h and overlap matrix s is invariant for qlat shift as q->q+G_0
  ! as long as we keep the set {G} and {G+G_0} invariant.
  !! input and output
  !!   qpg(ig) = tpiba * ( qin + matmul(qlat,igapwin(1:3,ig))) for h,o,hso
  !! internal
  !!   qpg(ig) = tpiba * ( q  + matmul(qlat,igapw(1:3,ig)))
  !!   NOTE: both qpg are the same for given ig.
  !!   qlat*igapw = qlat*igqwin + (qin-q) ---> igvapw = igvapwin + matmul(qlatinv,qin-q)
  
  !q=qin !when q is not shortned (but q=qin works well. So I don't know shortn3 is needed or not).
  qq= matmul(transpose(plat),qin) !qq:fractional on qlat, qin cartesian 
  call shortn3_qlat(qq)
  q=matmul(qlat,qq)
  inn=nint(matmul(transpose(plat),qin-q))
  do ig=1,napw
     igvapw(:,ig) = inn + igvapwin(:,ig) !a set {G} for qin is conver to {G} for q
  enddo
  h = 0d0 !Hamiltonian for the basis of MTO+APW
  s = 0d0 !Overlap matrix for the basis of MTO+APW
  call augmbl(isp,q, sv_p_osig,sv_p_otau,sv_p_oppi,ndimh, h,s)! Augmentation parts of h,s
  vavg = 0 ! vavg: optionally add average constant potential
  call smhsbl(vavg+vconst,       q,ndimh,napw,igvapw, h,s)
  call hsibl (k1,k2,k3,smpot,isp,q,ndimh,napw,igvapw, h)
  do i=1,ndimh
     h(i+1:ndimh,i)=dconjg(h(i,i+1:ndimh))
     s(i+1:ndimh,i)=dconjg(s(i,i+1:ndimh))
  enddo
  call tcx('hambl')
end subroutine hambl
