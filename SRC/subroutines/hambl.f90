subroutine hambl(isp,qin,smpot,vconst,osig,otau,oppi, h,s)! Make LDA/GGA Hamiltonian and overlap matrix for a k-point. No SOC added.
  use m_lmfinit,only: nbas , sspec=>v_sspec,nsp
  use m_igv2x,only: napw, igvapwin=>igv2x, ndimh
  use m_supot,only: k1,k2,k3
  use m_struc_def,only: s_rv1,s_cv1,s_rv4
  use m_lattic,only:plat=>lat_plat,qlat=>lat_qlat
  use m_augmbl,only: augmbl
  use m_hsibl,only:hsibl
  !i   isp   :spin index
  !i   qin    :Bloch vector (k-point).
  !i   k1,k2,k3 dimensions of smpot
  !i   smpot :smooth potential on uniform mesh (mkpot.f)
  !i   vconst: additional constant potential
  !i   osig,otau,oppi:  augmentation matrices
  !i   ndimh :dimension of hamiltonian and ovrlap h,s
  !o Outputs
  !o   h     :Hamiltonian matrix Bloch sum of T_ij + V_ij.
  !    See Eq.(C.2) and (C.3) in JPSJ84,034702(2015)
  !    ecalj/Document/PAPERandPRESENTATION/KotaniKinoAkai2015FormulationPMT.pdf [1]
  !o   s     :overlap matrix
  !    See Eq.(C.1)           in [1]
  !r  qpg(ig) = tpiba * ( qin + matmul(qlat,igapwin(1:3,ig))), ig=1,napw
  implicit none
  integer:: mode,isp,i
  type(s_cv1) :: oppi(3,nbas)
  type(s_rv1) :: otau(3,nbas)
  type(s_rv4) :: osig(3,nbas)
  real(8):: vconst, qin(3)
  complex(8):: smpot(k1,k2,k3,nsp), h(ndimh,ndimh),s(ndimh,ndimh)
  call tcn('hambl')
  h = 0d0 !Hamiltonian for the basis of MTO+APW
  s = 0d0 !Overlap matrix for the basis of MTO+APW
  call augmbl(isp,qin,osig,otau,oppi,ndimh, h,s)! Augmentation parts of h,s
  !    product sum f structure constant C_akL^i in Eq.(C.1)-(C.2) in Ref.[1].
  call smhsbl(vconst,qin,ndimh,napw,igvapwin,          h,s)!Smooth and Constant potential parts.
  call hsibl(k1,k2,k3,smpot,isp,qin,ndimh,napw,igvapwin, h)!Smooth potential part, 1st term of (C.3) in Ref.[1]
  do i=1,ndimh
     h(i+1:ndimh,i)=dconjg(h(i,i+1:ndimh)) !RU part = LD part*
     s(i+1:ndimh,i)=dconjg(s(i,i+1:ndimh))
  enddo
  call tcx('hambl')
end subroutine hambl
