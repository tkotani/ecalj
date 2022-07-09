subroutine gwcphi(sspec,isp,nsp,nlmax,ndham,nev,nbas,ipb,lmxax,nlindx,ndima,ppnl,aus,  cphi,cphin)
  use m_struc_def
  use m_density,only: pnzall
  use m_lmfinit,only: ispec
  !- Project (phi,phidot) onto MT sphere, Kotani's GW conventions
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: p pz lmxa
  !i     Stored:    *
  !i     Passed to: *
  !i   isp   :current spin channel (1 or 2)
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nlmax :leading dimension of aus
  !i   ndham :dimensions aus
  !i   nev   :number of eigenvectors to accumulate cphi
  !i   nbas  :number of sites in list
  !i   ipb   :index to true basis (excluding floating orbitals)
  !i         :given site index including those orbitals
  !i   lmxax :dimensions nlindx
  !i   nlindx:offset that set index in cphi for element (ipqn,l,ib)
  !i   ndima :leading dimension of cphi
  !i   ppnl  :NMTO potential parameters; see eg potpus.f
  !i         :This routine uses:
  !i         :(2) = s00 = <phi | phi>
  !i         :(7) = s11 = <phidot | phidot>
  !i         :(8)  = szz = <gz|gz> for local orbital
  !i         :(9)  = overlap <phi|gz> = s0z
  !i         :(10) = overlap <phidot|gz> = s1z
  !i   aus   :values of (phi,phidot) at MT sphere boundary; see makusq
  !o Outputs
  !o   cphi :coefficients to phi,phidot,phiz following Kotani conventions
  !o        :cphi(ichan,iv) :
  !o           ichan = orbital channel, i.e. one of (phi,phidot or phiz)
  !o           for a given site, l, and m; see nlindx.
  !o           iv = eigenvector
  !o   cphin:diagonal matrix elements, one for each eigenvector.
  !o        :cphin(1,iv) = <cphi(iv) | overlap | cphi(iv)>
  !o        :cphin(2,iv) = <cphi(iv) | overlap-from-phi-only | cphi(iv)>
  !r Remarks
  !r   gwcphi converts the projection of an eigenvector into coefficients
  !r   of (phi,phidot) pairs, or with local orbitals, (phi,phidot,phiz)
  !r   into Kotani's conventions.  The overlap matrix between
  !r   the original functions is (see ppnl above)
  !r
  !r             (s00   0      s0z  )
  !r      S =  = ( 0   s11     s1z  )
  !r             (s0z  s1z     szz  )
  !r
  !r *Linear transformation of augmentation orbitals.  This code was
  !r  originally designed to make a linear transformation of atomic
  !r  orbitals (phi,phidot,phiz) to an orthonormal set, and generate
  !r  cphi for the orthonormalized set.  Now it generates cphi for the
  !r  original (non-orthonormal) orbitals.  To compute the sphere
  !r  contribution to the normalization, the sphere contribution to matrix
  !r  elements between eigenvectors is returned in cphin.
  !r
  !r  Vector aus corresponds to coefficients of augmented functions
  !r  {phi,phidot,gz}, which are not orthnormal.  Let us choose a linear
  !r  transformation L that transforms new functions {f1~..fn~} back to
  !r  the original functions {f1...fn}.  L^-1 transforms {f1...fn} to
  !r  {f1~..fn~}.
  !r
  !r  A sum of functions sum_i C_i f_i can be expressed as
  !r  sum_i C~_i f~_i by provided C~ satisfy
  !r    sum_j C_j f_j = sum_i C~_i f~_i = sum_ij C~_i L^-1_ij f_j
  !r                  = sum_j (sum_i L+^-1_ji C~_i) f_j
  !r  Writing as vectors
  !r    C+ f> = C+ L L^-1 f> = (L+ C)+ (L^-1 f>)
  !r  Therefore
  !r    C~ = L+ C  and  f~> = (L^-1 f>) and then C+ f> = (C~)+ f~>
  !r
  !r  If S_ij is the overlap in the old basis, in the new basis it is
  !r    < f~_i f~_j> = sum_km L^-1_ik S_km L^-1_mj
  !r             S~  = L^-1 S L^-1+
  !r             S   = L S~ L+
  !r
  !r *Choice of L to orthormal basis.  Then S~ = 1 and
  !r    S = L L+
  !r
  !r  1. No local orbital; basis consists of orthogonal (phi,phidot)
  !r     See above for overlap matrix.
  !r
  !r          (sqrt(s00)  0          )          (1/sqrt(s00)  0          )
  !r     L  = (                      )   L^-1 = (                        )
  !r          ( 0          sqrt(s11) )          ( 0          1/sqrt(s11) )
  !r
  !r  2. Local orbital; basis consists of (phi,phidot,gz)
  !r     See above for overlap matrix.
  !r
  !r     Let D = sqrt(szz - s0z^2/s00 - s1z^2/s11)
  !r
  !r             (sqrt(s00)          0            0)
  !r             (                                 )
  !r      L    = ( 0              sqrt(s11)       0)
  !r             (                                 )
  !r             (s0z/sqrt(s00)  s1z/sqrt(s11)    D)
  !r
  !r
  !r             (1/sqrt(s00)         0          0 )
  !r             (                                 )
  !r      L^-1 = ( 0              1/sqrt(s11)    0 )
  !r             (                                 )
  !r             (-s0z/s00/D     -s1z/s11/D     1/D)
  !r
  !u Updates
  !u    5 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   25 Apr 02 Returns cphi as coefficients to original nonlocal
  !u             orbitals, and cphin as matrix elements of evecs.
  !u             Altered argument list.
  !u   19 Feb 02 Extended to local orbitals.
  !u   28 Mar 01 written by MvS.
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: isp,nsp,nlmax,ndham,nbas,nev,lmxax,ndima,ipb(nbas), nlindx(3,0:lmxax,nbas)
  integer :: n0,nppn
  parameter (n0=10,nppn=12)
  type(s_spec)::sspec(*)
  double precision :: ppnl(nppn,n0,nsp,*),cphin(2,nev)
  double complex aus(nlmax,ndham,3,nsp,*),cphi(ndima,nev)
  double precision :: lmat(3,3),pnz(n0,2)
  integer :: lmxa,ichan,ib,is,igetss,iv,ilm,l,im,k,ia,i,ibas
  double precision :: s00,s11,szz,s0z,s1z,D
  double complex au,as,az,sqrsz(3)
  call dpzero(lmat,9)
  call dpzero(cphin,2*nev)
  do  ib = 1, nbas
     is = ispec(ib) 
     ia = ipb(ib)
     pnz(:,1:nsp)=pnzall(:,1:nsp,ib)
     lmxa=sspec(is)%lmxa
     if (lmxa == -1) goto 10
     do  iv = 1, nev
        ilm = 0
        do  l = 0, lmxa
           k = l+1
           s00 = ppnl(2,k,isp,ib)
           s11 = ppnl(7,k,isp,ib)
           lmat(1,1) = sqrt(s00)
           lmat(2,2) = sqrt(s11)
           if (pnz(k,1) /= 0) then
              szz = ppnl(8,k,isp,ib)
              s0z = ppnl(9,k,isp,ib)
              s1z = ppnl(10,k,isp,ib)
              D = sqrt(szz - s0z**2/s00 - s1z**2/s11)
              lmat(3,1) = s0z/sqrt(s00)
              lmat(3,2) = s1z/sqrt(s11)
              lmat(3,3) = D
           else
              lmat(3,:) = 0
           endif
           if(.NOT.(pnz(k,1) == 0 .eqv. nlindx(3,l,ia) == -1)) call rx('gwcphi: nlindx mismatch')
           do  im = 1, 2*l+1
              ilm = ilm+1
              au = aus(ilm,iv,1,isp,ib)
              as = aus(ilm,iv,2,isp,ib)
              az = aus(ilm,iv,3,isp,ib)
              sqrsz = [lmat(1,1)*au + lmat(3,1)*az, lmat(2,2)*as + lmat(3,2)*az, lmat(3,3)*az]
              cphin(1,iv) = cphin(1,iv) +   sum(dconjg(sqrsz(:))*sqrsz(:))
              cphin(2,iv) = cphin(2,iv) +       dconjg(sqrsz(1))*sqrsz(1)
              ichan = nlindx(1,l,ia) + im
              cphi(ichan,iv) = au
              ichan = nlindx(2,l,ia) + im
              cphi(ichan,iv) = as
              if (nlindx(3,l,ia) >= 0) then
                 ichan = nlindx(3,l,ia) + im
                 cphi(ichan,iv) = az
              endif

           enddo
        enddo
     enddo
10   continue
  enddo
end subroutine gwcphi


