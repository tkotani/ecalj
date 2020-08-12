      subroutine gwcphi(ssite,sspec,isp,nsp,nlmax,ndham,nev,nbas,
     .ipb,lmxax,nlindx,ndima,ppnl,aus,cphi,cphin)

      use m_struc_def  !Cgetarg

C- Project (phi,phidot) onto MT sphere, Kotani's GW conventions
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: p pz lmxa
Ci     Stored:    *
Ci     Passed to: *
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmax :leading dimension of aus
Ci   ndham :dimensions aus
Ci   nev   :number of eigenvectors to accumulate cphi
Ci   nbas  :number of sites in list
Ci   ipb   :index to true basis (excluding floating orbitals)
Ci         :given site index including those orbitals
Ci   lmxax :dimensions nlindx
Ci   nlindx:offset that set index in cphi for element (ipqn,l,ib)
Ci   ndima :leading dimension of cphi
Ci   ppnl  :NMTO potential parameters; see eg potpus.f
Ci         :This routine uses:
Ci         :(2) = s00 = <phi | phi>
Ci         :(7) = s11 = <phidot | phidot>
Ci         :(8)  = szz = <gz|gz> for local orbital
Ci         :(9)  = overlap <phi|gz> = s0z
Ci         :(10) = overlap <phidot|gz> = s1z
Ci   aus   :values of (phi,phidot) at MT sphere boundary; see makusq
Co Outputs
Co   cphi :coefficients to phi,phidot,phiz following Kotani conventions
Co        :cphi(ichan,iv) :
Co           ichan = orbital channel, i.e. one of (phi,phidot or phiz)
Co           for a given site, l, and m; see nlindx.
Co           iv = eigenvector
Co   cphin:diagonal matrix elements, one for each eigenvector.
Co        :cphin(1,iv) = <cphi(iv) | overlap | cphi(iv)>
Co        :cphin(2,iv) = <cphi(iv) | overlap-from-phi-only | cphi(iv)>
Cr Remarks
Cr   gwcphi converts the projection of an eigenvector into coefficients
Cr   of (phi,phidot) pairs, or with local orbitals, (phi,phidot,phiz)
Cr   into Kotani's conventions.  The overlap matrix between
Cr   the original functions is (see ppnl above)
Cr
Cr             (s00   0      s0z  )
Cr      S =  = ( 0   s11     s1z  )
Cr             (s0z  s1z     szz  )
Cr
Cr *Linear transformation of augmentation orbitals.  This code was
Cr  originally designed to make a linear transformation of atomic
Cr  orbitals (phi,phidot,phiz) to an orthonormal set, and generate
Cr  cphi for the orthonormalized set.  Now it generates cphi for the
Cr  original (non-orthonormal) orbitals.  To compute the sphere
Cr  contribution to the normalization, the sphere contribution to matrix
Cr  elements between eigenvectors is returned in cphin.
Cr
Cr  Vector aus corresponds to coefficients of augmented functions
Cr  {phi,phidot,gz}, which are not orthnormal.  Let us choose a linear
Cr  transformation L that transforms new functions {f1~..fn~} back to
Cr  the original functions {f1...fn}.  L^-1 transforms {f1...fn} to
Cr  {f1~..fn~}.
Cr
Cr  A sum of functions sum_i C_i f_i can be expressed as
Cr  sum_i C~_i f~_i by provided C~ satisfy
Cr    sum_j C_j f_j = sum_i C~_i f~_i = sum_ij C~_i L^-1_ij f_j
Cr                  = sum_j (sum_i L+^-1_ji C~_i) f_j
Cr  Writing as vectors
Cr    C+ f> = C+ L L^-1 f> = (L+ C)+ (L^-1 f>)
Cr  Therefore
Cr    C~ = L+ C  and  f~> = (L^-1 f>) and then C+ f> = (C~)+ f~>
Cr
Cr  If S_ij is the overlap in the old basis, in the new basis it is
Cr    < f~_i f~_j> = sum_km L^-1_ik S_km L^-1_mj
Cr             S~  = L^-1 S L^-1+
Cr             S   = L S~ L+
Cr
Cr *Choice of L to orthormal basis.  Then S~ = 1 and
Cr    S = L L+
Cr
Cr  1. No local orbital; basis consists of orthogonal (phi,phidot)
Cr     See above for overlap matrix.
Cr
Cr          (sqrt(s00)  0          )          (1/sqrt(s00)  0          )
Cr     L  = (                      )   L^-1 = (                        )
Cr          ( 0          sqrt(s11) )          ( 0          1/sqrt(s11) )
Cr
Cr  2. Local orbital; basis consists of (phi,phidot,gz)
Cr     See above for overlap matrix.
Cr
Cr     Let D = sqrt(szz - s0z^2/s00 - s1z^2/s11)
Cr
Cr             (sqrt(s00)          0            0)
Cr             (                                 )
Cr      L    = ( 0              sqrt(s11)       0)
Cr             (                                 )
Cr             (s0z/sqrt(s00)  s1z/sqrt(s11)    D)
Cr
Cr
Cr             (1/sqrt(s00)         0          0 )
Cr             (                                 )
Cr      L^-1 = ( 0              1/sqrt(s11)    0 )
Cr             (                                 )
Cr             (-s0z/s00/D     -s1z/s11/D     1/D)
Cr
Cu Updates
Cu    5 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   25 Apr 02 Returns cphi as coefficients to original nonlocal
Cu             orbitals, and cphin as matrix elements of evecs.
Cu             Altered argument list.
Cu   19 Feb 02 Extended to local orbitals.
Cu   28 Mar 01 written by MvS.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer isp,nsp,nlmax,ndham,nbas,nev,lmxax,ndima,ipb(nbas),
     .nlindx(3,0:lmxax,nbas)
      integer n0,nppn
      parameter (n0=10,nppn=12)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

      double precision ppnl(nppn,n0,nsp,*),cphin(2,nev)
      double complex aus(nlmax,ndham,3,nsp,*),cphi(ndima,nev)
C ... Local parameters
      double precision lmat(3,3),pnu(n0,2),pnz(n0,2)
      integer lmxa,ichan,ib,is,igetss,iv,ilm,l,im,k,ia
      double precision s00,s11,szz,s0z,s1z,D
      double complex au,as,az,sqrsz(3)

      call dpzero(lmat,9)
      call dpzero(cphin,2*nev)
C     ichan0 = 0
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)

        ia = ipb(ib)

        i_copy_size=size(sspec(is)%p)
        call dcopy(i_copy_size,sspec(is)%p,1,pnu,1)
        i_copy_size=size(sspec(is)%pz)
        call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)
        lmxa=sspec(is)%lmxa

        if (lmxa .eq. -1) goto 10

        do  iv = 1, nev
C         ichan = ichan0
          ilm = 0
          do  l = 0, lmxa

            k = l+1

            s00 = ppnl(2,k,isp,ib)
            s11 = ppnl(7,k,isp,ib)
            lmat(1,1) = sqrt(s00)
            lmat(2,2) = sqrt(s11)
            if (pnz(k,1) .ne. 0) then
              szz = ppnl(8,k,isp,ib)
              s0z = ppnl(9,k,isp,ib)
              s1z = ppnl(10,k,isp,ib)
              D = sqrt(szz - s0z**2/s00 - s1z**2/s11)
              lmat(3,1) = s0z/sqrt(s00)
              lmat(3,2) = s1z/sqrt(s11)
              lmat(3,3) = D
            else
              lmat(3,1) = 0
              lmat(3,2) = 0
              lmat(3,3) = 0
            endif

            if (.not. (pnz(k,1) .eq. 0 .eqv. nlindx(3,l,ia) .eq. -1)) then
              call rx('gwcphi: nlindx mismatch')
            endif

            do  im = 1, 2*l+1
              ilm = ilm+1
              au = aus(ilm,iv,1,isp,ib)
              as = aus(ilm,iv,2,isp,ib)
              az = aus(ilm,iv,3,isp,ib)

              sqrsz(1) = lmat(1,1)*au + lmat(3,1)*az
              sqrsz(2) = lmat(2,2)*as + lmat(3,2)*az
              sqrsz(3) =                lmat(3,3)*az
              cphin(1,iv) = cphin(1,iv) +
     .        dconjg(sqrsz(1))*sqrsz(1) +
     .        dconjg(sqrsz(2))*sqrsz(2) +
     .        dconjg(sqrsz(3))*sqrsz(3)
              cphin(2,iv) = cphin(2,iv) +
     .        dconjg(sqrsz(1))*sqrsz(1)

C            cphin(2,iv) = cphin(2,iv) +
C     .                    dconjg(lmat(1,1)*au)*lmat(1,1)*au

C           ichan = ichan+1
              ichan = nlindx(1,l,ia) + im
              cphi(ichan,iv) = au
              ichan = nlindx(2,l,ia) + im
              cphi(ichan,iv) = as
              if (nlindx(3,l,ia) .ge. 0) then
                ichan = nlindx(3,l,ia) + im
                cphi(ichan,iv) = az
              endif

            enddo
          enddo
        enddo
C       ichan0 = ichan

   10   continue
      enddo

      end subroutine gwcphi


