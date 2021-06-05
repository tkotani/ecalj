      integer function roth(opt,nl,nbas,bas,iat,offH,indxsh,istab,g,ag,
     .q,ldha,ldhb,h0,h)
C- Transform hamiltonian h0(g^-1 q) into h(q)
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit:
Ci         :0, rotate entire h
Ci         :1  rotate site-diagonal h
Ci         :2  no longer supported
Ci         :3  rotate rows and columns connected to a particular site
Ci             ib only.  Site is specified by iat
Ci         :4  same as case 3, except that h is returned untouched;
Ci             and rotated h is returned in h0; see Remarks
Ci         :5  same as case 1, except that h consists only of
Ci             diagonal blocks.
Ci         10s digit
Ci         :0  h is in real spherical harmonics
Ci         :Add 2 if h is in true (complex) spherical harmonics
Ci         :    See cb2sph.f and s2sph.f for definition of true
Ci         :    spherical harmonics.
Ci         :Add 4 if phase convention phi = q * [(g R_j + a) - R_i]
Ci         :  should be scaled by -1
Ci        100s digit distinguishes how complex arithmetic is handled
Ci           0: h,h0 have real, imaginary separated
Ci              h = h(ldha,ldhb,2), with h(*,*,1..2) = real..imag
Ci           1: h,h0 are in complex*16 format (see Bugs)
Ci              h = h(2,ldha,ldhb), with s(1,*) = real, s(2,*) = imag
Ci           2: h,h0 have real, imaginary separated by columns
Ci              h = h(ldha,2,ldhb), with h(*,1..2,*) = real..imag
Ci       1000s digit
Ci         :For nonzero values, entire row dimension is rotated,
Ci         :but row dim of h may only consist of some subblock:
Ci         :0 row dimension consists only of lower block
Ci         :1 row dimension consists only of lower block (same as 0)
Ci         :2 row dimension consists only of intermediate block
Ci         :3 row dimension consists only of higher block
Ci         :4 row & col dimension consists lower+intermediate block
Ci      10000s digit
Ci         :Similar to 1000s digit, but applicable to both row
Ci         :and column dimensions.
Ci         :2 rotate both row and col dim in intermediate block only
Ci         :  NB: implies 1000 digit = 2.
Ci         :3 rotate both row and col dim in higher block only
Ci         :  NB: implies 1000 digit = 3.
Ci   nl    :(global maximum l) + 1
Ci   nbas  :number of atoms in the basis (input)
Ci   bas   :basis vectors (input)
Ci   iat   :if 1s digit opt=0 or opt=1, not used
Ci         :if 1s digit opt=3, h(1..nbas,iat) h(iat,1..nbas) are rotated
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   indxsh:permutations ordering orbitals in lower, intermediate,
Ci          and higher blocks (makidx.f)
Ci   istab :table of site permutations for each group op (symtbl.f)
Ci          Site istab(i,ig) is transformed into site i by grp op ig.
Ci          Obtain by calling symtbl with mode=1.
Ci          NB: istab(1)<0 flags that current symop is consists of a proper
Ci               symmetry operation * inversion; see Remarks
Ci          Special case
Ci   g,ag  :space group operation; see Remarks
Ci   q     :qp for which h is sought.
Ci   rmat  :a work array of dimension at least nl**4.  If the 10s digit
Ci         :of opt is set, rmat must be input.
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h
Cio Inputs/Outputs
Cio   rmat :on input rotation matrix for this g if 10s digit opt is zero
Cio        :on ouput rotation matrix for this g
Cio   h    :on input, unrotated hamiltonian
Cio        :on output, rotated hamiltonian, except when 1s digit opt=4
Cio   h0   :complex array, used as a workspace except for opt=4.
Cio        :h0 has the following dimensions:
Cio        :opt = 0,1 same dimensions as h
Cio        :opt = 3   nlma*ldhb + ldha*nlmb, with
Cio        :          nlma = row dimension of site iat
Cio        :          nlmb = column dimension of site iat
Cio        :opt = 4   2*dim for opt=3 case = 2*(nlma*ldhb + ldha*nlmb)
Cio        :          This is a special case: rotated h is returned in h0
Cio        :          For this case:
Cio        :h0(1..nlma*ldhb) holds rows of h corresponding to site iat
Cio        :h0 is an array dimensioned h0(nlma,ldhb), with
Cio        :          nlma = row dimension of site iat
Cio        :h0(1+off2..ldha*nlmb) holds cols of h b belonging to site iat
Cio        :h0(1+off2) is an array dimensioned h0(ldha,nlmb), with
Cio        :          nlmb = column dimension of site iat and
Cio        :          off2 = nlma*ldhb*2
Co Outputs
Co    roth :returns  0 if no transpose required
Co         :returns -1 transpose required but protht did not make it
Co         :returns  1 protht took transpose
Cl Local variables
Cl   ldhr  :leading dimension of h
Cl   ldmpa :offset to 1st orbital in downfolding block (row index)
Cl   ldmpb :offset to 1st orbital in downfolding block (col index)
Cl   ofhi  :offset to imaginary part of h
Cl   ... For opt0 gt 1:
Cl   ibla  :index to offh defining block to be rotated (row index)
Cl   iblb  :index to offh defining block to be rotated (col index)
Cl   nlma  :number of orbitals connected to site ib (row index)
Cl   nlma1 :number of orbitals in lower block (row index)
Cl   nlma2 :number of orbitals in upper block (row index)
Cl   nlmb  :number of orbitals connected to site ib (col index)
Cl   nlmb1 :number of orbitals in lower block (col index)
Cl   nlmb2 :number of orbitals in upper block (col index)
Cl   ofa1  :offset to lower block belonging to ib (row index)
Cl   ofa2  :offset to intermediate block belonging to ib (row index)
Cl   ofb1  :offset to lower block belonging to ib (col index)
Cl   ofb2  :offset to intermediate block belonging to ib (col index)
Cl   lroti :set to .true. if rotation consists of proper symmetry op
Cl         :* inversion.  (Flagged by istab(1)<0).  See Remarks.
Cr Remarks
Cr   Input hamiltonian h is rotated.  g must correspond to a group
Cr   operation that maps the qp corresponding to h0 to the input q, i.e.
Cr   h0 must correspond to q0 = g^-1 q.  roth applies the rotation R
Cr   making  h(q) = R h0 (g^-1 q) R^-1.
Cr
Cr   The phase shifts calculated here correspond to the fact that
Cr   R' = g R_j + ag can differ from the corresponding basis position
Cr   R_i by a lattice translation vector.  The phase shifts below
Cr   adjust for the unwonted shift in ( g R_j + ag - R_i ).
Cr
Cr  *Extended symmetry operations.  In certain cases the inversion
Cr   operation may be added to the space group (see mksym.f).
Cr   (This is flagged when istab(1)<0).
Cr   In that case, (-g,ag) should be a proper group operation.
Cr   roth simply uses (-g,ag) for the symmetry operation, and returns
Cr   h (or h0) for that symop.  However, the correct h (or h0) to use
Cr   is the transpose of that generated by initially by roth (see aginv).
Cr   roth calls a kernel, protht, afterwards to generate the transpose.
Cr   But because of its internal complexity, protht cannot handle all
Cr   cases; roth returns one of (0;1;-1) indicating one of (no transpose
Cr   needed; transpose taken; transpose needed but not taken).
Cr   See protht for a further discussion.
Cu Updates
Cu   17 Jun 02 Added option for phase convention (10s digit opt)
Cu   05 Feb 01 Added code to handle special symmetry operations:
Cu             (proper symmetry operations + inversion)
Cu             Turned roth into an integer function call
Cu   14 Mar 00 some bug fixes in downfolding cases.
Cu    2 Mar 00 Added 10000s digit options to opt
Cu   29 Oct 99 generalized in several respects; argument list changed.
Cu   10 Dec 99 handles h consisting of lower+intermediate blocks.
Cb Bugs
Cb   kcplx=1 complex storage makes an internal transformation to kcplx=2
Cb   and thus requires additional O(N^2) operations
Cb   roth does not work properly for 1s digit opt = 1 and sil case
Cb
Cb   This routine is too complicated!  It handles too many cases.
C ----------------------------------------------------------------------
      implicit none
      integer nkap0,n0H
      parameter (nkap0=3,n0H=5)
      integer opt,ldha,ldhb,nl,nbas,istab(nbas),offh(n0H,nkap0,1),
     .iat,indxsh(1)
      double precision bas(3,*),g(3,3),ag(3),q(3),h0(1),h(ldha,ldhb)
C Local variables
      integer nlx,n0
      parameter (nlx=7,n0=10)
      logical lsaveh,lunitg,lroti
      integer ibla,iblb,jb1,k1,k2,kcplx,ld01,ld02,ld0b,ib1,
     .ld0r,ldh1,ldh2,ldhr,ldima,ldimb,ldmpa,ldmpb,nl2,nlma,nlma1,
     .nlma2,nlmb,nlmb1,nlmb2,of0i,ofa,ofa1,ofa2,ofb,ofb1,ofb2,off1,
     .off2,ofhi,opt0,opt1,protht
      integer norb,offli(n0*nkap0)
      double precision wk(nlx**4),rmat(nlx**4),ddot,xx

      call tcn('roth')
C     call zprm('starting h',2,h,ldha,ldha,ldhb)

      if (nl .gt. nlx) call rxi('increase nlx in roth, need',nl)
      nl2 = nl*nl
      opt0 = mod(opt,10)
      if (opt0 .eq. 2)
     .call rx('roth: 1s digit mode 2 no longer supported')
C     lsaveh is true when h is preserved
      lsaveh = opt0 .eq. 4
      if (lsaveh) opt0 = 3
      if (opt0 .eq. 5) opt0 = 1
      kcplx = mod(opt/100,10)
      if (kcplx .eq. 1) then
        kcplx = 2
        call ztoyy(h,ldha,ldhb,ldha,ldhb,1,2)
      endif

C ... Get some dimensioning parameters
      call cplxdm(kcplx,ldha,ldhb,ldh1,ldh2,ldhr,ofhi)
      opt1 = mod(mod(opt/10,10),4)
      call pvrotd(opt,offH,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,
     .ldimb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2)


C     lunitg is true if unit rotation matrix
      lunitg = ddot(9,g,1,g,1)-3 .lt. 1d-10 .and.
     .abs(g(1,1)-1) .lt. 1d-10 .and.
     .abs(g(2,2)-1) .lt. 1d-10 .and.
     .abs(g(3,3)-1) .lt. 1d-10

C --- Hold on to rows and columns of unrotated h, to restore later ---
      if (lsaveh) then
        off1 = (nlma*ldhb + ldha*nlmb)*2
C       Just copy to final position if unit rotation matrix
        if (lunitg) off1 = 0
        off2 = off1 + nlma*ldhb*2
C       Save rows of h(iat) into h0
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
        call ymscop(0,nlma1,ldhb,ldhr,ld0r,ofa1,0,0,0,h,ofhi,
     .  h0(1+off1),of0i)
        call ymscop(0,nlma2,ldhb,ldhr,ld0r,ofa2,0,nlma1,0,h,ofhi,
     .  h0(1+off1),of0i)
C       Save columns of h(iat) into h0
        call cplxdm(kcplx,ldha,nlmb,ld01,ld02,ld0r,of0i)
        call ymscop(0,ldha,nlmb1,ldhr,ld0r,0,ofb1,0,0,h,ofhi,
     .  h0(1+off2),of0i)
        call ymscop(0,ldha,nlmb2,ldhr,ld0r,0,ofb2,0,nlmb1,h,ofhi,
     .  h0(1+off2),of0i)
C       call yprm('h0 (r)',kcplx+2,h0,nlma*ldhb,nlma,nlma,ldhb)
C       call yprm('h0 (c)',kcplx+2,h0(1+off2),ldha*nlmb,ldha,ldha,nlmb)
      endif

C     Nothing to rotate; just quit
      if (lunitg) goto 99

C --- Setup for pseudorotation (rotation + inversion) ---
      lroti = .false.
      if (istab(1) .lt. 0) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
        lroti = .true.
      endif

C --- Set up rotation matrix ---
      if (opt1 .eq. 0) then
        call ylmrtg(nl2,g,rmat)
      else
        call ylmrtg(nl2,g,wk)
        call s2sph(kcplx,nl,nl,wk,nl2,nl2,nl2,nl2,rmat)
C       call yprm('rmat',kcplx+2,rmat,nl2*nl2,nl2,nl2,nl2)
      endif
C     call yprm('rmat',1,rmat,nl2*nl2,nl2,nl2,nl2)

C --- Rotate rows of h into h0: jbas rotated into ibas ---
C     opt0=3 rotates rows and columns connected to a particular site.
C     This is accomplished in two function calls.
C     Let jb=istab(ib), and ib1=iat with jb1=istab(ib1)
C     call 1: h0(ib1,k) are rotated from h(jb1,k) for k=1..ldimb
C     call 2: h0(ib,k) are rotated from h(jb,k) for ib=1..nbas, k in jb1
      if (opt0 .eq. 3) then
        ib1 = iat
        k1 = 1
        k2 = ldimb-ldmpb
        off1 = 0
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
        call iprmsb(10,ib1,ib1,nbas,offH,indxsh,xx,offli,norb)
        call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,q,
     .  rmat,k1,k2,ib1,ib1,0,ld0r,offli,0,h0,of0i,ldhr,ldhb,h,ofhi)
C       call yprm('R H (r)',kcplx+2,h0,of0i,nlma,nlma,ldhb)

C   ... Rotate into h0(ib,k) for ib=1..nbas, k in jb1
C       offset = size needed to hold all rows of h above
        off2 = off1 + nlma*ldhb*2
        jb1 = istab(iat)
C       k1,k2 encompass iblb-waves for site jb1
        k1  = offH(iblb,1,jb1)+1
        k2  = offH(iblb,1,jb1+1)
C       Leading dimension encompasses iblb-waves for site jb1
        ld0b = k2-k1+1
C       Leading dimension encompasses l+i-waves for site jb1
        if (ibla .eq. 4) ld0b = offH(ibla,1,jb1+1) - offH(ibla,1,jb1)
C       ofb is col offset in h0 so that R(h(*,jb1)) placed at h0 start
        ofb = k1-1
        call cplxdm(kcplx,ldha,ld0b,ld01,ld02,ld0r,of0i)
        call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,q,
ckino  off0 is array     .  rmat,k1,k2,1,nbas,0,ld0r,-1,ofb,h0(1+off2),of0i,ldhr,
     .  rmat,k1,k2,1,nbas,0,ld0r,(/-1/),ofb,h0(1+off2),of0i,ldhr,
     .  ldhb,h,ofhi)
C       Rotate il block, when h is composite
        if (ibla .eq. 4) then
C         k1,k2 encompass iwaves for site jb1
          k1  = offH(1,1,nbas+1) + offH(2,1,jb1)+1
          k2  = offH(1,1,nbas+1) + offH(2,1,jb1+1)
C         ofb is col offset subtracted from k1<=k<=k2 when rotating
C         h(*,k) -> h0(P(*),k-ofb).  ofb = first iwave col in h,
C         + size of l-block so as to preserve contents of l-block
          ofb = k1-1 - (offH(1,1,jb1+1)-offH(1,1,jb1))
          call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,
ckino off0 is array     .    q,rmat,k1,k2,1,nbas,0,ld0r,-1,ofb,h0(1+off2),of0i,ldhr,
     .    q,rmat,k1,k2,1,nbas,0,ld0r,(/-1/),ofb,h0(1+off2),of0i,ldhr,
     .    ldhb,h,ofhi)
        endif
C       call yprm('R H (c)',kcplx+2,h0(1+off2),ldha*ld0b,ldha,ldha,ld0b)
      else
        if (opt0 .eq. 1 .and. (ibla .eq. 2 .or. ibla .eq. 3)) then
          if (iblb .eq. 1) call rxi('roth not implemented for mode',opt)
        endif
        k1 = 1
        k2 = ldimb-ldmpb
        call cplxdm(kcplx,ldha,ldhb,ld01,ld02,ld0r,of0i)
        call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,q,
ckino off0 is array     .  rmat,k1,k2,1,nbas,0,ld0r,-1,0,h0,of0i,ldhr,ldhb,h,ofhi)
     .  rmat,k1,k2,1,nbas,0,ld0r,(/-1/),0,h0,of0i,ldhr,ldhb,h,ofhi)
C       call yprm('R H',kcplx+2,h0,ldha*ldhb,ldha,ldha,ldha)
C       call yprm('R H',kcplx+2,h0,ldha*ldhb,ldha,ldha,k2-k1+1)
      endif

C --- Rotate columns of h0 into h ---
C     opt0=3 rotates rows and columns connected to a particular site.
C     This is accomplished in two passes.
C     Let jb=istab(ib), and ib1=iat with jb1=istab(ib1)
C     pass 1: h(k,ib1) are rotated from h0(k,jb1) for k=1..ldima
C     pass 2: h(k,ib) are rotated from h0(k,jb) for ib=1..nbas, k in ib1
      if (opt0 .eq. 3) then
C   ... Rotate h0 into h(k,ib1) for k = 1 .. ldima
        ib1 = iat
C       jb1 = istab(ib1)
        k1 = 1
        k2 = offH(ibla,1,nbas+1)
        call cplxdm(kcplx,ldha,ld0b,ld01,ld02,ld0r,of0i)
        call iprmsb(10,ib1,ib1,nbas,offH,indxsh,xx,offli,norb)
        call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,q,
     .  rmat,k1,k2,ib1,ib1,0,ld0r,offli,0,h0(1+off2),of0i,ldhr,ldhb,h,
     .  ofhi)

C       call yprm('R H R^-1 (a)',kcplx+2,h,ofhi,ldha,ldhb,ldhb)
C   ... Rotate h0 into h(k,ib) for ib=1..nbas, k = in ib1
C       k1,k2 encompass l-waves for site ib1
        k1  = offH(ibla,1,ib1)+1
        k2  = offH(ibla,1,ib1+1)
        if (ibla .eq. 4) then
          k1  = offH(1,1,ib1)+1
          k2  = offH(1,1,ib1+1)
        endif
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
C       ofa is col offset in h0 so that R(h(*,jb1)) placed at h0 start
        ofa = k1-1
        call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,q,
ckino off0 is array     .  rmat,k1,k2,1,nbas,ib1,ld0r,-1,ofa,h0(1+off1),of0i,ldhr,ldhb,h,
     .  rmat,k1,k2,1,nbas,ib1,ld0r,(/-1/),ofa,h0(1+off1),of0i,ldhr,ldhb,h,
     .  ofhi)
C       Rotate il block, when h is composite
        if (ibla .eq. 4) then
C         k1,k2 encompass iwaves for site ib1
          k1  = offH(1,1,nbas+1) + offH(2,1,ib1)+1
          k2  = offH(1,1,nbas+1) + offH(2,1,ib1+1)
C         ofa is row offset subtracted from k1<=k<=k2 when rotating
C         h0(k-ofa,P(*)) -> h(k,*).  ofa = first iwave col in h,
C         + size of l-block so as to preserve contents of l-block
          ofa = k1-1 - (offH(1,1,ib1+1)-offH(1,1,ib1))
          call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,
ckino off0 is array     .    q,rmat,k1,k2,1,nbas,ib1,ld0r,-1,ofa,h0(1+off1),of0i,ldhr,
     .    q,rmat,k1,k2,1,nbas,ib1,ld0r,(/-1/),ofa,h0(1+off1),of0i,ldhr,
     .    ldhb,h,ofhi)
        endif
      else
        k1 = 1
        k2 = offH(ibla,1,nbas+1)
        call cplxdm(kcplx,ldha,ldhb,ld01,ld02,ld0r,of0i)
        call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,q,
ckino off0 is array     .  rmat,k1,k2,1,nbas,0,ld0r,-1,0,h0,of0i,ldhr,ldhb,h,ofhi)
     .  rmat,k1,k2,1,nbas,0,ld0r,(/-1/),0,h0,of0i,ldhr,ldhb,h,ofhi)
      endif
C     call yprm('R H R^-1',kcplx+2,h,ofhi,ldha,ldha,ldhb)

C --- Load h0 with rotated rows and columns of h; restore h ---
      if (lsaveh) then
        off1 = 0
        off2 = off1 + nlma*ldhb*2
C       Copy rows,columns of rotated h(iat) into h0
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
        call ymscop(0,nlma1,ldhb,ldhr,ld0r,ofa1,0,0,0,h,ofhi,
     .  h0(1+off1),of0i)
        call ymscop(0,nlma2,ldhb,ldhr,ld0r,ofa2,0,nlma1,0,h,ofhi,
     .  h0(1+off1),of0i)
        call cplxdm(kcplx,ldha,nlmb,ld01,ld02,ld0r,of0i)
        call ymscop(0,ldha,nlmb1,ldhr,ld0r,0,ofb1,0,0,h,ofhi,
     .  h0(1+off2),of0i)
        call ymscop(0,ldha,nlmb2,ldhr,ld0r,0,ofb2,0,nlmb1,h,ofhi,
     .  h0(1+off2),of0i)
        off1 = (nlma*ldhb + ldha*nlmb)*2
        off2 = off1 + nlma*ldhb*2
C       Restore rows,columns of h(iat)
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
        call ymscop(0,nlma1,ldhb,ld0r,ldhr,0,0,ofa1,0,h0(1+off1),
     .  of0i,h,ofhi)
        call ymscop(0,nlma2,ldhb,ld0r,ldhr,nlma1,0,ofa2,0,h0(1+off1),
     .  of0i,h,ofhi)
        call cplxdm(kcplx,ldha,nlmb,ld01,ld02,ld0r,of0i)
        call ymscop(0,ldha,nlmb1,ld0r,ldhr,0,0,0,ofb1,h0(1+off2),
     .  of0i,h,ofhi)
        call ymscop(0,ldha,nlmb2,ld0r,ldhr,0,nlmb1,0,ofb2,h0(1+off2),
     .  of0i,h,ofhi)
      endif

C --- Cleanup for pseudorotation (rotation + inversion) ---
      if (lroti) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
      endif

   99 continue
      kcplx = mod(opt/100,10)
      if (kcplx .eq. 1) then
        call ztoyy(h,ldha,ldhb,ldha,ldhb,2,1)
        if (lsaveh) then
          call ztoyy(h0,nlma,ldhb,nlma,ldhb,2,1)
          off2 = nlma*ldhb*2
          call ztoyy(h0(1+off2),ldha,nlmb,ldha,nlmb,2,1)
        endif
      endif
      call tcx('roth')

      roth = protht(opt,nbas,iat,offH,indxsh,istab,ldha,ldhb,h0,h)

C     call zprm('ending h',2,h,ldha,ldha,ldhb)
      end
      subroutine pvrotd(opt,offH,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,
     .ldimb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2)
C- Sets up dimensions for roth.
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :see roth
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   nbas  :size of basis
Ci   iat   :see roth
Co Outputs
Co   ibla  :index to offH that define which block of hamiltonian
Co         :is to be rotated in row dimension.
Co         :ibla=(1,2,3) -> (lower,intermediate,higher) block
Co         :ibla=4 -> lower+intermediate blocks
Co   iblb  :like ibla, but for column dimension of hamiltonian
Co   ldmpa :last orbital in block prior to current block, row dimension
Co   ldima :last orbital in current block, row dimension
Co   ldmpb :last orbital in block prior to current block, col dimension
Co   ldimb :last orbital in current block, row dimension
Co   The next quantities are defined only for 1s digit opt>=2
Co   (when only rows, cols associated with site iat are rotated).
Co   ofa1  :offset to lower block belonging to ib (row index)
Co   ofa2  :(if compound l+i block) offset to corresponding intermed. block
Co   nlma  :number of orbitals connected to site ib (row index)
Co   nlma1 :number of orbitals in lower block (row index)
Co   nlma2 :number of orbitals in upper block (row index)
Co   ofb1  :offset to lower block belonging to ib (col index)
Co   ofb2  :(if compound l+i block) offset to corresponding intermed. block
Co   nlmb  :number of orbitals connected to site ib (col index)
Co   nlmb1 :number of orbitals in lower block (col index)
Co   nlmb2 :number of orbitals in upper block (col index)
Cu Updates
Cu   05 Feb 01
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,ldimb,nkap0,n0H,
     .ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2
      parameter (nkap0=3,n0H=5)
      integer offH(n0H,nkap0,1)
C ... Local parameters
      integer isum

      ibla = max(mod(opt/1000,10),1)
      iblb = max(mod(opt/10000,10),1)
      if (iblb .gt. 1 .and. (ibla .ne. iblb .or. iblb .gt. 3))
     .call rxi('roth: inconsistent opt,',opt)

C ... Row dimensions and offsets
      if (ibla .lt. 4) then
        ldima = isum(ibla,offH(1,1,nbas+1),1)
        ldmpa = ldima - offH(ibla,1,nbas+1)
        ldmpb = 0
        ldimb = offH(1,1,nbas+1)
        if (iblb .gt. 1) then
          ldmpb = ldmpa
          ldimb = ldima
        endif

C       Dimensions and offsets for site iat
C       If mod(opt,10)<=1, these quantities are undefined
        if (mod(opt,10) .gt. 1) then
          ofa1 = offH(ibla,1,iat)
          ofa2 = ofa1
          nlma = offH(ibla,1,iat+1) - offH(ibla,1,iat)
          nlma1 = nlma
          nlma2 = 0

          ofb1 = offH(iblb,1,iat)
          ofb2 = ofb1
          nlmb = offH(iblb,1,iat+1) - offH(iblb,1,iat)
          nlmb1 = nlmb
          nlmb2 = 0
        endif

C ... Dimensions and offsets when h a composite of downfolding blocks
      elseif (ibla .eq. 4) then
        ldmpa = 0
        ldima = offH(1,1,nbas+1) + offH(2,1,nbas+1)
        ldmpb = 0
        ldimb = ldima

C       Dimensions and offsets for site iat
        if (mod(opt,10) .gt. 1) then
          ofa1 = offH(1,1,iat)
          ofa2 = offH(1,1,nbas+1) + offH(2,1,iat)
          nlma = offH(ibla,1,iat+1) - offH(ibla,1,iat)
          nlma1 = offH(1,1,iat+1) - offH(1,1,iat)
          nlma2 = nlma - nlma1

          ofb1 = ofa1
          ofb2 = ofa2
          nlmb = nlma
          nlmb1 = nlma1
          nlmb2 = nlma2
        endif

      else
        call rx1('roth: illegal value for opt (%i)',opt)
      endif
      end
      subroutine prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,
     .q,rmat,k1,k2,ib1,ib2,iskip,ld0a,off0,ofb,h0,of0i,ldha,ldhb,
     .h,ofhi)
C- Rotate h0 <- rmat h for some block of h
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit:
Ci         :0  rotate entire h
Ci         :1  rotate site-diagonal h
Ci         :5  same as case 1, but h contains only blocks diagonal in l
Ci         :10s digit
Ci         :0  h is in real spherical harmonics
Ci         :Add 2 if h is in true (complex) spherical harmonics
Ci         :    See cb2sph.f and s2sph.f for definition of true
Ci         :    spherical harmonics.
Ci         :Add 4 if phase convention phi = q * [(g R_j + a) - R_i]
Ci         :  should be scaled by -1
Ci   nl2   :leading dimension of rmat
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   ldmpa :offset to first orbital in downfolding block of h to rotate
Ci   ldima :points to last orbital in  downfolding block of h to rotate
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operation
Ci   ag    :translation part of space group
Ci   q     :wave number corresponding to Bloch sum
Ci   rmat  :YLM rotation matrix for this group operation
Ci   k1    :columns h(*,k1..k2) are rotated into h0(P(*)),k1..k2) ;
Ci   k2    :see Remarks
Ci   ib1   :subblocks of h for sites ib1..ib2 are rotated; see Remarks
Ci   ib2   : -//-
Ci   iskip :site ib = iskip is excluded from rotation
Ci   ld0a  :leading dimension of h0
Ci   off0  :Special indexing for h0 when off0(1)>=0 .
Ci         :(For default, use off0(1) = -1; off0 is then not used).
Ci         :Special indexing is used when the caller needs to specify
Ci         :the row ordering of h0, e.g. when a subblock of h is needed,
Ci         :and the mapping of the subblock is not into contiguous rows
Ci         :as occurs when concatenating orbitals of different downfolding
Ci         :subblocks.
Ci   ofb   :subtracted from the column offset when storing into h0.
Ci         :When only a portion of h is rotated (columns belonging to
Ci         :only one site), it conserves memory to store into first
Ci         :columns of h0, rather than the columns corresponding to site
Ci         :Then ofb = k1-1.
Ci   of0i  :offset to imaginary part of h0
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h (not needed)
Ci   h     :matrix to be rotated
Ci   ofhi  :offset to imaginary part of h
Co Outputs
Co   h0    :contains rotated subblock of h
Cl Local variables
Cl   korb  :l orbital column index.  Rotations are for blocks of
Cl          contiguous columns.  korb is an index that marks current
Cl          group of columns are being rotated.
Cr Remarks
Cr   This routine rotates rmat h, storing result into h0 for a block of
Cr   rows corresponding to sites ib1..ib2 and columns k1..k2.  Each site
Cr   ib is permuted into site istab(ib), thus permuting the rows of h
Cr   in addition to rotating by YLM(g).  Finally if matrix h corresponds
Cr   to the periodic part of Bloch hamiltonian, there is a phase factor
Cr     exp(i q (g R_j + a) - R_i)) where  (g R_j + ag) - R_i
Cr   is a lattice vector between sites i and j, and q is the wave number
Cr
Cr   NB: in the site-diagonal case (1s digit opt=1), columns k1 and k2
Cr   are not fixed but depend on the rows.  In that case k1 and k2
Cr   are set internally, and the input values are ignored.
Cr
Cr   prothl can handle only a single downfolding block in the row
Cr   dimension, because it assumes the orbitals corresponding to one
Cr   site lie in contiguous rows.  The columns are not permuted and
Cr   can be of any number.
Cu Updates
Cu   17 Jun 02 Added option for phase convention (10s digit opt)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,ld0a,ldha,ldhb,nl2,nbas,istab(nbas),off0(1),
     .indxsh(1),ldmpa,ldima,iskip,k1,k2,ofb,of0i,ofhi
      double precision bas(3,*),g(3,3),ag(3),q(3),
     .rmat(nl2,nl2),h0(ld0a,*),h(ldha,ldhb)
C ... Local parameters
      logical nphase
      integer n0,nkap0,i,ib1,ib2,ibas,korb,j,jbas,k,li,ncut,
     .nk,nlm,nlmi,norb,offh00,offr,opt0,opt1,opt4,iorb
      parameter (n0=10,nkap0=3)
      integer offi,offj,ltab(n0*nkap0),ktab(n0*nkap0),
     .offli(n0*nkap0),offlj(n0*nkap0)
      double precision twopi,sp,cosP,sinP,tbas(3),xx
#if SGI
      parameter (ncut=9)
#else
      parameter (ncut=9)
#endif

      if (k2 .lt. k1 .and. mod(opt,10) .ne. 1) return

C     call tcn('rothl')

      opt0 = mod(opt,10)
      opt1 = mod(mod(opt/10,10),4)
      nphase = mod(opt/10,10) .ge. 4
      opt4 = mod(opt/10000,10)
      twopi = 8d0*datan(1d0)

      offh00 = 0
      do  10  ibas = ib1, ib2
        jbas = istab(ibas)

C   ... Skip over diagonal site second pass
        if (ibas .eq. iskip) goto 10

C   ... Setup for multiplication by phase
C       tbas = (g R_j + a) - R_i; should be a lattice vector
        do  34  i = 1, 3
          tbas(i) = ag(i) - bas(i,ibas)
          do  36  j = 1, 3
            tbas(i) = tbas(i) + g(i,j)*bas(j,jbas)
   36     continue
   34   continue
C       exp(i q tbas)
        sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
        if (nphase) sp = -sp
        cosP = dcos(sp)
        sinP = dsin(sp)

C   ... Row, column offsets to h, h0
        korb = 0
        call orbl(ibas,ldmpa,ldima,indxsh,norb,ltab,ktab,xx,offli,nlm)
        call orbl(jbas,ldmpa,ldima,indxsh,norb,ltab,ktab,xx,offlj,nlm)
   12   if (nlm .eq. 0) goto 10

C       Calculate k1,k2 for site-diagonal part of h (opt0 = 1)
        if (opt0 .eq. 1) then
C         We cannot assume orbitals are contiguous.  Count the number of
C         contiguous orbitals.  korb will be index to last contiguous
C         orbital, and k1..k2 will encompass the offsets in h.  Must
C         be repeated until all blocks are exhausted.
          korb = korb+1
          k1 = offlj(korb)+1
          k2 = offlj(korb)
   16     continue
          k2 = k2 + 2*ltab(korb)+1
          if (korb .lt. norb) then
C             If contiguous, increment korb and increase k2
            if (offlj(korb+1) .eq. k2) then
              korb = korb + 1
              goto 16
            endif
          endif
C         For ii (or hh) case, h, h0 consists only of ii (or hh) block
          if (opt4 .gt. 1) then
            k1 = k1-ldmpa
            k2 = k2-ldmpa
          endif
        endif

C   ... For each l, do rmat h(jbas,lm) into h0(ibas).  Only rows mix
        offi = offh00
        iorb = 0
C       When h is in l-diagonal form, column offset starts at zero
C       and increments by 2*l+1 for each orbital in the list.
        if (opt0 .eq. 5) k2 = 0
        do  20  iorb = 1, norb
          li = ltab(iorb)
          if (off0(1) .lt. 0) then
            offi = offli(iorb) - ldmpa
          else
            offi = off0(iorb)
          endif
          offj =  offlj(iorb) - ldmpa
          nlmi = 2*li + 1
          offr = li**2

C    ... Debugging ... copy rmat into h
C        do  j = 1, nlmi
C        do  i = 1, nlmi
C          h0(i+offi,j+offj) = rmat(i+offr,j+offr)
C          h0(i+offi+of0i,j+offj) = 0
C        enddo
C        enddo
C        goto 20

C         When h is in l-diagonal form,  see comment above.
          if (opt0 .eq. 5) then
            k1 = k2+1
            k2 = k2+nlmi
          endif
          nk = k2-k1+1
          if (opt1 .gt. 1) then
            call yygemm('N','N',nlmi,nk,nlmi,1d0,rmat(offr+1,offr+1),
     .      rmat(offr+1,offr+1+nl2),nl2,h(offj+1,k1),h(offj+ofhi+1,k1)
     .      ,ldha,0d0,h0(offi+1,k1-ofb),h0(offi+of0i+1,k1-ofb),ld0a)
          elseif (nlmi .gt. ncut) then
            call dgemm('N','N',nlmi,nk,nlmi,1d0,rmat(offr+1,offr+1),nl2,
     .      h(offj+1,k1),ldha,0d0,h0(offi+1,k1-ofb),ld0a)
            call dgemm('N','N',nlmi,nk,nlmi,1d0,rmat(offr+1,offr+1),nl2,
     .      h(offj+ofhi+1,k1),ldha,0d0,h0(offi+of0i+1,k1-ofb),ld0a)
          else
            do  222  k = k1, k2
            do  221  i = 1, nlmi
              h0(i+offi,k-ofb) = 0
              h0(i+offi+of0i,k-ofb) = 0
            do  22  j = 1, nlmi
              h0(i+offi,k-ofb) =      h0(i+offi,k-ofb) +
     .            rmat(i+offr,j+offr)*h(j+offj,k)
              h0(i+of0i+offi,k-ofb) = h0(i+of0i+offi,k-ofb) +
     .            rmat(i+offr,j+offr)*h(j+offj+ofhi,k)
   22       continue
 221        continue
 222        continue
          endif

C     --- Multiply rows by phase ---
          if (dabs(sp) .gt. 1d-8) then
            do  24  i = 1, nlmi
C             print 456, offi+i,cosP,sinP
C  456        format(' phase',i4,2f15.10)
              call yscal(nk,cosP,sinP,
     .        h0(offi+i,k1-ofb),h0(offi+of0i+i,k1-ofb),ld0a)
   24       continue
          endif
   20   continue
C       Case column blocks are not contiguous
        if (opt0 .eq. 1 .and. korb .lt. norb) goto 12
C       Increment offh00 by the number of orbitals in this site
        offh00 = offh00 + nlm
   10 continue

C     call yprm('h0',kcplx+2,h0,of0i,ld0a,ld0a,k2-k1+1)

C     call tcx('rothl')

      end

      subroutine prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,
     .q,rmat,k1,k2,ib1,ib2,iskip,ld0a,off0,ofa,h0,of0i,ldha,ldhb,
     .h,ofhi)
C- Rotate h <- h0 rmat for some block of h0
C ----------------------------------------------------------------------
Ci   opt   :1s digit:
Ci         :0, rotate entire h
Ci         :1  rotate site-diagonal part of h
Ci         :10s digit
Ci         :0  h is in real spherical harmonics
Ci         :Add 2 if h is in true (complex) spherical harmonics
Ci         :    See cb2sph.f and s2sph.f for definition of true
Ci         :    spherical harmonics.
Ci   nl2   :leading dimension of rmat
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   ldmpb :offset to first orbital in downfolding block of h to rotate
Ci   ldimb :points to last orbital in  downfolding block of h to rotate
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operation
Ci   ag    :translation part of space group
Ci   q     :wave number corresponding to Bloch sum
Ci   rmat  :YLM rotation matrix for this group operation
Ci   k1    :rows h0(k1..k2,*) are rotated into h(k1..k2,P(*))
Ci   k2    :see Remarks
Ci   ib1   :subblocks of h for sites ib1..ib2 are rotated; see Remarks
Ci   ib2   : -//-
Ci   iskip :site ib = iskip is excluded from rotation
Ci   ld0a  :leading dimension of h0
Ci   off0  :Special indexing for h0 when off0(1)>=0 .
Ci         :(For default, use off0(1) = -1; off0 is then not used).
Ci         :Special indexing is used when the caller needs to specify
Ci         :the row ordering of h0, e.g. when a subblock of h is needed,
Ci         :and the mapping of the subblock is not into contiguous rows
Ci         :as occurs when concatenating orbitals of different downfolding
Ci         :subblocks.
Ci
Ci   cth   :1 assemble the columns from the different l channels from
Ci         :contiguous columns in h0, starting from the first column.
Ci         :Needed if h is partitioned into lower+higher blocks, so that
Ci         :consecutive l channels are not in contiguous columns.
Ci   ofa   :subtracted from the row offset when storing into h0.
Ci         :When only a portion of h is rotated (rows belonging to
Ci         :only one site), it conserves memory to store into first
Ci         :rows of h0, rather than the rows corresponding to site
Ci         :Then ofb = k1-1.
Co   h0    :contains rotated subblock of h
Ci   of0i  :offset to imaginary part of h0
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h (not needed)
Co Outputs
Co   h     :contains rotated subblock of h0; see Remarks
Cl Local variables
Cl   korb  :l orbital row index.  Rotations are for blocks of
Cl          contiguous rows.  korb is an index that marks current
Cl          group of rows are being rotated.
Cr Remarks
Cr   This routine rotates h0 rmat, storing result into h for a block of
Cr   columns corresponding to sites ib1..ib2 and rows k1..k2.  See
Cr   prothl for description of a similar rotation.
Cu Updates
Cu   17 Jun 02 Added option for phase convention (10s digit opt)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,ld0a,ldha,ldhb,nl2,nbas,istab(nbas),indxsh(1),ldmpb,
     .ldimb,iskip,k1,k2,ofa,of0i,ofhi,off0(1)
      double precision bas(3,*),g(3,3),ag(3),q(3),rmat(nl2,nl2),
     .h0(ld0a,*),h(ldha,ldhb)
C ... Local parameters
      logical nphase
      integer i,ib1,ib2,ibas,korb,j,jbas,k,li,ncut,nk,nlm,nlmi,norb,
     .offi,offj,offr,opt0,opt1,opt4,iorb,n0,nkap0
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offli(n0*nkap0),
     .offlj(n0*nkap0)
      double precision twopi,sp,cosP,sinP,tbas(3),xx
#if SGI
      parameter (ncut=9)
#else
      parameter (ncut=9)
#endif

      if (k2 .lt. k1) return
C     call tcn('rothr')
      opt0 = mod(opt,10)
      opt1 = mod(mod(opt/10,10),4)
      nphase = mod(opt/10,10) .ge. 4
      opt4 = mod(opt/10000,10)
      twopi = 8d0*datan(1d0)

      do  10  ibas = ib1, ib2
        jbas = istab(ibas)

C   ... Skip over diagonal site second pass
        if (ibas .eq. iskip) goto 10

C   ... Setup for multiplication by phase
C       tbas = (g R_j + a) - R_i; should be a lattice vector
        do  34  i = 1, 3
          tbas(i) = ag(i) - bas(i,ibas)
          do  36  j = 1, 3
            tbas(i) = tbas(i) + g(i,j)*bas(j,jbas)
   36     continue
   34   continue
C       exp(i q tbas)
        sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
        if (nphase) sp = -sp
        cosP = dcos(sp)
        sinP = dsin(sp)

C   ... Offset to first orbital in ibas,jbas.  Uses offi,nlm
        korb = 0
        call orbl(ibas,ldmpb,ldimb,indxsh,norb,ltab,ktab,xx,offli,nlm)
        call orbl(jbas,ldmpb,ldimb,indxsh,norb,ltab,ktab,xx,offlj,nlm)
   12   if (nlm .eq. 0) goto 10

C       Calculate k1,k2 for diagonal part of h (opt0 = 1)
        if (opt0 .eq. 1) then
C         We cannot assume orbitals are contiguous.  Count the number of
C         contiguous orbitals.  korb will be index to last contiguous
C         orbital, and k1..k2 will encompass the offsets in h.  Must
C         be repeated until all blocks are exhausted.
          korb = korb+1
          k1 = offli(korb)+1
          k2 = offli(korb)
   16     continue
          k2 = k2 + 2*ltab(korb)+1
          if (korb .lt. norb) then
C             If contiguous, increment korb and increase k2
            if (offli(korb+1) .eq. k2) then
              korb = korb + 1
              goto 16
            endif
          endif
C         For ii (or hh) case, h, h0 consist only of ii (or hh) block
          if (opt4 .gt. 1) then
            k1 = k1-ldmpb
            k2 = k2-ldmpb
          endif
        endif

C   ... For each l, do h0(jbas,lm) rmat into h(ibas). Only columns mix
C       When h is in l-diagonal form,  Assume k1,k2 assume role of
C       offi and offi starts at zero and increments by 2*l+1 for
C       each orbital in the list.
        offi = 0
        do  20  iorb = 1, norb
          li = ltab(iorb)
          nlmi = 2*li + 1
          offr = li**2
C         When h is in l-diagonal form,  see comment above.
          if (opt0 .eq. 5) then
            k1 = 1    + offli(iorb) - ldmpb
            k2 = nlmi + offli(iorb) - ldmpb
            offj = offi
          else
            offi = offli(iorb) - ldmpb
            if (off0(1) .lt. 0) then
              offj =  offlj(iorb) - ldmpb
            else
              offj =  off0(iorb)
            endif
          endif
          nk = k2-k1+1
          if (opt1 .gt. 1) then
            call yygemm('N','C',nk,nlmi,nlmi,1d0,h0(k1-ofa,offj+1),
     .      h0(k1+of0i-ofa,offj+1),ld0a,rmat(offr+1,offr+1),rmat(offr+
     .      1,offr+1+nl2),nl2,0d0,h(k1,offi+1),h(k1+ofhi,offi+1),ldha)
          elseif (nlmi .gt. ncut) then
            call dgemm('N','T',nk,nlmi,nlmi,1d0,h0(k1-ofa,offj+1),
     .      ld0a,rmat(offr+1,offr+1),nl2,0d0,h(k1,offi+1),ldha)
            call dgemm('N','T',nk,nlmi,nlmi,1d0,h0(k1+of0i-ofa,offj+1),
     .      ld0a,rmat(offr+1,offr+1),nl2,0d0,h(k1+ofhi,offi+1),ldha)
          else
            do  21  i = 1, nlmi
            do  21  k = k1, k2
              h(k,i+offi) = 0
              h(k+ofhi,i+offi) = 0
   21       continue
            do  22  i = 1, nlmi
            do  22  j = 1, nlmi
            do  22  k = k1, k2
              h(k,i+offi) = h(k,i+offi) +
     .            h0(k-ofa,j+offj)*rmat(i+offr,j+offr)
              h(k+ofhi,i+offi) = h(k+ofhi,i+offi) +
     .            h0(k+of0i-ofa,j+offj)*rmat(i+offr,j+offr)
   22       continue
          endif

C     --- Multiply columns by phase ---
          if (dabs(sp) .gt. 1d-8) then
            do  24  i = 1, nlmi
              call yscal(nk,cosP,-sinP,h(k1,offi+i),h(k1+ofhi,offi+i),1)
   24       continue
          endif

C         Increment offi in diagonal case
          offi = offi + nlmi
   20   continue
C       Case column blocks are not contiguous
        if (opt0 .eq. 1 .and. korb .lt. norb) goto 12
   10 continue

C     call tcx('rothr')

      end
      integer function protht(opt,nbas,iat,offH,iprmb,istab,ldha,ldhb,
     .h0,h)
C- Make appropriate transposes for pseudorotation, depending on mode
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :see roth
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   h0    :see roth
Ci   h     :small h = cnu-enu+sqrdel*S^beta*sqrdel
Ci   iat   :if opt=0 or opt=1, not used (1s digit)
Ci         :if opt=3, h(1..nbas,iat) h(iat,1..nbas) are rotated
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nbas  :size of basis
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h
Co Outputs
Co   h0    :some portion or all of h0 is transposed, depending on opt
Co   h     :some portion or all of h is transposed, depending on opt
Co   protht:returns  0 if no transpose required
Co         :returns -1 transpose required but protht did not make it
Co         :returns  1 protht took transpose
Cl Local variables
Ci  nrh    :spacing between adjacent row elements
Ci  nch    :spacing between adjacent column elements
Cr Remarks
Cr   This subroutine is intended as a 'cleanup' to roth to handle cases
Cr   when the inversion operation is not part of the space group but
Cr   is included anyway because of time reversal symmetry.
Cr   (This is flagged when istab(1)<0).  In such a case, roth generates
Cr   the transpose of h.  For the cases it is able to handle, protht
Cr   overwrites h (or part of it) with its transpose.  For it to do so,
Cr   protht requires rows, columns to correspond to the same orbitals.
Cr   it cannot handle opt=2000, 2003, 2004.  protht could be made to
Cr   handle opt=4003 but it hasn't been implemented.
Cb Bugs
Cb   kcplx<>0 not checked for mode 5 (kcplx=1 will not work).
Cu Updates
Cu   05 Feb 01 Created and checked using troth.
C ----------------------------------------------------------------------

      implicit none
      integer opt,ldha,ldhb,istab(1),iat,nbas,iprmb(1)
      integer nkap0,n0H
      parameter (nkap0=3,n0H=5)
      integer offH(n0H,nkap0,1)
      double precision h(ldha,*),h0(1)
C
      integer nrh,nch,ofhi,ldh1,ldh2,w,opt0,kcplx,ibla,iblb,
     .ldmpa,ldmpb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,
     .nlmb2,ldima,ldimb,off2,ndim

      protht = 0
      if (istab(1) .gt. 0) return
      protht = -1

C ... Get some dimensioning parameters
      opt0 = mod(opt,10)
      call pvrotd(opt,offH,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,
     .ldimb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2)
C     protht requires rows, columns to correspond to same orbitals
      if (ibla .ne. 4 .and. iblb .ne. ibla) return
C     if (ldima .ne. ldimb) return

      kcplx = mod(opt/100,10)
      call cplxdm(kcplx,ldha,ldhb,ldh1,ldh2,nch,ofhi)
      nrh = 1
      if (kcplx .eq. 1) nrh = 2

C ... rotate entire h
      if (opt0 .eq. 0) then
        ndim = ldima-ldmpa
        call ymtrns(0,h,nch,nrh,ofhi,w,w,w,w,1,ndim,1,ndim)
        protht = 1
C ... site-diagonal h
      elseif (opt0 .eq. 1) then
        call pvrots(opt0,nbas,ldha,ldhb,iprmb,ldmpa,ldima,nch,nrh,ofhi,h)
        protht = 0
C ... rows and columns connected with one site
      elseif (opt0 .eq. 3) then
        if (nlma.ne.nlmb .or. ldima.ne.ldimb .or. ofa1.ne.ofb1) return
        if (nlma2 .ne. 0) return
        ndim = ldima-ldmpa
        call ymtrns(0,h,nch,nrh,ofhi,w,w,w,w,ofa1+1,ofa1+nlma,1,ndim)
        protht = 1
      elseif (opt0 .eq. 4) then
        if (nlma.ne.nlmb .or. ldima.ne.ldimb .or. ldmpa.ne.ldmpb) return
        ldima = ldima - ldmpa
        ldimb = ldimb - ldmpb
        off2 = nlma*ldimb*2
        if (kcplx .ne. 0) then
          call ztoyy(h0,nlma,ldimb,nlma,ldimb,kcplx,0)
          call ztoyy(h0(1+off2),ldima,nlmb,ldima,nlmb,kcplx,0)
        endif
        call pvrott(nlma,ldima,h0,h0(1+off2))
        if (kcplx .ne. 0) then
          call ztoyy(h0,nlma,ldimb,nlma,ldimb,0,kcplx)
          call ztoyy(h0(1+off2),ldima,nlmb,ldima,nlmb,0,kcplx)
        endif
        protht = 1
C ... hh block
      elseif (opt0 .eq. 5) then
        if (ibla .ne. 3 .or. iblb .ne. 3) return
        call pvrots(opt0,nbas,ldha,ldhb,iprmb,ldmpa,ldima,nch,nrh,ofhi,h)
        protht = 0

      else
        return
      endif

      end
      subroutine pvrott(nlm,ndim,sr,sc)
      implicit none
      integer nlm,ndim
      double precision sr(nlm,ndim,2),sc(ndim,nlm,2)
      integer i,j
      double precision wkr,wki

      do  i = 1, nlm
        do  j = 1, ndim
          wkr = sr(i,j,1)
          wki = sr(i,j,2)
          sr(i,j,1) = sc(j,i,1)
          sr(i,j,2) = sc(j,i,2)
          sc(j,i,1) = wkr
          sc(j,i,2) = wki
        enddo
      enddo
      end

      subroutine pvrots(opt0,nbas,ldha,ldhb,iprmb,ldmpa,ldima,nch,nrh,
     .ofhi,h)
      implicit none
      integer nbas,opt0,ldha,ldhb,iprmb(1),ldmpa,ldima,nch,nrh,ofhi
      double precision h(ldha,ldhb)
C Local variables
      integer n0,nkap0,jb,k1,k2,norb,nlm,korb
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offlj(n0*nkap0)
      double precision xx

      do  jb = 1, nbas
        korb = 0
        call orbl(jb,ldmpa,ldima,iprmb,norb,ltab,ktab,xx,offlj,nlm)
        if (nlm .ne. 0) then
   12     continue
C       We cannot assume orbitals are contiguous.  Count the number of
C       contiguous orbitals.  korb will be index to last contiguous
C       orbital, and k1..k2 will encompass the offsets in h.  Must
C       be repeated until all blocks are exhausted.
          korb = korb+1
          k1 = offlj(korb)+1
          k2 = offlj(korb)
   16     continue
          k2 = k2 + 2*ltab(korb)+1
C         Look for contiguous blocks (not higher blocks)
          if (korb .lt. norb .and. opt0 .ne. 5) then
            if (offlj(korb+1) .eq. k2) then
              korb = korb + 1
              goto 16
            endif
          endif

C       For higher blocks h consists of hh block;
          if (opt0 .eq. 5) then
C         k1-ldmpa = offset rel to start of hh block
C         k1-offlj(1) = offset for start of this atom
            k1 = k1-offlj(1)
            k2 = k2-offlj(1)
            call ymtrns(0,h(k1-ldmpa+offlj(1),k1),nch,nrh,ofhi,xx,xx,xx,
     .      xx,k1,k2,k1,k2)
          else
            k1 = k1-ldmpa
            k2 = k2-ldmpa
            call ymtrns(0,h,nch,nrh,ofhi,xx,xx,xx,xx,k1,k2,k1,k2)
          endif

          if (korb .lt. norb) goto 12
        endif
      enddo

      end


