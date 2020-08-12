      subroutine rotycs(mode,a,nbas,nsp,lmaxu,sspec,ssite,lldau)

      use m_struc_def  !Cgetarg

C- Rotate matrix a from real to spherical harmonics
C for LDA+U objects densmat and vorb
C-------------------------------------
Ci mode =1 from real to spherical
Ci      -1 from spherical to real
Ci a matrix to be transformed a(m,m,isp,iblu)  could be vorb or densmat
Ci nbas : number of sites
Ci nsp  : number of spins
Ci lmaxu: lmax for U
Ci sspec: species info
Ci ssite: sites info
Ci lldau  :lldau(ib)=0 => no U on this site otherwise
Ci        :U on site ib with dmat in dmats(*,lldau(ib))
Co a rotated in place
Cr Remarks
Cr order of cubic harmonics ls (l-1)s,ms...1s 0 1c mc... (l-1)c lc
Cr order of spherical harmonics -l:l
Cr Yl-m=(Ylmc-iYlms)/sqrt(2)  Ylm=(-1)**m*conjg(Yl-m)
Cu Updates
Cu   18 Jan 06 A. Chantis changed rotation matrices in accordance with
Cu             the definition of real harmonics used in the rest of
Cu             the code (Hund's rules satisfied as indicated by orb. moment)
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   30 Apr 05 Lambrecht first created
C----------------------------------------------------------------
C     implicit none
C passed variables
      integer nbas,lldau(nbas),mode,lmaxu,nsp
      double complex a(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,*)
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)

C local variables
      integer ib,m,l,idu(4),lmxa,is,igetss,i,j,k,ll,isp,iblu
      double complex rot(2,2)
      double complex b(2,2),c(2,2),add
      double precision s2

      s2 = 1/dsqrt(2d0)
      iblu = 0
      do  ib = 1, nbas
        if (lldau(ib) .ne. 0) then
          is = int(ssite(ib)%spec)
          lmxa=sspec(is)%lmxa
          i_copy_size=size(sspec(is)%idu)
          call icopy(i_copy_size,sspec(is)%idu,1,idu,1)

          do  l = 0, min(lmxa,3)
            if (idu(l+1) .ne. 0) then
              iblu = iblu+1
              do  isp = 1, 2
                do  m = 1, l
                  b(1,1) = a(m,m,isp,iblu)
                  b(1,2) = a(m,-m,isp,iblu)
                  b(2,1) = a(-m,m,isp,iblu)
                  b(2,2) = a(-m,-m,isp,iblu)
                  if (mode .eq. -1) then
C rotation from spherical  to cubic basis
                    rot(1,1) = dcmplx(s2,0d0)
                    rot(1,2) = dcmplx(s2*(-1d0)**m,0d0)
C                   rot(2,1) = dcmplx(0d0,-s2)
C                   rot(2,2) = dcmplx(0d0,s2*(-1d0)**m)
                    rot(2,1) = dcmplx(0d0,s2)
                    rot(2,2) = dcmplx(0d0,-s2*(-1d0)**m)
                  elseif (mode .eq. 1) then
C rotation from cubic  to spherical  basis
                    rot(1,1) = dcmplx(s2,0d0)
C                   rot(1,2) = dcmplx(0d0,s2)
                    rot(1,2) = dcmplx(0d0,-s2)
                    rot(2,1) = dcmplx(s2*(-1d0)**m,0d0)
C                   rot(2,2) = dcmplx(0d0,-s2*(-1d0)**m)
                    rot(2,2) = dcmplx(0d0,s2*(-1d0)**m)
                  else
                    call rx('ROTYCS: mode must be 1 or -1')
                  endif
C calculate matrix product c=rot*b*rot^+
                  do  i = 1, 2
                    do  j = 1, 2
                      add = dcmplx(0d0,0d0)
                      do  k = 1, 2
                        do  ll = 1, 2
                          add = add + rot(i,k)*b(k,ll)*dconjg(rot(j,ll))
                        enddo
                      enddo
                      c(i,j) = add
                    enddo
                  enddo
C place c in appropriate place
                  a(m,m,isp,iblu) = c(1,1)
                  a(m,-m,isp,iblu) = c(1,2)
                  a(-m,m,isp,iblu) = c(2,1)
                  a(-m,-m,isp,iblu) = c(2,2)
                enddo
              enddo
            endif
          enddo
        endif
      enddo
      end subroutine rotycs


