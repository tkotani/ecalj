      subroutine invbl(lasa,lov,nbas,nl,nsp,ldim,plat,nsite,iax,
     .indxsh,ik,nk,bk,strxk,hk,sk,strx,hrs,srs)
C- Inverse Bloch transform to extract real space hamiltonian and overlap
C ----------------------------------------------------------------------
Ci Inputs:  lasa: T for lmasa, F for TBE which tranposes the strux.
Ci          lov:  TBE only indicates non orthogonal TB.
Ci          nbas,nl,ldim,alat,pos,npr,iax,nk,hk,sk
Ci          strxk is actually sqr(del) S(k) sqr(del) --see hmltns.f
Co Outputs: strx; hrs, srs: H and O in real space formatted like strux.
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer niax
      parameter (niax=10)
      logical lasa,lov
      integer nbas,nl,ldim,nsp,nsite,
     .ik,nk,iax(niax,nsite),indxsh(1)
      double precision bk(3),plat(3,3),strxk(ldim,ldim,2),
     .hk(ldim,ldim,2),sk(ldim,ldim,2),
     .strx(nl**2,nl**2,nsite),hrs(nl**2,nl**2,nsite),
     .srs(nl**2,nl**2,nsite)
C Local Variables
      integer i,j,k,mxorb,isite,ib,ia,offb,offa,ipb,ipa,lmb,lma
      integer ifi,fopn,iprint,nglob
      double precision vec(3),TdotK,twopi,cosT,sinT,wt

C Extrinsic functions
      double precision dcos

C Intrinsic functions
      intrinsic datan, dcos, dsin
      mxorb = nglob('mxorb')
      twopi = 8*datan(1d0)

      call rxx(nsp.ne.1,'INVBL: not set up for spin')

      wt = 1d0 / nk
      ifi = fopn('RSH')
      if (bk(1) .eq. 0d0 .and. bk(2) .eq. 0d0 .and. bk(3) .eq. 0d0
     ..and. iprint() .ge. 30) then
        call awrit0('Hamiltonian at Gamma:',' ',120,ifi)
        do  i = 1, ldim
          write (ifi,1) (hk(i,j,1),j=1,ldim)
        enddo
        if (lov) then
          call awrit0('Overlap at Gamma:',' ',120,ifi)
          do  i = 1, ldim
            write (ifi,1) (sk(i,j,1),j=1,ldim)
          enddo
        endif
      endif
    1 format (128f10.6)

C --- loop over all pairs ---
      do  isite = 1, nsite
        TdotK = 0
        do  j = 1, 3
          do  k = 1, 3
            TdotK = TdotK + twopi*bk(j)*plat(j,k)*iax(2+k,isite)
          enddo
        enddo
        cosT = dcos(TdotK)
        sinT = dsin(TdotK)
        if (lasa) then
          ib = iax(1,isite)
          ia = iax(2,isite)
        else
          ia = iax(1,isite)
          ib = iax(2,isite)
        endif
        offb = mxorb * (ib-1)
        do  lmb = 1, mxorb
          offb = offb + 1
          ipb = indxsh(offb)
          offa = mxorb * (ia-1)
          do  lma = 1, mxorb
            offa = offa + 1
            ipa = indxsh(offa)
            if (ipa .le. ldim .and. ipb .le. ldim) then
              hrs(lma,lmb,isite) = hrs(lma,lmb,isite) +
     .        (hk(ipa,ipb,1)*cosT + hk(ipa,ipb,2)*sinT) * wt
              if (lov) then
                srs(lma,lmb,isite) = srs(lma,lmb,isite) +
     .          (sk(ipa,ipb,1)*cosT + sk(ipa,ipb,2)*sinT) * wt
              endif
              if (lasa) then
                strx(lma,lmb,isite) = strx(lma,lmb,isite) +
     .          (strxk(ipa,ipb,1)*cosT + strxk(ipa,ipb,2)*sinT) * wt
              endif
            endif
          enddo
        enddo
      enddo
      end

      subroutine wrirsh(lasa,lov,lgamma,nl,nsp,nclass,nbas,plat,
     .pos,alat,dclabl,ipc,pp,alpha,nsite,iax,srx,
     .hrs,srs,h1rs,h2rs,s2rs,Don)
C- Write out real space H, O and strux
C ----------------------------------------------------------------------
Ci Inputs: lasa: std ordering of harmonics, if F empirical TB ordering
Ci         lov:  non orthogonal (write O as well as H)
Ci         lgamma: gamma representation used in secmat
Ci         nl,nbas,plat,pos,alat,dclabl,ipc nsite,iax (pairc.f)
Ci         hrs, srs, srx: real space H, O and strux
Ci                        formatted as in asastr.f
Ci                        note that srx is actually sqr(del) S sqr(del)
Ci                        (see invbl)
Ci         Don: workspace of size nbas*nl**4
Co Outputs:
Co         h1rs, h2rs, s2rs: Hamiltonian and overlap to first
Co                                 and second order in the two centre
Co                                 approximation
Cr Remarks
Cr         This makes a table of the real space H and O matrix elements
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer niax
      parameter (niax=10)
      integer nl,nsp,nclass,nbas,nsite,iax(niax,nsite),ipc(nbas)
      double precision plat(3,3),pos(3,nbas),alat,dclabl(1),
     .hrs(nl**2,nl**2,nsite),srs(nl**2,nl**2,nsite),
     .srx(nl**2,nl**2,nsite),h1rs(nl**2,nl**2,nsite),
     .h2rs(nl**2,nl**2,nsite),s2rs(nl**2,nl**2,nsite),
     .pp(6,0:nl-1,nsp,1),alpha(1),Don(nl**2,nl**2,nbas)
      logical lasa,lgamma,lov
C Local Variables
      integer isite,ib,ia,icb,ica,ix,i,j,k,L,lma,lmb,ll,
     .ifi,pvec(3)
      integer fopn
      character*8 clabli,clablj,orb1(9),orb2(9)
      character*50 outs
      logical onsite
      double precision dx,vec(4),c,enu,del,p,gam,alp,ob,oa,cmenu,
     .F1a,F2a,F2b,G1a,G2a,G2b,G3a,G3b,
     .dsd,h1,h2,o2
C Extrinsic functions
      double precision oalpha

      data orb1 /'s','y','z','x','xy','yz','3z^2-1','zx','x^2-y^2'/
      data orb2 /'s','x','y','z','xy','yz','zx','x^2-y^2','3z^2-1'/


C Connecting vector ...
      dx(ix,ia,i,j,k) = pos(ix,ia) - pos(ix,ib) +
     .plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k
C Potential parameters
      enu(i,L)          = pp(1,ll(L),1,i)
      c(i,L)            = pp(2,ll(L),1,i)
      cmenu(i,L)        = pp(2,ll(L),1,i) - pp(1,ll(L),1,i)
      del(i,L)          = pp(3,ll(L),1,i)
      p(i,L)            = pp(4,ll(L),1,i)
      alp(i,L)          = pp(5,ll(L),1,i)
      gam(i,L)          = pp(6,ll(L),1,i)

      ifi = fopn('RSH')

      call wripp(nl,nsp,nclass,pp,dclabl,ifi)

      call clean0(nl**4*nsite,hrs)
      if (lov) call clean0(nl**4*nsite,srs)
      if (lasa) call clean0(nl**4*nsite,srx)
      if (lasa) then
        if (lgamma) then
          write (ifi,1)
        else
          write (ifi,2)
        endif
      else
        if (lov) then
          write (ifi,3)
        else
          write (ifi,4)
        endif
      endif
      if (lasa) then
C .. get onsite strux
        do  isite = 1, nsite
          ib = iax(1,isite)
          ia = iax(2,isite)
          if (ia .eq. ib) then
            do  lma = 1, nl**2
              do  lmb = 1, nl**2
                Don(lma,lmb,ia) = srx(lma,lma,isite)
              enddo
            enddo
          endif
        enddo
      endif
      do isite = 1, nsite
        if (lasa) then
          ib = iax(1,isite)
          ia = iax(2,isite)
        else
          ia = iax(1,isite)
          ib = iax(2,isite)
        endif
        pvec(1) = iax(3,isite)
        pvec(2) = iax(4,isite)
        pvec(3) = iax(5,isite)
        icb = ipc(ib)
        ica = ipc(ia)
        call r8tos8(dclabl(icb),clabli)
        call r8tos8(dclabl(ica),clablj)
        do  ix = 1, 3
          vec(ix) = dx(ix,ia,pvec(1),pvec(2),pvec(3))
        enddo
        vec(4) = alat*sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
        onsite = vec(4) .lt. 1d-6
        do  lmb = 1, nl**2
          do  lma = 1, nl**2
            if (dabs(hrs(lma,lmb,isite)) .gt. 1d-6) then
              outs = ' '
              call strcop(outs,clabli,8,' ',i)
              call awrit1('(%i)',outs(i+1:50),5,0,ib)
              j = i + 1
              call skp2bl(outs,50,j)
              if (lasa) then
                call strcop(outs(j+2:50),orb1(lma),8,' ',i)
              else
                call strcop(outs(j+2:50),orb2(lma),8,' ',i)
              endif
              j = j + 2 + i
              call strcop(outs(j:50),'--',2,' ',i)
              j = j + 2
              call strcop(outs(j+1:50),clablj,8,' ',i)
              i = i + j
              call awrit1('(%i)',outs(i+1:50),5,0,ia)
              j = i + 2
              call skp2bl(outs,50,j)
              i = j
              if (lasa) then
                call strcop(outs(i+2:50),orb1(lmb),8,' ',i)
              else
                call strcop(outs(i+2:50),orb2(lmb),8,' ',i)
              endif
              i = i + j + 2
              if (lasa) then
                dsd = srx(lma,lmb,isite)
                h1 = dsd
                if (onsite .and. lmb .eq. lma) then
                  h1 = c(ica,lma) + dsd
                endif
                oa = oalpha(enu(ica,lma),c(ica,lma),del(ica,lma),
     .          alp(ica,lma),gam(ica,lma))
                ob = oalpha(enu(icb,lmb),c(icb,lmb),del(icb,lmb),
     .          alp(icb,lmb),gam(icb,lmb))
                F1a = 2d0*cmenu(ica,lma)*oa +
     .          cmenu(ica,lma)*cmenu(ica,lma)*p(ica,lma)
C                F1b = 2d0*cmenu(icb,lmb)*ob +
C     .                    cmenu(icb,lmb)*cmenu(icb,lmb)*p(icb,lmb)
                F2a = oa + p(ica,lma)*cmenu(ica,lma)
                F2b = ob + p(icb,lmb)*cmenu(icb,lmb)
                G1a = cmenu(ica,lma)*cmenu(ica,lma)*
     .          (oa + p(ica,lma)*enu(ica,lma))
     .          + 2d0*cmenu(ica,lma)*oa*enu(ica,lma)
C                G1b = cmenu(icb,lmb)*cmenu(icb,lmb)*
C     .                (ob + p(icb,lmb)*enu(icb,lmb))
C     .              + 2d0*cmenu(icb,lmb)*ob*enu(icb,lmb)
                G2a = oa*enu(ica,lma)
     .          + (oa + p(ica,lma)*enu(ica,lma))*cmenu(ica,lma)
                G2b = ob*enu(icb,lmb)
     .          + (ob + p(icb,lmb)*enu(icb,lmb))*cmenu(icb,lmb)
                G3a = oa + p(ica,lma)*enu(ica,lma)
                G3b = ob + p(icb,lmb)*enu(icb,lmb)
                if (onsite) then
                  if (lma .eq. lmb) then
                    o2 = 1d0 + F1a + 2d0*dsd*F2a + dsd*dsd*p(ica,lma)
                    h2 = c(ica,lma) + dsd + G1a +
     .              2d0*dsd*G2a + dsd*dsd*G3a
                  else
                    o2 = 0d0
                    h2 = dsd
                  endif
                else
                  o2 = dsd*(F2a + F2b)
     .            + dsd*p(icb,lmb)*Don(lmb,lmb,ib)
     .            + dsd*p(ica,lma)*Don(lma,lma,ia)
                  h2 = dsd + dsd*(G2a + G2b)
     .            + dsd*G3b*Don(lmb,lmb,ib)
     .            + dsd*G3a*Don(lma,lma,ia)
                endif
                h1rs(lma,lmb,isite) = h1
                h2rs(lma,lmb,isite) = h2
                s2rs(lma,lmb,isite) = o2
                write (ifi,6) ib,ia,(pvec(k),k=1,3),(vec(k),k=1,3),
     .          vec(4),lma,lmb,hrs(lma,lmb,isite),h1,h2,
     .          srs(lma,lmb,isite),o2,outs(1:i+1)
              else
                if (lov) then
                  write (ifi,7) ib,ia,(pvec(k),k=1,3),(vec(k),k=1,3),
     .            vec(4),lma,lmb,hrs(lma,lmb,isite),
     .            srs(lma,lmb,isite),outs(1:i+1)
                else
                  write (ifi,8) ib,ia,(pvec(k),k=1,3),(vec(k),k=1,3),
     .            vec(4),lma,lmb,hrs(lma,lmb,isite),outs(1:i+1)
                endif
              endif
            endif
          enddo
        enddo
      enddo
      call fclose('RSH')
    1 format ('  R   R''',2x,'-- plat --',8x,'  vector',12x,
     .'dist',4x,'L',1x,'L''',5x,'H',6x,'H_1(gam)',2x,
     .'H_2(gam)',5x,'O',6x,'O_2(gam)')
    2 format ('  R   R''',2x,'-- plat --',8x,'  vector',12x,
     .'dist',4x,'L',1x,'L''',5x,'H',6x,'H_1(alp)',2x,
     .'H_2(alp)',5x,'O',6x,'O_2(alp)')
    3 format ('  R   R''',2x,'-- plat --',8x,'  vector',12x,
     .'dist',4x,'L',1x,'L''',5x,'H',10x,'O')
    4 format ('  R   R''',2x,'-- plat --',8x,'  vector',12x,
     .'dist',4x,'L',1x,'L''',5x,'H')
    6 format (i3,i4,2x,3i3,1x,3f8.4,1x,f10.6,1x,2i2,1x,5f10.6,2x,a)
    7 format (i3,i4,2x,3i3,1x,3f8.4,1x,f10.6,1x,2i2,1x,2f10.6,2x,a)
    8 format (i3,i4,2x,3i3,1x,3f8.4,1x,f10.6,1x,2i2,1x,f10.6,2x,a)
      end

      subroutine clean0(len,a)
C- Clean up "zero's"
C ----------------------------------------------------------------------
Ci Inputs:
Ci
Co Outputs:
Co
Cr Remarks
Cr   I don't know why you need to do this. I think it's the Intel
Cr   compiler, but numbers that are zero can appear as -0 and this
Cr   confuses awrite. d1mach(1) or d1mach(3) doesn't work
Cr   either. It's a nightmare. To coin a phrase.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer len
      double precision a(len)
C Local Variables
      integer i
      double precision dabs
      do  i = 1, len
        if (dabs(a(i)) .lt. 1d-7) a(i) = 0d0
      enddo
      end

      subroutine wripp(nl,nsp,nclass,pp,dclabl,ifi)
C- Write potential parameters to RSH file
C ----------------------------------------------------------------------
Ci Inputs:
Ci
Co Outputs:
Co
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nl,nsp,nclass
      double precision pp(6,nl,nsp,nclass),dclabl(nclass)
C Local Variables
      integer ifi,ic,l,k
      double precision sqd,d,oalp,pgam,gmina,dsign
      character*8 clabl

      do  ic = 1, nclass
        call r8tos8(dclabl(ic),clabl)
        write (ifi,1) clabl,nl
        do  l = 1, nl
          gmina = pp(6,l,1,ic) - pp(5,l,1,ic)
          oalp = gmina/
     .    (pp(3,l,1,ic)**2 + gmina*(pp(1,l,1,ic) - pp(2,l,1,ic)))
          sqd = pp(3,l,1,ic)
          d = sqd**2*dsign(1d0,sqd)
          pgam = pp(4,l,1,ic) - oalp**2
          write(ifi,2) l, (pp(k,l,1,ic), k=1,2),
     .    d,pp(4,l,1,ic),pgam,(pp(k,l,1,ic), k=5,6), oalp
        enddo
      enddo
      write (ifi,3)
    1 format(/'PPAR:  ',a,'  nl=',i1,/' l',5x,'e_nu',10x,'c',8x,
     .'+/-del',7x,'p^alp',7x,'p^gam',8x,'gam',9x,'alp',10x,'o')
    2 format(i2,3f12.8,f12.7,5f12.8)
    3 format (/)
      end

