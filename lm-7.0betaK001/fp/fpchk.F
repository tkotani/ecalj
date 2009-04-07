      subroutine fpchk(sspec,ssite)
C- Routines to check various quantities in FP code
C ----------------------------------------------------------------------
Ci Inputs
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: rmt rsma lmxa kmxt
Ci     Stored:    name
Ci     Passed to: chkxpn spacks uspecb
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    spec
Ci     Passed to: chkxpn
Cu Updates
Cu   23 Apr 02 Added option (mode=0) to find MT radii
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision sspec(1),ssite(1)

      call info(5,2,0,
     .   ' ----- Check accuracy of augmentation expansion -----',0,0)
      call chkxpn(ssite,sspec)
      end

      subroutine chkxpn(ssite,sspec)
C- Check accuracy of Pkl expansions of augmented functions
C     implicit none
C ... Passed parameters
      double precision ssite(1),sspec(1)
C ... Local parameters
      integer k0,nlm0,nrx,n0,nkap0,stdo,lgunit
      parameter (k0=20, nlm0=25, nrx=201)
      parameter (nkap0=3,n0=10)
      integer nglob,nbas,nspec,ib1,is,oips,
     .  iclbsj,nr,j1,nch,lmxa,kmax,lh(nkap0),nkapi,ik,l1,l2,l,
     .  lmax,nlm,k,i,ilm,ipr
      double precision c(0:k0,nlm0),
     .  rofi(nrx),rwgt(nrx),pkl(0:k0,0:n0),
     .  xi(0:n0),sumh(0:n0),difh(0:n0),errh(0:n0),errj(0:n0),
     .  cj(0:k0,0:n0),phi(0:n0),psi(0:n0),difj(0:n0),sumj(0:n0),
     .  rsmh(n0,nkap0),eh(n0,nkap0),
     .  a,rmt,rsma,e,rsm,aa,r,rl,fith,fitj,fcth,fctj,wgt
      character*8 spid
      character lch(0:4)
C ... Heap
      integer w(1)
      common /w/ w
      data lch /'s','p','d','f','g'/

C ... setup
      stdo = lgunit(1)
      nbas  = nglob('nbas')
      nspec = nglob('nspec')
      call getpr(ipr)
      write(stdo,400)
  400 format(/' RMS errors in expanded head and tail functions.'
     .   /' Only column ''tail'' matters for program accuracy.')

      call defi(oips,nbas)
      call spackv(10,'site spec',ssite,1,nbas,w(oips))

C --- Loop over species ---
      do  40  is = 1, nspec
        ib1 = iclbsj(is,w(oips),-nbas,1)
        if (ib1 .lt. 0) then
          if (ipr .ge. 20) write(stdo,
     .   '('' chkxpn (warning) no sites corresponding to species'',i3)')
     .      is
          goto 40
        endif

        call spacks(0,'spec name',sspec,spid,is,is)
        call upack('spec rmt rsma lmxa kmxt',sspec,is,rmt,rsma,lmxa,
     .    kmax)
C        call u_spec_p(sspec,is,!z,lmxl,lmxa,pnu,idmod,!a0,!nr0,rmt)
C        call u_spec_a(sspec,is,lmxa,lmxh,kmax,rsma)
        a = 0.02d0
        nr = 101

        call word(spid,1,j1,nch)
        write(stdo,600) spid(1:nch),rmt,rsma,lmxa,kmax
  600   format(/' Species ',a,':   rmt=',f7.4,'   rsma=',f7.4
     .     /' augment to  lmax=',i2,'   kmax=',i3)

        if (nr .gt. nrx) call rxi('chkxpn: need nrx ge',nr)
        call radmsh(rmt,a,nr,rofi)
        call radwgt(rmt,a,nr,rwgt)
        call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)

C  ...  Loop over orbitals for this species
        write(stdo,948)
        do  ik = 1, nkapi
        l2 = -1
        do  l1  = 0, lh(ik)
          e = eh(l1+1,ik)
          rsm = rsmh(l1+1,ik)
          if (rsm .ne. 0 .and. l1 .gt. l2) then
            l2 = l1-1
   12       l2 = l2+1
            if (l2 .lt. lh(ik)) then
            if (rsmh(l2+2,ik).eq.rsm .and. eh(l2+2,ik).eq.e) goto 12
            endif

          lmax = l2
          nlm = (lmax+1)**2

          if (nlm .gt. nlm0) call rxi('chkxpn: need nlm0 ge',nlm)
          if (kmax .gt. k0)  call rxi('chkxpn: need k0 ge',kmax)

C     ... Compare expanded and exact head,tail functions
          do  l = 0, lmax
            difh(l) = 0d0
            sumh(l) = 0d0
            difj(l) = 0d0
            sumj(l) = 0d0
          enddo

          call hxpos(rsmh,rsma,eh,kmax,nlm,k0,c)
C         call hxps(rsmh,rsma,e,kmax,nlm,k0,c)
          call jxpos(rsma,e,kmax,lmax,k0,cj)
          aa = 1d0/rsm

          do  i = 1, nr
            r = rofi(i)
            call radpkl(r,rsma,kmax,lmax,k0,pkl)
            call hansmr(r,e,aa,xi,lmax)
            call bessl(e*r*r,lmax,phi,psi)
            do  l = 0, lmax
              if (l .eq. 0) then
                rl = 1
              else
                rl = r**l
              endif
              ilm = l*l+1
              fith = 0d0
              fitj = 0d0
              do  k = 0, kmax
                fith = fith + c(k,ilm)*pkl(k,l)*rl
                fitj = fitj + cj(k,l)*pkl(k,l)*rl
              enddo
              fcth = xi(l)*rl
              fctj = phi(l)*rl
              wgt = rwgt(i)*r*r*rl
              difh(l) = difh(l) + wgt* (fith-fcth)**2
              sumh(l) = sumh(l) + wgt* fcth**2
              difj(l) = difj(l) + wgt* (fitj-fctj)**2
              sumj(l) = sumj(l) + wgt* fctj**2
            enddo
          enddo

          do  l = l1, l2
            errh(l) = dsqrt(difh(l)/sumh(l))
            errj(l) = dsqrt(difj(l)/sumj(l))
            write(stdo,947) l,rsm,e,lch(l),errh(l),errj(l)
  947       format(i4,f10.3,f9.3,2x,3x,a1,1p,6d12.2)
  948       format(/' block    rsm        e      l',6x,
     .         'head',8x,'tail')


          enddo
          endif
        enddo
        enddo
  40  continue
      end

