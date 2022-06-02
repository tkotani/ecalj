subroutine mkpdos(mode,ssite,sspec,isp,nsp,nlmax,ndham,nev, &
     nchan,npd,lsite,nsite,ppnl,aus,doswt)
  use m_struc_def  !Cgetarg
  !- Decomposition of norm from projection of w.f into MT sphere
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i         :0 or 3 DOS resolved by site R
  !i         :1 or 4 DOS resolved by R and l
  !i         :2 or 5 DOS resolved by R,l,m
  !i         :10s digit (not implemented)
  !i         :1 doswt = separate phi,phidot contributions to norm, i.e.
  !i         :          <phi | evec> ;  <phidot | evec>
  !i         :          In this case, ndp should be at least 2.
  !i         :100s digit
  !i         :0  aus = (phi,phidot) of w.f. at MT boundary
  !i         :1  aus = (val,slo) of w.f. at MT boundary (not implemented)
  !i   ssite :struct containing site-specific information
  !i          Elts read: spec
  !i   sspec :struct containing species-specific information
  !i          Elts read: lmxa
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nlmax :leading dimension of aus
  !i   ndham :dimensions aus
  !i   nev   :number of eigenvectors to accumulate DOS
  !i   nchan :total number of channels
  !i   npd   :second dimension of doswt.
  !i   lsite :list of sites for which to accumulate partial dos
  !i         :if lsite(1)=0, a site list 1..nsite is used instead
  !i   nsite :number of sites in list
  !i   ppnl  :NMTO potential parameters; see eg potpus.f
  !i   aus   :values of (phi,phidot) MT sphere boundary; see makusq
  !o Outputs
  !o   doswt :dos weights for this qp
  !u Updates
  !u   13 Feb 02 Extended to local orbitals.
  !u   28 Mar 01 written by MvS.
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,isp,nsp,nlmax,ndham,nchan,nsite,lsite(nsite),nev,npd
  integer :: n0,nppn
  parameter (n0=10,nppn=12)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)

  double precision :: ppnl(nppn,n0,nsp,*),doswt(nchan,npd,nev)
  double complex aus(nlmax,ndham,3,nsp,*)
  ! ... Local parameters
  integer :: lmxa,lmxax,lmdim,ichan,i,ib,is,igetss,iv,ilm,l,m,nlma,ll
  double precision :: summ(121),suml(11),s11,s22,s33,s12,s13,s23, &
       suma,sum,rmt
  double complex au,as,az
  logical:: l_dummy_isanrg,isanrg
  ! ino isanrg is logical function,       call isanrg(mode,0,5,' mkpdos','mode',.true.)
  l_dummy_isanrg=isanrg(mode,0,5,' mkpdos','mode',.true.)
  lmxax = ll(nlmax)
  ichan = 0
  do  i = 1, nsite
     if (lsite(1) <= 0) then
        ib = i
     else
        ib = lsite(i)
     endif
     is  = ssite(ib)%spec
     lmxa= sspec(is)%lmxa
     rmt = sspec(is)%rmt
     lmxa = min(lmxa,lmxax)
     if (lmxa == -1) goto 10

     nlma = (lmxa+1)**2
     !       Decompose by site R
     if (mode == 0 .OR. mode == 3) then
        lmdim = 1
        !       Decompose by Rl
     elseif (mode == 1 .OR. mode == 4) then
        lmdim = 1+lmxa
        !       Decompose by Rlm
     elseif (mode == 2 .OR. mode == 5) then
        lmdim = nlma
     endif
     !       call pp2hvs(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)
     !       call phvsfp(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)
     do  iv = 1, nev
        call dpzero(suml,1+lmxa)
        suma = 0
        ilm = 0

        do  l = 0, lmxa
           do  m = -l, l
              ilm = ilm+1
              au = aus(ilm,iv,1,isp,ib)
              as = aus(ilm,iv,2,isp,ib)
              az = aus(ilm,iv,3,isp,ib)

              !           If (au,as) are coefficients to (u,s), use this
              !           s11 = dconjg(au)*au*sab(1,l+1,isp)
              !           s12 = 2*dconjg(au)*as*sab(2,l+1,isp)
              !           s22 = dconjg(as)*as*sab(4,l+1,isp)

              !           If (au,as) are coefficients to (phi,phidot), use this
              s11 = dconjg(au)*au*ppnl(2,l+1,isp,ib)
              s22 = dconjg(as)*as*ppnl(7,l+1,isp,ib)
              s33 = dconjg(az)*az*ppnl(8,l+1,isp,ib)
              s12 = 0
              s13 = 2*dconjg(au)*az*ppnl(9,l+1,isp,ib)
              s23 = 2*dconjg(as)*az*ppnl(10,l+1,isp,ib)

              sum = s11+s22+s33 + s12+s13+s23

              !            if (iv .le. 5 .and. l .eq. 2 .and. s33 .ne. 0) then
              !              print 333, iv,m,s11+s22+s33,s12+s13+s23
              !  333         format('pdos',2i4,2f12.6)
              !            endif

              suma = suma + sum
              suml(l+1) = suml(l+1) + sum
              summ(ilm) = sum
           enddo
        enddo
        if (lmdim == 1) then
           call dcopy(1,suma,1,doswt(1+ichan,1,iv),1)
        else if (lmdim == 1+lmxa) then
           call dcopy(1+lmxa,suml,1,doswt(1+ichan,1,iv),1)
        else
           call dcopy(nlma,summ,1,doswt(1+ichan,1,iv),1)
        endif
     enddo
     ichan = ichan + lmdim
10   continue
  enddo

  ! ino isanrg is logical function,       call isanrg(ichan,1,nchan,' mkpdos','nchan',.true.)
  l_dummy_isanrg=isanrg(ichan,1,nchan,' mkpdos','nchan',.true.)

end subroutine mkpdos


