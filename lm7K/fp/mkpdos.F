Cgetarg...info...           structure ['mkpdos', 'site', 'ssite']
Cgetarg...info...           structure ['mkpdos', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['mkpdos', 'm_struc_def']
      subroutine mkpdos(mode,ssite,sspec,isp,nsp,nlmax,ndham,nev,
     .  nchan,npd,lsite,nsite,ppnl,aus,doswt)
       
       use m_struc_def  !Cgetarg

C- Decomposition of norm from projection of w.f into MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci         :0 or 3 DOS resolved by site R
Ci         :1 or 4 DOS resolved by R and l
Ci         :2 or 5 DOS resolved by R,l,m
Ci         :10s digit (not implemented)
Ci         :1 doswt = separate phi,phidot contributions to norm, i.e.
Ci         :          <phi | evec> ;  <phidot | evec>
Ci         :          In this case, ndp should be at least 2.
Ci         :100s digit
Ci         :0  aus = (phi,phidot) of w.f. at MT boundary
Ci         :1  aus = (val,slo) of w.f. at MT boundary (not implemented)
Ci   ssite :struct containing site-specific information
Ci          Elts read: spec
Ci   sspec :struct containing species-specific information
Ci          Elts read: lmxa
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmax :leading dimension of aus
Ci   ndham :dimensions aus
Ci   nev   :number of eigenvectors to accumulate DOS
Ci   nchan :total number of channels
Ci   npd   :second dimension of doswt.
Ci   lsite :list of sites for which to accumulate partial dos
Ci         :if lsite(1)=0, a site list 1..nsite is used instead
Ci   nsite :number of sites in list
Ci   ppnl  :NMTO potential parameters; see eg potpus.f
Ci   aus   :values of (phi,phidot) MT sphere boundary; see makusq
Co Outputs
Co   doswt :dos weights for this qp
Cu Updates
Cu   13 Feb 02 Extended to local orbitals.
Cu   28 Mar 01 written by MvS.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,isp,nsp,nlmax,ndham,nchan,nsite,lsite(nsite),nev,npd
      integer n0,nppn
      parameter (n0=10,nppn=12)
Cgetarg       double precision ssite(1),sspec(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

      double precision ppnl(nppn,n0,nsp,*),doswt(nchan,npd,nev)
      double complex aus(nlmax,ndham,3,nsp,*)
C ... Local parameters
      integer lmxa,lmxax,lmdim,ichan,i,ib,is,igetss,iv,ilm,l,m,nlma,ll
      double precision summ(121),suml(11),s11,s22,s33,s12,s13,s23,
     .  suma,sum,rmt
      double complex au,as,az

      call isanrg(mode,0,5,' mkpdos','mode',.true.)
      lmxax = ll(nlmax)

      ichan = 0
      do  i = 1, nsite
        if (lsite(1) .le. 0) then
          ib = i
        else
          ib = lsite(i)
        endif
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec lmxa rmt',sspec,is,lmxa,rmt,0,0)
         
         lmxa=sspec(is)%lmxa
         rmt=sspec(is)%rmt

        lmxa = min(lmxa,lmxax)
        if (lmxa .eq. -1) goto 10

        nlma = (lmxa+1)**2
C       Decompose by site R
        if (mode .eq. 0 .or. mode .eq. 3) then
          lmdim = 1
C       Decompose by Rl
        elseif (mode .eq. 1 .or. mode .eq. 4) then
          lmdim = 1+lmxa
C       Decompose by Rlm
        elseif (mode .eq. 2 .or. mode .eq. 5) then
          lmdim = nlma
        endif

C       call pp2hvs(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)
C       call phvsfp(1,nsp,lmxa,ppnl(1,1,1,ib),rmt,sab,sab,sab)

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

C           If (au,as) are coefficients to (u,s), use this
C           s11 = dconjg(au)*au*sab(1,l+1,isp)
C           s12 = 2*dconjg(au)*as*sab(2,l+1,isp)
C           s22 = dconjg(as)*as*sab(4,l+1,isp)

C           If (au,as) are coefficients to (phi,phidot), use this
            s11 = dconjg(au)*au*ppnl(2,l+1,isp,ib)
            s22 = dconjg(as)*as*ppnl(7,l+1,isp,ib)
            s33 = dconjg(az)*az*ppnl(8,l+1,isp,ib)
            s12 = 0
            s13 = 2*dconjg(au)*az*ppnl(9,l+1,isp,ib)
            s23 = 2*dconjg(as)*az*ppnl(10,l+1,isp,ib)

            sum = s11+s22+s33 + s12+s13+s23

C            if (iv .le. 5 .and. l .eq. 2 .and. s33 .ne. 0) then
C              print 333, iv,m,s11+s22+s33,s12+s13+s23
C  333         format('pdos',2i4,2f12.6)
C            endif

            suma = suma + sum
            suml(l+1) = suml(l+1) + sum
            summ(ilm) = sum

          enddo
          enddo

          if (lmdim .eq. 1) then
            call dcopy(1,suma,1,doswt(1+ichan,1,iv),1)
          else if (lmdim .eq. 1+lmxa) then
            call dcopy(1+lmxa,suml,1,doswt(1+ichan,1,iv),1)
          else
            call dcopy(nlma,summ,1,doswt(1+ichan,1,iv),1)
          endif

        enddo
        ichan = ichan + lmdim
   10   continue
      enddo

      call isanrg(ichan,1,nchan,' mkpdos','nchan',.true.)

Cgetarg       end
       end subroutine mkpdos 


