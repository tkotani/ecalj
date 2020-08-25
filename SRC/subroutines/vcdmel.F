      subroutine vcdmel(nl,ssite,sspec,slat,nlmax,ndham,ndimh,
     .nq,nsp,nspc,ef,evl,aus,nsite,isite,iclsl,iclsn,sbz,dosw)
      use m_work,only: iv_a_oidtet, rv_a_ocg , iv_a_ojcg , iv_a_oidxcg
      use m_struc_def
      use m_ext,only: sname !extention for file
C- Valence-core dipole matrix elements
C ----------------------------------------------------------------------
Ci Inputs:
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: nl
Ci   ssite :struct containing site-specific information
Ci     Elts read: spec pnu ov0 pz
Ci   sspec :struct containing species-specific information
Ci     Elts read: a nr rmt z lmxa
Ci   slat  :struct containing information about the lattice
Ci     Elts read: ocg ojcg oidxcg ocy
Ci   nlmax :first dimension of aus; largest augmentation (l+1)^2
Ci   ndham :second dimension of aus, at least as large as ndimh
Ci   ndimh :number of eigenvalues
Ci   nq    :number of k-points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ef    :Fermi energy
Ci   evl   :energy bands at the nq k-points
Ci   aus   :values and slopes of eigenstates at MT sphere surfaces
Ci          (makusq)
Ci   nsite,isite,iclsl,iclsn see suclst
Co Outputs:
Co   weights for each channel output in iomoms style
Cr Remarks
Cu Updates
Cu   08 Jul 08 Dimension aus separately from ndimh
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   19 Sep 03 (ATP) Bug fixes
Cu   28 Mar 01 (MvS) rearrangement of indices to accommodate new makusq.
Cu   20 Mar 01 (ATP) extended to handle multiple core levels
Cu   20 Feb 01 Written by ATP
C ----------------------------------------------------------------------
      implicit none
      type(s_bz)  ::  sbz
C Passed Parameters
      integer nlmax,ndham,ndimh,nq,nsp,nspc,nsite
      integer isite(nsite),iclsl(nsite),iclsn(nsite)
      real(8):: ef , evl(ndham,nsp,nq)
      type(s_ctrl)::sctrl
      type(s_site)::ssite(nsite)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      double complex aus(nlmax,ndham,3,nsp,nsite,nq)
C Local Variables
      integer n0,lmxax
      parameter (n0=10,lmxax=10)
      integer ifi,fopna,isp,ib,is,lcls,ncls,nl,i,j, 
     .nr,lmxa,iq,nlma,igets,igetss,i1mach,nfstg,nchan
      integer lh(10)
      real(8) ,allocatable :: rofi_rv(:)
      real(8) ,allocatable :: ul_rv(:)
      real(8) ,allocatable :: sl_rv(:)
      real(8) ,allocatable :: gz_rv(:)
      real(8) ,allocatable :: ruu_rv(:)
      real(8) ,allocatable :: rus_rv(:)
      real(8) ,allocatable :: rss_rv(:)
      real(8) ,allocatable :: g_rv(:)
      real(8) ,allocatable :: s_rv(:,:,:)
      double precision pnu(n0,2),pnz(n0,2),a,rmt,z,xx,rsml(n0),ehl(n0)
      double precision ume(0:lmxax,nsp,nsite),sme(0:lmxax,nsp,nsite)
      character clabl*8

      integer:: nbandx,nspx,npts,ifid,fopnn,ikp,ichib
      logical:: lidos
      real(8),allocatable:: wk_rv(:),dos_rv(:,:,:)
      real(8):: emin,dosw(2),emax,del
      call tcn ('vcdmel')
      rsml=0d0 !call dpzero(rsml,n0)
      ehl=0d0  !call dpzero(ehl,n0)
      do  i = 1, nsite
        ib = isite(i)
        ncls = iclsn(i)
        lcls = iclsl(i)
        is = ssite(ib)%spec
        call dcopy(n0*2,ssite(ib)%pnu,1,pnu,1)
        call dcopy(n0*2,ssite(ib)%pz,1,pnz,1)
        clabl=sspec(is)%name
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        z=sspec(is)%z
        lmxa=sspec(is)%lmxa
        if (lmxa .gt. lmxax) call rxi('vcdmel needs lmxax ',lmxa)
        if (lmxa .eq. -1) cycle !  if (lmxa .eq. -1) goto 10
        allocate(rofi_rv(nr))
        call radmsh ( rmt , a , nr , rofi_rv )
C   --- Augmented wave functions u,s
        allocate(ul_rv(nr*(lmxa+1)*nsp))
        allocate(sl_rv(nr*(lmxa+1)*nsp))
        allocate(gz_rv(nr*(lmxa+1)*nsp))
        allocate(ruu_rv(nr*(lmxa+1)*2*nsp))
        allocate(rus_rv(nr*(lmxa+1)*2*nsp))
        allocate(rss_rv(nr*(lmxa+1)*2*nsp))
        call uspecb(0,4,sspec,is,is,lh,rsml,ehl,j)
        call makusp ( n0 , z , nsp , rmt , lmxa , ssite(ib)%rv_a_ov0 , a , nr , 
     .   xx , xx , pnu , pnz , rsml , ehl , ul_rv , sl_rv , gz_rv 
     .   , ruu_rv , rus_rv , rss_rv )
C   --- Matrix elements of u,s with core
c        call skpblb(clabl,8,j)
c        call awrit3('%N CLS atom %i ('//clabl(1:j+1)//') n=%i l=%i',
c     .       ' ',120,i1mach(2),ib,ncls,lcls)
        write(6,"('CLS atom: ib name n l=',i5,' ',a,2i5)") ib, trim(clabl),ncls,lcls
        allocate(g_rv(nr*2))
        call pvcdm1 ( ncls , lcls , g_rv , z , lmxa , ssite(ib)%rv_a_ov0 , a 
     .   , nr , rofi_rv , ul_rv , sl_rv , nsp , lmxax , ume ( 
     .   0 , 1 , i ) , sme ( 0 , 1 , i ) )
        deallocate(g_rv,rss_rv,rus_rv,ruu_rv,gz_rv,sl_rv,ul_rv,rofi_rv)
      enddo
C --- Open CLS weights file and write first line
      allocate( s_rv(3*nsite*2*ndimh,nsp,nq) )
      nfstg = 11
      nchan = 3*nsite !this means only half of 3*nsite*2. S
c      ifi = fopna('cls',-1,4)
c      j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,0,1,0,0,0,0,0,
c     .0d0,0d0,0d0,0d0,0d0,0d0)
C --- For each qp, make <nk|x,y,z|core> at each site and save to disk in
C     iomoms-compatible format ---
      s_rv=0d0
      do   iq = 1, nq
        do  isp = 1, nsp
c          call dpzero ( s_rv , 3 * ndimh * nsite * 2 )
          do  i = 1, nsite
            lcls = iclsl(i)
            ib   = isite(i)
            is   = ssite(ib)%spec
            lmxa = sspec(is)%lmxa
            nlma = (lmxa+1)**2
            if (lmxa .gt. -1) then
              call pvcdm2 ( i , nsite , ndham , ndimh , nlma , nlmax , aus 
     .              ( 1 , 1 , 1 , isp , i , iq ) , ume ( 0 , isp , i ) , sme ( 0 , isp , i ) , lcls ,
     &         rv_a_ocg , iv_a_ojcg , iv_a_oidxcg, s_rv(1,isp,iq))
            endif
          enddo
C --- Scale weights arbitrarily by 100 for plotting etc ..
          call dscal ( 3 * ndimh * nsite * 2, 1d2 , s_rv(1,isp,iq) , 1 )
c          j = iomoms ( - ifi , nl , nsp , nspc , nq , ndimh , nfstg , 1 
c     .     , 1 , 1 , ndimh , ndimh , nchan , nchan , ndimh , evl ( 1 , isp 
c     .     , iq ) , 0d0 , s_rv(1,isp,iq) , 0d0 , 0d0 , 0d0 )   !only save s_rv(1:nchan,ndimh)
        enddo
      enddo
c      iq = nq*nsp + 1
c      j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,iq,1,ndimh,ndimh,
c     .nchan,nchan,ndimh,0d0,0d0,0d0,0d0,ef,0d0)
c      call fclose(ifi)
c$$$ccccccccccccccccccccccccc
      nbandx = ndimh*nspc       !q-dependent only for no PW case.
      nspx = nsp / nspc
      lidos=.false.
      npts=sbz%ndos
      allocate(wk_rv(npts))
c     print *,'ewindow'//trim(sname)//'xxxx'
c     open(newunit=ifid,file='ewindow.'//trim(sname) )
c     read(ifid,*) emin
c     close(ifid)
      emin = dosw(1)
      emax = sbz%dosmax + ef
!     ! cf. aaax in asados.F at aug2020
c     print *,'aaax', ndham , nsp , nspx , nbandx , nchan ,
c     &         sbz%nkabc(1),sbz%nkabc(2),sbz%nkabc(3), sbz%ntet, npts, emin, emax, lidos
c     do ib=1,nbandx
c     do isp=1,nsp
c     do ikp=1,nq
c     write(6,"('aaa ',3i5,f12.5)")ib,isp,ikp,evl(ib,isp,ikp)
c     enddo
c     enddo
c     enddo
c     do iq=1,nq
c     do isp=1,nsp
c     do ichib=1,3*ndimh*nsite
c     write(6,"('aaa ',3i5,f12.5)")ichib,isp,iq,s_rv(ichib,isp,iq)
c     enddo
c     enddo
c     enddo
      allocate(dos_rv(npts,nsp,nchan))
      call dostet ( ndham , nsp , nspx , nbandx , nchan ,
     &     sbz%nkabc(1),sbz%nkabc(2),sbz%nkabc(3), sbz%ntet , iv_a_oidtet, evl
     .     , s_rv(1:nchan*nbandx,1:nsp,1:nq), npts, emin, emax, lidos , wk_rv , dos_rv )
c     print *,'end of dostet xxxxx', npts,nsp,nchan
      del = 0d0
      call dosio ( dos_rv , npts * ( 1 + nsp - i ) , nchan , npts 
     .     , nchan , emin , emax , i , ef , del , 1 , - fopnn ('DOS-VCDMEL') )
c$$$ccccccccccccccccccccccccc
      if (allocated(s_rv)) deallocate(s_rv)
      call tcx ('vcdmel')
      end subroutine vcdmel

      
      subroutine pvcdm1(ncls,lcls,gcore,z,lmxa,v,a,nr,rofi,ul,sl,nsp,
     .lmxax,ume,sme)
C- Radial matrix elements < (u,s) | r | core >
C     implicit none
      integer ncls,lcls,lmxa,nr,nsp,lmxax
      double precision a,z,gcore(nr,2),rofi(nr),v(nr,nsp),
     .ul(nr,0:lmxa,nsp),sl(nr,0:lmxa,nsp),ume(0:lmxax,nsp),
     .sme(0:lmxax,nsp)
C Local
      integer nodes,l,nre,isp,ll,ir,i1mach
      double precision e1,e2,slo,val,rmax,b,ecore,tol,yyy,dlml,slo1,
     .r,wgt,uc,sc,ecor0,sum

      do  isp = 1, nsp
        if (nsp .eq. 2) then
          call info2(30,0,0,' Spin %i ..',isp,0)
        endif
C   --- gcore <- core level wave function * r ---
        tol = 1.d-8
        e1 = -2.5d0*z*z - 5
        e2 = 20.d0
        val = 1.d-30
        slo = -val
        l = lcls
        rmax = rofi(nr)
        b = rmax/(dexp(a*nr-a)-1.d0)
        nodes = ncls - (l+1)
        ecore = (e1+e2)/2
        call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,v(1,isp),gcore,sum,
     .  a,b,rofi,nr,nre)
        ecor0 = ecore
C   ... Correct core energy by using hankel bc's
        yyy = ecore - v(nr,isp) + 2*z/rmax
        if(nre .eq. nr .and. yyy .lt. 0.d0) then
          dlml = -1.d0-dsqrt(-yyy)*rmax
          do  ll = 1, l
            dlml = -yyy*rmax*rmax/dlml - (2*ll+1)
          enddo
          slo1 = val*(dlml+l+1)/rmax
          call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo1,v(1,isp),gcore,
     .    sum,a,b,rofi,nr,nre)
        endif
c        call awrit2(' vcdmel: ecor0=%,8d ecore=%,8d%N'
c     .  //' (not including electrostatic potential shift)',
c     .  ' ',120,i1mach(2),ecor0,ecore)
        write(6,"('vcdmel: ecor0=',f15.8,' ecore=',f15.8)")ecor0,ecore
        write(6,"('(not including electrostatic potential shift)')")
C --- Matrix elements < (u,s) | r | core > ---
        print 332
  332   format( '   l',3x,'<u|core>',5x,'<s|core>',4x,
     .  '<u|r|core>',2x,'<s|r|core>')
        do  l = 0, lmxa
          ume(l,isp) = 0
          sme(l,isp) = 0
          uc = 0
          sc = 0
          do  ir = 2, nre-1
            r = rofi(ir)
            wgt = (mod(ir+1,2)+1) * (r+b)
            uc     = uc + wgt * ul(ir,l,isp) * gcore(ir,1)
            sc     = sc + wgt * sl(ir,l,isp) * gcore(ir,1)
            ume(l,isp) = ume(l,isp)
     .      + wgt * ul(ir,l,isp) * r * gcore(ir,1)
            sme(l,isp) = sme(l,isp)
     .      + wgt * sl(ir,l,isp) * r * gcore(ir,1)
          enddo
          ir = nre
          r = rofi(ir)
          wgt = .5d0 * (r+b)
          uc     = uc + wgt * ul(ir,l,isp) * gcore(ir,1)
          sc     = sc + wgt * sl(ir,l,isp) * gcore(ir,1)
          ume(l,isp) = ume(l,isp) + wgt * ul(ir,l,isp) * r * gcore(ir,1)
          sme(l,isp) = sme(l,isp) + wgt * sl(ir,l,isp) * r * gcore(ir,1)
          uc = uc*2d0*a/3d0
          sc = sc*2d0*a/3d0
          ume(l,isp) = ume(l,isp)*2d0*a/3d0
          sme(l,isp) = sme(l,isp)*2d0*a/3d0
          print 335, l,uc,sc,ume(l,isp),sme(l,isp)
  335     format(i4,4f12.6)
        enddo
      enddo
      end subroutine pvcdm1

      subroutine pvcdm2(isite,nsite,ndham,ndimh,nlma,nlmax,aus,ume,sme,
     .lcls,cg,jcg,indxcg,s)

C- Kernel called by vcmdel
C ----------------------------------------------------------------------
Ci Inputs
Ci   isite :
Ci   nsite :
Ci   ndimh :
Ci   nlma  :
Ci   nlmax :
Ci   aus   :
Ci   ume   :
Ci   sme   :
Ci   lcls  :
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Co Outputs
Co   s     :Matrix elements
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer isite,lcls,ndham,ndimh,nlma,nlmax,nsite,indxcg(*),jcg(*)
      double precision cg(*),ume(0:*),sme(0:*),s(3,nsite,ndimh,2)
      double complex aus(nlmax,ndham,2)
C ... Local parameters
      integer kk(4),mlm,lm,ll,klm,ii,indx,icg1,icg2,icg,llm,ib,k
      double complex cxx
C     Transposes (y,z,x) to (x,y,z)
      data kk /0,2,3,1/
C ... Loop over lm of (u,s)
      do  11  mlm = 1, nlma
        lm = ll(mlm)
C       Selection rule would be handled by CG anyway:
        if (lm .eq. lcls-1 .or. lm .eq. lcls+1) then
C     ... Loop over m index of r
          do  14  klm = 2, 4
            ii = max0(mlm,klm)
            indx = (ii*(ii-1))/2 + min0(mlm,klm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            do  15  icg = icg1, icg2
C             lm of core
              llm  = jcg(icg)
              if (ll(llm) .eq. lcls) then
                do  10  ib = 1, ndimh
                  cxx =  cg(icg)*
     .            (dconjg(aus(mlm,ib,1))*ume(lm) +
     .            dconjg(aus(mlm,ib,2))*sme(lm))
C                  write (*,200) kk(klm),mlm,klm,llm,aus(mlm,ib,1),
C     .              aus(mlm,ib,2),cg(icg),cxx
                  s(kk(klm),isite,ib,1) = s(kk(klm),isite,ib,1)
     .            + dble(cxx)
                  s(kk(klm),isite,ib,2) = s(kk(klm),isite,ib,2)
     .            + dimag(cxx)
   10           continue
              endif
   15       continue
   14     continue
        endif
   11 continue
      do  20  k = 1, 3
      do  20  ib = 1, ndimh
        s(k,isite,ib,1) = s(k,isite,ib,1)*s(k,isite,ib,1)
     .    + s(k,isite,ib,2)*s(k,isite,ib,2)
   20 continue
      end subroutine pvcdm2


