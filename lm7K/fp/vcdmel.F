Cdelw1 var ['orofi', '1', 'orofi', 'nr', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['oul', '1', 'oul', 'nr*(lmxa+1)*nsp', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['osl', '1', 'osl', 'nr*(lmxa+1)*nsp', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['ogz', '1', 'ogz', 'nr*(lmxa+1)*nsp', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['oruu', '1', 'oruu', 'nr*(lmxa+1)*2*nsp', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['orus', '1', 'orus', 'nr*(lmxa+1)*2*nsp', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['orss', '1', 'orss', 'nr*(lmxa+1)*2*nsp', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['ov0', '-1', 'undef', '0', 'undef', 'x', 'vcdmel']
Cdelw1 var ['og', '1', 'og', 'nr*2', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['os', '0', 'os', '3*ndimh*nsite*2', 'real(8)', 'rv_w_', 'vcdmel']
Cdelw1 var ['ocg', '-1', 'undef', '0', 'undef', 'x', 'vcdmel']
Cdelw1 var ['ojcg', '-1', 'undef', '0', 'undef', 'x', 'vcdmel']
Cdelw1 var ['oidxcg', '-1', 'undef', '0', 'undef', 'x', 'vcdmel']
Cdelw1 undel [['undef', 'wvar', 'wref'], 'ov0', 'vcdmel']
Cdelw1 undel [['undef', 'wvar', 'wref'], 'ocg', 'vcdmel']
Cdelw1 undel [['undef', 'wvar', 'wref'], 'ojcg', 'vcdmel']
Cdelw1 undel [['undef', 'wvar', 'wref'], 'oidxcg', 'vcdmel']
      subroutine vcdmel(sctrl,ssite,sspec,slat,nlmax,ndham,ndimh,
     .nq,nsp,nspc,ef,evl,aus,nsite,isite,iclsl,iclsn)
      use m_struc_def  
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
C Passed Parameters
      integer nlmax,ndham,ndimh,nq,nsp,nspc,nsite
      integer isite(nsite),iclsl(nsite),iclsn(nsite)
      real(8):: ef , evl(ndham,nsp,nq)
      type(s_ctrl)::sctrl
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex aus(nlmax,ndham,3,nsp,nsite,nq)
C Local Variables
      integer n0,lmxax
      parameter (n0=10,lmxax=10)
      integer ifi,fopna,isp,ib,is,lcls,ncls,nl,i,j,iomoms,
     .nr,lmxa,iq,nlma,igets,igetss,i1mach,nfstg,nchan
      integer lh(10)
Cdelw1 do not change  [['undef', 'wvar', 'wref'], 'ocg', 'vcdmel']
Cdelw1 do not change  [['undef', 'wvar', 'wref'], 'ojcg', 'vcdmel']
Cdelw1 do not change  [['undef', 'wvar', 'wref'], 'oidxcg', 'vcdmel']
Cdelw1 do not change  [['undef', 'wvar', 'wref'], 'ov0', 'vcdmel']
Cdelw1       integer ocg,ojcg,oidxcg,ocy,ov0,orofi,oul,osl,ogz,oruu,orus,orss,
Cdelw1      .og,os
       integer:: ocg , ojcg , oidxcg , ocy , ov0 
       real(8) ,allocatable :: rv_w_orofi(:)
       real(8) ,allocatable :: rv_w_oul(:)
       real(8) ,allocatable :: rv_w_osl(:)
       real(8) ,allocatable :: rv_w_ogz(:)
       real(8) ,allocatable :: rv_w_oruu(:)
       real(8) ,allocatable :: rv_w_orus(:)
       real(8) ,allocatable :: rv_w_orss(:)
       real(8) ,allocatable :: rv_w_og(:)
       real(8) ,allocatable :: rv_w_os(:)

      double precision pnu(n0,2),pnz(n0,2),a,rmt,z,xx,rsml(n0),ehl(n0)
      double precision ume(0:lmxax,nsp,nsite),sme(0:lmxax,nsp,nsite)
      character clabl*8
C ... Heap
      integer w(1)
      common /w/ w

      call tcn ('vcdmel')
      rsml=0d0 !call dpzero(rsml,n0)
      ehl=0d0  !call dpzero(ehl,n0)
Cdelw1 warning(1) , probably  ocg  is not defined yet at linenumber= 66
Cdelw1 do not change ocg because of ['undef', 'wvar', 'wref']
      ocg=slat%ocg
Cdelw1 warning(1) , probably  ojcg  is not defined yet at linenumber= 67
Cdelw1 do not change ojcg because of ['undef', 'wvar', 'wref']
      ojcg=slat%ojcg
Cdelw1 warning(1) , probably  oidxcg  is not defined yet at linenumber= 68
Cdelw1 do not change oidxcg because of ['undef', 'wvar', 'wref']
      oidxcg=slat%oidxcg
      ocy=slat%ocy

      do  i = 1, nsite
        ib = isite(i)
        ncls = iclsn(i)
        lcls = iclsl(i)
        is = ssite(ib)%spec
        !i_copy_size=size(ssite(ib)%pnu)
        !call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
        !pnu(1:n0,1:2)= ssite(ib)%pnu(1:n0,1:2)
        call dcopy(n0*2,ssite(ib)%pnu,1,pnu,1)
        ov0=ssite(ib)%ov0
        !i_copy_size=size(ssite(ib)%pz)
        !call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
        !pnz(1:n0,1:2)= ssite(ib)%pz(1:n0,1:2)
        !do i_spacks=is,is
        !  call spacks_copy('u',sspec(i_spacks)%name,is,is,clabl,i_spacks)
        !enddo
        call dcopy(n0*2,ssite(ib)%pz,1,pnz,1)
        call r8tos8(sspec(is)%name,clabl)
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        z=sspec(is)%z
        lmxa=sspec(is)%lmxa
        if (lmxa .gt. lmxax) call rxi('vcdmel needs lmxax ',lmxa)
        if (lmxa .eq. -1) cycle !  if (lmxa .eq. -1) goto 10
Cdelw1         call defrr (orofi, nr)
        allocate(rv_w_orofi(nr))
Cdelw1         call radmsh (rmt,a,nr,w(orofi))
        call radmsh ( rmt , a , nr , rv_w_orofi ) 
C   --- Augmented wave functions u,s
Cdelw1         call defrr(oul,   nr*(lmxa+1)*nsp)
        allocate(rv_w_oul(nr*(lmxa+1)*nsp))
Cdelw1         call defrr(osl,   nr*(lmxa+1)*nsp)
        allocate(rv_w_osl(nr*(lmxa+1)*nsp))
Cdelw1         call defrr(ogz,   nr*(lmxa+1)*nsp)
        allocate(rv_w_ogz(nr*(lmxa+1)*nsp))
Cdelw1         call defrr(oruu,  nr*(lmxa+1)*2*nsp)
        allocate(rv_w_oruu(nr*(lmxa+1)*2*nsp))
Cdelw1         call defrr(orus,  nr*(lmxa+1)*2*nsp)
        allocate(rv_w_orus(nr*(lmxa+1)*2*nsp))
Cdelw1         call defrr(orss,  nr*(lmxa+1)*2*nsp)
        allocate(rv_w_orss(nr*(lmxa+1)*2*nsp))

        call uspecb(0,4,sspec,is,is,lh,rsml,ehl,j)
Cdelw1 warning(1) , probably  ov0  is not defined yet at linenumber= 105
Cdelw1 warning(1) , probably  ov0  is not defined yet at linenumber= 105
Cdelw1         call makusp(n0,z,nsp,rmt,lmxa,w(ov0),a,nr,xx,xx,pnu,pnz,rsml,
Cdelw1      .  ehl,w(oul),w(osl),w(ogz),w(oruu),w(orus),w(orss))
        call makusp ( n0 , z , nsp , rmt , lmxa , w ( ov0 ) , a , nr 
     .  , xx , xx , pnu , pnz , rsml , ehl , rv_w_oul , rv_w_osl , rv_w_ogz 
     .  , rv_w_oruu , rv_w_orus , rv_w_orss ) 

C   --- Matrix elements of u,s with core
        call skpblb(clabl,8,j)
        call awrit3('%N CLS atom %i ('//clabl(1:j+1)//') n=%i l=%i',
     .  ' ',120,i1mach(2),ib,ncls,lcls)
Cdelw1         call defrr(og,   nr*2)
         allocate(rv_w_og(nr*2))
         if (nr*2<0) rv_w_og(:)=0.0d0

Cdelw1 warning(1) , probably  ov0  is not defined yet at linenumber= 112
Cdelw1 warning(1) , probably  ov0  is not defined yet at linenumber= 112
Cdelw1         call pvcdm1(ncls,lcls,w(og),z,lmxa,w(ov0),a,nr,w(orofi),
Cdelw1      .  w(oul),w(osl),nsp,lmxax,ume(0,1,i),sme(0,1,i))
        call pvcdm1 ( ncls , lcls , rv_w_og , z , lmxa , w ( ov0 ) , 
     .  a , nr , rv_w_orofi , rv_w_oul , rv_w_osl , nsp , lmxax , ume 
     .  ( 0 , 1 , i ) , sme ( 0 , 1 , i ) ) 

Cdelw1 rlse name= orofi old_list= orofi oul osl ogz oruu orus orss og 
Cdelw1 rlse name= orofi new_list= (None)
Cdelw1         call rlse (orofi)
         deallocate(rv_w_og,rv_w_orss,rv_w_orus,rv_w_oruu,rv_w_ogz,rv_w_osl,rv_w_oul,rv_w_orofi)
c   10   continue
      enddo
C --- Open CLS weights file and write first line
Cdelw1       call defrr(os,3*ndimh*nsite*2)
       allocate(rv_w_os(3*ndimh*nsite*2))
       if (3*ndimh*nsite*2<0) rv_w_os(:)=0.0d0

      nl = sctrl%nl
      nfstg = 11
      nchan = 3*nsite
      ifi = fopna('cls',-1,4)
      j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,0,1,0,0,0,0,0,
     .0d0,0d0,0d0,0d0,0d0,0d0)
C --- For each qp, make <nk|x,y,z|core> at each site and save to disk in
C     iomoms-compatible format ---
      do   iq = 1, nq
        do  isp = 1, nsp
Cdelw1           call dpzero(w(os),3*ndimh*nsite*2)
           call dpzero ( rv_w_os , 3 * ndimh * nsite * 2 ) 

          do  i = 1, nsite
            lcls = iclsl(i)
            ib = isite(i)
            is = ssite(ib)%spec
            lmxa=sspec(is)%lmxa
            nlma = (lmxa+1)**2
            if (lmxa .gt. -1) then
Cdelw1 warning(1) , probably  ocg  is not defined yet at linenumber= 137
Cdelw1 warning(1) , probably  ojcg  is not defined yet at linenumber= 137
Cdelw1 warning(1) , probably  oidxcg  is not defined yet at linenumber= 137
Cdelw1 warning(1) , probably  ocg  is not defined yet at linenumber= 137
Cdelw1 warning(1) , probably  ojcg  is not defined yet at linenumber= 137
Cdelw1 warning(1) , probably  oidxcg  is not defined yet at linenumber= 137
Cdelw1               call pvcdm2(i,nsite,ndham,ndimh,nlma,nlmax,
Cdelw1      .        aus(1,1,1,isp,i,iq),ume(0,isp,i),sme(0,isp,i),
Cdelw1      .        lcls,w(ocg),w(ojcg),w(oidxcg),w(os))
               call pvcdm2 ( i , nsite , ndham , ndimh , nlma , nlmax , aus 
     .         ( 1 , 1 , 1 , isp , i , iq ) , ume ( 0 , isp , i ) , sme ( 0 
     .         , isp , i ) , lcls , w ( ocg ) , w ( ojcg ) , w ( oidxcg ) , 
     .         rv_w_os ) 

            endif
          enddo
C --- Scale weights arbitrarily by 100 for plotting etc ..
Cdelw1           call dscal(3*ndimh*nsite*2,1d2,w(os),1)
           call dscal ( 3 * ndimh * nsite * 2 , 1d2 , rv_w_os , 1 ) 

Cdelw1           j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,1,1,ndimh,ndimh,
Cdelw1      .    nchan,nchan,ndimh,evl(1,isp,iq),0d0,w(os),0d0,0d0,0d0)
           j = iomoms ( - ifi , nl , nsp , nspc , nq , ndimh , nfstg , 1 
     .     , 1 , 1 , ndimh , ndimh , nchan , nchan , ndimh , evl ( 1 , isp 
     .     , iq ) , 0d0 , rv_w_os , 0d0 , 0d0 , 0d0 ) 

        enddo
      enddo
      iq = nq*nsp + 1
      j = iomoms(-ifi,nl,nsp,nspc,nq,ndimh,nfstg,1,iq,1,ndimh,ndimh,
     .nchan,nchan,ndimh,0d0,0d0,0d0,0d0,ef,0d0)
Cdelw1 rlse name= os old_list= os 
Cdelw1 rlse name= os new_list= (None)
Cdelw1       call rlse(os)
       if (allocated(rv_w_os)) deallocate(rv_w_os)

      call fclose(ifi)
      call tcx ('vcdmel')
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ov0 ocg ojcg oidxcg]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine vcdmel


      subroutine pvcdm1(ncls,lcls,gcore,z,lmxa,v,a,nr,rofi,ul,sl,nsp,
     .lmxax,ume,sme)
C- Radial matrix elements < (u,s) | r | core >
C     implicit none
      integer ncls,lcls,lmxa,nr,nsp,lmxax
      double precision a,z,gcore(nr,2),rofi(1),v(nr,nsp),
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
        call awrit2(' vcdmel: ecor0=%,8d ecore=%,8d%N'
     .  //' (not including electrostatic potential shift)',
     .  ' ',120,i1mach(2),ecor0,ecore)
C      call prrmsh('gcore',rofi,gcore,nr,nr,1)
C      call prrmsh('ul',rofi,ul,nr,nr,1+lmxa)
C      call prrmsh('sl',rofi,sl,nr,nr,1+lmxa)

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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
      integer isite,lcls,ndham,ndimh,nlma,nlmax,nsite,indxcg(1),jcg(1)
      double precision cg(1),ume(0:1),sme(0:1),s(3,nsite,ndimh,2)
      double complex aus(nlmax,ndham,2)
C ... Local parameters
      integer kk(4),mlm,lm,ll,klm,ii,indx,icg1,icg2,icg,llm,ib,k
      double complex cxx
C     Transposes (y,z,x) to (x,y,z)
      data kk /0,2,3,1/

C      write (*,100)
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
C  100 format (2x,'n',2x,'L(val)',1x,'L(r)',1x,'L(core)',11x,'au',25x,
C     .        'as',14x,'CG',20x,'cxx')
C  200 format (1x,i3,3i5,2x,'(',g12.5,',',g12.5,')',
C     .  '(',g12.5,',',g12.5,')',g12.5,'(',g12.5,',',g12.5,')')

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine pvcdm2


