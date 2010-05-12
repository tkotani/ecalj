Cdelw1 var ['olst1', '0', 'olst1', 'nbas', 'integer', 'iv_w_', 'shoshl']
Cdelw1 var ['olst2', '0', 'olst2', 'nbas', 'integer', 'iv_w_', 'shoshl']
Cdelw1 var ['odist', '0', 'odist', 'nbas*3', 'real(8)', 'rv_w_', 'shoshl']
Cdelw1 var ['omark', '0', 'omark', 'nbas', 'real(8)', 'rv_w_', 'shoshl']
Cdelw1 var ['ortab', '1', 'ortab', '3*nttab', 'real(8)', 'rv_w_', 'shoshl']
Cdelw1 var ['oiax', '-1', 'undef', '0', 'undef', 'x', 'shoshl']
Cdelw1 var ['owk', '0', 'owk', 'mxnbr', 'real(8)', 'rv_w_', 'shoshl']
Cdelw1 var ['olstc', '-1', 'undef', '0', 'undef', 'x', 'shoshl']
Cdelw1 var ['oslist', '0', 'oslist', 'nbas', 'integer', 'iv_w_', 'shoang']
Cdelw1 var ['oblist', '0', 'oblist', 'nbas', 'integer', 'iv_w_', 'shoang']
Cdelw1 var ['oiax', '0', 'oiax', 'niax*mxnbr', 'integer', 'iv_w_', 'shoang']
Cdelw1 var ['owk', '0', 'owk', 'mxnbr', 'real(8)', 'rv_w_', 'shoang']
Cdelw1 var ['onum', '1', 'onum', '-nspec*npr(1)', 'integer', 'iv_w_', 'shoang']
Cdelw1 var ['oang', '1', 'oang', 'nspec**2*nshell**2*nmxshl**2', 'real(8)', 'rv_w_', 'shoang']
Cdelw1 var ['od', '1', 'od', 'nshell', 'real(8)', 'rv_w_', 'shoang']
Cdelw1 undel [['undef', 'wref', 'redef'], 'oiax', 'shoshl']
Cdelw1 undel [['undef', 'wref'], 'olstc', 'shoshl']
Cdelw1 undel [['wref'], 'oslist', 'shoang']
      subroutine shoshl(sopts,nbas,bas,plat,mxnbr0,z,slabl,dclabl,ips,
     .  ipc,ves,eula,nclass)

C- Print nearest-neighbor shells
C ----------------------------------------------------------------
Ci Inputs
Ci   Everything is input
Ci   sopts: a set of modifiers, with the syntax
Ci          [:v][:e][:r=#][:sites:site-list][:pairs:pair-list] ..
Ci          [:tab[=#]][:disp=fnam][:nn][:fn=fnam]
Ci          :
Ci          :r=# sets range for shells
Ci          :v prints out electrostatic potential
Ci          :e prints out Euler angles
Ci          :r=# restricts neighbor table to range #
Ci          :sites:list print table only for sites in site-list
Ci          :pairs:pair-list print table only for pairs in pair-list
Ci          :tab prints out neighbor table
Ci          :tab=# prints out neighbor table in style #
Ci          :  ...The following only apply to :tab
Ci          :  fn=fnam write table to file fnam
Ci          :  disp=fnam : read another site positions file; neighbor
Ci          :              table for both original and displaced
Ci          :              positions is written.
Ci          :  nn only print first entry for a given pair site
Ci          :     in neighbor table
Ci  OLD: doesn't work
Ci          :i[=style-#]:list  restricts neighbors in shell to list.
Ci                             This must be the last modifier.
Cu Updates
Cu   19 Apr 03 Can write displaced neigbhor table
Cu   12 Apr 03 Can write neigbhor table
Cr   24 Nov 97 changed modifier list
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nclass,mxnbr0,ips(nbas),ipc(nbas),nttab
      double precision bas(3,nbas),plat(3,3),dclabl(nclass),ves(1),
     .  eula(1),z(1)
      character sopts*(*),slabl(nbas)*8
C Local parameters
      integer niax
      parameter (niax=10)
      logical lves,leula,a2bin
      double precision avwsr,avw,range,xx,stdo
Cdelw1 do not change  [['undef', 'wref'], 'olstc', 'shoshl']
Cdelw1       integer owk,ontab,ortab,olstc,olst1,olst2,odist,omark
      integer:: ontab , olstc
      real(8) ,allocatable :: rv_w_owk(:)
      real(8) ,allocatable :: rv_w_ortab(:)
      integer ,allocatable :: iv_w_olst1(:)
      integer ,allocatable :: iv_w_olst2(:)
      real(8) ,allocatable :: rv_w_odist(:)
      real(8) ,allocatable :: rv_w_omark(:)

Cdelw1 do not change  [['undef', 'wref', 'redef'], 'oiax', 'shoshl']
      integer npr(2),mxnbr,ib,ic,oiax,iclbsj,j,j1,j2,lstyle,ival,
     .  scrwid,ltab,nlstc,nlst1,nlst2,mxcsiz,iprint,ifi,fopn,lgunit,
     .  ldisp,lnn
      character*8  dc*1
      parameter (scrwid=120)
C heap:
      integer w(1)
      common /w/ w

C --- Parse modifiers ---
      stdo = lgunit(1)
      ifi = stdo
      lves  = .false.
      leula = .false.
      lstyle = 1
      range = 2
      nlstc = 0
Cdelw1 warning(1) , probably  olstc  is not defined yet at linenumber= 64
Cdelw1 do not change olstc because of ['undef', 'wref']
      olstc = 1
      ltab = 0
      nlst1 = 0
      nlst2 = 0
Cdelw1       call defi(olst1,nbas)
      allocate(iv_w_olst1(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_olst1(:)=0

Cdelw1       call defi(olst2,nbas)
      allocate(iv_w_olst2(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_olst2(:)=0

      lnn = 0
      ldisp = 0
Cdelw1       call defrr(odist,nbas*3)
      allocate(rv_w_odist(nbas*3))
Ctakao_ZeroClear_NotRequiered        if (nbas*3<0) rv_w_odist(:)=0.0d0

Cdelw1       call defrr(omark,nbas)
      allocate(rv_w_omark(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) rv_w_omark(:)=0.0d0

Cdelw1       call dcopy(nbas*3,bas,1,w(odist),1)
      call dcopy ( nbas * 3 , bas , 1 , rv_w_odist , 1 )


      if (sopts .ne. ' ') then
C       ls = len(sopts)
        j1 = 1
        dc = sopts(j1:j1)
        j1 = j1+1

C   ... Return here to resume parsing for arguments
   40   continue
        call nwordg(sopts,0,dc//' ',1,j1,j2)

C   ... Parse special arguments
        if (j2 .ge. j1) then
C         print *, sopts(j1:j2)
          if (sopts(j1:j2) .eq. 'v')  then
            lves = .true.

          elseif (sopts(j1:j2) .eq. 'e')  then
            leula = .true.

          elseif (sopts(j1:j1+1) .eq. 'r=') then
            j = 0
            if (.not. a2bin(sopts(j1+2:),range,4,0,' ',j,j2-j1-2))
     .        goto 999

          elseif (sopts(j1:j1+2) .eq. 'fn=')  then
            if (j1+3 .le. j2) ifi = fopn(sopts(j1+3:j2))

          elseif (sopts(j1:j1+3) .eq. 'tab=')  then
            j = 0
            if (.not. a2bin(sopts(j1+4:),ltab,2,0,' ',j,j2-j1-4))
     .        goto 999

          elseif (sopts(j1:j1+4) .eq. 'disp=')  then
Cdelw1             call iopos(.false.,-1,sopts(j1+5:j2),nbas,w(odist))
            call iopos ( .false. , - 1 , sopts ( j1 + 5:j2 ) , nbas , rv_w_odist 
     .       )

            ldisp = 1

          elseif (sopts(j1:j1+1) .eq. 'nn')  then
            lnn = 1

          elseif (sopts(j1:j1+2) .eq. 'tab')  then
            ltab = 1

          elseif (sopts(j1:j1+4) .eq. 'sites')  then
Cdelw1             call baslst(0,11,sopts(j1+5:),j,ips,nbas,slabl,z,0,
Cdelw1      .        ' ',w,nlst1,w(olst1))
            call baslst ( 0 , 11 , sopts ( j1 + 5: ) , j , ips , nbas , slabl 
     .       , z , 0 , ' ' , w , nlst1 , iv_w_olst1 )

            j2 = j1+5+j-2

          elseif (sopts(j1:j1+4) .eq. 'pairs')  then
Cdelw1             call baslst(0,11,sopts(j1+5:),j,ips,nbas,slabl,z,0,
Cdelw1      .        ' ',w,nlst2,w(olst2))
            call baslst ( 0 , 11 , sopts ( j1 + 5: ) , j , ips , nbas , slabl 
     .       , z , 0 , ' ' , w , nlst2 , iv_w_olst2 )

            j2 = j1+5+j-2

          else
            call rxs('shoshl: failed to parse --shell switch: ',
     .        sopts(j1:j2))
            goto 999
          endif
          j1 = j2+2
          goto 40
        endif
      endif

C --- Print neighbor table for each site ---
      if (ltab .ne. 0) then
        mxcsiz = 0
        call pshpr(iprint()-20)
Cdelw1 warning(1) , probably  oiax  is not defined yet at linenumber= 142
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
        call pairs(nbas,nbas,1d0,plat,range/2*(1+1d-6),bas,
     .    -1,3,-1,w,nttab,ontab,oiax,mxcsiz)
        call poppr
Cdelw1         call defrr(ortab,3*nttab)
        allocate(rv_w_ortab(3*nttab))
Ctakao_ZeroClear_NotRequiered          if (3*nttab<0) rv_w_ortab(:)=0.0d0

Cdelw1 warning(1) , probably  oiax  is not defined yet at linenumber= 146
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  oiax  is not defined yet at linenumber= 146
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
Cdelw1         call mkrtab(000,1d0,plat,bas,w(oiax),nttab,bas,w(ortab))
        call mkrtab ( 000 , 1d0 , plat , bas , w ( oiax ) , nttab , bas 
     .   , rv_w_ortab )

C       if (ldisp .eq. 1) ltab = 2
Cdelw1 warning(1) , probably  oiax  is not defined yet at linenumber= 148
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  oiax  is not defined yet at linenumber= 148
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
Cdelw1         call psho1(ldisp*10+ltab,lnn,w(oiax),nbas,nttab,nlst1,w(olst1),
Cdelw1      .    nlst2,w(olst2),w(ortab),plat,bas,w(odist),w(omark),ipc,dclabl,
Cdelw1      .    ifi)
        call psho1 ( ldisp * 10 + ltab , lnn , w ( oiax ) , nbas , nttab 
     .   , nlst1 , iv_w_olst1 , nlst2 , iv_w_olst2 , rv_w_ortab , plat 
     .   , bas , rv_w_odist , rv_w_omark , ipc , dclabl , ifi )

Cdelw1 rlse name= ontab old_list= olst1 olst2 odist omark ortab

Cki#error ERROR, try to release name= ontab ,but list does not have  ontab at linenumber= 151 list= olst1 olst2 odist omark ortab

Cdelw1 rlse name= ontab new_list= olst1 olst2 odist omark ortab
Cdelw1         call rlse(ontab)

Cki#error, have return with len(w_varlist)>0 at line 152
Cdelw1 w_varlist remains: olst1 olst2 odist omark ortab
Cdelw1 w_varlistundel: [oiax olstc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: [olst1 olst2 odist omark ortab]
Cdelw1         return
        if (allocated(rv_w_ortab)) deallocate(rv_w_ortab)
        if (allocated(rv_w_omark)) deallocate(rv_w_omark)
        if (allocated(rv_w_odist)) deallocate(rv_w_odist)
        if (allocated(iv_w_olst2)) deallocate(iv_w_olst2)
        if (allocated(iv_w_olst1)) deallocate(iv_w_olst1)
        return

      endif

C --- Show shells for each class ---
      if (mxnbr0 .eq. 0) then
        mxnbr = 100*range**3
      else
        mxnbr = mxnbr0
      endif
Cdelw1 w_varlistall_append_uniq ERROR, define again with another type: oiax undef and integer
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
      call defi(oiax,niax*mxnbr)
Cdelw1       call defdr(owk,mxnbr)
      allocate(rv_w_owk(mxnbr))
Ctakao_ZeroClear_NotRequiered        if (mxnbr<0) rv_w_owk(:)=0.0d0

      avw = avwsr(plat,1d0,xx,nbas)
      avw = 1
c     call pshprt(50)
      do  20  ic = 1, nclass
        ib = iclbsj(ic,ipc,nbas,1)
        if (nlst1 .gt. 0) then
Cdelw1           call hunti(w(olst1),nlst1,ib,0,j)
          call hunti ( iv_w_olst1 , nlst1 , ib , 0 , j )

Cdelw1           if (ival(w(olst1),j+1) .ne. ib) goto 20
          if ( ival ( iv_w_olst1 , j + 1 ) .ne. ib ) goto 20

        endif
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
Cdelw1         call nghbor(nbas,plat,bas,range,range,ib,
Cdelw1      .              mxnbr,npr,w(oiax),w(owk))
        call nghbor ( nbas , plat , bas , range , range , ib , mxnbr 
     .   , npr , w ( oiax ) , rv_w_owk )

Cdelw1 do not change oiax because of ['undef', 'wref', 'redef']
Cdelw1 warning(1) , probably  olstc  is not defined yet at linenumber= 174
Cdelw1 do not change olstc because of ['undef', 'wref']
        call xxsho(npr(1),nbas,plat,bas,w(oiax),ipc,dclabl,nlstc,
     .    w(olstc),lves,ves,leula,eula,z)
   20 continue
Cdelw1 rlse name= oiax old_list= ortab omark odist olst2 olst1 oiax owk
Cdelw1 rlse name= oiax new_list= ortab omark odist olst2 olst1
Cdelw1 not deallocate oiax because of [undef wref redef]
Cdelw1       call rlse(oiax)
      if (allocated(rv_w_owk)) deallocate(rv_w_owk)
      call rlse(oiax)

Cki#error, have return with len(w_varlist)>0 at line 178
Cdelw1 w_varlist remains: ortab omark odist olst2 olst1
Cdelw1 w_varlistundel: [oiax olstc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: [ortab omark odist olst2 olst1]
Cdelw1       return
      if (allocated(iv_w_olst1)) deallocate(iv_w_olst1)
      if (allocated(iv_w_olst2)) deallocate(iv_w_olst2)
      if (allocated(rv_w_odist)) deallocate(rv_w_odist)
      if (allocated(rv_w_omark)) deallocate(rv_w_omark)
      if (allocated(rv_w_ortab)) deallocate(rv_w_ortab)
      return


  999 call rxs('shoshl: failed to parse ',sopts)
Cgetarg       end
Cdelw1 w_varlist remains: olst1 olst2 odist omark ortab
Cdelw1 w_varlistundel: [oiax olstc]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: [olst1 olst2 odist omark ortab]
      if (allocated(rv_w_ortab)) deallocate(rv_w_ortab)
      if (allocated(rv_w_omark)) deallocate(rv_w_omark)
      if (allocated(rv_w_odist)) deallocate(rv_w_odist)
      if (allocated(iv_w_olst2)) deallocate(iv_w_olst2)
      if (allocated(iv_w_olst1)) deallocate(iv_w_olst1)

      end subroutine shoshl

      subroutine xxsho(npr,nbas,plat,bas,iax,ipc,dclabl,nlstc,lstc,lves,
     .  ves,leul,eula,z)

C- Kernel called by shoshl
C  nlstc,lstc:  a list of classes to include as pairs (nlstc>0)
C     implicit none
      logical lves,leul
      integer npr,nbas,niax,ipc(1),nlstc,lstc(nlstc)
      parameter (niax=10)
      integer iax(niax,1)
      double precision plat(3,3),bas(3,1),dclabl(1),ves(1),eula(nbas,3),
     .  z(32)
      integer ih(2,120),scrwid
      parameter (scrwid=120)
      integer i,l,ishell,nshell,j,k,ii,kk,ic,jc,i1,lgunit,awrite,iclbsj,
     .  ib
      double precision dr(3),d,drr2,dshell,fuzz,z1(3),z2(3),alfa,beta,
     .  angle,pi,ddot
      character*8 clabl,outs1*25,outs2*(scrwid),outsv*(scrwid),
     .  outse*(scrwid)

      pi = 4*datan(1d0)
      fuzz = 1d-3
      dshell = 0
      nshell = 0
      ishell = 1
      if (leul) then
        alfa = eula(iax(1,1),1)
        beta = eula(iax(1,1),2)
        z1(1) = dcos(alfa)*dsin(beta)
        z1(2) = dsin(alfa)*dsin(beta)
        z1(3) = dcos(beta)
      endif
      ic = ipc(iax(1,1))
      call r8tos8(dclabl(ic),clabl)
      print 302, clabl, ic, nint(z(ic))
  302 format(/' Shell decomposition for class ',a,
     .        '  class',i4,'  z=',i2/
     .        ' shell   d     nsh csiz  class ...')

      do  10  i = 1, npr
        d = dsqrt(drr2(plat,bas(1,iax(1,i)),bas(1,iax(2,i)),
     .    iax(3,i),iax(4,i),iax(5,i),dr))
C   ... new shell, or close of last shell
        if (dabs(d-dshell) .gt. fuzz .or. i .eq. npr) then
          i1 = i-1
          if (i .eq. npr) i1 = i
          nshell = nshell+1
          write (outs1,301) nshell, dshell, i1+1-ishell, i1
  301     format(i4,f10.6,i4,i5,2x)
          call iinit(ih,2*(i-ishell))
C     ... ii is the number of different classes in this shell
          ii = 0
          do  12  j = ishell, i1
            ic = ipc(iax(2,j))
C       ... See whether already found one of these or if not in list
            kk = 0
            if (nlstc .gt. 0) then
              kk = -1
              do  15  jc = 1, nlstc
                if (lstc(jc) .gt. ic) goto 17
                if (lstc(jc) .eq. ic) kk = 0
   15         continue
   17         continue
            endif
            if (kk .eq. 0) then
              do  14  k = 1, ii
                if (ih(2,k) .ne. ic) goto 14
                kk = k
   14         continue
            endif
C       ... We haven't --- increment ii and add this one
            if (kk .eq. 0) then
              ii = ii+1
              kk = ii
              ih(2,kk) = ic
            endif
C       ... Increment number of occurrences of this species
            ih(1,kk) = ih(1,kk)+1
   12     continue

C     ... Setup for printout
          outs2 = ' '
          outsv = ' '
          outse = ' '
          kk = 0
          do  16  k = 1, ii
            kk = kk+1
            call r8tos8(dclabl(ih(2,k)),clabl)
            l = awrite('%a  %np%i:'//clabl//
     .        '%a%?;n>1;(%i);%j;',outs2,len(outs2),0,
     .        (kk-1)*14,ih(2,k),ih(1,k),ih(1,k),0,0,0,0)
            if (lves) call awrit2('%np%d',outsv,len(outsv),0,
     .        (kk-1)*14,ves(ih(2,k)))
            ib = iclbsj(ih(2,k),ipc,-nbas,1)
            if (leul .and. ib .gt. 0) then
              alfa = eula(ib,1)
              beta = eula(ib,2)
              z2(1) = dcos(alfa)*dsin(beta)
              z2(2) = dsin(alfa)*dsin(beta)
              z2(3) = dcos(beta)
              angle = dacos(max(-1d0,min(1d0,ddot(3,z1,1,z2,1))))
              if (angle .gt. pi) angle = angle - 2*pi
              call awrit2('%np%d',outse,len(outse),0,(kk-1)*14,angle)
            endif
            if (l .gt. scrwid-35) then
              call awrit0(outs1//outs2,' ',-scrwid,lgunit(1))
              if(lves) call awrit0('v%26f'//outsv,' ',-scrwid,lgunit(1))
              if(leul) call awrit0('e%26f'//outse,' ',-scrwid,lgunit(1))
              kk = 0
              outs2 = ' '
              outsv = ' '
            endif
   16     continue
          if (outs2 .ne. ' ') then
            call awrit0(outs1//outs2,' ',-scrwid,lgunit(1))
            if (lves) call awrit0('v%26f'//outsv,' ',-scrwid,lgunit(1))
            if (leul) call awrit0('e%26f'//outse,' ',-scrwid,lgunit(1))
          endif
          outs1 = ' '

          ishell = i
          dshell = d
        endif
   10 continue
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine xxsho

      subroutine shoang(sopts,nbas,bas,plat,mxnbr0,slabl,ips)

C- Print bond angles
C ----------------------------------------------------------------
Ci Inputs
Ci   sopts :a set of modifiers, with the syntax
Ci         :  [:r=#][:spec=spec-list][:spec=spec-list]
Ci         :  :r=# sets range for shells
Ci         :  :sites=site-list collects angles only for sites within list
Ci         :  :spec=list       prints angles for bonds connecting to
Ci                             species in list
Ci   bas   :basis vectors, in units of alat
Ci   plat  :primitive lattice vectors, in units of alat
Ci   mxnbr0
Ci   slabl :list of species labels
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   nbas  :size of basis
Cu Updates
Cu   13 Sep 01 Added options sopts.  Altered argument list.
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,mxnbr0
      double precision bas(3,nbas),plat(3,3),slabl(*)
      integer ips(nbas)
      character sopts*(*)
C ... Local parameters
      double precision avwsr,avw,range,xx,zz
Cdelw1 do not change  [['wref'], 'oslist', 'shoang']
Cdelw1       integer npr(2),mxnbr,ib,nshell,nmxshl,niax,j,j1,j2,
Cdelw1      .  nspec,mxint,oiax,owk,onum,oang,od,nsites,oslist,nbonds,oblist,
Cdelw1      .  ilst,parg,m,iv(10),ilast
      integer:: npr(2) , mxnbr , ib , nshell , nmxshl , niax , j , 
     . j1 , j2 , nspec , mxint , nsites , oslist , nbonds , ilst , parg 
     . , m , iv(10) , ilast
      integer ,allocatable :: iv_w_oiax(:)
      real(8) ,allocatable :: rv_w_owk(:)
      integer ,allocatable :: iv_w_onum(:)
      real(8) ,allocatable :: rv_w_oang(:)
      real(8) ,allocatable :: rv_w_od(:)
      integer ,allocatable :: iv_w_oblist(:)

      parameter (niax=10)
      character dc*1
C ... Heap
      integer w(1)
      common /w/ w
C ... External calls
      external defdr,defi,nghbor,nwordg,pvang1,pvang2,rlse,rxs

C ... Setup
      range = 2.5d0
      nspec = mxint(nbas,ips)
      nsites = 0
      nbonds = 0
Cdelw1 do not change oslist because of ['wref']
      call defi(oslist,nbas)
Cdelw1       call defi(oblist,nbas)
      allocate(iv_w_oblist(nbas))
Ctakao_ZeroClear_NotRequiered        if (nbas<0) iv_w_oblist(:)=0


C ... Switches
      dc = sopts(1:1)
      if (dc .ne. ' ') then
        j2 = 0
C   ... Return here to resume parsing for arguments
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) .eq. dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 .ge. j1) then
          if (.false.) then

C         range
          elseif (sopts(j1:j1+1) .eq. 'r=')  then
            m = 0
            j = parg('r=',4,sopts(j1:),m,len(sopts(j1:)),
     .        dc//' ',1,1,iv,range)
            if (j .le. 0) goto 999

C         Site list
          elseif (sopts(j1:j1+4) .eq. 'sites') then
            if (sopts(j1+5:j1+5) .eq. '=') sopts(j1+5:j1+5) = dc
Cdelw1 do not change oslist because of ['wref']
            call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zz,0,' ',
     .        xx,nsites,w(oslist))

C         Bond list
          elseif (sopts(j1:j1+4) .eq. 'bonds') then
            if (sopts(j1+5:j1+5) .eq. '=') sopts(j1+5:j1+5) = dc
Cdelw1             call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zz,0,' ',
Cdelw1      .        xx,nbonds,w(oblist))
            call baslst ( 0 , 10 , sopts ( j1 + 5: ) , ilast , ips , nbas 
     .       , slabl , zz , 0 , ' ' , xx , nbonds , iv_w_oblist )

          endif
          goto 10

        endif
      endif

      if (mxnbr0 .eq. 0) then
        mxnbr = 2*range**3
      else
        mxnbr = mxnbr0
      endif
Cdelw1       call defi(oiax,niax*mxnbr)
      allocate(iv_w_oiax(niax*mxnbr))
Ctakao_ZeroClear_NotRequiered        if (niax*mxnbr<0) iv_w_oiax(:)=0

Cdelw1       call defdr(owk,mxnbr)
      allocate(rv_w_owk(mxnbr))
Ctakao_ZeroClear_NotRequiered        if (mxnbr<0) rv_w_owk(:)=0.0d0

      avw = avwsr(plat,1d0,xx,nbas)

C --- For every site in list, generate tables of bond angles ---
      ilst = 0
      do  ib = 1, nbas
        if (nsites .gt. 0) then
Cdelw1 do not change oslist because of ['wref']
          if (w(oslist+ilst) .ne. ib) goto 31
        endif
        ilst = ilst+1

C   ... Get neighbor lists

Cdelw1         call nghbor(nbas,plat,bas,range*avw,range*avw,ib,
Cdelw1      .              mxnbr,npr,w(oiax),w(owk))
        call nghbor ( nbas , plat , bas , range * avw , range * avw , 
     .   ib , mxnbr , npr , iv_w_oiax , rv_w_owk )


C   ... Get shell dimensions
Cdelw1         call defi(onum,-nspec*npr(1))
        allocate(iv_w_onum(abs(-nspec*npr(1))))
        if (-nspec*npr(1)<0) iv_w_onum(:)=0

Cdelw1         call pvang1(npr(1),nbas,plat,bas,w(oiax),ips,w(onum),nshell,
Cdelw1      .    nmxshl)
        call pvang1 ( npr ( 1 ) , nbas , plat , bas , iv_w_oiax , ips 
     .   , iv_w_onum , nshell , nmxshl )

Cdelw1 rlse name= onum old_list= oslist oblist oiax owk onum
Cdelw1 rlse name= onum new_list= oslist oblist oiax owk
Cdelw1         call rlse(onum)
        if (allocated(iv_w_onum)) deallocate(iv_w_onum)

Cdelw1         call defi(onum, -nspec**2*nshell**2)
        allocate(iv_w_onum(abs(-nspec**2*nshell**2)))
        if (-nspec**2*nshell**2<0) iv_w_onum(:)=0

Cdelw1         call defdr(oang, nspec**2*nshell**2*nmxshl**2)
        allocate(rv_w_oang(nspec**2*nshell**2*nmxshl**2))
Ctakao_ZeroClear_NotRequiered          if (nspec**2*nshell**2*nmxshl**2<0) rv_w_oang(:)=0.0d0

Cdelw1         call defdr(od,   nshell)
        allocate(rv_w_od(nshell))
Ctakao_ZeroClear_NotRequiered          if (nshell<0) rv_w_od(:)=0.0d0


C   ... Print bond angles
Cdelw1         call pvang2(npr(1),nbas,nspec,nshell,nmxshl,plat,bas,w(oiax),
Cdelw1      .    ips,slabl,nbonds,w(oblist),w(onum),w(oang),w(od))
        call pvang2 ( npr ( 1 ) , nbas , nspec , nshell , nmxshl , plat 
     .   , bas , iv_w_oiax , ips , slabl , nbonds , iv_w_oblist , iv_w_onum 
     .   , rv_w_oang , rv_w_od )

Cdelw1 rlse name= onum old_list= oslist oblist oiax owk onum oang od
Cdelw1 rlse name= onum new_list= oslist oblist oiax owk
Cdelw1         call rlse(onum)
        if (allocated(rv_w_od)) deallocate(rv_w_od)
        if (allocated(rv_w_oang)) deallocate(rv_w_oang)
        if (allocated(iv_w_onum)) deallocate(iv_w_onum)

   31   continue
      enddo
Cdelw1 rlse name= oslist old_list= oslist oblist oiax owk
Cdelw1 rlse name= oslist new_list= (None)
Cdelw1 not deallocate oslist because of [wref]
Cdelw1       call rlse(oslist)
      if (allocated(rv_w_owk)) deallocate(rv_w_owk)
      if (allocated(iv_w_oiax)) deallocate(iv_w_oiax)
      if (allocated(iv_w_oblist)) deallocate(iv_w_oblist)
      call rlse(oslist)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oslist]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      return

  999 call rxs('shoang: failed to parse ',sopts)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oslist]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine shoang


      subroutine pvang1(npr,nbas,plat,bas,iax,ips,num,nshell,nmxshl)

C- Help routine for shoang
C ----------------------------------------------------------------------
Ci Inputs
Ci   npr   :number of pairs in neighbor table
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   iax   :neighbor table containing pair information for one site.
Ci         :Table must be sorted by increasing distance from iax(1)
Ci   ips   :species table: site ib belongs to species ips(ib)
Co Outputs
Co   num   :(ishell,is) number of pairs in shell ishell of species is
Co   nshell:number of shells
Co   nmxshl:max value of num
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npr,nbas,nshell,nmxshl,niax,ips(1),num(npr,1)
      parameter (niax=10)
      integer iax(niax,1)
      double precision plat(3,3),bas(3,1)
C ... Local parameters
      integer i,is
      double precision d,wk(0:3),tol

      tol = 1d-6
C --- Get no. of shells and max no. of atoms in 1 shell and 1 class ---
      nshell = 0
      nmxshl = 0
      d = 0d0
      do  10  i = 2, npr
        is = ips(iax(2,i))
        call dlmn(nbas,plat,bas,iax(1,i),wk)
C       distance changed by more than tol ... new shell
        if (dabs(wk(0)-d) .gt. tol) then
          nshell = nshell + 1
          d = wk(0)
        endif
        num(nshell,is) = num(nshell,is) + 1
        nmxshl = max0(nmxshl,num(nshell,is))
   10 continue

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine pvang1


      subroutine pvang2(npr,nbas,nspec,nshell,nmxshl,plat,bas,iax,
     .  ips,slabl,nbonds,blist,num,ang,d)

C- Kernel called by shoang
C ----------------------------------------------------------------------
Ci Inputs
Ci   npr   :number of neighbors connecting site ib=iax(1,1)
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   nshell:number of shells
Ci   nmxshl:dimensions ang
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   slabl :struct containing global strings
Ci   num
Co Outputs
Co   ang   :table of angles
Co   d     :table of distances for each shell
Co   Angles and distances are printed out
Cu Updates
Cu   13 Sep 01
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npr,nbas,nspec,nshell,nmxshl,niax,ips(1),
     .  num(nspec,nspec,nshell,nshell),nbonds,blist(nbonds)
      parameter (niax=10)
      integer iax(niax,1)
      character*8 slabl(*)
C ... Local parameters
      double precision plat(3,3),bas(3,1),d(nshell),
     .  ang(nmxshl**2,nspec,nspec,nshell,nshell)
      integer i,j,n,is,js,nsh1,nsh2,nmx2,k
      double precision rdtodg,d1,d2,dp,ddot,wk1(0:3),wk2(0:3)
C ... External calls
      external dlmn,rxx

      nmx2 = nmxshl**2
      rdtodg = 45d0 / datan(1.d0)

C --- Accumulate bond angles by shell and class ---
      nsh1 = 0
      d1 = 0d0
      do  20  i = 2, npr
        is = ips(iax(2,i))
        call dlmn(nbas,plat,bas,iax(1,i),wk1)
        if (dabs(wk1(0)-d1) .gt. 1d-6) then
          nsh1 = nsh1 + 1
          d1 = wk1(0)
          d(nsh1) = d1
        endif
        nsh2 = nsh1
        d2 = d1
        if (nbonds .gt. 0) then
          k = 0
          call hunti(blist,nbonds,iax(2,i),0,k)
          if (k .ge. nbonds) goto 20
          if (blist(k+1) .ne. iax(2,i)) goto 20
        endif
        do  10  j = i+1, npr
          js = ips(iax(2,j))
          call dlmn(nbas,plat,bas,iax(1,j),wk2)
          if (dabs(wk2(0)-d2) .gt. 1d-6) then
            nsh2 = nsh2 + 1
            d2 = wk2(0)
          endif
          if (nbonds .gt. 0) then
            k = 0
            call hunti(blist,nbonds,iax(2,j),0,k)
            if (k .ge. nbonds) goto 10
            if (blist(k+1) .ne. iax(2,j)) goto 10
          endif
          dp = ddot(3,wk1(1),1,wk2(1),1)
          if (dp .gt.  1d0) dp =  1d0
          if (dp .lt. -1d0) dp = -1d0
          if (nsh1 .eq. nsh2 .and. js .lt. is) then
            num(js,is,nsh1,nsh2) = num(js,is,nsh1,nsh2) + 1
            n = num(js,is,nsh1,nsh2)
            ang(n,js,is,nsh1,nsh2) = rdtodg*dacos(dp)
          else
            num(is,js,nsh1,nsh2) = num(is,js,nsh1,nsh2) + 1
            n = num(is,js,nsh1,nsh2)
            ang(n,is,js,nsh1,nsh2) = rdtodg*dacos(dp)
          endif
          call rxx(n .gt. nmx2,'PVANG2: num gt nmx2')
   10   continue
   20 continue
      call rxx(nsh1 .ne. nshell,'PVANG2: nsh1 ne nshell')

C --- Printout ---
      print 400, iax(1,1), slabl(ips(iax(1,1)))
  400 format(/' Bond angles for site',i4,', species ',a/
     .  ' shl1    d1    shl2    d2     cl1      cl2       angle(s) ...')

      do  60  nsh1 = 1, nshell
        do  50  nsh2 = nsh1, nshell
          do  40  is = 1, nspec
            do  30  js = 1, nspec
              n = num(is,js,nsh1,nsh2)
              if (n .ne. 0) then
                print 401, nsh1, d(nsh1), nsh2, d(nsh2), slabl(is),
     .            slabl(js), (ang(i,is,js,nsh1,nsh2), i = 1, n)
              endif
   30       continue
   40     continue
   50   continue
   60 continue

  401 format(2(1x,i3,1x,f9.6),1x,2(1x,a8),20(4(1x,f7.2):/47x))

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine pvang2

      subroutine psho1(ltab,lnn,iax,nbas,nttab,nlst1,lst1,nlst2,lst2,
     .  rtab,plat,pos,pos2,mark,ipc,dclabl,ifi)

C- Kernel called by supcel to displace pairs radially
C ----------------------------------------------------------------------
Ci Inputs
Ci   ltab  :style which to print out neighbor table
Ci         :1s digit
Ci         :0 do nothing ;return.  Else print table as:
Ci         :1 (standard mode)
Ci         :    ib jb dpos(1,jb) dpos(2,jb) dpos(3,jb)
Ci         :2  (just the positions)
Ci         :   dpos(1,jb) dpos(2,jb) dpos(3,jb)
Ci         :3 (sparse matrix format)
Ci         :     1 jb dpos(1,jb)
Ci         :     2 jb dpos(2,jb)
Ci         :     3 jb dpos(3,jb)
Ci         :10s digit
Ci         :1 print out neighbor table for pos and
Ci         :  displaced pos2 as well
Ci         :
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nbas  :size of basis
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   nlst1 :number of sites of "center" type
Ci   lst1  :list of sites of "center" type
Ci   nlst2 :number of sites of "neighbor" type
Ci   lst2  :list of sites of "neighbor" type
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   pos2  :displaced basis vectors (ltab >= 10)
Ci   mark  :work array of dimension nbas
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   dclabl:class name, packed as a real number
Ci   ifi   :file handle
Co Inputs/Outputs
Cu Updates
Cu   08 Aug 07 case ltab=12: allow for numerical imprecision in vector
Cu   19 Apr 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ltab,nlst1,lst1(nlst1),nlst2,lst2(nlst2),ifi,nbas,lnn
      integer niax,nttab,ipc(*),mark(nbas)
      double precision pos(3,*),pos2(3,*),rtab(3,*),dclabl(*),plat(9)
      parameter (niax=10)
      integer iax(niax,1)
C ... Local parameters
      integer iat,low1,low2,i1,i2,ipr,stdo,lgunit,i,oldi1,oldiat,
     .  ic,jc
      character clabl*8, clabl2*8
      double precision d,ddot,dsqrt,dpos(3),qlat(9),tol
C     logical latvec

      tol = 1d-6

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (ltab .eq. 0) return
      call mkqlat(plat,qlat,d)

      call info5(20,1,1,' ... shoshl: making neighbor list for'//
     .  ' %?#n#%-1j%i%j#%i# atom(s), style %i',nlst1,nbas,ltab,0,0)

      if (ltab .eq. 3)
     .  call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)

      call getpr(ipr)
      ipr = 100
      stdo = lgunit(1)
      call isanrg(mod(ltab,10),1,3,' shoshl:','tab',.true.)

      if (ltab .eq. 3)
     .  call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)
      if (ltab .eq. 13)
     .  call awrit1('%% rows 6 cols %i sparse',' ',80,ifi,nbas)

      low1 = 0
      low2 = 0
      oldiat = 0
      oldi1 = 0
      do  iat = 1, nttab

        i1 = iax(1,iat)
        i2 = iax(2,iat)
C   ... If site i1 isn't in the supplied list, skip this pair
        if (nlst1 .ne. 0) then
          call hunti(lst1,nlst1,i1,0,low1)
          if (i1 .ne. lst1(low1+1)) goto 10
        endif

C   ... New shell
        if (i1 .ne. oldi1) then
          call iinit(mark,nbas)
          oldi1 = i1
          oldiat = iat
          ic = ipc(i1)
          call r8tos8(dclabl(ic),clabl)
          call awrit1('# neighbor list for site %i, class '//clabl,
     .      ' ',80,ifi,i1)
C         Don't print out on-site entry
C         goto 10
        endif

C   ... If site i2 isn't in the supplied lst2, skip this pair
        if (nlst2 .ne. 0) then
          call hunti(lst2,nlst2,i2,0,low2)
          if (i2 .ne. lst2(low2+1)) goto 10
        endif

C   ... If i2 already marked, skip this pair
        if (mark(i2) .ne. 0) goto 10

        if (mod(ltab,10) .eq. 1) then
          jc = ipc(i2)
          call r8tos8(dclabl(jc),clabl2)
          d = dsqrt(ddot(3,rtab(1,iat),1,rtab(1,iat),1))
          if (ltab .lt. 10) then
            write(ifi,331) i1,i2,(rtab(i,iat),i=1,3),d,clabl,clabl2
          else
            do  i = 1, 3
              dpos(i) = pos2(i,i2)-pos(i,i2)
            enddo
            call shorbz(dpos,dpos,plat,qlat)
C           Only print out if displacement nonzero
            if (ddot(3,dpos,1,dpos,1) .ne. 0) then
              write(ifi,1331) i1,i2,(rtab(i,iat),i=1,3),
     .          (dpos(i),i=1,3),d,clabl,clabl2
              mark(i2) = lnn
            endif
          endif
  331     format(1x,2i4,3f12.7,2x,f12.7,2x,a,1x,a)
 1331     format(1x,2i4,3f12.7,2x,3f12.7,2x,f12.7,2x,a,1x,a)
        elseif (ltab .eq. 2) then
          write(ifi,332) (rtab(i,iat),i=1,3)
          mark(i2) = lnn
  332     format(1x,3f12.7)
        elseif (ltab .eq. 12) then
          do  i = 1, 3
            dpos(i) = pos2(i,i2)-pos(i,i2)
          enddo
C         Only print out if displacement nonzero
          call shorbz(dpos,dpos,plat,qlat)
          if (ddot(3,dpos,1,dpos,1) .gt. tol*tol) then
C         if (.not. latvec(1,tol,qlat,dpos)) then
            write(ifi,1332)
     .        (rtab(i,iat),i=1,3),(dpos(i),i=1,3)
            mark(i2) = lnn
          endif
 1332     format(1x,3f12.7,2x,3f12.7)
        elseif (ltab .eq. 3) then
          do  i = 1, 3
            write(ifi,333) i, i2, rtab(i,iat)
            mark(i2) = lnn
          enddo
  333     format(1x,i3,i5,3f12.7)
        endif

   10   continue
      enddo

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end subroutine psho1



