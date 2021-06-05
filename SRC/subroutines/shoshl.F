
#define USEOPTIONAL

      module m_xxsho
      contains
      subroutine xxsho(npr,nbas,plat,bas,iax,ipc,dclabl,nlstc,lstc,lves,
     .ves,leul,eula,z)

C- Kernel called by shoshl
C  nlstc,lstc:  a list of classes to include as pairs (nlstc>0)
C     implicit none
      logical lves,leul
#ifdef USEOPTIONAL
      integer npr,nbas,niax,ipc(*),nlstc
      integer,optional::lstc(nlstc)
#else
      integer npr,nbas,niax,ipc(*),nlstc,lstc(nlstc)
#endif
      parameter (niax=10)
      integer iax(niax,*)
      double precision plat(3,3),bas(3,*),ves(*),eula(nbas,3),z(32)
      character(8):: dclabl(*)
      integer ih(2,120),scrwid
      parameter (scrwid=120)
      integer i,l,ishell,nshell,j,k,ii,kk,ic,jc,i1,lgunit,awrite,iclbsj,
     .ib
      double precision dr(3),d,drr2,dshell,fuzz,z1(3),z2(3),alfa,beta,
     .angle,pi,ddot
      character*8 clabl,outs1*25,outs2*(scrwid),outsv*(scrwid),
     .outse*(scrwid)

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
c      call r8tos8(dclabl(ic),clabl)
      clabl=dclabl(ic)
      print 302, clabl, ic, nint(z(ic))
  302 format(/' Shell decomposition for class ',a,
     .'  class',i4,'  z=',i2/
     .' shell   d     nsh csiz  class ...')

      do  10  i = 1, npr
        d = dsqrt(drr2(plat,bas(1,iax(1,i)),bas(1,iax(2,i)),
     .  iax(3,i),iax(4,i),iax(5,i),dr))
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
c            call r8tos8(dclabl(ih(2,k)),clabl)
            clabl=dclabl(ih(2,k))
            l = awrite('%a  %np%i:'//clabl//
     .      '%a%?;n>1;(%i);%j;',outs2,len(outs2),0,
     .      (kk-1)*14,ih(2,k),ih(1,k),ih(1,k),0,0,0,0)
            if (lves) call awrit2('%np%d',outsv,len(outsv),0,
     .      (kk-1)*14,ves(ih(2,k)))
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
      end subroutine xxsho
      end module m_xxsho


c$$$      subroutine shoshl(sopts,nbas,bas,plat,mxnbr0,z,slabl,dclabl,ips,
c$$$     .ipc,ves,eula,nclass)
c$$$      use m_pairs
c$$$      use m_xxsho
c$$$C- Print nearest-neighbor shells
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   Everything is input
c$$$Ci   sopts: a set of modifiers, with the syntax
c$$$Ci          [:v][:e][:r=#][:sites:site-list][:pairs:pair-list] ..
c$$$Ci          [:tab[=#]][:disp=fnam][:nn][:fn=fnam]
c$$$Ci          :
c$$$Ci          :r=# sets range for shells
c$$$Ci          :v prints out electrostatic potential
c$$$Ci          :e prints out Euler angles
c$$$Ci          :r=# restricts neighbor table to range #
c$$$Ci          :sites:list print table only for sites in site-list
c$$$Ci          :pairs:pair-list print table only for pairs in pair-list
c$$$Ci          :tab prints out neighbor table
c$$$Ci          :tab=# prints out neighbor table in style #
c$$$Ci          :  ...The following only apply to :tab
c$$$Ci          :  fn=fnam write table to file fnam
c$$$Ci          :  disp=fnam : read another site positions file; neighbor
c$$$Ci          :              table for both original and displaced
c$$$Ci          :              positions is written.
c$$$Ci          :  nn only print first entry for a given pair site
c$$$Ci          :     in neighbor table
c$$$Ci  OLD: doesn't work
c$$$Ci          :i[=style-#]:list  restricts neighbors in shell to list.
c$$$Ci                             This must be the last modifier.
c$$$Cu Updates
c$$$Cu   19 Apr 03 Can write displaced neigbhor table
c$$$Cu   12 Apr 03 Can write neigbhor table
c$$$Cr   24 Nov 97 changed modifier list
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$C Passed parameters
c$$$      integer nbas,nclass,mxnbr0,ips(nbas),ipc(nbas),nttab
c$$$      double precision bas(3,nbas),plat(3,3),ves(1),
c$$$     .eula(1),z(32)
c$$$      character(8):: dclabl(*)
c$$$ckino z(32?)     .eula(1),z(1)
c$$$      character sopts*(*),slabl(nbas)*8
c$$$C Local parameters
c$$$      integer niax
c$$$      parameter (niax=10)
c$$$      logical lves,leula,a2bin
c$$$      double precision avwsr,avw,range,xx,stdo,wdummy
c$$$      integer:: olstc
c$$$ckino Dec.14.2011:         integer,pointer :: iv_p_ontab(:) =>NULL()
c$$$      integer,allocatable :: iv_a_ontab(:)
c$$$
c$$$      real(8) ,allocatable :: wk_rv(:)
c$$$      real(8) ,allocatable :: rtab_rv(:)
c$$$      integer ,allocatable :: lst1_iv(:)
c$$$      integer ,allocatable :: lst2_iv(:), w_oiax(:)
c$$$      real(8) ,allocatable :: dist_rv(:)
c$$$      real(8) ,allocatable :: mark_rv(:)
c$$$
c$$$      integer:: npr(2) , mxnbr , ib , ic , iclbsj , j , j1 , j2 , lstyle 
c$$$     . , ival , scrwid , ltab , nlstc , nlst1 , nlst2 , mxcsiz , iprint 
c$$$     . , ifi , fopn , lgunit , ldisp , lnn
c$$$ckino Dec.14.2011:         integer,pointer :: iv_p_oiax(:) =>NULL()
c$$$      integer,allocatable :: iv_a_oiax(:)
c$$$
c$$$      character*8  dc*1
c$$$      parameter (scrwid=120)
c$$$C heap:
c$$$      integer ::iwdummy
c$$$
c$$$
c$$$C --- Parse modifiers ---
c$$$      stdo = lgunit(1)
c$$$      ifi = stdo
c$$$      lves  = .false.
c$$$      leula = .false.
c$$$      lstyle = 1
c$$$      range = 2
c$$$      nlstc = 0
c$$$c      olstc = 1
c$$$      ltab = 0
c$$$      nlst1 = 0
c$$$      nlst2 = 0
c$$$      allocate(lst1_iv(nbas))
c$$$
c$$$      allocate(lst2_iv(nbas))
c$$$
c$$$      lnn = 0
c$$$      ldisp = 0
c$$$      allocate(dist_rv(nbas*3))
c$$$
c$$$c      allocate(mark_rv(nbas))
c$$$
c$$$      call dcopy ( nbas * 3 , bas , 1 , dist_rv , 1 )
c$$$
c$$$
c$$$      if (sopts .ne. ' ') then
c$$$C       ls = len(sopts)
c$$$        j1 = 1
c$$$        dc = sopts(j1:j1)
c$$$        j1 = j1+1
c$$$
c$$$C   ... Return here to resume parsing for arguments
c$$$   40   continue
c$$$        call nwordg(sopts,0,dc//' ',1,j1,j2)
c$$$
c$$$C   ... Parse special arguments
c$$$        if (j2 .ge. j1) then
c$$$C         print *, sopts(j1:j2)
c$$$          if (sopts(j1:j2) .eq. 'v')  then
c$$$            lves = .true.
c$$$
c$$$          elseif (sopts(j1:j2) .eq. 'e')  then
c$$$            leula = .true.
c$$$
c$$$          elseif (sopts(j1:j1+1) .eq. 'r=') then
c$$$            j = 0
c$$$            if (.not. a2bin(sopts(j1+2:),range,4,0,' ',j,j2-j1-2))
c$$$     .      goto 999
c$$$
c$$$          elseif (sopts(j1:j1+2) .eq. 'fn=')  then
c$$$            if (j1+3 .le. j2) ifi = fopn(sopts(j1+3:j2))
c$$$
c$$$          elseif (sopts(j1:j1+3) .eq. 'tab=')  then
c$$$            j = 0
c$$$            if (.not. a2bin(sopts(j1+4:),ltab,2,0,' ',j,j2-j1-4))
c$$$     .      goto 999
c$$$
c$$$          elseif (sopts(j1:j1+4) .eq. 'disp=')  then
c$$$            call iopos ( .false. , - 1 , sopts ( j1 + 5:j2 ) , nbas , dist_rv 
c$$$     .      )
c$$$
c$$$            ldisp = 1
c$$$
c$$$          elseif (sopts(j1:j1+1) .eq. 'nn')  then
c$$$            lnn = 1
c$$$
c$$$          elseif (sopts(j1:j1+2) .eq. 'tab')  then
c$$$            ltab = 1
c$$$
c$$$          elseif (sopts(j1:j1+4) .eq. 'sites')  then
c$$$            call baslst ( 0 , 11 , sopts ( j1 + 5: ) , j , ips , nbas , slabl 
c$$$     .       , z , 0 , ' ' , iwdummy , nlst1 , lst1_iv )
c$$$
c$$$
c$$$            j2 = j1+5+j-2
c$$$
c$$$          elseif (sopts(j1:j1+4) .eq. 'pairs')  then
c$$$            call baslst ( 0 , 11 , sopts ( j1 + 5: ) , j , ips , nbas , slabl 
c$$$     .       , z , 0 , ' ' , iwdummy , nlst2 , lst2_iv )
c$$$
c$$$
c$$$            j2 = j1+5+j-2
c$$$
c$$$          else
c$$$            call rxs('shoshl: failed to parse --shell switch: ',
c$$$     .      sopts(j1:j2))
c$$$            goto 999
c$$$          endif
c$$$          j1 = j2+2
c$$$          goto 40
c$$$        endif
c$$$      endif
c$$$
c$$$C --- Print neighbor table for each site ---
c$$$      if (ltab .ne. 0) then
c$$$        mxcsiz = 0
c$$$        call pshpr(iprint()-20)
c$$$        call pairs ( nbas , nbas , 1d0 , plat , ( / range / 2 * ( 1 + 
c$$$     .   1d - 6 ) / ) , bas , ( / - 1 / ) , 3 , - 1 , (/iwdummy/) , nttab
c$$$ckino Dec.14.2011:       .   , iv_p_ontab , iv_p_oiax , mxcsiz )
c$$$ckino Dec.14.2011:       .   , iv_p_ontab , iv_a_oiax , mxcsiz ) 
c$$$     .   , iv_a_ontab , iv_a_oiax , mxcsiz )
c$$$
c$$$        call poppr
c$$$        allocate(rtab_rv(3*nttab))
c$$$
c$$$ckino Dec.14.2011:           call mkrtab ( 000 , 1d0 , plat , bas , iv_p_oiax , nttab , bas
c$$$        call mkrtab ( 000 , 1d0 , plat , bas , iv_a_oiax , nttab , bas 
c$$$     .   , rtab_rv )
c$$$
c$$$
c$$$C       if (ldisp .eq. 1) ltab = 2
c$$$ckino Dec.14.2011:           call psho1 ( ldisp * 10 + ltab , lnn , iv_p_oiax , nbas , nttab
c$$$        call psho1 ( ldisp * 10 + ltab , lnn , iv_a_oiax , nbas , nttab 
c$$$     .   , nlst1 , lst1_iv , nlst2 , lst2_iv , rtab_rv , plat , bas , 
c$$$     .   dist_rv , ipc , dclabl , ifi )
c$$$
c$$$
c$$$
c$$$Cki#error ERROR, try to release name= ontab ,but list does not have  ontab at linenumber= 151 list= olst1 olst2 odist omark ortab
c$$$
c$$$
c$$$Cki#error, have return with len(w_varlist)>0 at line 152
c$$$        if (allocated(rtab_rv)) deallocate(rtab_rv)
c$$$c        if (allocated(mark_rv)) deallocate(mark_rv)
c$$$        if (allocated(dist_rv)) deallocate(dist_rv)
c$$$        if (allocated(lst2_iv)) deallocate(lst2_iv)
c$$$        if (allocated(lst1_iv)) deallocate(lst1_iv)
c$$$        return
c$$$
c$$$      endif
c$$$
c$$$C --- Show shells for each class ---
c$$$      if (mxnbr0 .eq. 0) then
c$$$        mxnbr = 100*range**3
c$$$      else
c$$$        mxnbr = mxnbr0
c$$$      endif
c$$$c      call defi(oiax,niax*mxnbr)
c$$$      allocate(wk_rv(mxnbr),w_oiax(niax*mxnbr))
c$$$
c$$$      avw = avwsr(plat,1d0,xx,nbas)
c$$$      avw = 1
c$$$c     call pshprt(50)
c$$$      do  20  ic = 1, nclass
c$$$        ib = iclbsj(ic,ipc,nbas,1)
c$$$        if (nlst1 .gt. 0) then
c$$$          call hunti ( lst1_iv , nlst1 , ib , 0 , j )
c$$$
c$$$          if ( ival ( lst1_iv , j + 1 ) .ne. ib ) goto 20
c$$$
c$$$        endif
c$$$        call nghbor ( nbas , plat , bas , range , range , ib , mxnbr 
c$$$     .  , npr , w_oiax , wk_rv )
c$$$
c$$$#ifdef USEOPTIONAL
c$$$        call xxsho(npr=npr(1),nbas=nbas,plat=plat,bas=bas,iax=w_oiax,ipc=ipc,
c$$$     .  dclabl=dclabl,nlstc=nlstc,
c$$$     .  lves=lves,ves=ves,leul=leula,eula=eula,z=z)
c$$$#else
c$$$        call xxsho(npr(1),nbas,plat,bas,w_oiax,ipc,dclabl,nlstc,
c$$$     .  wdummy,lves,ves,leula,eula,z)
c$$$#endif
c$$$   20 continue
c$$$      deallocate(wk_rv,w_oiax) !call rlse(oiax)
c$$$
c$$$Cki#error, have return with len(w_varlist)>0 at line 178
c$$$      if (allocated(lst1_iv)) deallocate(lst1_iv)
c$$$      if (allocated(lst2_iv)) deallocate(lst2_iv)
c$$$      if (allocated(dist_rv)) deallocate(dist_rv)
c$$$      if (allocated(mark_rv)) deallocate(mark_rv)
c$$$      if (allocated(rtab_rv)) deallocate(rtab_rv)
c$$$      return
c$$$
c$$$
c$$$  999 call rxs('shoshl: failed to parse ',sopts)
c$$$      if (allocated(rtab_rv)) deallocate(rtab_rv)
c$$$      if (allocated(mark_rv)) deallocate(mark_rv)
c$$$      if (allocated(dist_rv)) deallocate(dist_rv)
c$$$      if (allocated(lst2_iv)) deallocate(lst2_iv)
c$$$      if (allocated(lst1_iv)) deallocate(lst1_iv)
c$$$
c$$$      end subroutine shoshl


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
      integer:: npr(2) , mxnbr , ib , nshell , nmxshl , niax , j , 
     .j1 , j2 , nspec , mxint , nsites ,  nbonds ,  parg  !oslist ,ilst ,
     ., m , iv(10) , ilast
      integer ,allocatable :: iax_iv(:)
      real(8) ,allocatable :: wk_rv(:)
      integer ,allocatable :: num_iv(:)
      real(8) ,allocatable :: ang_rv(:)
      real(8) ,allocatable :: d_rv(:)
      integer ,allocatable :: blist_iv(:),slist(:)

      parameter (niax=10)
      character dc*1
C ... Heap
c      integer w(1)
c      common /w/ w
C ... External calls
ckino defdr,defi,rlse
      external nghbor,nwordg,pvang1,pvang2,rxs

C ... Setup
      range = 2.5d0
      nspec = mxint(nbas,ips)
      nsites = 0
      nbonds = 0
c      call defi(oslist,nbas)
      allocate(blist_iv(nbas),slist(nbas))


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
     .      dc//' ',1,1,iv,range)
            if (j .le. 0) goto 999

C         Site list
          elseif (sopts(j1:j1+4) .eq. 'sites') then
            if (sopts(j1+5:j1+5) .eq. '=') sopts(j1+5:j1+5) = dc
            call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zz,0,' ',
     .      xx,nsites,slist)

C         Bond list
          elseif (sopts(j1:j1+4) .eq. 'bonds') then
            if (sopts(j1+5:j1+5) .eq. '=') sopts(j1+5:j1+5) = dc
            call baslst ( 0 , 10 , sopts ( j1 + 5: ) , ilast , ips , nbas 
     .      , slabl , zz , 0 , ' ' , xx , nbonds , blist_iv )

          endif
          goto 10

        endif
      endif

      if (mxnbr0 .eq. 0) then
        mxnbr = 2*range**3
      else
        mxnbr = mxnbr0
      endif
      allocate(iax_iv(niax*mxnbr))

      allocate(wk_rv(mxnbr))

      avw = avwsr(plat,1d0,xx,nbas)

C --- For every site in list, generate tables of bond angles ---
c      ilst = 0
      do  ib = 1, nbas
        if (nsites .gt. 0) then
c          if (w(oslist+ilst) .ne. ib) goto 31
          if (slist(ib) .ne. ib) goto 31
        endif
c        ilst = ilst+1

C   ... Get neighbor lists

        call nghbor ( nbas , plat , bas , range * avw , range * avw , 
     .  ib , mxnbr , npr , iax_iv , wk_rv )


C   ... Get shell dimensions
        allocate(num_iv(abs(-nspec*npr(1))))
        if (-nspec*npr(1)<0) num_iv(:)=0

        call pvang1 ( npr ( 1 ) , nbas , plat , bas , iax_iv , ips 
     .  , num_iv , nshell , nmxshl )

        if (allocated(num_iv)) deallocate(num_iv)

        allocate(num_iv(abs(-nspec**2*nshell**2)))
        if (-nspec**2*nshell**2<0) num_iv(:)=0

        allocate(ang_rv(nspec**2*nshell**2*nmxshl**2))

        allocate(d_rv(nshell))


C   ... Print bond angles
        call pvang2 ( npr ( 1 ) , nbas , nspec , nshell , nmxshl , plat 
     .  , bas , iax_iv , ips , slabl , nbonds , blist_iv , num_iv 
     .  , ang_rv , d_rv )

        if (allocated(d_rv)) deallocate(d_rv)
        if (allocated(ang_rv)) deallocate(ang_rv)
        if (allocated(num_iv)) deallocate(num_iv)

   31   continue
      enddo
      if (allocated(wk_rv)) deallocate(wk_rv)
      if (allocated(iax_iv)) deallocate(iax_iv)
      deallocate(blist_iv)
c      call rlse(oslist)
      deallocate(slist)
      return

  999 call rxs('shoang: failed to parse ',sopts)
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
      integer iax(niax,*)
      double precision plat(3,3),bas(3,1)
C ... Local parameters
      integer i,is
      double precision d,wk(0:3),tol

      tol = 1d-6
C --- Get no. of shells and max no. of atoms in 1 shell and 1 class ---
      nshell = 0
      nmxshl = 0
      d = 0d0
      do  i = 2, npr
        is = ips(iax(2,i))
        call dlmn(nbas,plat,bas,iax(1,i),wk)
C       distance changed by more than tol ... new shell
        if (dabs(wk(0)-d) .gt. tol) then
          nshell = nshell + 1
          d = wk(0)
        endif
        num(nshell,is) = num(nshell,is) + 1
        nmxshl = max0(nmxshl,num(nshell,is))
      enddo
      end subroutine pvang1


      subroutine pvang2(npr,nbas,nspec,nshell,nmxshl,plat,bas,iax,
     .ips,slabl,nbonds,blist,num,ang,d)

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
     .num(nspec,nspec,nshell,nshell),nbonds,blist(nbonds)
      parameter (niax=10)
      integer iax(niax,*)
      character*8 slabl(*)
C ... Local parameters
      double precision plat(3,3),bas(3,1),d(nshell),
     .ang(nmxshl**2,nspec,nspec,nshell,nshell)
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
     .' shl1    d1    shl2    d2     cl1      cl2       angle(s) ...')

      do  60  nsh1 = 1, nshell
        do  50  nsh2 = nsh1, nshell
          do  40  is = 1, nspec
            do  30  js = 1, nspec
              n = num(is,js,nsh1,nsh2)
              if (n .ne. 0) then
                print 401, nsh1, d(nsh1), nsh2, d(nsh2), slabl(is),
     .          slabl(js), (ang(i,is,js,nsh1,nsh2), i = 1, n)
              endif
   30       continue
   40     continue
   50   continue
   60 continue

  401 format(2(1x,i3,1x,f9.6),1x,2(1x,a8),20(4(1x,f7.2):/47x))

      end subroutine pvang2

      subroutine psho1(ltab,lnn,iax,nbas,nttab,nlst1,lst1,nlst2,lst2,
     .rtab,plat,pos,pos2,ipc,dclabl,ifi)

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
      double precision pos(3,*),pos2(3,*),rtab(3,*),plat(9)
      character(8):: dclabl(*)
      parameter (niax=10)
      integer iax(niax,*)
C ... Local parameters
      integer iat,low1,low2,i1,i2,ipr,stdo,lgunit,i,oldi1,oldiat,
     .ic,jc
      character clabl*8, clabl2*8
      double precision d,ddot,dsqrt,dpos(3),qlat(9),tol
C     logical latvec
      logical::isanrg, l_dummy_isanrg

      tol = 1d-6

      if (ltab .eq. 0) return
      call mkqlat(plat,qlat,d)

      call info5(20,1,1,' ... shoshl: making neighbor list for'//
     .' %?#n#%-1j%i%j#%i# atom(s), style %i',nlst1,nbas,ltab,0,0)

      if (ltab .eq. 3)
     .call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)

      call getpr(ipr)
      ipr = 100
      stdo = lgunit(1)
Ckino isanrg is logical function,       call isanrg(mod(ltab,10),1,3,' shoshl:','tab',.true.)
      l_dummy_isanrg=isanrg(mod(ltab,10),1,3,' shoshl:','tab',.true.)

      if (ltab .eq. 3)
     .call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)
      if (ltab .eq. 13)
     .call awrit1('%% rows 6 cols %i sparse',' ',80,ifi,nbas)

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
c          call r8tos8(dclabl(ic),clabl)
          clabl=dclabl(ic)
          call awrit1('# neighbor list for site %i, class '//clabl,
     .    ' ',80,ifi,i1)
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
c          call r8tos8(dclabl(jc),clabl2)
          clabl2=dclabl(jc)
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
     .        (dpos(i),i=1,3),d,clabl,clabl2
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
     .      (rtab(i,iat),i=1,3),(dpos(i),i=1,3)
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

      end subroutine psho1



