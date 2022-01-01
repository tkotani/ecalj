      module m_xxsho_unused
      private
      contains
      subroutine xxsho(npr,nbas,plat,bas,iax,ipc,dclabl,nlstc,lstc,lves,
     .ves,leul,eula,z)
C- Kernel called by shoshl
C  nlstc,lstc:  a list of classes to include as pairs (nlstc>0)
C     implicit none
      logical lves,leul
      integer npr,nbas,niax,ipc(*),nlstc
      integer,optional::lstc(nlstc)
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
      end module m_xxsho_unused


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


