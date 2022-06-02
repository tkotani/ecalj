      subroutine splcls(nosplt,bas,nbas,ng,istab,nspec,slabl,nclass,ipc,
     .ics,nrclas)
      use m_lgunit,only:stdo
C- Splits species into classes
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nosplt:   T copy class and species
Ci   bas,nbas: dimensionless basis vectors, and number
Ci   nspec:    number of species
Ci   ipc:      on input, site j belongs to species ipc(j)
Ci   slabl:    on input, slabl is species label
Ci   ng:       number of group operations
Ci   istab:    site ib is transformed into istab(ib,ig) by grp op ig
Co  Outputs:
Co   slabl:    class labels
Co   ipc:      on output, site j belongs to class ipc(j)
Co   ics:      class i belongs to species ics(i)
Co   nclass:   number of classes
Co   nrclas:   number of classes per each species
Cu Updates
Cu   04 Apr 03 Search for equivalent classes more thorough
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters:
      logical nosplt
      integer nbas,nspec,nclass,ng,istab(nbas,ng),ipc(nbas),
c     .ics(nspec),nrclas(nspec)
     .ics(*),nrclas(nspec)
      double precision bas(3,*)
      character*8 slabl(*)
C Local parameters:
      integer ib,ic,icn,iclbsj,ig,jb,m,i,is,ipr,idx,ispec,j
      logical lyetno
      character*80 outs,clabl*8

c      stdo = lgunit(1)
ccccccccccccccccccccccccccc
ctakao
      call getpr(ipr)
c      ipr=40
c     print * ,'nosplt=',nosplt
ccccccccccccccccccccccccccc
      call icopy(nspec,1,0,nrclas,1)
      nclass = nspec
      do  5  i = 1, nspec
        ics(i) = i
    5 continue

      if (nosplt) then
      else
C --- For each species, split to make equivalent classes ---
        ic = 1
   10   if (ic .le. nclass) then
          is = ics(ic)
          ib = iclbsj(ic,ipc,-nbas,1)
C   ... No sites of this class ... skip
          if (ib .lt. 0) goto 11
          lyetno = .true.
C   ... For each basis atom in this class, do
          do  20  jb = 1, nbas
C         print *, 'jb=',ib,jb,nclass
            if (ipc(jb) .eq. ic) then
C      ... If there is a g mapping ib->jb, sites are equivalent
              do  22  ig = 1, ng
                if (istab(ib,ig) .eq. jb) goto 20
   22         continue
C      ... If there is a g mapping jb->ib, sites are equivalent
              do  23  ig = 1, ng
                if (istab(jb,ig) .eq. ib) goto 20
   23         continue
C 04 Apr 03 consider these other possibilities
C      ... If there is a g mapping ib->kb,jb, sites are equivalent
              do  24  ig = 1, ng
                if (istab(istab(ib,ig),ig) .eq. jb) goto 20
   24         continue
C      ... If there is a g mapping jb->kb,ib, sites are equivalent
              do  25  ig = 1, ng
                if (istab(istab(jb,ig),ig) .eq. ib) goto 20
   25         continue
C      ... There wasn't one
              if (ipr .ge. 70) then
                write(stdo,400) slabl(is),ib,(bas(m,ib),m = 1,3),
     .          jb,(bas(m,jb),m = 1,3)
              endif
C      ... If the classes haven't been split yet, do so
              if (lyetno) then
                nclass = nclass+1
                icn  =  nclass
                ics(icn) = is
                nrclas(is) = nrclas(is)+1
                lyetno = .false.
              endif
              if (nclass .gt. nbas) then
                call rx('splcls:  problem with istab')
              endif
              icn  =  nclass
              ipc(jb)=  icn
            endif
   20     continue
   11     ic = ic + 1
          goto 10
        endif
      endif
      if(ipr>=30) then
c      write(stdo,"(a,i0,x,i0)")'SPLCLS: ',nspec,nclass
      write(stdo,"(a)")        ' SPLCLS: ibas iclass ispec label(ispec)'
      do j=1,nbas
         ic   =ipc(j)
         ispec=ics(ic)
         write(stdo,"(a,3i5,a)")" SPLCLS ",j,ic,ispec,'     '//trim(slabl(ispec))
      enddo
      endif
ccccccccccccccccccccccccc
ctakao
c      print *,' ipr=',ipr,nclass,nspec
c      print *,' ics=',ics(1:nspec)
ccccccccccccccccccccccccccc


C --- Printout ---
      if (nclass .eq. nspec .or. ipr .lt. 20) return
c      call awrit2('%N SPLCLS:  %i species split into %i classes'//
c     .'%N Species  Class      Sites...',' ',80,stdo,nspec,nclass)
c      write(stdo,"(a,i0,x,i0)")'SPLCLS: species ==> classes',nspec,nclass
c      write(stdo,"(a)")' Species  Class      Sites...'
c      if (ipr .le. 30) return
c$$$      do  40  is = 1, nspec
c$$$        if (nrclas(is) .eq. 1 .and. ipr .lt. 40) goto 40
c$$$        outs = ' '//slabl(is)
c$$$        do  42  idx = 1, nbas
c$$$          ic = iclbsj(is,ics,-nclass,idx)
c$$$          if (ic .gt. 0) then
c$$$            call clabel(slabl,is,idx,clabl)
c$$$            write(stdo,*) trim(slabl(is))//' '//trim(clabl),ic
c$$$c$$$            outs = ' '
c$$$c$$$            if (idx .eq. 1) outs = ' '//slabl(is)
c$$$c$$$            outs= '         '//trim(i2char(ic))//' '//clabl
c$$$c$$$!     call awrit1('%(n>9?9:10)p%i:'//clabl,outs,80,0,ic)
c$$$c$$$            do  44  ib = 1, nbas
c$$$c$$$              if (ipc(ib) .eq. ic)
c$$$c$$$              outs=trim(outs)//trim(i2char(ib) 
c$$$c$$$c     .        call awrit1('%a%(p>20?p:20)p %i',outs,80,0,ib)
c$$$c$$$   44       continue
c$$$c$$$            write(stdo,*) trim(outs) !call awrit0(outs,' ',-80,stdo)
c$$$          else
c$$$            goto 43
c$$$          endif
c$$$   42   continue
c$$$   43   continue
c$$$   40 continue

  400 format(' SPLCLS: species: ',a,'has inequivalent positions:'/
     .'  IB: ',i3,',  POS=',3f10.5/
     .'  JB: ',i3,',  POS=',3f10.5)

      end

