subroutine splcls(nosplt,bas,nbas,ng,istab,nspec,slabl,nclass,ipc, &
     ics,nrclas)
  use m_lgunit,only:stdo
  !- Splits species into classes
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   nosplt:   T copy class and species
  !i   bas,nbas: dimensionless basis vectors, and number
  !i   nspec:    number of species
  !i   ipc:      on input, site j belongs to species ipc(j)
  !i   slabl:    on input, slabl is species label
  !i   ng:       number of group operations
  !i   istab:    site ib is transformed into istab(ib,ig) by grp op ig
  !o  Outputs:
  !o   slabl:    class labels
  !o   ipc:      on output, site j belongs to class ipc(j)
  !o   ics:      class i belongs to species ics(i)
  !o   nclass:   number of classes
  !o   nrclas:   number of classes per each species
  !u Updates
  !u   04 Apr 03 Search for equivalent classes more thorough
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters:
  logical :: nosplt
  integer :: nbas,nspec,nclass,ng,istab(nbas,ng),ipc(nbas), &
       !     .ics(nspec),nrclas(nspec)
       ics(*),nrclas(nspec)
  double precision :: bas(3,*)
  character(8) :: slabl(*)
  ! Local parameters:
  integer :: ib,ic,icn,iclbsj,ig,jb,m,i,is,ipr,idx,ispec,j
  logical :: lyetno
  character(80) :: outs,clabl(8)

  !      stdo = lgunit(1)
  ! ccccccccccccccccccccccccc
  ! akao
  call getpr(ipr)
  !      ipr=40
  !     print * ,'nosplt=',nosplt
  ! ccccccccccccccccccccccccc
  call icopy(nspec,1,0,nrclas,1)
  nclass = nspec
  do  5  i = 1, nspec
     ics(i) = i
5 enddo

  if (nosplt) then
  else
     ! --- For each species, split to make equivalent classes ---
     ic = 1
10   if (ic <= nclass) then
        is = ics(ic)
        ib = iclbsj(ic,ipc,-nbas,1)
        !   ... No sites of this class ... skip
        if (ib < 0) goto 11
        lyetno = .true.
        !   ... For each basis atom in this class, do
        do  20  jb = 1, nbas
           !         print *, 'jb=',ib,jb,nclass
           if (ipc(jb) == ic) then
              !      ... If there is a g mapping ib->jb, sites are equivalent
              do  22  ig = 1, ng
                 if (istab(ib,ig) == jb) goto 20
22            enddo
              !      ... If there is a g mapping jb->ib, sites are equivalent
              do  23  ig = 1, ng
                 if (istab(jb,ig) == ib) goto 20
23            enddo
              ! 04 Apr 03 consider these other possibilities
              !      ... If there is a g mapping ib->kb,jb, sites are equivalent
              do  24  ig = 1, ng
                 if (istab(istab(ib,ig),ig) == jb) goto 20
24            enddo
              !      ... If there is a g mapping jb->kb,ib, sites are equivalent
              do  25  ig = 1, ng
                 if (istab(istab(jb,ig),ig) == ib) goto 20
25            enddo
              !      ... There wasn't one
              if (ipr >= 70) then
                 write(stdo,400) slabl(is),ib,(bas(m,ib),m = 1,3), &
                      jb,(bas(m,jb),m = 1,3)
              endif
              !      ... If the classes haven't been split yet, do so
              if (lyetno) then
                 nclass = nclass+1
                 icn  =  nclass
                 ics(icn) = is
                 nrclas(is) = nrclas(is)+1
                 lyetno = .false.
              endif
              if (nclass > nbas) then
                 call rx('splcls:  problem with istab')
              endif
              icn  =  nclass
              ipc(jb)=  icn
           endif
20      enddo
11      ic = ic + 1
        goto 10
     endif
  endif
  if(ipr>=30) then
     !      write(stdo,"(a,i0,x,i0)")'SPLCLS: ',nspec,nclass
     write(stdo,"(a)")        ' SPLCLS: ibas iclass ispec label(ispec)'
     do j=1,nbas
        ic   =ipc(j)
        ispec=ics(ic)
        write(stdo,"(a,3i5,a)")" SPLCLS ",j,ic,ispec,'     '//trim(slabl(ispec))
     enddo
  endif
  ! ccccccccccccccccccccccc
  ! akao
  !      print *,' ipr=',ipr,nclass,nspec
  !      print *,' ics=',ics(1:nspec)
  ! ccccccccccccccccccccccccc


  ! --- Printout ---
  if (nclass == nspec .OR. ipr < 20) return
  !      call awrit2('%N SPLCLS:  %i species split into %i classes'//
  !     .'%N Species  Class      Sites...',' ',80,stdo,nspec,nclass)
  !      write(stdo,"(a,i0,x,i0)")'SPLCLS: species ==> classes',nspec,nclass
  !      write(stdo,"(a)")' Species  Class      Sites...'
  !      if (ipr .le. 30) return
  !$$$      do  40  is = 1, nspec
  !$$$        if (nrclas(is) .eq. 1 .and. ipr .lt. 40) goto 40
  !$$$        outs = ' '//slabl(is)
  !$$$        do  42  idx = 1, nbas
  !$$$          ic = iclbsj(is,ics,-nclass,idx)
  !$$$          if (ic .gt. 0) then
  !$$$            call clabel(slabl,is,idx,clabl)
  !$$$            write(stdo,*) trim(slabl(is))//' '//trim(clabl),ic
  !$$$c$$$            outs = ' '
  !$$$c$$$            if (idx .eq. 1) outs = ' '//slabl(is)
  !$$$c$$$            outs= '         '//trim(i2char(ic))//' '//clabl
  !$$$c$$$!     call awrit1('%(n>9?9:10)p%i:'//clabl,outs,80,0,ic)
  !$$$c$$$            do  44  ib = 1, nbas
  !$$$c$$$              if (ipc(ib) .eq. ic)
  !$$$c$$$              outs=trim(outs)//trim(i2char(ib)
  !$$$c$$$c     .        call awrit1('%a%(p>20?p:20)p %i',outs,80,0,ib)
  !$$$c$$$   44       continue
  !$$$c$$$            write(stdo,*) trim(outs) !call awrit0(outs,' ',-80,stdo)
  !$$$          else
  !$$$            goto 43
  !$$$          endif
  !$$$   42   continue
  !$$$   43   continue
  !$$$   40 continue

400 format(' SPLCLS: species: ',a,'has inequivalent positions:'/ &
       '  IB: ',i3,',  POS=',3f10.5/ &
       '  JB: ',i3,',  POS=',3f10.5)

end subroutine splcls
