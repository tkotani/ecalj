      subroutine packs(snamein,sstrn,packingw)
C- Pack words=packingw into sstrn. Use upacks to take out the words below.
Ci snamein: sname defined in sstrn_sname.
Ci sstrn: long string. sstrings.
Ci packingw]
      use m_sstrnsname
      implicit none
      character*(*):: snamein,sstrn,packingw
      integer:: i,lll
      if(initialized==0) then
        call initialize()
        initialized=1
      endif
      do i=1,nsname
        lll=len_trim(packingw)
        if(lll>nx) then
          call rxs('packs: too large packing word.'//
     &    'Need to modify code',packingw)
        endif
        if(trim(sname(i))==trim(snamein)) then
          iend(i) = init(i) +lll -1
          sstrn(init(i):iend(i))= packingw
          return
        endif
      enddo
      call rxs('packs: no sname is found in sstrn_sname. '//
     &'bug or need to modify codes',snamein)
      end

      subroutine upacks(snamein,init_,iend_)
C- Return init and iend for sstring(init:iend) for given snamein. Corresponds to packs.
Ci snamein
Co init_
Co iend_
      use m_sstrnsname
      implicit none
      integer::init_,iend_,i
      character*(*):: snamein
      if(initialized==0) then
        call initialize()
        initialized=1
      endif
      do i=1,nsname
        if(trim(sname(i))==trim(snamein)) then
          init_=init(i)
          iend_=iend(i)
          return
        endif
      enddo
      call rxs('upacks: no sname is found for snamein. '//
     &'bug or need to modify codes',snamein)
      end
