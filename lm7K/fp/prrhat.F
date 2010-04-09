Cgetarg...info...           structure ['prrhat', 'site', 'ssite']
Cgetarg...info...           structure ['prrhat', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['prrhat', 'm_struc_def']
      subroutine prrhat(nbas,ssite,sspec,orhoat)
       
       use m_struc_def  !Cgetarg

C- Print info about the atomic charge densities
C     implicit none
C ... Passed parameters
      integer nbas,orhoat(3,nbas)
Cgetarg       double precision ssite(1),sspec(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

C ... Local parameters
      integer stdo,lgunit,ib,lmxl,nr,orofi,orwgt,nlml,is,igetss,nsp,
     .  nglob
      double precision a,rmt
C ... Heap
      integer w(1)
      common /w/ w

      nsp  = nglob('nsp')
      stdo = lgunit(1)
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec lmxl a nr rmt',sspec,is,lmxl,a,nr,rmt)
         
         lmxl=sspec(is)%lmxl
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt

        if (lmxl .eq. -1) goto 10
        call defrr(orofi, nr)
        call defrr(orwgt, nr)
        call radmsh(rmt,a,nr,w(orofi))
        call radwgt(rmt,a,nr,w(orwgt))

        nlml = (lmxl+1)**2
        write(stdo,200) ib,rmt,nr,nlml
  200   format(/' Density at site',i3,'   rmt=',f8.4,
     .     '   nr=',i5,'   nlml=',i3)
        call prlrho('true density',nr,nlml,nsp,w(orofi),w(orwgt),
     .     w(orhoat(1,ib)))

        call prlrho('smooth density',nr,nlml,nsp,w(orofi),w(orwgt),
     .     w(orhoat(2,ib)))

        call rlse(orofi)

   10   continue
      enddo

Cgetarg       end
       end subroutine prrhat 


      subroutine prlrho(str,nr,nlm,nsp,rofi,rwgt,rho)

C- Print info about site density
C     implicit none
C ... Passed parameters
      integer nr,nlm,nsp
      double precision rho(nr,nlm,nsp),rofi(nr),rwgt(nr)
      character*(*) str
C ... Local parameters
      integer stdo,lgunit,ilm,l,ll,itop,ibot,i
      double precision xx,pi,srfpi,top,bot,sum

      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)

      write(stdo,601) str
      do  ilm = 1, nlm
        l = ll(ilm)
        top = -1d10
        bot = 1d10
        sum = 0d0
        itop = 0
        ibot = 0
        do  i = 1, nr
          xx = (rho(i,ilm,1)+rho(i,ilm,nsp))/(3-nsp)
          if (xx .gt. top) then
            itop = i
            top = xx
          endif
          if (xx .lt. bot) then
            ibot = i
            bot = xx
          endif
          sum = sum + rwgt(i)*xx*(rofi(i)+1d-32)**l
        enddo
        xx = dmax1(dabs(top),dabs(bot))
        if (xx .gt. 1d-6) then
          if (ilm .eq. 1) then
            write(stdo,600) ilm,bot,ibot,top,itop,sum,srfpi*sum
          else
            write(stdo,600) ilm,bot,ibot,top,itop,sum
          endif
        endif
  600   format(i4,f15.8,i5,f15.8,i5,f15.8,f12.6)
  601   format(/' prlrho: ',a/'  ilm',6x,'rhomin    pnt',
     .     7x,'rhomax    pnt       moment')
      enddo
Cgetarg       end
       end subroutine prlrho 


