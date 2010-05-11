Cdelw1 var ['orofi', '1', 'orofi', 'nr', 'real(8)', 'rv_w_', 'prrhat']
Cdelw1 var ['orwgt', '1', 'orwgt', 'nr', 'real(8)', 'rv_w_', 'prrhat']
Cdelw1 undel [['call'], 'orhoat', 'prrhat']
Cgetarg...info...           structure ['prrhat', 'site', 'ssite']
Cgetarg...info...           structure ['prrhat', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['prrhat', 'm_struc_def']
      subroutine prrhat(nbas,ssite,sspec,orhoat)
       
       use m_struc_def  !Cgetarg

C- Print info about the atomic charge densities
C     implicit none
C ... Passed parameters
Cdelw1 do not change  [['call'], 'orhoat', 'prrhat']
      integer nbas,orhoat(3,nbas)
Cgetarg       double precision ssite(1),sspec(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

C ... Local parameters
Cdelw1       integer stdo,lgunit,ib,lmxl,nr,orofi,orwgt,nlml,is,igetss,nsp,
Cdelw1      .  nglob
       integer:: stdo , lgunit , ib , lmxl , nr , nlml , is , igetss 
     . , nsp , nglob 
       real(8) ,allocatable :: rv_w_orofi(:)
       real(8) ,allocatable :: rv_w_orwgt(:)

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
Cdelw1         call defrr(orofi, nr)
         allocate(rv_w_orofi(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orofi(:)=0.0d0

Cdelw1         call defrr(orwgt, nr)
         allocate(rv_w_orwgt(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orwgt(:)=0.0d0

Cdelw1         call radmsh(rmt,a,nr,w(orofi))
         call radmsh ( rmt , a , nr , rv_w_orofi ) 

Cdelw1         call radwgt(rmt,a,nr,w(orwgt))
         call radwgt ( rmt , a , nr , rv_w_orwgt ) 


        nlml = (lmxl+1)**2
        write(stdo,200) ib,rmt,nr,nlml
  200   format(/' Density at site',i3,'   rmt=',f8.4,
     .     '   nr=',i5,'   nlml=',i3)
Cdelw1 1 is digit, dropped
Cdelw1         call prlrho('true density',nr,nlml,nsp,w(orofi),w(orwgt),
Cdelw1      .     w(orhoat(1,ib)))
         call prlrho ( 'true density' , nr , nlml , nsp , rv_w_orofi , 
     .   rv_w_orwgt , w ( orhoat ( 1 , ib ) ) ) 


Cdelw1 2 is digit, dropped
Cdelw1         call prlrho('smooth density',nr,nlml,nsp,w(orofi),w(orwgt),
Cdelw1      .     w(orhoat(2,ib)))
         call prlrho ( 'smooth density' , nr , nlml , nsp , rv_w_orofi 
     .   , rv_w_orwgt , w ( orhoat ( 2 , ib ) ) ) 


Cdelw1 rlse name= orofi old_list= orofi orwgt 
Cdelw1 rlse name= orofi new_list= (None)
Cdelw1         call rlse(orofi)
         if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)
         if (allocated(rv_w_orofi)) deallocate(rv_w_orofi)


   10   continue
      enddo

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orhoat]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine prlrho 


