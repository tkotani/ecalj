subroutine xlgen(plat,rmax,rmax2,nvmax,opts,mode,nv,vecs)
  !- Generate a list of lattice vectors, subject to constraints
  ! ----------------------------------------------------------------
  !i Inputs
  !i   plat  :dimensionless primitive lattice vectors
  !i   rmax  :largest radius for vector
  !i   rmax2 :if nonzero, largest radius for vector after
  !i         :multiples of plat are added
  !i   nvmax :maximum number of vectors allowed
  !i   opts  :1s digit:
  !i           1 add +/- any plat(j) to all other lattice vectors
  !i             found if plat(j) not in original list. (Ewald sums)
  !i          10s digit:
  !i           1 sort lattice vectors by increasing length
  !i           2 return nv only (or upper limit if 1s digit of opts is 1)
  !i           4 for padding, add a single vector and its reciprocal
  !i             instead of replicating entire original set.
  !i         100s digit:
  !i           1 return in vecs the multiples of plat(j) that
  !i             make up the lattice vector
  !i   mode  :vector of length 3 governing shifts along selected axes.
  !i         :0 suppresses shifts along plat(j)
  !i         :1 same as 0
  !i         :2 shifts to minimize length of pos
  !o Outputs
  !o   nv    :number of vectors found
  !o   vecs  :list of vectors
  !u Updates
  !u  17 Mar 04 Bug fix for rmax2
  !u   2 Mar 04 New rmax2: truncate radius of lattice vectors to rmax2
  !u            when list has to be padded in order to include at
  !r            least one lattice vector.
  ! ----------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nvmax,opts,mode(3)
  double precision :: plat(3,3),vecs(3,1),rmax,rmax2
  ! ... Local parameters
  double precision :: rsqr,v2,vj(3)
  integer :: i,j,k,imx(3),nv,m,ivck(3),iprint,lgunit,iv,jv,oiwk,owk
  !      integer w(1)
  !      common /w/ w
  integer,allocatable:: w_oiwk(:)
  real(8),allocatable:: w_owk(:)
  !     call setpr(110)

  ! --- Setup ---
  call latlim(plat,rmax,imx(1),imx(2),imx(3))
  do  10  i = 1, 3
     !   ... Switches flagging whether this plat in lattice vectors
     ivck(i) = 0
     if (mod(opts,10) == 1) ivck(i) = 1
     imx(i) = max(imx(i),1)
     if (mode(i) == 0) then
        imx(i) = 0
        ivck(i) = 0
     endif
10 enddo
  rsqr = rmax*rmax
  nv = 0

  ! --- Loop over all triples, paring out those within rmax ---
  do  202  i = -imx(1), imx(1)
     do  201  j = -imx(2), imx(2)
        do  20  k = -imx(3), imx(3)
           v2 = 0
           do  21  m = 1, 3
              vj(m) = i*plat(m,1) + j*plat(m,2) + k*plat(m,3)
              v2 = v2 + vj(m)**2
21         enddo

           !   --- A lattice vector found ---
           if (v2 > rsqr) goto 20

           !   ... Flag any plat in this vec as being present
           if (iabs(i) + iabs(j) + iabs(k) == 1) then
              if (i == 1) ivck(1) = 0
              if (j == 1) ivck(2) = 0
              if (k == 1) ivck(3) = 0
           endif

           !   ... Increment nv and copy to vec(nv)
           nv = nv+1
           if (nv > nvmax .AND. mod(opts/10,10) /= 2) &
                call fexit(-1,111,' xlgen: too many vectors, n=%i',nv)
           if (mod(opts/10,10) == 2) then
           elseif (mod(opts/100,10) == 1) then
              vecs(1,nv) = i
              vecs(2,nv) = j
              vecs(3,nv) = k
           else
              vecs(1,nv) = vj(1)
              vecs(2,nv) = vj(2)
              vecs(3,nv) = vj(3)
           endif
20      enddo
201  enddo
202 enddo

  ! --- Add plat if ivck ne 0 ---
  if (ivck(1)+ivck(2)+ivck(3) /= 0) then
     if (mod(opts/10,10) == 2) then
        nv = 3*nv
     elseif (mod(opts/10,10) == 4) then
        do  33 i = 1, 3
           if (ivck(i) == 1) then
              call dcopy(3,plat(1,i),1,vj,1)
              if (mod(opts/100,10) == 1) then
                 vj(1) = 0
                 vj(2) = 0
                 vj(3) = 0
                 vj(i) = ivck(i)
              endif
              vecs(1,nv+1) =  vj(1)
              vecs(2,nv+1) =  vj(2)
              vecs(3,nv+1) =  vj(3)
              vecs(1,nv+2) = -vj(1)
              vecs(2,nv+2) = -vj(2)
              vecs(3,nv+2) = -vj(3)
              nv = nv+2
           endif
33      enddo
     else
        if (iprint() >= 20) print 333, ivck, rmax, rmax2
333     format(/' xlgen: added missing plat: ivck=',3i2, &
             '  rmax=',f8.3,'  rpad*rmax=',f8.3)
        if (3*nv > nvmax) &
             call fexit(-1,111,' xlgen: too many vectors, n=%i',3*nv)
        if (ivck(1)+ivck(2)+ivck(3) /= 1) &
             call rx('lgen: more than 1 missing plat')
        do  31  m = 1, 3
           v2 = ivck(1)*plat(m,1)+ivck(2)*plat(m,2)+ivck(3)*plat(m,3)
           if (mod(opts/100,10) == 1) v2 = ivck(m)
           call dcopy(nv,vecs(m,1),3,vecs(m,nv+1),3)
           call dcopy(nv,vecs(m,1),3,vecs(m,2*nv+1),3)
           call daxpy(nv, 1d0,v2,0,vecs(m,nv+1),3)
           call daxpy(nv,-1d0,v2,0,vecs(m,2*nv+1),3)
31      enddo
        nv = 3*nv

        !   ... Find and eliminate any replicas
        !          call defi(oiwk,nv)
        !          call dvshel(3,nv,vecs,w(oiwk),1)
        allocate(w_oiwk(nv))
        call dvshel(3,nv,vecs,w_oiwk,1)
        !       call awrit2('%n:1i',' ',80,6,nv,w(oiwk))
        k = 0
        !   ... Mark any replica iv by iwk(i) -> -iv
        do  32  i = nv-1, 1, -1
           !            iv = w(oiwk+i) + 1
           !            jv = w(oiwk+i-1) + 1
           iv = w_oiwk(i+1) + 1
           jv = w_oiwk(i) + 1
           v2 = (vecs(1,iv)-vecs(1,jv))**2 + &
                (vecs(2,iv)-vecs(2,jv))**2 + &
                (vecs(3,iv)-vecs(3,jv))**2
           !            if (v2 .lt. 1d-10) w(oiwk+i) = -iv
           if (v2 < 1d-10) w_oiwk(i+1) = -iv
32      enddo

        !   ... Flag vectors with radius > rmax2
        if (rmax2 > 0) then
           rsqr = rmax2*rmax2
           k = 0
           do  37  i = 0, nv-1
              !              if (w(oiwk+i) .ge. 0) then
              !                iv = w(oiwk+i) + 1
              if (w_oiwk(i+1) >= 0) then
                 iv = w_oiwk(i+1) + 1
                 v2 = vecs(1,iv)**2 + vecs(2,iv)**2 + vecs(3,iv)**2
                 if (v2 > rsqr) then
                    if (iv < 0) call rx('bug in xlgen')
                    !                  w(oiwk+i) = -iv
                    w_oiwk(i+1) = -iv
                    k = k+1
                 endif
              endif
37         enddo
           !         print *, 'rpad reduced nv by',k
        endif
        !       call prmx('unsorted vecs',vecs,3,3,nv)

        !   ... Make a sorted list of replicas (any of iwk < 0)
        call ishell(nv,w_oiwk)
        !          call ishell(nv,w(oiwk))
        !   ... For each replica, put lastmost vec into replica's place
        k = nv
        do  34  i = 0, nv-1
           !            iv = -w(oiwk+i)
           iv = -w_oiwk(i+1)
           if (iv <= 0) goto 35
           call dpcopy(vecs(1,k),vecs(1,iv),1,3,1d0)
           k = k-1
34      enddo
35      continue
        nv = k
        !          call rlse(oiwk)
        deallocate(w_oiwk)
     endif
  endif
  !     call prmx('after purging vecs',vecs,3,3,nv)


  if (mod(opts/10,10) == 2) return

  ! --- Sort vectors by increasing length ---
  if (mod(opts/10,10) == 1) then
     !        call defi(oiwk,nv)
     !        call dvshel(3,nv,vecs,w(oiwk),11)
     !        call defi(owk,nv*3)
     !        call dvperm(3,nv,vecs,w(owk),w(oiwk),.true.)
     !        call rlse(oiwk)
     allocate(w_oiwk(nv))
     call dvshel(3,nv,vecs,w_oiwk,11)
     allocate(w_owk(nv*3))
     call dvperm(3,nv,vecs,w_owk,w_oiwk,.true.)
     deallocate(w_owk,w_oiwk)
  endif

  ! --- Printout ---
  !$$$      if (iprint() .le. 70) return
  !$$$      call awrit5(' xlgen: opts=%i  mode=%3:1i  rmax=%;4d  plat='//
  !$$$     .'%9:1;4d  nv=%i',' ',80,lgunit(1),opts,mode,rmax,plat,nv)
  !$$$      if (iprint() .lt. 110) return
  !$$$      print 345
  !$$$  345 format('  iv',6x,'px',8x,'py',8x,'pz',8x,'l')
  !$$$      do  40  i = 1, nv
  !$$$        v2 = vecs(1,i)**2 + vecs(2,i)**2 + vecs(3,i)**2
  !$$$        print 346, i, (vecs(m,i), m=1,3), dsqrt(v2)
  !$$$  346   format(i4,3f10.4,f10.3)
  !$$$   40 continue
end subroutine xlgen
!      subroutine fmain

!      implicit none
!      integer wksize,nvmx,nv,opts,mode(3),ovecs
!      double precision plat(9),rmax
!      parameter(wksize=500000)
!      integer w(wksize)
!      common /w/ w

!      data plat /.5d0,.5d0,0d0, .5d0,-.5d0,0d0, 0d0,2d0,2d0/


!      call wkinit(wksize)
!      nvmx = 500
!      rmax = 2

!      mode(1) = 0
!      mode(2) = 2
!      mode(3) = 2
!      opts = 0
!      call initqu(.true.)
!      call query('opts=?',2,opts)
!      call defrr(ovecs,3*nvmx)
!      call lgen(plat,rmax,nv,nvmx,w(ovecs))
!      print *, 'old lgen found nv=', nv
!      print *, 'call xlgen, opt=',opts
!      call xlgen(plat,rmax,nvmx,opts,mode,nv,w(ovecs))
!      end

subroutine ishell(n,iarray)
  integer :: n
  integer :: iarray(1)
  integer :: lognb2,i,j,k,l,m,nn,it
  if (n <= 1) return
  lognb2 = int(log(float(n+1))*1.4426950)
  m = n
  do  12  nn = 1, lognb2
     m = m/2
     k = n - m
     do  11  j = 1, k
        i = j
3       continue
        l = i + m
        if (iarray(l) < iarray(i)) then
           it = iarray(i)
           iarray(i) = iarray(l)
           iarray(l) = it
           i = i - m
           if (i >= 1) goto 3
        endif
11   enddo
12 enddo
  return
end subroutine ishell

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine lgen(bas,bmax,nv,nvmax,vecs,work)
  !  generates lattice vectors.
  implicit real*8 (a-h,p-z), integer(o)
  implicit integer(i-n)
  dimension bas(3,3),v(3),vecs(3,*),work(*) ! MIZUHO-IR
  call latlim(bas,bmax,imax,jmax,kmax)
  bmax2=bmax*bmax
  nv=0
  do 202 i=-imax,imax
     do 201 j=-jmax,jmax
        do 20 k=-kmax,kmax
           do 21 m=1,3
              v(m)=i*bas(m,1)+j*bas(m,2)+k*bas(m,3)
21         enddo
           v2=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
           if(v2 > bmax2) goto 20
           nv=nv+1
           if(nv > nvmax) write(6,633) nvmax,i,imax
           if(nv > nvmax) call rx( '')
633        format(/' --- nv=',i6,'  exceeded,   i=',i3,'  imax=',i3)
           do 22 m=1,3
              vecs(m,nv)=v(m)
22         enddo
           vsm=dabs(v(1))+dabs(v(2))+dabs(v(3))
           work(nv)=v2+vsm/1000.
20      enddo
201  enddo
202 enddo
  ! --- sort by length -----------
  do 30 iv=1,nv
     ilow=iv
     alow=work(iv)
     do 31 jv=iv,nv
        if(work(jv) < alow) then
           alow=work(jv)
           ilow=jv
        endif
31   enddo
     if(ilow == iv) goto 30
     do 32 m=1,3
        xx=vecs(m,iv)
        vecs(m,iv)=vecs(m,ilow)
        vecs(m,ilow)=xx
32   enddo
     work(ilow)=work(iv)
     xx=work(ilow)
30 enddo
  ! ---- add neighbor layers if basis vec 3 is not in list ------
  do 41 iv=1,nv
     ddd=(bas(1,3)-vecs(1,iv))**2+(bas(2,3)-vecs(2,iv))**2 &
          +(bas(3,3)-vecs(3,iv))**2
     if(ddd < 1.d-8) return
41 enddo
  write(6,650)
650 format(/' basis vec 3 not in list - include 2 more planes')
  if(3*nv > nvmax) write(6,643) nvmax
  if(3*nv > nvmax) call rx( '')
643 format( '--- lgen needs nvmax at least',i7)
  do iv=1,nv
     do m=1,3
        vecs(m,iv+nv)=vecs(m,iv)+bas(m,3)
        vecs(m,iv+2*nv)=vecs(m,iv)-bas(m,3)
     enddo
  enddo
  nv=3*nv
  return
end subroutine lgen