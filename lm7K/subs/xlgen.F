      subroutine xlgen(plat,rmax,rmax2,nvmax,opts,mode,nv,vecs)
C- Generate a list of lattice vectors, subject to constraints
C ----------------------------------------------------------------
Ci Inputs
Ci   plat  :dimensionless primitive lattice vectors
Ci   rmax  :largest radius for vector
Ci   rmax2 :if nonzero, largest radius for vector after
Ci         :multiples of plat are added
Ci   nvmax :maximum number of vectors allowed
Ci   opts  :1s digit:
Ci           1 add +/- any plat(j) to all other lattice vectors
Ci             found if plat(j) not in original list. (Ewald sums)
Ci          10s digit:
Ci           1 sort lattice vectors by increasing length
Ci           2 return nv only (or upper limit if 1s digit of opts is 1)
Ci           4 for padding, add a single vector and its reciprocal
Ci             instead of replicating entire original set.
Ci         100s digit:
Ci           1 return in vecs the multiples of plat(j) that
Ci             make up the lattice vector
Ci   mode  :vector of length 3 governing shifts along selected axes.
Ci         :0 suppresses shifts along plat(j)
Ci         :1 same as 0
Ci         :2 shifts to minimize length of pos
Co Outputs
Co   nv    :number of vectors found
Co   vecs  :list of vectors
Cu Updates
Cu  17 Mar 04 Bug fix for rmax2
Cu   2 Mar 04 New rmax2: truncate radius of lattice vectors to rmax2
Cu            when list has to be padded in order to include at
Cr            least one lattice vector.
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nvmax,opts,mode(3)
      double precision plat(3,3),vecs(3,1),rmax,rmax2
C ... Local parameters
      double precision rsqr,v2,vj(3)
      integer i,j,k,imx(3),nv,m,ivck(3),iprint,lgunit,iv,jv,oiwk,owk
c      integer w(1)
c      common /w/ w
      integer,allocatable:: w_oiwk(:)
      real(8),allocatable:: w_owk(:)
C     call setpr(110)

C --- Setup ---
      call latlim(plat,rmax,imx(1),imx(2),imx(3))
      do  10  i = 1, 3
C   ... Switches flagging whether this plat in lattice vectors
        ivck(i) = 0
        if (mod(opts,10) .eq. 1) ivck(i) = 1
        imx(i) = max(imx(i),1)
        if (mode(i) .eq. 0) then
          imx(i) = 0
          ivck(i) = 0
        endif
   10 continue
      rsqr = rmax*rmax
      nv = 0

C --- Loop over all triples, paring out those within rmax ---
      do  20  i = -imx(1), imx(1)
      do  20  j = -imx(2), imx(2)
      do  20  k = -imx(3), imx(3)
        v2 = 0
        do  21  m = 1, 3
          vj(m) = i*plat(m,1) + j*plat(m,2) + k*plat(m,3)
          v2 = v2 + vj(m)**2
   21   continue

C   --- A lattice vector found ---
        if (v2 .gt. rsqr) goto 20

C   ... Flag any plat in this vec as being present
        if (iabs(i) + iabs(j) + iabs(k) .eq. 1) then
          if (i .eq. 1) ivck(1) = 0
          if (j .eq. 1) ivck(2) = 0
          if (k .eq. 1) ivck(3) = 0
        endif

C   ... Increment nv and copy to vec(nv)
        nv = nv+1
        if (nv .gt. nvmax .and. mod(opts/10,10) .ne. 2) 
     .      call fexit(-1,111,' xlgen: too many vectors, n=%i',nv)
        if (mod(opts/10,10) .eq. 2) then
        elseif (mod(opts/100,10) .eq. 1) then
          vecs(1,nv) = i
          vecs(2,nv) = j
          vecs(3,nv) = k
        else
          vecs(1,nv) = vj(1)
          vecs(2,nv) = vj(2)
          vecs(3,nv) = vj(3)
        endif
   20 continue

C --- Add plat if ivck ne 0 ---
      if (ivck(1)+ivck(2)+ivck(3) .ne. 0) then
        if (mod(opts/10,10) .eq. 2) then
          nv = 3*nv
        elseif (mod(opts/10,10) .eq. 4) then
          do  33 i = 1, 3
            if (ivck(i) .eq. 1) then
              call dcopy(3,plat(1,i),1,vj,1)
              if (mod(opts/100,10) .eq. 1) then
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
   33     continue
        else
          if (iprint() .ge. 20) print 333, ivck, rmax, rmax2
  333     format(/' xlgen: added missing plat: ivck=',3i2,
     .    '  rmax=',f8.3,'  rpad*rmax=',f8.3)
          if (3*nv .gt. nvmax)
     .    call fexit(-1,111,' xlgen: too many vectors, n=%i',3*nv)
          if (ivck(1)+ivck(2)+ivck(3) .ne. 1)
     .    call rx('lgen: more than 1 missing plat')
          do  31  m = 1, 3
            v2 = ivck(1)*plat(m,1)+ivck(2)*plat(m,2)+ivck(3)*plat(m,3)
            if (mod(opts/100,10) .eq. 1) v2 = ivck(m)
            call dcopy(nv,vecs(m,1),3,vecs(m,nv+1),3)
            call dcopy(nv,vecs(m,1),3,vecs(m,2*nv+1),3)
            call daxpy(nv, 1d0,v2,0,vecs(m,nv+1),3)
            call daxpy(nv,-1d0,v2,0,vecs(m,2*nv+1),3)
   31     continue
          nv = 3*nv

C   ... Find and eliminate any replicas
c          call defi(oiwk,nv)
c          call dvshel(3,nv,vecs,w(oiwk),1)
          allocate(w_oiwk(nv))
          call dvshel(3,nv,vecs,w_oiwk,1)
*       call awrit2('%n:1i',' ',80,6,nv,w(oiwk))
          k = 0
C   ... Mark any replica iv by iwk(i) -> -iv
          do  32  i = nv-1, 1, -1
c            iv = w(oiwk+i) + 1
c            jv = w(oiwk+i-1) + 1
            iv = w_oiwk(i+1) + 1
            jv = w_oiwk(i) + 1
            v2 = (vecs(1,iv)-vecs(1,jv))**2 +
     .      (vecs(2,iv)-vecs(2,jv))**2 + 
     .      (vecs(3,iv)-vecs(3,jv))**2
c            if (v2 .lt. 1d-10) w(oiwk+i) = -iv
            if (v2 .lt. 1d-10) w_oiwk(i+1) = -iv
   32     continue

C   ... Flag vectors with radius > rmax2
          if (rmax2 .gt. 0) then
            rsqr = rmax2*rmax2
            k = 0
            do  37  i = 0, nv-1
c              if (w(oiwk+i) .ge. 0) then
c                iv = w(oiwk+i) + 1
              if (w_oiwk(i+1) .ge. 0) then
                iv = w_oiwk(i+1) + 1
                v2 = vecs(1,iv)**2 + vecs(2,iv)**2 + vecs(3,iv)**2
                if (v2 .gt. rsqr) then
                  if (iv .lt. 0) call rx('bug in xlgen')
c                  w(oiwk+i) = -iv
                  w_oiwk(i+1) = -iv
                  k = k+1
                endif
              endif
   37       continue
C         print *, 'rpad reduced nv by',k
          endif
C       call prmx('unsorted vecs',vecs,3,3,nv)

C   ... Make a sorted list of replicas (any of iwk < 0)
          call ishell(nv,w_oiwk)
c          call ishell(nv,w(oiwk))
C   ... For each replica, put lastmost vec into replica's place
          k = nv
          do  34  i = 0, nv-1
c            iv = -w(oiwk+i)
            iv = -w_oiwk(i+1)
            if (iv .le. 0) goto 35
            call dpcopy(vecs(1,k),vecs(1,iv),1,3,1d0)
            k = k-1
   34     continue
   35     continue
          nv = k
c          call rlse(oiwk)
          deallocate(w_oiwk)
        endif
      endif
C     call prmx('after purging vecs',vecs,3,3,nv)


      if (mod(opts/10,10) .eq. 2) return

C --- Sort vectors by increasing length ---
      if (mod(opts/10,10) .eq. 1) then
c        call defi(oiwk,nv)
c        call dvshel(3,nv,vecs,w(oiwk),11)
c        call defi(owk,nv*3)
c        call dvperm(3,nv,vecs,w(owk),w(oiwk),.true.)
c        call rlse(oiwk)
        allocate(w_oiwk(nv))
        call dvshel(3,nv,vecs,w_oiwk,11)
        allocate(w_owk(nv*3))
        call dvperm(3,nv,vecs,w_owk,w_oiwk,.true.)
        deallocate(w_owk,w_oiwk)
      endif

C --- Printout ---
      if (iprint() .le. 70) return
      call awrit5(' xlgen: opts=%i  mode=%3:1i  rmax=%;4d  plat='//
     .'%9:1;4d  nv=%i',' ',80,lgunit(1),opts,mode,rmax,plat,nv)
      if (iprint() .lt. 110) return
      print 345
  345 format('  iv',6x,'px',8x,'py',8x,'pz',8x,'l')
      do  40  i = 1, nv
        v2 = vecs(1,i)**2 + vecs(2,i)**2 + vecs(3,i)**2
        print 346, i, (vecs(m,i), m=1,3), dsqrt(v2)
  346   format(i4,3f10.4,f10.3)
   40 continue
      end
C      subroutine fmain
C
C      implicit none
C      integer wksize,nvmx,nv,opts,mode(3),ovecs
C      double precision plat(9),rmax
C      parameter(wksize=500000)
C      integer w(wksize)
C      common /w/ w
C
C      data plat /.5d0,.5d0,0d0, .5d0,-.5d0,0d0, 0d0,2d0,2d0/
C
C
C      call wkinit(wksize)
C      nvmx = 500
C      rmax = 2
C
C      mode(1) = 0
C      mode(2) = 2
C      mode(3) = 2
C      opts = 0
C      call initqu(.true.)
C      call query('opts=?',2,opts)
C      call defrr(ovecs,3*nvmx)
C      call lgen(plat,rmax,nv,nvmx,w(ovecs))
C      print *, 'old lgen found nv=', nv
C      print *, 'call xlgen, opt=',opts
C      call xlgen(plat,rmax,nvmx,opts,mode,nv,w(ovecs))
C      end

