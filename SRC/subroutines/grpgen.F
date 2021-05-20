      subroutine grpgen(gen,ngen,symops,ng,ngmx)
      use m_lmfinit,only: stdo


C- Generate all point symmetry operations from the generation group
C ----------------------------------------------------------------
Ci Inputs
Ci   gen,ngen,ngmx
Ci   ng  (12 Sep 96): if>0 , add symops to the ng already in list.
Co Outputs
Co   symops,ng
Cr Remarks
Cr   This works for point groups only and is set up for integer
Cr   generators.
C ----------------------------------------------------------------
C     implicit none
      integer ngen,ng,ngmx
      double precision gen(9,ngen),symops(9,ngmx)
      double precision h(9),hh(9),e(9),sig(9),ae(3)
      integer igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,n,m
      integer ipr
      logical grpeql
      character*80 sout
      data e /1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/, ae/0d0,0d0,0d0/

cChangenglob      stdo = nglob('stdo')
c      stdo = globalvariables%stdo
      call getpr(ipr)
      sout = ' '
      call grpcop(e,symops)
      if (ng .lt. 1) ng = 1
      do  80  igen = 1, ngen
        call grpcop(gen(1,igen),sig)
C ---   Extend the group by all products with sig ---
        do  9  ig = 1, ng
          if (grpeql(symops(1,ig),sig) .and. ipr .gt. 30)
     .    call awrit2(' Generator %i already in group as element %i',
     .    ' ',80,stdo,igen,ig)
          if (grpeql(symops(1,ig),sig)) goto 80
    9   continue

C ---   Determine order ---
        call grpcop(sig,h)
        do  1  itry = 1, 100
          iord = itry
          if (grpeql(h,e)) goto 2
          call grpprd(sig,h,h)
    1   continue
C --- Products of type  g1 sig**p g2 ---
    2   nnow = ng
        if(ipr .ge. 40) call awrit2('%a  %i is %i,',sout,80,0,igen,iord)
        do  8  j = 1, ng
          call grpcop(symops(1,j),h)
          do  10  ip = 1, iord-1
C ... h = sig**ip
            call grpprd(sig,h,h)
            do  11  i = 1, ng
C ... hh = symops_i sig**ip
              call grpprd(symops(1,i),h,hh)
              do  12  k = 1, nnow
                if ( grpeql(symops(1,k),hh) ) goto 11
   12         continue
              nnow = nnow+1
              if (nnow .gt. ngmx) goto 99
              call grpcop(hh,symops(1,nnow))
C              print 333, (symops(k,nnow), k=1,9), nnow
C  333         format(9f12.6,i3)
   11       continue
   10     continue
          if (j .eq. 1) n2 = nnow
    8   continue

C --- Products with more than one sandwiched sigma-factor ---
        m1 = ng+1
        m2 = nnow
        do  20  i = 2, 50
          do  21  n = ng+1, n2
          do  21  m = m1, m2
            call grpprd(symops(1,n),symops(1,m),h)
            do  22  k = 1, nnow
              if (grpeql(symops(1,k),h)) goto 21
   22       continue
            nnow = nnow+1
            if (nnow .gt. ngmx) goto 99
            call grpcop(h,symops(1,nnow))
   21     continue
          if (m2 .eq. nnow) goto 25
          m1 = m2 + 1
          m2 = nnow
   20   continue
   25   continue
        ng = nnow
   80 continue

C --- Printout ---
      if (ipr .ge. 30) then
        if (sout .ne. ' ' .and. ipr .ge. 60) call awrit0
     .  (' Order of generator'//sout//'%a%b',' ',80,stdo)
        call awrit2(' GRPGEN: %i symmetry operations from %i '//
     .  'generator(s)',' ',80,stdo,ng,ngen)
      endif
      if (ipr .ge. 80 .and. ng .gt. 1) then
        write(stdo,'('' ig  group op'')')
        do  60  ig = 1, ng
          call asymop(symops(1,ig),ae,' ',sout)
          write(stdo,'(i4,2x,a)') ig,sout(1:35)
   60   continue
      endif

C      if (ipr .ge. 110) then
C        print *, 'group operations:'
C        call ywrm(0,' ',1,i1mach(2),'(5f12.6)',symops,1,9,9,ng)
C      endif
      return
   99 call rx('GRPGEN: too many elements')
      end
      subroutine grpcop(g,h)
C- Copy matrix
C     implicit none
      double precision h(9),g(9)
      integer i
      do  10  i = 1, 9
        h(i) = g(i)
   10 continue
      end
      logical function grpeql(g1,g2)
C- Checks if G1 is equal to G2
C     implicit none
      double precision g1(9),g2(9),dabs,tol,x1,x2
      parameter (tol = 1d-8)
      logical ddif
      integer i
      ddif(x1,x2) = dabs(x1-x2) .gt. tol
      grpeql = .false.
      do  10  i = 1, 9
        if (ddif(g1(i),g2(i))) return
   10 continue
      grpeql = .true.
      end
      subroutine grpprd(g1,g2,g1xg2)
C- Returns the product of two point group operations
C     implicit none
      double precision g1(3,3),g2(3,3),g1xg2(3,3),h(3,3),sum
      integer i,j,k

      do  10  i = 1, 3
      do  10  j = 1, 3
        sum = 0d0
        do  11  k = 1, 3
          sum = sum + g1(i,k)*g2(k,j)
   11   continue
        h(i,j) = sum
   10 continue
      do  12  j = 1, 3
      do  12  i = 1, 3
        g1xg2(i,j) = h(i,j)
   12 continue


      end

