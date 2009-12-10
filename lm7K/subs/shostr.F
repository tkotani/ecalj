      subroutine shostr(nds,nsite,nbas,plat,bas,lio,alpha,
     .                  iax,ntab,s,nkaps,nkapn,ekap,io,scale)
C- Prints out structure constants
C ----------------------------------------------------------------------
Ci Inputs
Ci   nds   :leading dimension of s, alpha
Ci   nsite :number of pairs in neighbor table
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   lio   :specifies data contents; patterned after iostr.f
Ci         :digit  specification
Ci         :1-10s   bits in 1-10s digit control program flow.  Not used here
Ci         :100-1000 bits in these digits contain info about structure of s
Ci                  1 s is complex (otherwise s is real)
Ci                  2 s is stored in slices according to site.
Ci                    row ordering that of the cluster determined
Ci                    by iax table.
Ci                    Default storage (2's bit zero):
Ci                    s is stored by (nds,nds) blocks
Ci                  4 s is the structure matrix
Ci                    Default: s is the matrix of expansion coffs
Ci                  8 file contains energy derivative of s
Ci         :10000-100000 info about the functions s corresponds to
Ci                  0 2nd generation Hankel functions with Andersen's
Ci                    conventions
Ci                  1 NMTO Hankel functions with Tank and Andersen
Ci                    conventions: alpha=alpha(2,2,nl**2,nbas,nkapn)
Ci                    and s=s(nl**2,nl**2,nsite,nkapn)
Ci                  2 Conventional Hankel functions, Methfessel defs
Ci   alpha :tight-binding screening parameters
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   s     :real-space structure constants, in order of iax table
Ci   nkaps :number of energies a strux contains
Ci   nkapn :number sets of nkaps-strux
Ci   ekap  :muffin-tin zero plus kap2
Ci   io:    0, print out all structure constants
Ci          1, print out only sc's whose ss element differs from last
Ci   scale: Multiply structure constants by scale in printout
Co Outputs
Co   strx are printed to stdo
Cv Verbosity
Cv  >=10: prints out L=0 S for each site
Cv  >=30: prints out and diagonal strux S_LL for each site
Cv  >=40: prints out and full strux S_LL for each site
Cv  
Cv  >30: prints out iax for each site
Cv   40: prints out structure constant matrix for any site whose
Cv       ss element differs from previous
Cv  >40: prints out full structure constant matrix
Cv   50: prints out alpha's too
Cb Bugs
Cb   alpha dimensions should the like dimensions likes s
Cr Remarks
Cu Updates
Cu   06 Aug 06 Redesigned to printout 2-kappa strux
Cu   08 Jun 04 iostr can treat alpha as tral matrix.  New argument lalph
Cu    7 Sep 00 changed conventions for lio
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nds,nsite,nbas,io,nkaps,nkapn,lio,niax,ntab(1)
      parameter (niax=10)
      integer iax(niax,nsite)
      double precision s(nds,nds,nkaps,nkaps,nsite,nkapn),scale,
     .                 plat(3,3),bas(3,*),ekap(2),
     .                 alpha(nds,nbas,nkaps,nkapn)
C Local parameters
      double precision dx
      integer i,j,k,il1,il2,ib,jb,iat2,ix,lm,ipr,getdig,
     .  lgunit,ii,jj,iic,jjc,nlma,nlmb,ip,i0,stdo,isw,
     .  lio23,lio45,nmto,ik,ikap,jkap,ll
      character outs*120,out2*120,stra*7
      logical lcmplx

      dx(ix,iat2,i,j,k) = bas(ix,iat2) - bas(ix,ib) +
     .                    plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k

      call getpr(ipr)
      stdo = lgunit(1)
      lio23 = mod(lio/100,100)
      lio45 = mod(lio/10000,100)
      lcmplx = getdig(lio23,0,2) .eq. 1
      nmto = getdig(lio45,0,2)

      call info8(10,1,0,
     .  ' SHOSTR:  '//
     .  '%?#(n==0)#2nd generation ##%-1j'//
     .  '%?#(n==1)#NMTO ##%-1j'//
     .  '%?#(n==2)#Screened ##'//
     .  '%?#n#sdot#strux#'//
     .  '%?#(n<2)#, OKA defs##'//
     .  '%?#(n==0)#, 1-center##'//
     .  '%?#(n>1)#, %-1j%i energies##'//
     .  '%?#n#, scaled by %d##'//
     .  ' ',
     .  lio45,
     .  getdig(lio23,3,2),
     .  lio45,
     .  getdig(lio23,2,2),
     .  max(nkaps,nkapn),
     .  isw(scale.ne.1),
     .  scale,
     .  0)

      outs = ' '
      call info5(10,0,0,
     .  '%10f%i pairs'
     .  //'%?#n# (%i inequivalent)#%j#'
     .  //'  nkap=%i%-1j, E =%n:1;4,4d',
     .  ntab(nbas+1),ntab(nbas+1)-nsite,nsite,max(nkaps,nkapn),ekap)

      do  10  ik = 1, nkapn
        if (nmto .eq. 1) call info2(10,1,0,
     .    ' --- Energy %i, e*avw**2 = %d ---',ik,ekap(ik))

C     do  i = 1, min(nsite,ntab(nbas+1))
      do  i = 1, ntab(nbas+1)

          ip = 3

        if (i .eq. 1 .or. iax(1,i) .ne. iax(1,i-1)) then
C          if (i .ne. 1 .and. ipr .le. 30 .and. ipr .ge. 10
C     .        .and. out2 .ne. ' ') print '(a)', out2(1:ip)
          i0 = i-1
          out2 = 'S00'
          ip = 3
          ib  = iax(1,i)
          nlma = iax(9,i)
          jb = 0
          if (iax(8,i) .ne. i) then
            j = iax(1,iax(8,i))
            call info5(10,1,0,
     .        ' Site  %i  has  %i  neighbors'//
     .        ' (mapped to pair %i, in cluster of site %i)',
     .        ib, ntab(ib+1)-ntab(ib),iax(8,i),j,0)
          else
          call info5(0,1,0,' Site  %i  has  %i  neighbors',
     .      ib, ntab(ib+1)-ntab(ib), ipr, 0, 0)
C          call info5(0,1,0,' Site  %i  has  %i  neighbors'//
C     .      '%?#(n>40)#:  alpha, B_R''R:##%-1j'//
C     .      '%?#(40>n&n>30)#    diagonal strux B_LL##%-1j'//
C     .      '%?#(30>=n&n>=10)#  L=0 strux B_00##%-1j'//
C     .      '%?#(40>n)# (inequivalent blocks only)##%-1j',
C     .      ib, ntab(ib+1)-ntab(ib), ipr, 0, 0)
C         Handle NMTO case
          stra = ' alpha:'
          if (ipr .gt. 40 .and. nkaps .eq. 1) then
            write(stdo,885) stra,(alpha(lm**2,ib,ik,1),lm=1,1+ll(nlma))
          elseif (ipr .gt. 40) then
            do  ikap = 1, nkaps
            do  jkap = 1, nkaps
              write(stdo,885) stra, (alpha(lm**2,ib,ikap,jkap),
     .          lm=1,1+ll(nlma))
              stra = ' '
            enddo
            enddo
  885       format(a,5f12.6)
          endif
          endif
        endif

        jb = jb+1
        ii = iax(8,i)
        if (ii .eq. 0) ii = i
        if (ipr .ge. 30) then
          write(outs,300) i, i-i0,iax(2,i), iax(1,i), (iax(j,i),j=3,5),
     .      (dx(ix,iax(2,i),iax(3,i),iax(4,i),iax(5,i)), ix=1,3)
  300     format(' pair',i5,'(',i2,')',' R''=',i3,
     .           ' R=',i3,' plat=',3i3,' tau=',3f9.4)
          if (iax(8,i) .ne. i) call awrit1('%a(%i)',outs,80,0,ii)
          call awrit0('%a',outs,80,-stdo)
        endif

        jj = iax(8,max(i-1,1))
        if (jj .eq. 0) jj = max(i-1,1)
        iic = ii
        jjc = jj
        if (lcmplx) then
          iic = 2*iic-1
          jjc = 2*jjc-1
        endif
        if (io .eq. 0 .or. i .eq. 1 .or. ipr .gt. 40 .or.
     .      dabs(s(1,1,1,1,iic,ik) - s(1,1,1,1,jjc,ik)) .gt. 1.e-8) then
  110     continue
          nlmb = iax(9,iic)
          if (ipr .ge. 40) then
            do  ikap = 1, nkaps
            do  jkap = 1, nkaps
            do  il1 = 1, nlmb
C       ... most digits, but ugly
C            if (ipr .ge. 60) then
C              outs = ' '
C              ip = 0
C              do  117  il2 = 1, nlma
C                call bin2a('F11',1,0,scale*s(il1,il2,iic,ik),4,0,
C     .            len(outs),outs,ip)
C  117         continue
C              print '(a)', outs(1:ip)
              if (ipr .gt. 50) then
                print 302, (scale*s(il1,il2,ikap,jkap,iic,ik),
     .            il2=1,nlma)
              elseif (ipr .ge. 40) then
                print 301, (scale*s(il1,il2,ikap,jkap,iic,ik),
     .            il2=1,nlma)
              endif
            enddo
            if (nkaps .gt. 1) write(stdo,'('' ------'')')
            enddo
            enddo
          elseif (ipr .ge. 30) then
            do  ikap = 1, nkaps
            do  jkap = 1, nkaps
              outs = 'SLL '
              ip = 3
              do  il2 = 1, min(nlmb,nlma)
                call bin2a('F6',1,0,scale*
     .            s(il2,il2,ikap,jkap,iic,ik),4,0,len(outs),outs,ip)
              enddo
              print '(a)', outs(1:ip)
            enddo
            enddo
          elseif (ipr .ge. 10) then
            do  ikap = 1, nkaps
            do  jkap = 1, nkaps
              call bin2a('F8',1,0,scale*s(1,1,ikap,jkap,iic,ik),4,0,
     .          len(outs),out2,ip)
              if (ip .ge. 72) then
                print '(a)', out2(1:ip)
                out2 = ' '
                ip = 3
              endif
            enddo
            enddo
          endif
          if (lcmplx .and. mod(iic,2) .eq. 1) then
            iic = iic+1
            print '(1x)'
            goto 110
          endif
        endif

      enddo
      if (i .ne. 1 .and. ipr .le. 30 .and. ipr .ge. 10
     .    .and. out2 .ne. ' ') print '(a)', out2(1:ip)

   10 continue
  301 format(1x,9f8.4)
  302 format(1x,9f12.8)
C 303 format(' SLL ',9f8.4)
      end

