      subroutine mixpq(nclas,nl,lmx,nsp,pnu,qnu,pold,qold,ido,tjmax,tj,
     .fnam,beta,nmix,mmix,a,namx,rms2,rmsdel,lq)
C- Mix to make a new vector of P and Q.
C ----------------------------------------------------------------
Ci Inputs
Ci   nclas,nl,lmx,nsp
Ci   a:    Workspace of four vectors
Ci   beta: mixing parameter: mixes beta * new  +  (1-beta) * old
Ci         (see function amix, below)
Ci   mmix: number of previous iterations for mixing.
Ci   namx: upper limit to number of elements in the vector
Ci   pold,qold: input for last iteration
Ci   pnu,qnu:   output of last iteration
Ci   abet: automatically adjust beta when true;
Ci   lq:   mix only 0th moments when true;
Ci   ido:  passed to amix to govern control of tj (see amix)
Co Outputs
Co   pnu,qnu are mixed
Co   rmsdel is the rms change in parameters
Cr Remarks
Cr   sign of beta is passed as sign of npmix<0 in function amix
Cr   Improvments: suggestion of PYB, take lc of moms 0 and 2
C ----------------------------------------------------------------
C     implicit none
      logical lq
      integer nclas,nl,nsp,mmix,namx,lmx(0:nclas-1),nmix,ido
      character*(*) fnam
      double precision pnu(0:nl-1,0:nsp-1,0:nclas-1),
     .                 pold(0:nl-1,0:nsp-1,0:nclas-1),
     .                 qnu(3,0:nl-1,0:nsp-1,0:nclas-1),
     .                 qold(3,0:nl-1,0:nsp-1,0:nclas-1),
     .                 a(0:namx-1,0:mmix+1,2),beta,rms2,rmsdel(2),
     .                 tj(1),tjmax
C Local variables:
      logical lmix
      integer na,ic,isp,l,i,j,k,n1,imix,nelts,onorm,okpvt,awrite,
     .        fopn,ipr,amix,jmix,lgunit
      external fopn
      double precision ddot,dsqrt
      character*80 outs

C heap:
      integer w(1)
      common /w/ w

      call getpr(ipr)

C --- Copy P's and Q's into mixing array ---
      na = 0
      do  10  ic = 0, nclas-1
        do  10  isp = 0, nsp-1
        do  10  l = 0, lmx(ic)
          if (lq) then
            a(na,0,1) = qnu(1,l,isp,ic)
            a(na,0,2) = qold(1,l,isp,ic)
            na = na+1
          else
            a(na,0,1) = pnu(l,isp,ic)
            a(na,0,2) = pold(l,isp,ic)
            do  20  i = 1, 3
              a(na+i,0,1) = qnu(i,l,isp,ic)
              a(na+i,0,2) = qold(i,l,isp,ic)
   20       continue
            na = na + 4
          endif
   10 continue
      if (na .ne. namx) call rx('MIXPQ: elements count mismatch')

C --- Anderson mixing setup ---
      imix = 0
      if (mmix .gt. 0) then
        n1 = fopn(fnam)
        read(n1,err=8,end=8) lmix, imix, nelts
  103   format(2i5)
        if (lmix) goto 8
        if (nelts .eq. namx .and. imix .le. mmix) goto 12
        goto 9
    8   continue
        print *, 'MIXPQ:  read error in file ', fnam,' ... restart'
    9   imix = 0
        goto 40
C When reached this point, have found previous iter w/ correct # elts
   12   continue
        read(n1) (((a(i,j,k), i=0,namx-1), j=1, imix), k=1,2)
  104   format(1p,4d18.11)
      endif
C ... Now have read in imix previous iterations, if there were any
   40 continue

C --- Calculate new rmsdel ---
      rms2 =  dsqrt(dabs(ddot(namx,a,1,a,1)
     .        -2*ddot(namx,a,1,a(0,0,2),1)
     .        + ddot(namx,a(0,0,2),1,a(0,0,2),1))/namx)
C ... Historical (JEK) : adjust beta by
C ... min(max((1-newrms/oldrms)/beta,beta/bdelm,betmin),1d0,beta*bdelp)
      if (ipr .ge. 30) then
        print *, ' '
        j = awrite(' MIXPQ:  %i iter in file',
     .    outs,len(outs),0,imix,0,0,0,0,0,0,0)
        outs(j+2:len(outs)) = fnam
        call awrit1('%a.  RMS delta=%1,3;3e',outs,80,0,rms2)
        if (rmsdel(1) .ne. 0)
     .  call awrit1('%a  last it=%1,3;3e',outs,80,0,rmsdel)
        do  15  j = 1, 2
   15   call awrit0('%a',outs,-len(outs),-lgunit(j))
        call query('beta',4,beta)
      endif
      rmsdel(1) = rms2
      rmsdel(2) = beta

C --- do the mixing ---
      call defdr(onorm,mmix**2)
      call defi(okpvt,mmix)
      jmix = imix
      if (beta .lt. 0) jmix = -imix
      nmix = amix(namx,jmix,mmix,ido,dabs(beta),ipr,tjmax,
     .            w(onorm),w(okpvt),a,tj,rms2)

      if (ipr .ge. 40) call query('continue',-1,0)

C --- Save this iteration into mixing file ---
      imix = min(imix+1,mmix)
      if (mmix .gt. 0) then
        rewind n1
        write(n1) .false., imix, namx
        write(n1) (((a(i,j,k), i=0,namx-1), j=1, imix), k=1,2)
        call fclose(n1)
      endif

C --- RMS change in Q alone ---
      na = 0
      rms2 = 0
      do  11  ic = 0, nclas-1
        do  11  isp = 0, nsp-1
        do  11  l = 0, lmx(ic)
          if (lq) then
            qnu(1,l,isp,ic) = a(na,0,2)
            rms2 = rms2 + (qold(1,l,isp,ic) - a(na,0,2))**2
            na = na+1
            pnu(l,isp,ic) = pold(l,isp,ic)
            qnu(2,l,isp,ic) = qold(2,l,isp,ic)
   22       qnu(3,l,isp,ic) = qold(3,l,isp,ic)
          else
            pnu(l,isp,ic) = a(na,0,2)
            rms2 = rms2 + (qold(1,l,isp,ic) - a(na+1,0,2))**2
            do  21  i = 1, 3
              qnu(i,l,isp,ic) = a(na+i,0,2)
   21       continue
            na = na+4
          endif
   11 continue
      rms2 = dsqrt((rms2*4)/na)
      if (lq) rms2 = rms2/2
      if (ipr .ge. 30) print *, ' '

      end

