      subroutine mkppar(dclabl,nl,nsp,nc,lmx,z,rmax,avw,amsh,nrmsh,
     .  pnu,idmod,ves,qnu,pp)
C- Make pp's for all atoms for which potential is available.
C     implicit none
      integer nl,nsp,lmx(1),nrmsh(1),nc,idmod(nl,nc)
      double precision z(1),rmax(1),avw,amsh(1),pnu(nl,nsp,1),
     .  qnu(3,nl,nsp,nc),ves(nc),pp(6,nl,nsp,nc),dclabl(1)
C Heap allocation
      integer W(1)
      common /w/ w
      logical aiopot,sw,aiopar
      integer nrmx,nr,orofi,ov,ic,og,ogp,ir,ifi,fopna,iprint,i1mach
      integer lrel,lgrad,idu(5)
      parameter (nrmx=1501)
      double precision a,b,ea,thrpv,exc(2),thrpvl(10),rpb,rmx
      double precision ekap,sop
      logical lso
      character*8 clabl

      call iinit(idu,5)
      lso = .false.
      ekap = 0d0

      do  10  ic = 1, nc
        call r8tos8(dclabl(ic),clabl)
        rmx = rmax(ic)
        nr = nrmsh(ic)
        a  = amsh(ic)
        call rx('lrel,lgrad not set')
c        call rmesh(z(ic),rmx,lrel,lgrad,nrmx,a,b,nr)
        call defdr(ov,   nr*nsp)
        call defdr(orofi,nr)
        ifi = fopna(clabl,30,0)
        if (aiopot(nr,nsp,a,rmx,-99d0,w(ov),ifi)) then
          ea = dexp(a)
          b = rmx/(dexp(a*(nr-1)) - 1d0)
          rpb = b
          do  11  ir = 1, nr
          call dvset(w(orofi),ir,ir,rpb-b)
   11     rpb = rpb*ea
          call defdr(og,   nr*2)
          call defdr(ogp,  nr*2*4)
          if (iprint() .ge. 60) then
            print *, 'mkppar: potential parms before calling potpar'
            sw = aiopar(clabl,0,pp(1,1,1,ic),w,ves(ic),nl,lmx(ic),nsp,
     .        -i1mach(2))
          endif
          call rx('need mpolp')
          call potpar(nl,nsp,lmx(ic),z(ic),rmx,avw,ekap,lso,.false.,
     .      .false.,a,nr,w(orofi),w(ov),pnu(1,1,ic),idmod(1,ic),ves(ic),
     .      exc,qnu(1,1,1,ic),idu,w,w,thrpv,thrpvl,w(og),w(ogp),
     .      pp(1,1,1,ic),w,sop,w,w)

C     ... Shift enu and c by ves
          call daxpy(nl*nsp,1d0,ves(ic),0,pp(1,1,1,ic),6)
          call daxpy(nl*nsp,1d0,ves(ic),0,pp(2,1,1,ic),6)

          if (iprint() .ge. 60) then
            print *, 'mkppar: potential parms after calling potpar'
            sw = aiopar(clabl,0,pp(1,1,1,ic),w,ves(ic),nl,lmx(ic),nsp,
     .        -i1mach(2))
          endif
        else
          if (iprint() .gt. 20) print *,
     .      'mkppar: missing potential for class ',clabl
        endif
        call fclose(ifi)
        call rlse(ov)
   10 continue
      end

