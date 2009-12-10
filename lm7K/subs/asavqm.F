      subroutine asavqm(mode,sctrl,spot,slat,sarray,sspec,nlmf,vh,qmp,
     .  vval)
C- Electrostatic Multipole moments from wave function products
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  1s digit:  
Ci   
Ci   sarray->qpp:
Ci      qpp(*,1,ib) phi-phi   products for (L<=L')
Ci      qpp(*,2,ib) phi-dphi  products for all (L,L')
Ci      qpp(*,4,ib) dphi-dphi products for (L<=L')
Co Outputs
Co   Electrostatic potential
Co   vh
Co   qmp
Co   vval
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nlmf
      double precision sctrl(1),spot(1),slat(1),sarray(1),sspec(1),
     .  vval(nlmf,*),qmp(nlmf,*),vh(*)
C heap:
      integer w(1)
      common /w/ w
C Local variables
      double precision alat,plat(9),qlat(9),awald,vol,tau(3)
      integer ib,ipr,iprint,is,jb,jc,ll,lmax,lmxb,lmxl,nbas,nclasp,
     .  nkd,nkq,nl,nlm0,nlmb,nlmbx,nlml,nlmp,npow,nrx,
     .  ic,ocf,ocg,odlv,oics,oidxcg,oikl,
     .  oip,oipc,oips,ojcg,ojj,ojkl(0:20),olmxa,olmxf,olmxl,
     .  opos,opot0,oqlv,ormax,os,owk,oyl,oqc,oz,onrc
      parameter (nlm0=121)
      double precision dl(2,nlm0),hl(nlm0),rmax,yy,hh,xi(0:20)
      double precision valj(nlm0),pot0(nlm0),q(3)
      double precision dval
      character*4 cib

C ... Setup
      call upack('array opos oipc oips oics',sarray,opos,oipc,oips,oics,
     .  0)
      call upack('array nclasp onrc ormax',sarray,nclasp,onrc,ormax,0,0)
      call upack1('pot oqc',spot,oqc)
C     call upack1('array omad',sarray,omad)

      call upack('ctrl nbas nl',sctrl,nbas,nl,0,0,0)

      call upack('lat alat plat qlat vol',slat,alat,plat,qlat,vol,0)
      call upack('lat ocg ojcg oidxcg',slat,ocg,ojcg,oidxcg,0,0)
      call upack('lat awald nkd nkq odlv oqlv',slat,awald,nkd,nkq,odlv,
     .  oqlv)
      ipr = iprint()
      if (mode .ne. 0) call rxi('asavqm: not ready for mode',mode)
C     y0 = 1/sqrt(4*pi)

      call sp2cls('spec lmxa',sspec,w(oics),1,1,nclasp,olmxa)
      call sp2cls('spec lmxf',sspec,w(oics),1,1,nclasp,olmxf)
      call sp2cls('spec lmxl',sspec,w(oics),1,1,nclasp,olmxl)
      call sp2cls('spec z',sspec,w(oics),1,1,nclasp,oz)

C     call iinit(w(olmxl),nclasp)

C --- Multipole moments and cofficients to Hankel expansion ---
      call defrr(opot0,nlmf*nbas)
      call pshpr(ipr-10)
      call asaqmp(0,sctrl,spot,slat,sarray,w(olmxa),w(olmxf),nlmf,qmp)
      call poppr
      call pvvqm1(nbas,w(oipc),w(olmxf),w(oz),w(oqc),nlmf,qmp,
     .  w(opot0))

C --- For each site, calculate electrostatic at MT boundary ---
C ... Expand about ib from jb
      do  10  ib = 1, nbas
      write(cib,'(i4)') ib
      is = w(oips+ib-1)
      ic = w(oipc+ib-1)
      call upack('spec lmxl',sspec,is,lmxl,rmax,0,0)
      rmax = dval(w(ormax),ic)
      lmxl = w(olmxl+ic-1)
      nlml = (lmxl+1)**2
      
C ... Setup for reduced strux
      lmax = 2*ll(nlmf)
      nrx = max(nkd,nkq)
      call defrr(owk, nrx*(2*lmax+10))
      call defrr(oyl, nrx*(lmax+1)**2)
C ... Setup for strux
      call strxsu(nlml,-99,2*nl-2,1,0,nbas,w(oips),w(ocg),w(ojcg),
     .  w(oidxcg),nlmbx,nlmp,npow,ocf,oip,oikl,ojkl)
      if (nlmp .gt. nlm0) call rxi('asavqm: increase nlm0 to',nlmp)
      call defrr(os, nlml*nlmbx)

C ... Loop over all connecting vectors to ib
      do  20  jb = 1, nbas
C ... lmxb is for now the l-cutoff for the head
      jc = w(oipc+ib-1)
      lmxb = w(olmxf+jc-1)
      nlmb = (lmxb+1)**2
C ... Connecting vector is pos(jb)-pos(ib) for expansion at ib
      call dpscop(w(opos),tau,3,3*jb-2,1,1d0)
      call dpsadd(tau,w(opos),3,1,3*ib-2,-1d0)
      call shorbz(tau,tau,plat,qlat)
C ... Real reduced strux for this connecting vector
      call dpzero(q,3)
      call hsmqe0(lmxb+lmxl,0d0,110,q,tau,nrx,nlm0,w(owk),w(oyl),
     .  awald,alat,w(oqlv),nkq,w(odlv),nkd,vol,dl)
      call dcopy(nlmp,dl,2,hl,1)
C ... Structure constant matrix
      ojj = ojkl(lmxb)
      call hstrux(0d0,nlml,nlmb,nlmp,npow,1,1,w(oikl),w(ojj),w(oip),
     .  w(ocf),hl,w(os))
C ... ropbes avoids OKA conventions for Bessel functions
      call ropbes(rmax,0d0,lmxl,yy,hh,xi,1)
      call dpscop(w(opot0),pot0,nlmb,1+nlmf*(jb-1),1,1d0)
      call dpzero(valj,nlml)
      call sumtl1(nlml,nlmb,w(os),1,1,pot0,xi,hh,rmax,valj)
      call dpsadd(vval,valj,nlml,1+nlmf*(ib-1),1,1d0)
   20 continue
      call rlse(owk)
   10 continue
C ... Make vh and Printout
      call pvvqm2(1,nbas,w(oipc),w(olmxl),w(oz),w(ormax),w(oqc),nlmf,
     .  qmp,w(opot0),vval,vh)
      call rlse(olmxa)
      end

      subroutine sumtl1(nlml,nlmb,s,nr,nr1,poti,xi,h,rofi,vl)
C- Add bessel tails for potential
C     implicit none
      integer nlmb,nlml,nr,nr1
      double precision vl(nr1,nlml),rofi(1),xi(nr,0:1),s(nlml,nlmb),
     .  poti(nlmb),h(1)
C Local
      integer i,ilmb,ilml,l,ll,lmxl,m
      double precision sum

      lmxl = ll(nlml)
      do  1  i = 1, nr1
    1 h(i) = 1d0
      ilml = 0
      do  16  l = 0, lmxl
        do 14 m = -l,l
          ilml = ilml+1
          sum = 0d0
          do  15  ilmb = 1, nlmb
   15     sum = sum + s(ilml,ilmb)*poti(ilmb)
          do  3  i = 1, nr1
    3     vl(i,ilml) = vl(i,ilml) + sum*h(i)*xi(i,l)
   14   continue
        do  2  i = 1, nr1
    2   h(i) = h(i)*rofi(i)
   16 continue
      end

      subroutine pvvqm1(nbas,ipc,lmxf,z,qc,nlmf,qmp,pot0)
C- Make coffs pot0 from qmom
C     implicit none
      integer nbas,nlmf,ipc(*),lmxf(*)
      double precision qmp(nlmf,nbas),pot0(nlmf,nbas),z(*),qc(1)
      integer ib,ilm,k,l,ll,nlm,ic
C     integer ipr,iclbsj,iprint
      double precision df,pi,y0,xx
      character*4 cib

C     ipr = iprint()
      pi = 4*datan(1d0)
      y0 = 1/sqrt(4*pi)

C     if (ipr .ge. 30) print 649
      do  10  ib = 1, nbas
        ic = ipc(ib)
        write(cib,'(i4)') ib

        nlm = (lmxf(ic)+1)**2
        do  11  ilm = 1, nlm
          l = ll(ilm)
          df = 1d0
          do  3  k = 0, l
    3     df = df*(2*k+1)
          xx = qmp(ilm,ib)
          if (ilm .eq. 1) xx = xx + (qc(ic)-z(ic))*y0
          pot0(ilm,ib) = 2d0*xx*4d0*pi/df
C          if (ipr .ge. 30 .and. dabs(qmp(ilm,ib)) .gt. 1d-5) then
C            if (iclbsj(ic,ipc,nbas,1) .eq. ib .or. ipr .ge. 50) then
C            if (ilm .eq. 1) then
C              print 650, cib,ilm,qmp(ilm,ib),z(ic)-qc(ic),pot0(ilm,ib)
C            else
C              print 651, cib,ilm,qmp(ilm,ib),pot0(ilm,ib)
C            endif
C            endif
C          endif
C         cib = ' '
   11   continue
   10 continue

C  649 format(/'  ib','  ilm',7x,'Qmp',6x,'Z-Qc',7x,'Pot0')
C  650 format(a4,i4,2x,f12.6,f8.3,f12.6)
C  651 format(a4,i4,2x,f12.6,8x,f12.6)

      end

      subroutine pvvqm2(ib1,ib2,ipc,lmxl,z,rmax,qc,nlmf,qmp,pot0,vval,
     .  vh)
C- Printout
C     implicit none
      integer ib1,ib2,nlmf,ipc(*),lmxl(*)
      double precision qmp(nlmf,ib2),pot0(nlmf,ib2),vval(nlmf,ib2),
     .  z(*),rmax(*),qc(*),vh(*)
      logical lpr
      integer ib,ilm,ipr,iprint,nlm,ic,iclbsj
      double precision pi,y0,xx,qt
      character*8 cib

      ipr = iprint()
      pi = 4*datan(1d0)
      y0 = 1/sqrt(4*pi)

      if (ipr .ge. 30) print 649
      do  10  ib = ib1, ib2
        ic = ipc(ib)
        write(cib,'(2i4)') ib,ic

        nlm = (lmxl(ic)+1)**2
        do  11  ilm = 1, nlm
          xx = qmp(ilm,ib)
          if (ilm .eq. 1) xx = xx + (qc(ic)-z(ic))*y0
          lpr = iclbsj(ic,ipc,ib2,1) .eq. ib .and. ipr .ge. 30
     .          .or. ipr .ge. 50
          if (ilm .eq. 1) then
            qt = qmp(ilm,ib)/y0 + (qc(ic)-z(ic))
            vh(ib) = vval(1,ib)*y0 + 2*qt/rmax(ic)
            if (lpr)
     .        print 650, cib,1,qmp(1,ib),z(ic)-qc(ic),pot0(1,ib),
     .        vval(1,ib),vval(1,ib)*y0,vh(ib)
          else
            if (lpr .and. dabs(qmp(ilm,ib)) .gt. 1d-6)
     .        print 651, cib,ilm,qmp(ilm,ib),pot0(ilm,ib),vval(ilm,ib)
          endif
          cib = ' '
   11   continue
   10 continue
      
  649 format('  ib  ic  ilm',7x,'Qmp',6x,'Z-Qc',7x,'Pot0',8x,'Vval',
     .  8x,'Vmad',6x,'Vh(rmt)')
  650 format(a8,i4,2x,f12.6,f8.3,4f12.6)
  651 format(a8,i4,2x,f12.6,8x,2f12.6)

      end

