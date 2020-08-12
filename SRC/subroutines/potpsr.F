      subroutine potpsr2(z,rmax,lmxa,v,vdif,a,nr,rofi,g,gp,gpp,
     .  pnu,enu,c,srdel,qpar,ppar,hab,vab,sab,nsp,
     &   gx,gpx)

C- Makes potential parameters for potential v. Boundary condition
c  for phi is pnu=.5-atan(dnu)/pi+(princ.quant.number). scal-rel.
cm  02.09.94  Potential vdif is included by perturbation.
cs  02.01.95  spin polarized
c  make sure g,gp,gpp are big enough... two components.
      implicit real*8 (a-h,p-z), integer(o)
      dimension rofi(1),v(nr,1),pnu(0:8,nsp),konfig(0:10),dnu(0:10),
     .  hab(4,9,nsp),sab(4,9,nsp),vab(4,9,nsp),g(1),gp(1),gpp(1),
     .  vdif(nr,nsp),c(0:8,1),enu(0:8,1),srdel(0:8,1),qpar(0:8,1),
     .  ppar(0:8,1)
      real(8):: gx(2*nr,0:lmxa,nsp),gpx(2*nr,0:lmxa,nsp) ! takao

      ipr = iprint()
      tol = 1.d-8
      tolval = 1d-12
      eb1 = -20d0
      eb2 = 20d0
      pi = 4d0*datan(1d0)
      b = rmax/(dexp(a*nr-a)-1d0)
      rpb = b
      do  5  ir = 1, nr
        rofi(ir) = rpb-b
        rpb = rpb*dexp(a)
    5 continue

c|      write(6,916) v(nr,1),vdif(nr,1)
c|  916 format(' at rmax: v,vdif=',f10.5,f10.5)

C --- Loop over spins, lmxa ---
      do  80  i = 1, nsp
        do  1  l = 0, lmxa
          konfig(l) = pnu(l,i)
          dnu(l) = dtan(pi*(0.5d0-pnu(l,i)))
    1   continue
        if (nsp .eq. 2 .and. ipr .ge. 40) then
          write(6,550) i,(pnu(l,i),l=0,lmxa)
          write(6,552) (dnu(l),l=0,lmxa)
        elseif (ipr .ge. 40) then
          write(6,551) (pnu(l,i),l=0,lmxa)
          write(6,552) (dnu(l),l=0,lmxa)
        endif
  550   format(/' potpsr:  spin: ',i1,
     .  '   input log derivatives are'/(' pnu=',7f11.5))
  551   format(/' potpsr:  input log derivatives are'/(' pnu=',7f11.5))
  552   format(' dnu=',7f11.5)
        if (ipr .ge. 40) write(6,698)
  698   format(/'  l',7x,'enu',9x,'v',11x,'c',10x,'srdel',
     .   8x,'qpar',8x,'ppar')

C --- Loop over l ---
C     rk = -1d0
C     rj = 1d0/rmax
        do  10  l = 0, lmxa
          lp1 = l+1
C     rk = rk*(2*l-1)/rmax
C     rj = rj*rmax/(2*l+1)
          nn = konfig(l)-l-1
          e = -0.5d0
          val = rmax
          slo = dnu(l)+1d0

          call rseqsr(eb1,eb2,e,tolval,z,l,nn,val,slo,v(1,i),
     .   g,sum,a,b,rofi,nr,nre,000)
ccccccccccccc
          print *,' xxx sum=',sum
ccccccccccccc
          val = val/dsqrt(sum)
          slo = slo/dsqrt(sum)
          call phdfsr(z,l,v(1,i),e,a,b,rofi,nr,g,val,slo,gp,gpp,
     .  phi,dphi,phip,dphip,p,tolval,nn)
c
          gx (1:2*nr,l,i) = g (1:2*nr)
          gpx(1:2*nr,l,i) = gp(1:2*nr)
c
          enu(l,i) = e
          dlphi = rmax*dphi/phi
          dlphip = rmax*dphip/phip
          umegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
          umegap = -(phi/phip)*(l-dlphi)/(l-dlphip)
          phplus = phi+umegap*phip
          phmins = phi+umegam*phip
          c(l,i) = e+umegam
          vl = e+umegap
          srdel(l,i) = phmins*dsqrt(0.5d0*rmax)
          q = phmins/(2*(2*l+1)*phplus)
          qpar(l,i) = 1d0/q
          ppar(l,i) = 1d0/dsqrt(p)
          if(ipr.ge.40)
     .  write(6,699) l,e,vl,c(l,i),srdel(l,i),qpar(l,i),ppar(l,i)
  699     format(i3,4f12.5,2f12.4)

c ------ make integrals of products of phi,phidot times potential ----
          fllp1 = l*(l+1)
          cc = 274.074d0
          v00 = 0d0
          v10 = 0d0
          v11 = 0d0
          d00 = 0d0
          d10 = 0d0
          d11 = 0d0

          do  15  ir = 2, nr
            r = rofi(ir)
            drdi = a*(r+b)
            wgt = 2d0*(mod(ir+1,2)+1)/3d0
            if (ir.eq.nr) wgt = 1d0/3d0
            tmc = cc-(v(ir,i)-2d0*z/r-e)/cc
            jr = ir+nr
            xxx = drdi*wgt*(v(ir,i)-2d0*z/r)
            yyy = drdi*wgt*vdif(ir,i)
            gfac = 1d0+fllp1/(tmc*r)**2
            d00 = d00+yyy*(gfac*g(ir)*g(ir)+g(jr)*g(jr))
            d10 = d10+yyy*(gfac*gp(ir)*g(ir)+gp(jr)*g(jr))
            d11 = d11+yyy*(gfac*gp(ir)*gp(ir)+gp(jr)*gp(jr))
            v00 = v00+xxx*(gfac*g(ir)*g(ir)+g(jr)*g(jr))
            v10 = v10+xxx*(gfac*gp(ir)*g(ir)+gp(jr)*g(jr))
            v11 = v11+xxx*(gfac*gp(ir)*gp(ir)+gp(jr)*gp(jr))

   15     continue
c|    write(6,290) l,v00,v10,v11
c|290 format(' l=',i4,'   vij=',3f12.6)
c|    write(6,291) l,d00,d10,d11
c|291 format(' l=',i4,'   dij=',3f12.6)
c|    write(6,292) l,v00+d00,v10+d10,v11+d11
c|292 format(' l=',i4,'   sum=',3f12.6)
          v00 = v00 + d00
          v10 = v10 + d10
          v11 = v11 + d11
          v01 = v10
          s00 = 1d0
          s10 = 0d0
          s01 = 0d0
          s11 = p
          h00 = e   + d00
          h01 = 1d0 + d10
          h10 = 0d0 + d10
          h11 = e*p + d11


C --- Integrals of u-s products from phi,phidot products ---
          det = phi*dphip - dphi*phip
          au = dphip/det
          bu = -dphi/det
          as = -phip/det
          bs = phi/det
          k  = l+1

          hab(1,k,i) = au*h00*au + au*h01*bu + bu*h10*au + bu*h11*bu
          hab(2,k,i) = au*h00*as + au*h01*bs + bu*h10*as + bu*h11*bs
          hab(3,k,i) = as*h00*au + as*h01*bu + bs*h10*au + bs*h11*bu
          hab(4,k,i) = as*h00*as + as*h01*bs + bs*h10*as + bs*h11*bs
C ... next 2 lines put in Wronskian explicitly
          hab(2,k,i) = 0.5d0*(hab(2,k,i)+hab(3,k,i)-rmax**2)
          hab(3,k,i) = hab(2,k,i)+rmax**2

          sab(1,k,i) = au*s00*au + au*s01*bu + bu*s10*au + bu*s11*bu
          sab(2,k,i) = au*s00*as + au*s01*bs + bu*s10*as + bu*s11*bs
          sab(3,k,i) = as*s00*au + as*s01*bu + bs*s10*au + bs*s11*bu
          sab(4,k,i) = as*s00*as + as*s01*bs + bs*s10*as + bs*s11*bs

          vab(1,k,i) = au*v00*au + au*v01*bu + bu*v10*au + bu*v11*bu
          vab(2,k,i) = au*v00*as + au*v01*bs + bu*v10*as + bu*v11*bs
          vab(3,k,i) = as*v00*au + as*v01*bu + bs*v10*au + bs*v11*bu
          vab(4,k,i) = as*v00*as + as*v01*bs + bs*v10*as + bs*v11*bs

   10   continue

C --- Printout ---
        if (ipr .ge. 60) then
          write(6,810)
  810     format(/'  l',10x,'val-val',5x,'val-slo',5x,'slo-val',
     .    5x,'slo-slo')
          do  60  lp1 = 1, lmxa+1
            write(6,811) lp1-1,(sab(k,lp1,i),k=1,4)
            write(6,812) (hab(k,lp1,i),k=1,4)
            write(6,813) (vab(k,lp1,i),k=1,4)
   60     continue
  811     format(i3,'   s=',4f12.5)
  812     format(3x,'   h=',4f12.5)
  813     format(3x,'   v=',4f12.5)
        endif

   80 continue
      end
