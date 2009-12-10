      subroutine shftpq(nc,nrclas,nsp,nl,lmx,rmax,avw,pp,amom,
     .                  idmod,swmod,pmin,pnu,qnu,sumevm)
C- Shift p and q according to idmod
C ----------------------------------------------------------------
Ci Inputs
Ci   nc    :number of classes
Ci   nrclas:nrclas(i) = number of atoms in the ith class
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :(global maximum l) + 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   pp    :potential parameters (atomsr.f)
Ci   amom  :input values of moments qnu
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
Ci   swmod :see Remarks
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   nc,nsp,lmx,nl,rmax,pp,amom,idmod,swmod
Co Outputs
Co   pnu,qnu:pnu possibly floated, and qnu computed from amom but
Co          :1st and second energy moments shifted according to
Co          :shifted pnu
Co   sumevm :
Cv Verbosity
Cv   >=30: print out moments and shift in enu
Cr Remarks
Cr   when swmod=.true. prompt for each idmod=2 for shift in
Cr   amom and qnu can point to the same address space
Cr   command argument -enu=val sets each idmod=2 and each enu to val
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      logical swmod
      integer nl,nsp,nc,idmod(0:nl-1,nc),lmx(nc),nrclas(nc)
      double precision avw,sumevm
      double precision amom(3,0:nl-1,nsp,nc),pp(6,nl,nsp,nc),rmax(nc)
      double precision pmin(0:*),pnu(0:nl-1,nsp,nc),qnu(3,0:nl-1,nsp,nc)
C ... Local parameters
      logical cmdopt,a2bin
      double precision eb,pold,pminl
      integer ic,isp,l,i,j,iprint,lgunit,idm
      external iprint,lgunit
      character*72 outs

      if (iprint() .ge. 32) then
        do  10  j = 1, 2
   10   write(lgunit(j),334)
  334   format(/' CLASS L    Q0',9x,'Q1',9x,'Q2',9x,'EB',
     .    9x,'POLD',7x,'PNU')
      endif
      do  13  ic = 1, nc
        do  13  isp = 1, nsp
        do  13  l = 0, lmx(ic)
        if (iprint() .ge. 32) then
          do  11  j = 1, 2
   11     write(lgunit(j),331) ic,l,(amom(i,l,isp,ic), i=1,3)
        endif
        eb = 0
        idm = idmod(l,ic)
        if (cmdopt('-enu=',5,0,outs)) then
          j = 5
          if (.not. a2bin(outs,eb,4,0,' ',j,72))
     .      call rx('shftpq: bad value -enu')
          eb = eb - pp(1,l+1,isp,ic)
          idm = 2
          idmod(l,ic) = 2 + (idmod(l,ic)-mod(idmod(l,ic),10))
        endif
        if (swmod .and. idm .eq. 2) call query('eb=',4,eb)
        pold = pnu(l,isp,ic)
        pminl = pmin(l)
        call enutod(l,rmax(ic),avw,pp(1,l+1,isp,ic),amom(1,l,isp,ic),
     .    0d0,idm,pminl,pnu(l,isp,ic),qnu(1,l,isp,ic),eb)
C .. I think ...
        if (eb.ne.0 .and. idm.eq.2)
     .    pp(1,l+1,isp,ic) = pp(1,l+1,isp,ic) + eb
        if (iprint() .ge. 32 .and. eb .ne. 0) then
          do  12  j = 1, 2
            if (pminl .eq. pnu(l,isp,ic)) then
              write(lgunit(j),332) (qnu(i,l,isp,ic), i=1,3), eb,
     .                              pold, pnu(l,isp,ic), ' *'
            else
              write(lgunit(j),332) (qnu(i,l,isp,ic), i=1,3), eb,
     .                              pold, pnu(l,isp,ic)
            endif
   12     continue
         endif
   13 continue
      sumevm = 0d0
      do  1  ic = 1, nc
      do  1  isp = 1, nsp
      do  1  l = 0, lmx(ic)
        sumevm = sumevm + (qnu(2,l,isp,ic) +
     .                     qnu(1,l,isp,ic)*pp(1,l+1,isp,ic))*nrclas(ic)
    1 continue
  331 format(2i4,f10.6,5f11.6)
  332 format(8x,f10.6,5f11.6:a)
      end

