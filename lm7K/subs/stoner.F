      subroutine stoner(nl,nbas,nclass,clabl,iclass,idxdn,ef,dr,npts,
     .zos,index,stni,ammax,mnpts,io,Nbar,emag,
     .mwk,amom,e)
C- Do generalised self consistent Stoner model for the d-band
C ----------------------------------------------------------------------
Ci Inputs: {zos,npts,dr} is the number-of-states function from dostet
Ci  or dosspl called in bndtb or bndasa; ammax for each class
Ci  is the upper limit on the local moment (see remarks); mnpts (see
Ci  remarks); io (see remarks); Nbar, emag, mwk (work arrays); stni
Ci  is the Stoner parameter input from CTRL for each class; index
Ci  points to the d-electron channels in zos.
Co Outputs: amom(nclass),e(nclass) are local moment and magnetic energy.
Cr Remarks: Uses the nonmagnetic DOS to determine the moment using the
Cr  generalised Stoner criterion of Andersen (Varenna notes and J. de
Cr  Physique Coll. C8-17, 1988). The scheme is essentially exact within
Cr  the rigid band approximation. The magnetic energy is calculated on
Cr  a mesh of mnpts points as a function of fixed moment up to a maximum
Cr  moment. The lowest value is then returned as the magnetic energy and
Cr  the corresponding moment as the local magnetic moment. If switch io
Cr  is set then the moment, N^bar and magnetic energy are tabulated in
Cr  the STNR file for plotting. (These quantities are defined in eqs.
Cr  (11) and (13) of the J. de Phys paper.)
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      logical io
      integer nl,nbas,nclass,iclass(nbas),idxdn(0:nl-1,0:1),npts,mnpts,
     .index(0:1)
      double precision ef,dr(2),zos(0:npts-1,0:1),stni(0:1),
     .ammax(0:1),Nbar(mnpts),emag(mnpts),
     .mwk(mnpts),amom(0:nclass-1),e(0:nclass-1)
      character*8 clabl(0:nclass-1)
C Local Variables
      integer j,ic,n1,n2,m,ifi,fopn,iprint,i1mach
      double precision mom,astep,e1,e2,ekin,I,df,d1mach

      if (nl .lt. 3) return

      if (io) ifi = fopn('STNR')
      write (*,30)
      do  3  ic = 0, nclass-1
        if (idxdn(2,ic) .ne. 1) goto 3
        astep = ammax(ic) / (mnpts - 1)
        I = stni(ic)
        if (io) call awrit2('%% rows %i cols %i',' ',80,ifi,mnpts,4)
        if (io) call awrit5('#class %i: '//clabl(ic)//' I=%d npts=%i,'//
     .  'moment=(%d,%d)',
     .  ' ',80,ifi,ic,I,mnpts,0d0,ammax(ic))
        mom =  - astep
        do  1  m = 1, mnpts
          mom = mom + astep
          mwk(m) = mom
          call getde(mom,zos(0,index(ic)),npts,dr(1),dr(2),ef,n1,n2,
     .    e1,e2,df)
          if (dabs(e2 - e1) .lt. d1mach(3)) then
            Nbar(m) = df
            emag(m) = 0
            goto 1
          endif
          if (e2 .lt. e1) call fexit(-1,111,
     .    'STONER bug: e2 < e1 at mom=%d',mom)
          Nbar(m) = mom / (e2 - e1)
          call intnb(m,mwk,Nbar,astep,ekin)
          emag(m) = ekin - I*mom**2/4
          if (ekin .eq. 0) emag(m) = 0
    1   continue
        if (io) write (ifi,10)
        if (io) write (ifi,20) (mwk(m),Nbar(m),emag(m),
     .  emag(m)+I*mwk(m)**2/4,m=1,mnpts)
C --- find global minimum in emag and get moment ---
        e(ic) = 0
        j = 1
        do  2  m = 1, mnpts
          if (emag(m) .lt. e(ic)) then
            e(ic) = emag(m)
            j = m
          endif
    2   continue
        if(j .eq. mnpts .and. iprint() .gt. 10) call awrit3
     .  ('STONER ***warning*** class %i: '//clabl(ic)//
     .  'moment(%d) = mmax(%d). Increase mmax in CTRL',' ',120,
     .  i1mach(2),ic,mwk(j),ammax(ic))
        amom(ic) = mwk(j)
        write (*,40) clabl(ic),amom(ic),e(ic),amom(ic)/Nbar(j),I
    3 continue
      if (io) call fclose(ifi)
   10 format ('#',6x,'moment',9x,'N^bar',10x,'E_mag',10x,'E_kin')
   20 format (4f15.9)
   30 format ('STONER: ',
     .'CLASS   local moment      magnetic energy     E(+)-E(-)'
     .,7x,'I')
   40 format (9x,a,3x,f10.6,10x,f10.6,7x,f10.6,4x,f6.3)
      end

      subroutine getde(mom,zos,npts,emin,emax,ef,ndown,nup,edown,eup,df)
C- Get band splitting e_down and e_up from number-of-states function
C ----------------------------------------------------------------------
Ci Inputs: moment,{zos,npts,emin,emax},ef
Co Outputs: eup, edown gives the band splitting; nup(down) is the lower
Co          pointer to the bin containing eup(down); df DOS at E_F
C ----------------------------------------------------------------------
C     implicit none
      integer npts,nup,ndown
      double precision mom,zos(0:npts-1),emin,emax,ef,edown,eup,df
      integer nf,n1,n2,n
      double precision step,qf,q0,q1,q2,q3,d1,d2,q,e1,e2

      step = (emax - emin) / (npts - 1)

C --- find bin containing Fermi level ---
      nf = (ef - emin) / step

C --- interpolate to get N(E_F) ---
      e1 = emin + nf*step
      e2 = e1 + step
      q0 = zos(nf-1)
      q1 = zos(nf)
      q2 = zos(nf+1)
      q3 = zos(nf+2)
      qf = q1 + ((ef - e1)/(e2 - e1))*(q2 - q1)

C --- get DOS(E_F) ---
      d1 = (q2 - q0) / (2d0*step)
      d2 = (q3 - q1) / (2d0*step)
      if (d2 .gt. d1) then
        df = d1 + ((ef - e1)/(e2 - e1))*(d2 - d1)
      else
        df = d2 + ((ef - e1)/(e2 - e1))*(d1 - d2)
      endif

C --- find e-limits by brute force and interpolation ---
      q = qf - mom/2d0
      do  1  n = nf, 0, -1
        if (zos(n) .lt. q) goto 2
    1 continue
      if (zos(n) .gt. 0) then
        call fexit
     .  (-1,111,'GETDE: moment m=%d too large for NOS range',mom)
      else
        call fexit
     .  (-1,111,'GETDE: moment m=%d exceeds maximum possible',mom)
      endif
    2 continue
      n1 = n
      e1 = emin + n1*step
      e2 = e1 + step
      q1 = zos(n1)
      q2 = zos(n1+1)
      edown = e1 + ((q - q1)/(q2 - q1))*(e2 - e1)
      ndown = n1

      q = qf + mom/2d0
      do  3  n = nf+1, npts-1
        if (zos(n) .gt. q) goto 4
    3 continue
      call fexit(-1,111,
     .'GETDE: moment m=%d at top of band, MMAX too large.',mom)
    4 continue
      n1 = n - 1
      n2 = n
      e1 = emin + n1*step
      e2 = emin + n2*step
      q1 = zos(n1)
      q2 = zos(n1+1)
      eup = e1 + ((q - q1)/(q2 - q1))*(e2 - e1)
      nup = n1
      end

      subroutine intnb(npts,m,Nbar,step,ekin)
C- trapezium integration of (1/2) int_0^m mdm/Nbar
C     implicit none
      integer npts,i
      double precision m(npts),Nbar(npts),step,ekin

      ekin = 0
      if (npts .lt. 4) return
      ekin = 0.5d0*(m(1)/Nbar(1) + m(npts)/Nbar(npts))
      do  1  i = 2, npts-1
        ekin = ekin + m(i)/Nbar(i)
    1 continue
      ekin = ekin*step/2d0
      end

      subroutine getkin(zos,npts,emin,emax,ndown,nup,edown,eup,ef,ekin)
C- Get magnetic band energy by partial integration of zos between limits
C ----------------------------------------------------------------------
Ci Inputs: zos,npts,emin,emax (number-of-states function; eup(down)
Ci   are the band split and nup(down) the lower mesh point pointing
Ci   to the bin containing eup(down) in zos (from GETDE); efermi
Co Outputs: magnetic energy from bands: ekin
Cr Remarks
Cr   ekin = \int_ef^eup de e n(e) - \int_edown^ef de e n(e)
Cr   done by partial integration of zos between limits
C ----------------------------------------------------------------------

C     implicit none
      integer npts,ndown,nup
      double precision zos(0:npts-1),emin,emax,edown,eup,ef,ekin
      integer nf,i
      double precision step,sum1,sum2,qf,q1,q2,q,e1,e2,A,A1,A2,qdown,
     .qup,add1,add2


      step = (emax - emin)/(npts - 1)

C --- find bin containing Fermi level ---
      nf = (ef - emin) / step

      if (nf-ndown .lt. 5 .or. nup-nf .lt. 6) then
        ekin = 0
        return
      endif

C --- contribution from bin at edown ---
      e1 = emin + ndown*step
      e2 = e1 + step
      q1 = zos(ndown)
      q2 = zos(ndown+1)
      q = q1 + ((edown - e1)/(e2 - e1))*(q2 - q1)
      add1 = 0.5d0*(e2 - edown)*(q2 + q)
      qdown = q

C --- contribution from bin at eup ---
      e1 = emin + nup*step
      e2 = e1 + step
      q1 = zos(nup)
      q2 = zos(nup+1)
      q = q1 + ((eup - e1)/(e2 - e1))*(q2 - q1)
      add2 = 0.5d0*(eup - e1)*(q1 + q)
      qup = q

C --- contributions from bin at ef ---
C A is area of trapezium at ef, A1(2) are contributions above(below) ef
      e1 = emin + nf*step
      e2 = e1 + step
      q1 = zos(nf)
      q2 = zos(nf+1)
      qf = q1 + ((ef - e1)/(e2 - e1))*(q2 - q1)
      A  = 0.5d0*(e2 - e1)*(q1 + q2)
      A1 = 0.5d0*(e2 - ef)*(q1 + qf)
      A2 = A - A1

C --- trapezium integration for remaining filled bins ---
      sum1 = 0d0
      do  1  i = ndown+2, nf-1
        sum1 = sum1 + zos(i)
    1 continue
      sum1 = sum1 + 0.5d0*(zos(ndown+1) + zos(nf))
      sum1 = sum1*step
      sum2 = 0d0
      do  2  i = nf+2, nup-2
        sum2 = sum2 + zos(i)
    2 continue
      sum2 = sum2 + 0.5d0*(zos(nf+1) + zos(nup-1))
      sum2 = sum2*step

C --- end trapezium rule; add contributions at nup, ndown and ef; and
C     1st term partial integration ---
      sum1 = sum1 + A2 + add1
      sum1 = ef*qf - edown*qdown - sum1
      sum2 = sum2 + A1 + add2
      sum2 = eup*qup - ef*qf - sum2

      ekin = (sum2 - sum1)
      end

