      subroutine iorbtm(sspec,ics,nl,nlo,nclass,nsp,orbtm)

      use m_struc_def  !Cgetarg

C- Printout of orbital moments
C ----------------------------------------------------------------------
Ci Inputs
Ci   sspec :struct containing species-specific information
Ci     Elts read: name
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   nl    :(global maximum l) + 1
Ci   nlo   :number of l or lm
Ci   nclass:number of inequivalent classes
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   orbtm :orbital moments
Co Outputs
Cu Updates
Cu   09 Aug 04 (A. Chantis) Correct sign of orbl
Cu   08 Dec 00 First implementation
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlo,nl,nsp,nclass,ics(*)
      real(8):: orbtm(nlo,nsp,nclass)
      type(s_spec)::sspec(*)

C ... Local parameters
      integer ic,isp,l,im,lm,m,l1,ipr,stdo,lgunit,is
      double precision amom,orbl(10)
      character*8 slabl
C ... External calls
ckino      external dpzero,getpr,spacks
      external dpzero,getpr

      call getpr(ipr)
      if (ipr .lt. 20) return
      stdo = lgunit(1)
      write(stdo,332)
  332 format(/' IORBTM:  orbital moments :'/
     .' Class  Spec        spin   Moment decomposed by l ...')

C     totmom = 0
      do  ic = 1, nclass
        is = ics(ic)
c        do i_spacks=is,is
c          call spacks_copy('u',sspec(i_spacks)%name,is,is,slabl,i_spacks)
c        enddo
        slabl=sspec(is)%name

        amom = 0
        do  isp = 1, nsp
          call dpzero(orbl,nl)
          lm = 0
          do  l = 0, nl-1
            l1 = l+1
            im = l
            if (nl .eq. nlo) im = 0
            do  m = -im, im
              lm = lm+1
C              print *, l,m,isp,ic,orbtm(lm,isp,ic)
              orbl(l1) = orbl(l1) + orbtm(lm,isp,ic)
              amom = amom + orbtm(lm,isp,ic)
            enddo
          enddo
          write(stdo,333) ic,slabl,isp,(orbl(l1),l1=1,nl)
  333     format(i5,4x,a8,i6,8f12.6)
        enddo
        write(stdo,334) ic, amom
  334   format(' total moment for class',i4,':',f12.6)
      enddo

      end subroutine iorbtm


