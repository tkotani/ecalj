      subroutine gtpcor(sspec,is,kcore,lcore,qcore)
C- Unpacks parameters related to partial core occpation
C ----------------------------------------------------------------------
Ci Inputs
Ci   sspec :struct containing species-specific information
Ci   is    :species for which to unpack kcore,lcore,qcore
Co Outputs
Co   kcore  :p.q.n for occupation
Co   lcore  :l quantum for occupation
Co   qcore  :core charge and magnetic moment
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer is,kcore,lcore
      double precision sspec(1),qcore(2)
C ... Local parameters
      character*8 ch

      kcore = 0
      lcore = -1
      qcore(1) = 0
      qcore(2) = 0
      call spacks(0,'spec coreh',sspec,ch,is,is)
      if (ch .eq. ' ') return
      call upack2('spec coreq',sspec,is,qcore)
      read (ch,'(i1)') kcore
      if (ch(2:2).eq.'s' .or. ch(2:2).eq.'S') lcore = 0
      if (ch(2:2).eq.'p' .or. ch(2:2).eq.'P') lcore = 1
      if (ch(2:2).eq.'d' .or. ch(2:2).eq.'D') lcore = 2
      if (ch(2:2).eq.'f' .or. ch(2:2).eq.'F') lcore = 3
      end

