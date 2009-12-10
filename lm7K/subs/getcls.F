      subroutine getcls(recrd,clabl,nclass,ifreqd,i)
C- return index to a class
C ----------------------------------------------------------------
Ci Inputs
Ci   clabl,nclass
Ci   ifreqd: true if class is required
Ci   search starts at common variable iend
Co Outputs
Co   i: index to class (-1 if none found)
Cr Remarks
Cr
C ----------------------------------------------------------------
C     implicit none
      character*8 clabl(0:*)
      character*1 recrd(0:*)
      integer i,nclass
      logical ifreqd

C local variables
      integer i0,partok,iprint
      character*8 alabl
C --- Common blocks ---
      integer recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .        iend,ichoos,nchoos,optio
      logical noerr

      common /iolib/ recoff,reclen,nrecs,maxlen,catbeg,catsiz,subsiz,
     .               iend,ichoos,nchoos,noerr,optio

      if (optio .eq. 2) alabl = clabl(i)
      subsiz = catsiz

      i = -1
      i0 = partok(recrd(catbeg),'ATOM=','=',8,alabl,1,1,iend,ifreqd)
      if (.not. noerr) return
      call tokmat(alabl,clabl,nclass,8,' ',i,i0,ifreqd)
      if (i .eq. -1 .and. iprint() .ge. 30)
     .  print *, 'getcls: (warning) class ',alabl,' not in class table'
      call subcat(recrd(catbeg),'ATOM=','=',iend)

      end

