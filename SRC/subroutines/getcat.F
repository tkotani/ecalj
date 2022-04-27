      module m_getcat
      public scat
      private
      contains
      logical function scat(ifi,categ,term,lrewnd)
C- Scan file for a category
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit
Ci   categ :seek category 'categ'
Ci   term  :string that terminates category
Ci   lrewnd:true, rewind file before searching; false, do not
Co Outputs
Co   scat  :true if category found, false if not
Co         :file unit is placed after record for which category was found
Cr Remarks
Cr   This is a file version of subroutine getcat
C ----------------------------------------------------------------
      integer ifi
      logical lrewnd
      character categ,term
      character*72 a
C Local variables
      integer catl0,recl0
      parameter (catl0=7,recl0=72)
      integer i
      logical rdstrn,lsequ
      external rdstrn,lsequ
      if (lrewnd) rewind ifi
      scat = .false.
   10 if (.not. rdstrn(ifi,a,recl0,.false.)) return
      if (.not. lsequ(categ,a,catl0,term,i)) goto 10
      scat = .true.
      end
      end module
