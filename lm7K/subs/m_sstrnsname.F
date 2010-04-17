C this is original but not need to care ordering in out scheme
C
Cr --- explanation of  ustrn ---
Cr   off offe  name    purpose
Cr         1   amix    amix mode
Cr         2   gemb    file names for embedding GF
Cr         3   gfopt   GF options
Cr         4   jobid   job id
Cr         5   map     mapping
Cr         6   mix     mixing mode
Cr         7   mmham   micromagnetics hamiltonian
Cr         8   sxopt   SX options
Cr         9   symg    symmetry group
C ----------------------------------------------------------------
      module m_sstrnsname
      integer,parameter:: nsname=8,nx=2000
      integer:: init(nx),iend(nx)
      character*20 ::sname(10)
      integer,save :: initialized=0
      contains
      subroutine initialize()
      integer::ix
      iend=-1
      sname='xxxxxxxxxxxxxxxxxxxx'
      sname(1)='strn mix'   
      sname(2)='strn symg'  
      sname(3)='strn jobid'
      sname(4)='strn map'
      sname(5)='strn amix'
      sname(6)='strn mmham' 
      sname(7)='strn sxopt' 
      sname(8)='strn gfopt'
      do ix=1,nsname
        init(ix)=nx*(ix-1)+1
      enddo
      end subroutine
      end module
