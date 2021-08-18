      subroutine wrhomt(filnam,descr,ib,rhol,rofi,nr,nlml,nsp)
C- Writes augmented charge density or potential to file for 1 sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   filnam:file name
Ci   descr :string describing rhol (only for information)
Ci   ib    :site index
Ci   rhol  :sphere density, tabulated on a radial mesh
Ci         :rl = full charge density * r**2, written as:
Ci         :rl = sum_ilm rhol(ilm) Y_L(ilm)
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co   Radial mesh rofi and sphere density rhol are written to filnam.ib
      integer ib,nr,nlml,nsp
      double precision rofi(nr), rhol(nr,nlml,nsp)
      character :: descr*(*),filnam*(*)
      integer   jfi,iprint
      character*8 xtxx
      character(256):: fnam
      fnam=trim(filnam)//xtxx(ib)
      if(iprint()>30) write(6,"(a)")' Writing Sphere '//trim(descr)//' to '//trim(fnam)
      open(newunit=jfi,file=trim(fnam),form='unformatted')
      write (jfi) nr,1,1,0
      write (jfi) rofi
      write (jfi) nr,nlml*nsp,1,0,nsp
      write (jfi) rhol
      close(jfi)
      end
      
      character*8 function xtxx(num) !taken from xt in extension.F
      integer(4) :: num
      xtxx=''
      if(num>0)     xtxx = char(48+mod(num,10))
      if(num>9)     xtxx = char(48+mod(num/10,10))//xtxx
      if(num>99)    xtxx = char(48+mod(num/100,10))//xtxx
      if(num>999)   xtxx = char(48+mod(num/1000,10))//xtxx
      if(num>9999)  xtxx = char(48+mod(num/10000,10))//xtxx
      if(num>99999) xtxx = char(48+mod(num/100000,10))//xtxx
      if(num>999999) call rx( ' xtxx:can not produce')
      end

