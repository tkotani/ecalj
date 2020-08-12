      subroutine wrhomt(filnam,descr,ib,rhol,rofi,nr,nlml,nsp)
C- Writes augmented charge density or potential to file for 1 sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   filnam:file name
Ci   descr :string describing rhol (informational, for output)
Ci   ib    :site index.  If ib>0, append as extension to file name
Ci   rhol  :sphere density, tabulated on a radial mesh
Ci         :rl = full charge density * r**2, written as:
Ci         :rl = sum_ilm rhol(ilm) Y_L(ilm)
Ci   rofi  :radial mesh points
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   nsp   :2 for spin-polarized case, otherwise 1
Co Outputs
Co   Radial mesh rofi and sphere density rhol are written
Co   to binary file rhoMT.ib in mc-compatible format
Cl Local variables
Cr Remarks
Cr   File is written in format mc program can read.
Cr   To read, use, e.g.
Cr     mc -r:br:open rhoMT.1 -r:br rhoMT.1 -ccat
Cr   To integrate the l=0 density for the sphere charge:
Cr     mc -r:br:open rhoMT.1 -av:nr,1 rmt -a ri -r:br rhoMT.1 -a rhoin
Cr        ri rhoin -ccat -coll 1,2 -int 0 rmt -s'sqrt(4*pi)'
Cr
Cr   A more sophisticated example: integrate rho*(v-2*z/r)
Cr     set rfile =
Cr     '-r:br:open rhoMT.1 -av:nr,1 rmt -a ri -r:br rhoMT.1 -a rin
Cr      -r:br:s=2 vtrue.1 -a v'
Cr     set v1 = 'ri v -ccat -rowl 2:nr -e1 x2-44/x1/0.282094791
Cr               -sub 0,nr,1,1 v -coll 2:nc -ccat -a v1'
Cr     set int = "ri -tog -ccat -int 0 rmt -coll 2:nc -csum"
Cr     mc $rfile rin $int
Cu Updates
Cu   11 Jun 05 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ib,nr,nlml,nsp
      double precision rofi(nr), rhol(nr,nlml,nsp)
      character filnam*(*),descr*(*)
C ... Local parameters
      integer fopng,ip,jfi,flen
      parameter (flen=60)
      character fnam*(flen),strn*(flen)

C ... get file name
      call words(filnam,ip)
      if (ip .ne. 1) call rxs('wrhomt : illegal file name: ',filnam)
      call word(filnam,1,jfi,ip)
      fnam = filnam
      if (ib .gt. 0) then
        call bin2a(' ',0,0,ib,2,0,flen,fnam,ip)
      endif
      strn = descr
      call info0(30,0,0,
     .' writing sphere '//strn//'%a to file '//fnam//'%a ...')
      if (ip .ge. flen) call rxs('wrhomt : illegal file name: ',filnam)

C ... dump results
C     Open binary, preserving case
C     call pshpr(120)
      jfi = fopng(fnam,-1,8+4)
C     call poppr
      write (jfi) nr,1,1,0
      write (jfi) rofi
      write (jfi) nr,nlml*nsp,1,0,nsp
      write (jfi) rhol
      call fclose(jfi)

      end

