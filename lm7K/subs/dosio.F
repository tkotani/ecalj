      subroutine dosio(dos,nemx,ndmx,ndos,nld,emin,emax,nspin,
     .eferm,del,ifmt,ifile)
C- I/O for dos, MSM's format
C     implicit none
C     nemx: leading dimension of dos array (must be >ndos)
C     ndos: number of energy mesh points
C     ndmx: max number of dos allowed for input
C     nld:  number dos for i/o for each spin channel
C Passed parameters
      integer ifile,ifmt,ndmx,ndos,nemx,nld,nspin
      double precision del,eferm,emax,emin
      double precision dos(nemx,1)
C Local parameters
      integer ie,ild,iedum

C --- Write branch ---
      if (ifile .lt. 0) then
        write(-ifile,760) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ifmt .eq. 0) then
          do  10  ie = 1, ndos
            write(-ifile,761) ie,(dos(ie,ild),ild=1,nld*nspin)
   10     continue
        else
          do  11  ild = 1, nld*nspin
            write(-ifile,762) (dos(ie,ild),ie=1,ndos)
   11     continue
        endif
  761   format(i5,6f12.5/(5x,6f12.5))
  760   format(2f10.5,3i5,2f10.5,i5)
  762   format(5f14.6)
      endif
C --- Read branch ---
      if (ifile .gt. 0) then
        read(ifile,760) emin,emax,ndos,nld,nspin,eferm,del,ifmt
        if (ndos .gt. nemx) call rx('dosio: ndos.gt.nemx')
        if (nld*nspin .gt. ndmx) call rx('dosio: nld.gt.ndmx')
        if (ifmt .eq. 0) then
          do  20  ie = 1, ndos
            read(ifile,761) iedum,(dos(ie,ild),ild=1,nld*nspin)
   20     continue
        elseif (ifmt .eq. 1) then
          do  21  ild = 1, nld*nspin
            read(ifile,762) (dos(ie,ild),ie=1,ndos)
   21     continue
        else
          call rx('dosio: bad fmt')
        endif
      endif
      end

