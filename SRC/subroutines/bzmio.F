      subroutine bzmio(bzm,ndim,npln,nwmx,nqmx,nw,e,nq,vx,vy,x,y,z,
     .nqp,nsp,metal,n,w,wtbzm,efermi,ifile)
C- I/O for BZ maps
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   bzm : integrated DOS for each q-point, DOS channel,
Ci         energy window, and BZ plane
Ci   ndim : number of DOS channels
Ci   npln : total number of BZ planes
Ci   nwmx : maximum number of energy windows
Ci   nqmx : maximum number of k-points in one BZ plane
Ci   nw : number of energy windows for each BZ plane
Ci   e(1), e(2) : lower and upper limits of energy windows
Ci   nq,vx,vy,x,y,z : see subroutine getbzp
Ci   nqp : total number of q-points
Ci   nsp,metal,n,w,efermi
Ci   wtbzm : weight per q-point
Ci   ifile : file handle for BZMP, pos. for read, neg. for write
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim,npln,nwmx,nqmx,nqp,nsp,n,ifile
      integer nw(npln),nq(2,npln)
      double precision w,wtbzm,efermi
      double precision bzm(nqmx,ndim,nwmx,npln),e(2,nwmx,npln),
     .vx(3,npln),vy(3,npln),x(2,npln),y(2,npln),z(npln)
      logical metal

C Local parameters
      integer ip,nq0,iw,i,ich,iq,ndim1,npln1,nwmx1,nqmx1

      if (ifile .gt. 0) goto 100

C --- Write branch ---
      write (-ifile,500) ndim,npln,nwmx,nqmx,nqp,nsp
      write (-ifile,510) metal,n,w,wtbzm,efermi
  500 format(6i6)
  510 format(l3,i3,3f18.10)

C --- Loop over BZ planes, energy windows, and DOS channels ---
      do  20  ip = 1, npln
        write (-ifile,520) (nq(i,ip),i=1,2),nw(ip),(vx(i,ip),i=1,3),
     .  (x(i,ip),i=1,2),(vy(i,ip),i=1,3),(y(i,ip),i=1,2),z(ip)
        nq0 = nq(1,ip)*nq(2,ip)

        do  10  iw = 1, nw(ip)
          write (-ifile,530) (e(i,iw,ip),i=1,2)
          write (-ifile,540) ((bzm(iq,ich,iw,ip),iq=1,nq0),ich=1,ndim)
   10   continue

   20 continue

  520 format(3i5,2(/3f18.10/2f18.10)/f18.10)
  530 format(2f18.10)
  540 format(4f18.10)

      return

C --- Read branch ---
  100 continue
      read (ifile,500) ndim1,npln1,nwmx1,nqmx1,nqp,nsp
      call rxx(ndim1 .ne. ndim,'BZMIO: ndim1 ne ndim')
      call rxx(npln1 .ne. npln,'BZMIO: npln1 ne npln')
      call rxx(nwmx1 .ne. nwmx,'BZMIO: nwmx1 ne nwmx')
      call rxx(nqmx1 .ne. nqmx,'BZMIO: nqmx1 ne nqmx')
      read (ifile,510) metal,n,w,wtbzm,efermi

C --- Loop over BZ planes, energy windows, and DOS channels ---
      do  120  ip = 1, npln
        read (ifile,520) (nq(i,ip),i=1,2),nw(ip),(vx(i,ip),i=1,3),
     .  (x(i,ip),i=1,2),(vy(i,ip),i=1,3),(y(i,ip),i=1,2),z(ip)
        nq0 = nq(1,ip)*nq(2,ip)

        do  110  iw = 1, nw(ip)
          read (ifile,530) (e(i,iw,ip),i=1,2)
          read (ifile,540) ((bzm(iq,ich,iw,ip),iq=1,nq0),ich=1,ndim)
  110   continue

  120 continue

      end

