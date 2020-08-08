c calculate mean Tc
c Use QIBZ and jlist -------
      integer(4),parameter:: nx=100,nxx=1000
      integer(4):: nw,nwx,ifx=0
      real(8):: q(3),jzero,  freq(nx), dw=0.2,rrr,iii,polinta !meV
      real(8):: rchipmm(nx),ichipmm(nx),omg,rydberg,mmom,jval,www
     & ,drchi(nx),freqm(nx),czero,omgw,tc,omgwmean,qibz(3,nxx),wibz(nxx)
      character(30) fname,aaa
      character(5) ax
      real(8) ::bolz=11640.4

c--- Get mmom
c--- Get mmom
      ifchipm=1015
      open(ifchipm, file='QIBZ')
      iw=0
      omgwmean=0d0
      www=0d0
      read(ifchipm,*)
      do 
        read(ifchipm,*,end=2999) qibz(:,iw+1), wibz(iw+1)
c        write(6,"(3f19.15,2x,d23.15,f15.10,f16.8)")
c     &    q, jval, mmom, omgw
        www=www+ wibz(iw+1)
        iw=iw+1
      enddo
 2999 continue
      nwx=iw
      write(6,"('readin number of iqbz; www=',i5,d23.15)")nwx,www
      close(ifchipm)

      ifchipm=1015
      open(ifchipm, file='jlist')
      iw=0
      omgwmean=0d0
      do 
        read(ifchipm,*,end=1999) q, jval, mmom, omgw
c        write(6,"(3f19.15,2x,d23.15,f15.10,f16.8)")
c     &    q, jval, mmom, omgw
        do iwx=1,nwx
          if(sum(abs(q-qibz(:,iwx)))<1d-6) goto 888
        enddo
        do iwx=1,nwx
          write(6,"(3f20.15,3x 3f20.15,2x,3d13.5)") 
     &      q,qibz(:,iwx),sum(abs(q-qibz(:,iwx))) 
        enddo
        stop 'We can not find q=qibz '
 888    continue
        omgwmean = omgwmean + omgw * wibz(iwx)
        iw=iw+1
      enddo
 1999 continue
      nw = iw
      omgwmean = omgwmean/www
      tc= 2d0/3d0 * mmom/4d0* omgwmean
      write(6,"(' nw omgmean mean Tc=',i4,2f13.5)")
     &  nw,omgwmean, tc/1000d0*bolz
      end
