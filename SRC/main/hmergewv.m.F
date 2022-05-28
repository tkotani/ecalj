      program mergewv
      implicit integer (i-n)
      character(16) filename
      complex(8),allocatable:: zw(:,:)
c
c      call headver('hmergewv',0)
      open(newunit=ifwd,file='WV.d')
      read (ifwd,*) nprecx,mrecl,nblochpmx,nw,niw, nqnum
      close(ifwd)
      write(6,*) ' WV.d =', nprecx,mrecl,nblochpmx,nw,niw, nqnum
      allocate( zw(nblochpmx,nblochpmx) )
c
      ifWVRLIST =ifile_handle()
      ifWVILIST =ifile_handle()
      ifx    = ifile_handle()
c      ifrcw  = iopen('WVR',0,-1,mrecl)
c      ifrcwi = iopen('WVI',0,-1,mrecl)
      ifrcw  = ifile_handle()
      ifrcwi  = ifile_handle()
      open(ifrcw, file='WVR',form='unformatted',access='direct',recl=mrecl)
      open(ifrcwi,file='WVI',form='unformatted',access='direct',recl=mrecl)

      open (ifWVRLIST, file ='WVRLIST')
      open (ifWVILIST, file ='WVILIST')
c
      do
        read(ifWVRLIST,*,end=1012) filename,iqini, iqend
c        ifrcwin  = iopen(filename,0,-1,mrecl)
        ifrcwin  = ifile_handle()
        open(ifrcwin, file=trim(filename),form='unformatted',access='direct',recl=mrecl)
        write(6,*) filename,' is merging into WVR ...'
        do iq = iqini,iqend
          do iw  = 1,nw
            read (ifrcwin , rec=((iq-iqini)*nw +iw) ) zw
            write(ifrcw   , rec=((iq-2)*nw +iw)     ) zw
          enddo
        enddo
        close(ifrcwin,status='delete')
      enddo
 1012 continue
      do
        read(ifWVILIST,*,end=1013) filename, iqini, iqend
c        ifrcwiin  = iopen(filename,0,-1,mrecl)
        ifrcwiin  = ifile_handle()
        open(ifrcwiin,file=trim(filename),form='unformatted',access='direct',recl=mrecl)
        write(6,*) filename,' is merging into WVI ...'
        do iq  =  iqini,iqend
          do iw  = 1,niw
            read (ifrcwiin, rec=((iq-iqini)*niw +iw) ) zw
            write(ifrcwi  , rec=((iq-2)*niw +iw)     ) zw
          enddo
        enddo
        close(ifrcwiin,status='delete')
      enddo
 1013 continue
      call rx0( ' OK! mergewv WVI and WVR')
      end
