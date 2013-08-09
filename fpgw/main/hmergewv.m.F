      program mergewv
      character(16) filename
      complex(8),allocatable:: zw(:,:)
c
      call headver('hmergewv',0)
      ifwd = iopen('WV.d',1,-1,0)
      read (ifwd,*) nprecx,mrecl,nblochpmx,nw,niw, nqnum
      write(6,*) ' WV.d =', nprecx,mrecl,nblochpmx,nw,niw, nqnum
      allocate( zw(nblochpmx,nblochpmx) )
c
      ifWVRLIST =101
      ifWVILIST =102
      ifx    = 103
      ifrcw  = iopen('WVR',0,-1,mrecl)
      ifrcwi = iopen('WVI',0,-1,mrecl)
      open (ifWVRLIST, file ='WVRLIST')
      open (ifWVILIST, file ='WVILIST')
c
      do
        read(ifWVRLIST,*,end=1012) filename,iqini, iqend
        ifrcwin  = iopen(filename,0,-1,mrecl)
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
        ifrcwiin  = iopen(filename,0,-1,mrecl)
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
Cstop2rx 2013.08.09 kino      stop ' OK! mergewv WVI and WVR'
      call rx0( ' OK! mergewv WVI and WVR')
      end
