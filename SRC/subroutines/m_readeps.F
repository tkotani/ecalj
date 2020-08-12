      module m_readeps
!! Subroutines: read_eps,deallocate_eps
      implicit none

      complex(8),allocatable,protected :: epsinv(:,:),w_mu(:,:),llmat(:,:)
c      complex(8),allocatable,protected :: epsinv(:,:),w_mu(:,:)
c      complex(8),protected :: llmat(3,3)

      contains

      subroutine check_eps()
!! read info for inverse dielectric function. 
      integer:: iq,iq_read,ngb_read
      real(8):: q_read(3),tolq=1d-8
      integer::ifw4p,ifw4pchk,ngbq0,ix
      character(5)::charnum5
      logical::ext
      open(newunit=ifw4pchk,file='checkW4PHONON',status='replace')
      iq=0
      do 
         iq=iq+1
         inquire(file='W4PHONON.'//charnum5(iq),exist=ext)
         if(ext) then
            open(newunit=ifw4p,file='W4PHONON.'//charnum5(iq),form='unformatted',status='old')
            read(ifw4p) iq_read,q_read,ngb_read,ix !ix=0, or ix=1 for q=0 (iq=1)
            close(ifw4p)
            write(ifw4pchk,*) iq_read,q_read,ngb_read,ix !ix=0, or ix=1 for q=0 (iq=1)
         else
            exit
         endif
      end do
      close(ifw4pchk)
      end subroutine check_eps

      subroutine read_eps(q,iqx,ngb)
      intent(in)::        q,iqx,ngb
      integer:: iqx,ngb
      real(8):: q(3)
      integer:: iq,iq_read,ngb_read
      real(8):: q_read(3),tolq=1d-8
      integer::ifile_handle
      integer::ifw4p,ngbq0,ix
      character(5)::charnum5
      logical::ext
      iq=0
      do
         iq=iq+1
         inquire(file='W4PHONON.'//charnum5(iq),exist=ext)
         if(ext) then
            open(newunit=ifw4p,file='W4PHONON.'//charnum5(iq),form='unformatted',status='old')
            read(ifw4p) iq_read,q_read,ngb_read,ix !ix=0, or ix=1 for q=0 (iq=1)
c            if(iqx /=iq_read) call rx('m_readeps: iqx/=iqf')
            if(sum(abs(q-q_read)) <tolq) then
               if(ngb /=ngb_read) call rx('m_readeps: ngb/=ngb_read.')
c               write(*,*) 'm_readeps: ngb,ngb_read',ngb,ngb_read
               allocate(epsinv(ngb,ngb))
c               write(*,*) 'm_readeps: epsinv size',size(epsinv,1),size(epsinv,2)
               if(ix==1) then
                  epsinv(1,1:ngb)=0d0
                  epsinv(2:ngb,1)=0d0
               endif
               read(ifw4p) epsinv(ix+1:ngb,ix+1:ngb)
               exit
            endif
            close(ifw4p)
         else
            write(*,*) 'm_readeps err1: iq,q',iq,q
            call rx('m_readeps: W4PHONON'//charnum5(iq)//' not found for given q.')
         endif
      enddo
      inquire(file='W4PHONON.HeadWing',exist=ext)
      if(ext) then
         open(newunit=ifw4p,file='W4PHONON.HeadWing',form='unformatted',status='old')
         allocate(llmat(3,3))
         read(ifw4p) llmat(1:3,1:3),ngbq0 !for q~0
         write(*,*) 'm_readeps: llmat r1', llmat(1,:) 
         write(*,*) 'm_readeps: llmat r2', llmat(2,:)
         write(*,*) 'm_readeps: llmat r3', llmat(3,:)
         
         if(iqx ==1) then
            if(sum(abs(q))>tolq) call rx('m_readeps: q/=0 for iqx=1.')
            if(ngb /=ngbq0) call rx('m_readeps: ngb/=ngbq0.')
            allocate(w_mu(2:ngbq0,3))
            read(ifw4p) w_mu(2:ngbq0,3) !for q~0
            write(*,*) 'meps: w_mu 2,1', w_mu(2,1)
            write(*,*) 'meps: w_mu 2,2', w_mu(2,2)
            write(*,*) 'meps: w_mu 3,3', w_mu(3,3)
            write(*,*) 'meps: w_mu 4,1', w_mu(4,1)
            write(*,*) 'meps: w_mu ngb-1,1', w_mu(ngb-1,1)
            write(*,*) 'meps: w_mu ngb,1', w_mu(ngb,1)
            write(*,*) 'meps: w_mu ngb,2', w_mu(ngb,2)
            write(*,*) 'meps: w_mu ngb,3', w_mu(ngb,3)
         endif

         close(ifw4p)
      else
         call rx('m_readeps: W4PHONON.HeadWing not found.')
      endif
      end subroutine read_eps

      subroutine deallocate_eps(iqx)
      implicit none
      integer,intent(in):: iqx
      deallocate(epsinv)
      deallocate(llmat)
      if(iqx==1) deallocate(w_mu)
      end subroutine deallocate_eps

      end module
