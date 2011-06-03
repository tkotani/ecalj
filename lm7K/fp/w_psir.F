      module m_unusedfid
      contains
ckino  return: unusedfid>=start if success
ckino                   =0  if failure 
      integer function unusedfid(start,end)
      implicit none
      integer, optional:: start,end

      integer:: i,istart=10,iend=999
      logical::L


      if (present(start)) istart=start
      if (present(end)) iend=end

      do i=istart,iend
        inquire(i,opened=L)
        if (.not.L) then
          unusedfid=i
          return
        endif
      enddo
      unusedfid=0

      end function unusedfid
      end module m_unusedfid


      module m_w_psir
      contains

      subroutine psir_minmax(n1,n2,n3, k1,k2,k3, 
     i   f, 
     o   emin,emax,norm,imagmax,imaxi)
      implicit none
      integer,intent(in):: n1,n2,n3,k1,k2,k3
      complex(8),intent(in):: f(k1,k2,k3)
      real(8),intent(out):: emin(2),emax(2)
      real(8),intent(out),optional:: norm,imagmax
      integer,intent(out),optional:: imaxi(3)

      integer:: i1,i2,i3
ckino find min max
          i1=1;i2=1;i3=1
          emin(1)=dreal(f(i1,i2,i3))
          emin(2)=dimag(f(i1,i2,i3))
          emax(1)=dreal(f(i1,i2,i3))
          emax(2)=dimag(f(i1,i2,i3))
          if (present(imagmax)) imagmax=abs(dimag(f(i1,i2,i3)))
          if (present(imaxi)) imaxi=(/i1,i2,i3/)
          if (present(norm)) norm=0.0d0
          do i3=1,n3
          do i2=1,n2
          do i1=1,n1
             if (emin(1)> dreal(f(i1,i2,i3))) emin(1)=dreal(f(i1,i2,i3))
             if (emax(1)< dreal(f(i1,i2,i3))) emax(1)=dreal(f(i1,i2,i3))
             if (emin(2)> dimag(f(i1,i2,i3))) emin(2)=dimag(f(i1,i2,i3))
             if (emax(2)< dimag(f(i1,i2,i3))) emax(2)=dimag(f(i1,i2,i3))
             if (present(imagmax) .and. present(imaxi) ) then
               if (imagmax < abs(dimag(f(i1,i2,i3))) ) then
                 imaxi=(/i1,i2,i3/)
                 imagmax=abs(dimag(f(i1,i2,i3)))
               endif
             endif
             if (present(norm)) norm=norm+dreal(f(i1,i2,i3))**2+dimag(f(i1,i2,i3))**2
          enddo
          enddo
          enddo
       end subroutine psir_minmax

Ckino this part is from ioden2
      subroutine w_psir_xcrysden(ifi,plat,alat,nbas,pos,z
     i  ,n1,n2,n3,k1,k2,k3,psi  )
      implicit none
      real(8),intent(in):: plat(3,3),alat,pos(3,nbas),z(nbas)
      integer,intent(in):: nbas ,ifi
      integer,intent(in):: n1,n2,n3,k1,k2,k3
      complex(8),intent(in):: psi(k1,k2,k3)
ckino ,f0(k1,k2,k3)
      integer:: isp,i,i1,i2,i3

      character(10):: psiform='(6E15.5)'
      
          write(ifi,'("CRYSTAL")')
          write(ifi,'("PRIMVEC")')
          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
     .    ,i2=1,3)
          write(ifi,'("PRIMCOORD")')
          write(ifi,'(2i5)') nbas,1
          do i = 1, nbas
            write(ifi,'(i4,2x,3f10.5)') int(z(i)),
     .      (pos(i2,i)*alat*0.529177208,i2=1,3)
          enddo

          do isp=1,3
          write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
          write(ifi,'(" wavefunction_complex_",i1)') isp
          write(ifi,'(" BEGIN_DATAGRID_3D_complex_",i1)') isp
          write(ifi,'(3i4)') n1,n2,n3
          write(ifi,'(3F10.5)') 0.,0.,0.
          write(ifi,'(3F10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
     .    ,i2=1,3)
          if (isp.eq.1) then
            do i3=1,n3
            do i2=1,n2
            write(ifi,psiform) (dble(psi(i1,i2,i3)),i1=1,n1)
            enddo
            enddo
          elseif (isp.eq.2) then
            do i3=1,n3
            do i2=1,n2
            write(ifi,psiform)  
     .      (dimag(psi(i1,i2,i3)),i1=1,n1)
            enddo
            enddo
          elseif (isp.eq.3) then
            do i3=1,n3
            do i2=1,n2
            write(ifi,psiform)
     .      (dimag(psi(i1,i2,i3))**2+dreal(psi(i1,i2,i3))**2,i1=1,n1)
            enddo
            enddo
          endif
          write(ifi,'(" END_DATAGRID_3D")')
          write(ifi,'("END_BLOCK_DATAGRID_3D")')
          enddo



      end subroutine w_psir_xcrysden

Ckino from rsibl2
      subroutine w_psir(ng,nspc,nev,psi,n1,n2,n3,k1,k2,k3,kv,isp,q,iq
     i  ,n_eiglist, eiglist
     i  ,plat,alat,nbas,pos,z
     w ,f)
      use m_unusedfid
C- FT wave function to real space and add square into mesh density
C  and optionally make smpot * psi
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng    :number of G-vectors
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   nev   :number of wave functions
Ci   psi   :wave function in reciprocal space
Ci   n1..3 :size of FT mesh
Ci   k1..3 :dimensions smpot,smrho
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   isp  : spin index
Ci   n_eiglist: number of  eiglist
Ci   eiglist: eig list to print psir
Ci   plat: lattice  vector
Ci   alat: unit of plat 
Ci   nbas: # of atoms
Ci   pos: postion of atoms
Ci   z:   atomic number of atoms
Co Work
Co   f     :psi in real space
Cr Remarks
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer,intent(in):: k1,k2,k3,n1,n2,n3,ng,nspc,nev,kv(ng,3)
      integer,intent(in):: isp  ! spin ?

      complex(8),intent(in):: psi(ng,nspc,nev)
      complex(8)::f(k1,k2,k3)
ckino      integer,intent(in):: n_qlist
      integer,intent(in)::  n_eiglist
ckino      real(8),intent(in):: qlist(3,n_qlist)
      integer,intent(in):: eiglist(n_eiglist)
      integer,intent(in):: nbas
      real(8),intent(in):: plat(3,3),alat,pos(3,nbas),z(nbas),q(3)
ckino      complex(8),intent(in):: f0(k1,k2,k3,nspc,nev)

C ... Local parameters
      character(6):: thisfunc='w_psir' 
      integer i1,i2,i3,i,ispc,ieig,iq1,iq
      logical:: l_execute
      character(40):: fnamer,fnamei
      integer:: ifiler,ifilei,stat
      real(8):: emin(2),emax(2),norm,imagmax
      complex(8):: zphase
      integer:: imaxi(3)

      logical:: l_rot=.false.
      real(8),parameter:: zero=0.0d0

      write(*,*) thisfunc,':start '
 

      if (q(1).eq.zero .and. q(2).eq.zero .and. q(3).eq.zero ) then
         l_rot=.true.
         write(*,*)thisfunc,'rotation is done for gamma point'
      endif

#if 0
      write(*,'(a,9F15.5)')'plat=',plat
      write(*,'(a,F15.5)')'alat=',alat
      write(*,'(a,i5)')'nbas=',nbas
      write(*,*)'pos='
      write(*,'(20F15.5)')pos
      write(*,*)'z='
      write(*,'(20F15.5)')z
#endif

ckino get a free file descriptor
      ifiler=unusedfid(100) 
      if (ifiler.eq.0) then
          write(*,*) 'Error: failed to get a file descriptor in ',thisfunc,
     .      ' but continue.'
          return
      endif


      do  ispc = 1, nspc
        do ieig=1,n_eiglist
          i=eiglist(ieig)

          write(fnamer,'(a,I1,a,I1,a,I4.4,a,I4.4,a)') 'psi_s',isp,'_',ispc,'_q',iq,'_e',i,'.xsf'
          write(*,*)thisfunc,':open ',trim(fnamer),' for psi'
          open(ifiler,file=fnamer,status='unknown',iostat=stat,action='write')
          if (stat.ne.0) then
             write(*,*)thisfunc,':Error: failed to open ',fnamer,'. in ',thisfunc,
     .         ' but continue.'
             return
          endif

          call gvputf(ng,1,kv,k1,k2,k3,psi(1,ispc,i),f)
          call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,1)

ckino find min max
          call psir_minmax(n1,n2,n3, k1,k2,k3,
     i   f, 
     o   emin,emax,norm,imagmax,imaxi)
          write(*,'(a,a,4i4)') thisfunc,':ispc,ispc,iq,ie=',isp,ispc,iq,i
          write(*,'(a,a,3F15.7)') thisfunc,':real part, min,max and norm=',emin(1),emax(1),norm
          write(*,'(a,a,3F15.7)')thisfunc, ':imag part, min,max,|max|   =',emin(2),emax(2),imagmax

         zphase=f( imaxi(1),imaxi(2),imaxi(3) )
         if ( abs(dimag(zphase)) > abs(dreal(zphase))*1d-3 ) then
            zphase=dconjg(zphase)
            zphase=zphase / abs(zphase) 
ckino            write(*,'(a,2F15.5,a,2F15.5)') 'imagmax f=',
ckino     .   f( imaxi(1),imaxi(2),imaxi(3) ), 
ckino     .   'rot f=',f( imaxi(1),imaxi(2),imaxi(3) )*zphase
            do i3=1,n3
            do i2=1,n2
            do i1=1,n1
               f(i1,i2,i3)=f(i1,i2,i3)*zphase
            enddo
            enddo
            enddo

          call psir_minmax(n1,n2,n3, k1,k2,k3,
     i   f, 
     o   emin,emax,norm,imagmax,imaxi)
          write(*,'(a,a,4i4)') thisfunc,':isp,ispc,iq,ie=',isp,ispc,iq,i
          write(*,'(a,a,3F15.7)') thisfunc,':rot real part, min,max and norm=',emin(1),emax(1),norm
          write(*,'(a,a,3F15.7)') thisfunc,':rot imag part, min,max,|max|   =',emin(2),emax(2),imagmax

         endif

         write(ifiler,'(a,4i5,3F15.5)') '#isp,ispc,iq,ie,q=',isp,ispc,iq,i,q
         call w_psir_xcrysden(ifiler,plat,alat,nbas,pos,z
     i  ,n1,n2,n3,k1,k2,k3, f 
     i  )


         close(ifiler)

        enddo
      enddo
      write(*,*) 'end ', thisfunc
ckino      stop 'debug end: kino'
      end subroutine w_psir

      subroutine w_psir_cmdopt(sopts,m_listwf,n_listwf,listwf)
      implicit none
      character*(*),intent(in):: sopts
      integer,intent(in):: m_listwf
      integer,intent(out):: n_listwf,listwf(m_listwf)

      character(12)::thisfunc='w_psi_cmdopt'

      character:: dc
      integer:: j1,j2

      n_listwf=0
      j2=0
        j2=j2+1
        if (sopts(j2:j2).eq.dc) return
        j1=min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2.ge.j1) then
c           write(*,*)thisfunc,':j1,j2,sopts=',j1,j2,sopts(j1:j1+5)
           if (sopts(j1:j1+5).eq.':band=' ) then
c              write(*,*)thisfunc,':sopts=',sopts(j1+6:j2),m_listwf
              call mkils0(sopts(j1+6:j2),n_listwf,listwf(1))
              if (n_listwf>m_listwf) then
c                 write(*,*)thisfunc,': parse error, opt=',trim(sopts)
                 n_listwf=0
                 return
              endif
              call mkilst(sopts(j1+6:j2),n_listwf,listwf)
           endif
        endif

      end subroutine w_psir_cmdopt

      end module m_w_psir
