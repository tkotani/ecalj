      logical function aiomom(alabl,pl,ql,idmod,nl,lmax,nsp,rh,vrmax,
     .ifi)
C- File I/O for moments.  IFI>0 for read, <0 for write
C ----------------------------------------------------------------
Ci Inputs
Ci   alabl,lmax,nsp
Ci   nl:   dimensions pl,ql
Ci   ifi:  logical unit: positive for read, negative for write
Ci   pl:   related to log derivatives (see atom)
Ci   ql:   moments q (see atom) (ifi > 0)
Co   rh,vrmax density and Vxc at rmax (ifi>0)
Co Outputs
Co   pl,ql:   moments q (see atom) (ifi < 0)
Co   aiomom:true unless read error or category not found
Co   rh,vrmax density and Vxc at rmax (ifi<0)
Cr Remarks
Cr   Input for second spin defaults to that for first spin if not
Cr   available for second spin.
Cr   Bug in reading vrmax for nsp=2
Cu Updates
Cu   26 Apr 03 Added MPI calls
Cu   16 Apr 03 Read P,Q using a2vec
Cu   10 Apr 02 Set 'fill' pnu to free-electron value
Cu   16 May 01 On file read, do not rewind file
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,lmax,nsp,ifi
      character alabl*8
      integer idmod(0:nl-1)
      double precision pl(nl,nsp),ql(3,nl,nsp),rh,vrmax(2)
C Local parameters
      double precision pi,dl,xx(6),scl
      integer i,l,ii,j,k,ipr,ix(6),nl0,nspf,a2vec,lmx
      integer mpipid,procid
      logical scat,pars1v,parstr,rdstrn
      character s*100
CC#ifdefC AIX
CC      logical rdstrn
CC#endif
      call getpr(ipr)
      aiomom = .false.
      pi = 4*datan(1.d0)
      procid = mpipid(1)

C --- File READ ---
      if (ifi .gt. 0) then
        if (procid .eq. 0) then
          if (.not. scat(ifi,'MOMNTS:',':',.false.)) return
C   ... read file nl and nsp
          backspace ifi
          read(ifi,'(a72)') s
          nl0 = nl
          nspf = nsp
          scl = 1d0
          if (.not. pars1v(s,len(s),'nl=','=',2,nl0))
     .    call rx('AIOMOM: failed to find nl')
          if (lmax .ne. nl0-1 .and. ipr .ge. 10)
     .    print *, 'aiomom (warning) mismatch in nl, class '//alabl
          if (.not. pars1v(s,len(s),'nsp=','=',2,nspf))
     .    call rx('AIOMOM: failed to find nsp')
          if (nsp .ne. nspf) then
            if (nsp .gt. nspf) scl = dble(nspf)/dble(nsp)
            if (ipr .ge. 10)
     .      print *, 'aiomom (warning) mismatch in nsp, class '//alabl
          endif
          j = 0
          if (.not. parstr(s,'vrmax=',72,6,'=',j,i))
     .    call rx('AIOMOM: failed to find rho,vrmax')
          if (a2vec(s,72,i,4,' ',1,j,3,ix,xx) .ne. 3) 
     .    call rx('AIOMOM: failed to find rho,vrmax')
          rh = xx(1)
          vrmax(1) = xx(2)
          vrmax(2) = xx(3)
          read(ifi,'(a72)') s
          lmx = min(nl0-1,lmax)
          call dpzero(pl,nl*nsp)
          call dpzero(ql,nl*nsp*3)
          do  10  i = 1, nspf
            do  14  l = lmx+1, lmax
              pl(l+1,i) = l+1 + .5d0 - datan(dble(l))/pi
   14       continue
            do  16  l = 0, nl0-1
              if (.not. rdstrn(ifi,s,len(s),.false.)) goto 12
              ii = 0
              if (a2vec(s,len(s),ii,4,' ',1,2,6,ix,xx) .ne. 6) goto 12
              ii = nint(xx(1))
              k = nint(xx(6))
C         read(ifi,100,err=12) ii, (xx(j), j=2,5), k
              if (ii .ne. l) call rx('AIOMOM: bad l quantum number')
              if (l .le. lmx) then
                idmod(l) = k
                if (i .le. nsp) then
                  pl(l+1,i) = xx(2)
                  pl(l+1,nsp) = xx(2)
                else
                  pl(l+1,1) = (pl(l+1,1) + xx(2))/2
                endif
                do  15  ii = 1, 3
                  if (i .le. nsp) then
                    ql(ii,l+1,nsp) = xx(2+ii)*scl
                    ql(ii,l+1,i) = xx(2+ii)*scl
                  else
                    ql(ii,l+1,1) = xx(2+ii) + ql(ii,l+1,1)
                  endif
   15           continue
              endif
   16       continue
C --- Patch for bug in AIX err= ---
#if AIX
            if (i .eq. 1 .and. nsp .eq. 2) then
              if (.not. rdstrn(ifi,s,len(s),.false.)) goto 19
              if (s(1:5) .ne. '   0 ') goto 19
              backspace ifi
            endif
#endif
   10     continue
        endif
        aiomom = .true.
        goto 99

C --- File WRITE ---
      else

        if (procid .eq. 0) then
          write(-ifi,333) alabl, rh, vrmax, lmax+1, nsp
  333     format('MOMNTS:  ',a4,
     .    '  rho,vrmax=',3f10.6,' nl=',i1,' nsp=',i1)
          write(-ifi,891)
  891     format(
     .    '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
          do  110  i = 1, nsp
          do  11  l = 0, lmax
            dl = dtan(pi*(.5d0 - pl(l+1,i)))
            if (dabs(dl) .gt. 9999) dl = 0
            write(-ifi,100) l,pl(l+1,i),(ql(ii,l+1,i),ii=1,3),idmod(l),dl
   11     continue
 110      continue
        endif
      endif
      return

  100 format(i4,4f13.7,i4,f13.7)

C --- handle read exception ---
#if AIX
   19 l = 0
      i = 2
#endif
   12 continue
      if (l .eq. 0 .and. i .eq. 2) then
        if (ipr .ge. 20)
     .  print *, 'AIOMOM, ATOM=',alabl,
     .  ':  spin 2 input missing; spin 1 moments split'
        do  400  l = 0, lmax
        do  40  ii = 1, 3
          ql(ii,l+1,1) = ql(ii,l+1,2)
   40   continue
  400   continue
      else
        aiomom = .false.
      endif

   99 continue
C      call mpibc1(pl,nl*nsp,4,.false.,'aiomom','pl')
C      call mpibc1(ql,3*nl*nsp,4,.false.,'aiomom','ql')
C      call mpibc1(rh,1,4,.false.,'aiomom','ql')
C      call mpibc1(vrmax,2,4,.false.,'aiomom','ql')
      end

