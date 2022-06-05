subroutine atwf(mode,a,lmxa,nr,nsp,pnu,pnz,rsml,ehl,rmt,z,v0, &
     nphimx,ncore,konfig,ecore,gcore,gval,nmcore)
  use m_ftox
  !- Make properties related to core for one sphere
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 return ncore, and konfig, and nphimx only;
  !i         :  see description below for contents of nphimx
  !i         :1s digit
  !i         :1 return valence wave functions
  !i         :2 return core wave functions
  !i         :3 combination of 1+2
  !i         :10s digit concerns orthogonalization
  !i         :0 do not orthogonalize
  !i         :1 return orthogonalized to valence orbitals
  !i         :2 return orthogonalized to valence orbitals
  !i         :  using large component only
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i   lmxa  :augmentation l-cutoff
  !i   nr    :number of radial mesh points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
  !i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
  !i   pnz   :pnu por local orbitals
  !i   rmt   :MT boundary
  !i   z     :nuclear charge      (not used if mode=0)
  !i   v0    :spherical potential (not used if mode=0)
  !i   ehl   :energy of smoothed Hankel tail for extended local orbital
  !i   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
  ! o Inputs/Outputs
  ! o  nphimx:dimensions gval.  Must be at least as large as the
  ! o        :number of valence wave functions
  ! o        :For mode=0, nphimx is output and is assigned to
  !i         :maximum number radial wave functions for any l channel.
  !o Outputs
  !o   ncore :number of core levels
  !o   konfig:1s digit contains core configuration
  !o         :10s digit:
  !o         : 0 -> no local orbitals
  !o         : 1 -> local orbital with p.q.n. < pnu
  !o         : 2 -> local orbital with p.q.n. > pnu
  !o   ... The following are not used if mode=0
  !o   ecore :core eigenvalues
  !o   gcore :core wave functions
  !o   gval  :valence wave functions
  !o          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
  !o            i=0 : phi
  !o            i=1 : phidot
  !o            i=2 : local orbital
  !r Remarks
  !u Updates
  !u    4 Sep 04 Adapted to extended local orbitals
  !u   22 Dec 01 Adjustments to accomodate changes in phidx
  !u   22 Apr 01 Created by MvS
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,nr,nsp,lmxa,ncore,konfig(1+lmxa),n0,nrmx,nphimx
  parameter (n0=10,nrmx=1501)
  double precision :: rmt,z,a,v0(nr,nsp),pnu(n0,nsp),pnz(n0,nsp), &
       gval(nr*2,0:lmxa,nphimx,nsp),ecore(*),gcore(nr,2,*), &
       rsml(n0),ehl(n0)
  ! ... Local parameters
  logical :: lpz
  integer :: l,isp,konf,konfz,k,mode0,mode1,  nmcore
  double precision :: sumtc,sumec,e,ez,xx
  !     double precision hcrl,val(5),slo(5),pi,tol
  !     parameter (tol=1d-12)
  double precision :: rofi(nrmx),rwgt(nrmx),rhoc(nrmx,2),gp(2*nrmx*4)
  double precision :: phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp

  logical:: isanrg, l_dummy_isanrg

  mode0 = mod(mode,10)
  mode1 = mod(mode/10,10)

  ! --- Count number of core states ---
  lpz = .false.
  ncore = 0
  do  l = 0, lmxa
     k = l+1
     konfig(k) = pnu(k,1)
     konfz = mod(pnz(k,1),10d0)
     if (konfz == 0) konfz = konfig(k)
     l_dummy_isanrg=isanrg(konfz,konfig(k)-1,konfig(k)+1,'atwf:','pnuz',.true.)
     !       lpz = konfz .ne. konfig(k)
     do  konf = l+1, min(konfz,konfig(k))-1
        ncore = ncore+nsp
     enddo
     if (konfz < konfig(k)) then
        konfig(k) = konfz + 10
        lpz = .true.
     elseif (konfz > konfig(k)) then
        konfig(k) = konfig(k) + 20
        lpz = .true.
     endif
  enddo

  if (mode0 == 0) then
     nphimx = 2
     if (lpz) nphimx = 3
     return
  endif

  if (nr > nrmx) call rx('increase nrmx in atwf')
  call radmsh(rmt,a,nr,rofi)
  call radwgt(rmt,a,nr,rwgt)

  ! --- Valence wave functions ---
  if (mod(mode0,2) == 1) then
     do  l = 0, lmxa
        k = l+1
        do  isp = 1, nsp
           konf = pnu(k,1)

           !    ...  Make phi and phidot
           !         NB: Write gdot to gp, with extra space for higher derivatives
           !         nn  = konf-l-1
           !         pi = 4d0*datan(1d0)
           !         hcrl = 0
           !         val(1) = rofi(nr)
           !         slo(1) = 1 + dtan(pi*(0.5d0 - pnu(k,isp)))
           !         call phidx(0,z,l,v0(1,isp),hcrl,0d0,rofi,nr,4,tol,e,val,slo,
           !    .      nn,gval(1,l,1,isp),gp,xx,xx,xx,xx,pgam,xx,xx,xx,xx)
           call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnu(1,isp),4, &
                gval(1,l,1,isp),gp,e,phi,dphi,phip,dphip,p)
           ! cccccccccccccccccccc
           !            print *,' eeeeeeeeee=',l,isp,e
           ! cccccccccccccccccccc
           !         Copy 1st derivative to passed array
           call dcopy(2*nr,gp,1,gval(1,l,2,isp),1)
           !         phi,phidot already orthogonal if mode1=1
           if (mode1 == 2) &
                call ortrwf(10*(mode1-1)+2,z,l,v0(1,isp),nr,nr,nr,rofi,rwgt, &
                e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp),xx)

           !     ... Make local orbital
           if (konf /= konfig(k)) then
              ! ino isanrg is logical function,               call isanrg(nphimx,3,3,'atwf:','nphimx',.true.)
              l_dummy_isanrg=isanrg(nphimx,3,3,'atwf:','nphimx',.true.)
              call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnz(1,isp),2, &
                   gval(1,l,3,isp),gp,ez,phz,dphz,phzp,dphzp,p)

              ! ino isanrg is logical function,               call isanrg(mode1,0,2,'atwf:','10s digit mode',.true.)
              l_dummy_isanrg=isanrg(mode1,0,2,'atwf:','10s digit mode',.true.)
              if (mode1 == 0) then

                 !             Extra scaling
                 !              call ortrwf(0,z,l,v0(1,isp),nr,nr,nr,rofi,rwgt,e,e,ez,
                 !     .          gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp),xx)
                 !             call dscal(nr*2,1/xx,gval(1,l,3,isp),1)
                 !             phz = phz/xx
                 !             dphz = dphz/xx

                 call wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz, &
                      phzp,dphzp,pnz(1,isp),rsml,ehl, &
                      gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp))
              elseif (pnz(l+1,isp) < 10) then
                 call ortrwf(10*(mode1-1)+1,z,l,v0(1,isp),nr,nr,nr,rofi, &
                      rwgt,e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp), &
                      gval(1,l,3,isp),xx)
              endif
              !           call prrmsh('gz',rofi,gval(1,l,3,isp),nr,nr,2)
           endif

        enddo
     enddo

     !       call prrmsh('gval',rofi,gval,nr,nr,2*(1+lmxa))

  endif

  ! --- Core eigenfunctions and eigenvalues ---
  if (mode0 >= 2) then
     call getcor(1,z,a,pnu,pnz,nr,lmxa,rofi,v0,0,0,0d0,sumec,sumtc, &
          rhoc,ncore,ecore,gcore,nmcore) !nmcore jun2012
  endif

end subroutine atwf

subroutine atwf2l(ifi,jfi,iclass,a,lmxa,nr,nsp,pnz,rmt, &
     nphimx,konfig,ecore,gcore,gval)
  use m_lmfinit,only: stdo
  use m_ftox
  !- Translate radial wave functions from shifted to standard log mesh
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ifi   :file logical unit for valence w.f.
  !i   jfi   :file logical unit for core w.f.
  !i   iclass:for printout
  !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  !i         :See Remarks
  !i   lmxa  :augmentation l-cutoff
  !i   nr    :number of radial mesh points
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   pnz   :p.q.n for local orbitals
  !i   rmt   :MT boundary
  !i   nphimx:dimensions gval.  Should be 2, or 3 if any local orbitals
  !i   konfig:1s digit contains core configuration
  !i         :10s digit:
  !i         : 0 -> no local orbitals
  !i         : 1 -> local orbital with p.q.n. < pnu
  !i         : 2 -> local orbital with p.q.n. > pnu
  !i   ecore :core eigenvalues
  !i   gcore :core wave functions
  !i   gval  :valence wave functions
  !i          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
  !i            i=0 : phi
  !i            i=1 : phidot
  !i            i=2 : local orbital
  !r Remarks
  !r   This routine translates w.f. to a standard log mesh and writes
  !r   them in a spex-readable format.
  !r
  !r   Debugging:
  !r   set sqr = '-p -p -xe -tog -coll 1 -tog -coll 2:nc -ccat'
  !r   mc -qr out.gas $sqr -int 0 2.361911
  !u Updates
  !u   15 Jul 09 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ifi,jfi,iclass,nr,nsp,lmxa,konfig(0:lmxa),n0,nphimx
  parameter (n0=10)
  double precision :: rmt,a,pnz(n0),gval(nr,2,0:lmxa,nphimx,nsp), &
       ecore(*),gcore(nr,2,*)
  ! ... Local parameters
  integer :: ic,icore,iphi,iprint,ir,isp,jcore,jx,konf,l,nphi, &
       npoly,nrmx,PRTV
  double precision :: xxo,xx,fac,tolg,dot3
  double precision :: x,y,dy,dytop,dytopa  ! For interpolation
  parameter (npoly=6,nrmx=1501,PRTV=60)
  parameter (tolg=1d-8)
  real(8):: gwk(nrmx),gwk2(nrmx,2),rofi(nrmx),rwgt(nrmx), rwgtl(nrmx),rofil(nrmx), &
       onorm(2,0:lmxa,nphimx,nsp), norm(2,0:lmxa,nphimx,nsp)  ! norm on log mesh

  ! ... Setup
  ! angenglob      stdo = nglob('stdo')
  !      stdo = globalvariables%stdo
  if (nr > nrmx) call rx('increase nrmx in atwf')

  ! ... Shifted and standard log mesh
  call radmsh(rmt,a,nr,rofi)
  call radwgt(rmt,a,nr,rwgt)
  xx = rmt
  fac = exp(-a)
  do  ir = nr, 1, -1
     rofil(ir) = xx
     rwgtl(ir) = 2d0*a*xx/3d0   ! Simpson's rule for standard log mesh
     !       rwgtl(ir) = 1d0*a*xx
     xx = xx*fac
  enddo
  !     Complete Simpson's rule
  do  ir = 2, nr-1, 2
     rwgtl(ir) = 2d0*rwgtl(ir)
  enddo
  rwgtl(1) = rwgtl(1)/2
  rwgtl(nr) = rwgtl(nr)/2

  !     call prrmsh('standard log mesh',rofi,rofil,nr,nr,1)

  !      print *, '!! testing ... int x dx',
  !     .  nr,sngl(a),sngl(rmt)
  !      sum1 = 0
  !      sum2 = 0
  !      do  ir  = 1, nr
  !        sum1 = sum1 + rofi(ir)*rwgt(ir)
  !        sum2 = sum2 + rofil(ir)*rwgtl(ir)
  !      enddo
  !      print 333, rmt**2/2,sum1-rmt**2/2,sum2-rmt**2/2
  !  333 format('exact',f15.10,'  shifted mesh error',1pe10.2,
  !     .  '  standard mesh error',1pe10.2)
  !      stop

  ! --- Valence wave functions ---
  !     ic = 2*(lmxa+1)*nphimx*nsp
  !     call prrmsh('gval on shifted log mesh',rofi,gval,nr,nr,ic)

  write(ifi) iclass,lmxa,nsp ! 1st record in class: dimensioning
  if(iprint()>=PRTV) write(stdo,ftox)' Valence wave function normalizations for ' &
       //'class ',iclass,' phi  l  spin  sqrt(norm)',isp
  dytopa = 0
  do  isp = 1, nsp
     do  l = 0, lmxa
        nphi = 2
        if (pnz(l+1) > 0) nphi = 3
        do  iphi = 1, nphi
           write(ifi) l,iphi,isp   ! indices to wave function being written
           do  ic = 1, 2  ! large, small components
              onorm(ic,l,iphi,isp) = &
                   dot3(nr,gval(1,ic,l,iphi,isp),gval(1,ic,l,iphi,isp),rwgt)

              !         g(r) = phi(r)*r
              !         For small r, g(r)~r^(l+1) ... so fit g(r)/r^(l+1)
              do  ir = 2, nr
                 gwk(ir) = gval(ir,ic,l,iphi,isp)/rofi(ir)**(l+1)
              enddo
              jx = 0
              dytop = 0
              do  ir = 1, nr
                 x = rofil(ir)
                 call polint(rofi(2),gwk(2),nr-1,npoly,x,tolg,0,jx,y,dy)
                 dytop = max(dytop,abs(dy))
                 dytopa = max(dytopa,abs(dy))
                 gwk2(ir,ic) = y*rofil(ir)**(l+1)
                 !           gval(ir,ic,l,iphi,isp) = gwk2(ir,ic) ! Overwrite, for debuggging
              enddo
              norm(ic,l,iphi,isp) = dot3(nr,gwk2(1,ic),gwk2(1,ic),rwgtl)
              !          Debugging printout
              !          call info8(10,0,0,'interp ic=%i l=%i iphi=%i isp=%i:'//
              !     .      ' err ~ %;g  onorm = %;6d  nnorm = %;6d',
              !     .      ic,l,iphi,isp,dytop,onorm(ic,l,iphi,isp),
              !     .      norm(ic,l,iphi,isp),0)

           enddo
           write(ifi) gwk2(1:nr,1),gwk2(1:nr,2) ! Large component, followed by small component
           if (iprint() >= PRTV) write(stdo,345) iphi, l, isp, &
                dsqrt(norm(1,l,iphi,isp)+norm(2,l,iphi,isp))
345        format(2i4,i6,f15.8,f15.6)
        enddo
     enddo
  enddo
  write(ifi) -1,-1,-1 ! Flags last record for this class
  !     Debugging; need to comment overwrite above
  !      ic = 2*(lmxa+1)*nphimx*nsp
  !      call prrmsh('gval on standard log mesh',rofil,gval,nr,nr,ic)

  ! --- Core eigenfunctions and eigenvalues ---
  write(jfi) iclass,lmxa,nsp ! 1st record in class: dimensioning
  dytopa = 0
  if(iprint()>=PRTV) write(stdo,ftox)' Core wave function normalizations for ' &
       //'class ',iclass,'   n   l   spin    sqrt(norm)        ecore',isp 
  icore = 0
  do  l = 0, lmxa
     do  isp = 1, nsp
        jcore = 0
        do  konf = l+1, mod(konfig(l),10)-1
           icore = icore+1
           jcore = jcore+1
           xxo = dot3(nr,gcore(1,1,jcore),gcore(1,1,jcore),rwgt) + &
                 dot3(nr,gcore(1,2,jcore),gcore(1,2,jcore),rwgt)
           dytop = 0d0
           do  ic = 1, 2 ! large, small components
              !           g(r) = phi(r)*r
              !           For small r, g(r)~r^(l+1) ... so fit g(r)/r^(l+1)
              do  ir = 2, nr
                 gwk(ir) = gcore(ir,ic,jcore)/rofi(ir)**(l+1)
              enddo
              jx = 0
              do  ir = 1, nr
                 x = rofil(ir)
                 call polint(rofi(2),gwk(2),nr-1,npoly,x,tolg,0,jx,y,dy)
                 if (ic == 1) then
                    dy = dy*rofil(ir)**(l+1)
                    dytop  = max(dytop,dy)
                    dytopa = max(dytopa,dy)
                 endif
                 gwk2(ir,ic) = y*rofil(ir)**(l+1)! gcore(ir,ic,jcore) = gwk2(ir,ic) ! Overwrite, for debuggging
              enddo
           enddo
           xx = dot3(nr,gwk2(1,1),gwk2(1,1),rwgtl) + dot3(nr,gwk2(1,2),gwk2(1,2),rwgtl)
           if(iprint()>99)write(stdo,ftox)' atwf2l interp l=',l,' konf=',konf,'isp=',isp,  &
                'err~',ftod(dytop),'onorm=',ftof(xxo),'nnorm-onorm=',xx-xxo
           write(jfi) jcore, l, isp, konf, ecore(icore)
           write(jfi) gwk2(1:nr,1), gwk2(1:nr,2)
           if (iprint()>=PRTV) write(stdo,345) konf, l, isp, dsqrt(xx), ecore(icore)
        enddo
     enddo
  enddo
  write(jfi) -1,-1,-1,-1,-1d0 ! Flags last record for this class
end subroutine atwf2l

