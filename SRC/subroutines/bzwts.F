      subroutine bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,
     .metal,tetra,norder,npts,width,rnge,wtkp,eb,efmax,efermi,sumev,
     .wtkb,dosef,qval,ent,lfill)
      use m_lmfinit,only: globalvariables
C- BZ integration for fermi level, band sum and qp weights
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbmx  :leading dimension of eb
Ci   nevx  :leading dimension of wtkb and max number of evals calculated
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   n1..n3:number of divisions for the k-point mesh
Ci   nkp   :number of inequivalent k-points (bzmesh.f)
Ci   ntet  :number of inequivalent tetrahedra (tetirr.f)
Ci   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
Ci         :corners of tetrahedron;
Ci         :idtet(0,i) number of tetrahedra of the i'th kind
Ci   zval  :valence charge
Ci   metal :T => metal, F => nonmetal
Ci   tetra :T => tetrahedron integration
Ci   norder:(sampling) polynomial order in Methfessel-Paxton integration
Ci         :100s digit norder flags that metals treatment should apply
Ci         :     regardless of whether a gap is present or not
Ci   width :(sampling) gaussian width in Methfessel-Paxton integration
Ci   npts  :(sampling) number of points in DOS mesh
Ci   rnge  :(sampling) range over which sampling delta function is assumed
Ci         :to vanish, in units of width
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci   eb    :energy bands; alias eband
Ci   efmax :largest eigenvalue for which eigenvectors are known
Ci   eb    :energy bands
Co Outputs
Co   efermi:Fermi energy
Co   sumev :sum of eigenvalues
Co   wtkb  :integration weights (not generated for nonmetal case)
Co   dosef :DOS at Fermi level
Co   qval  :qval(1) = total charge; qval(2) = magnetic moment
Co   ent   :entropy term (actually TS)
Cl Local variables
Cl   lfill :true => insulator
Cu Updates
Cu   12 Jul 08 (ATP) bzwts now returns entropy term (actually kTS)
Cu   04 Jun 08 (ATP) Handles metal case when nkp=1
Cu    4 Aug 07 bzwts can make and tabulate specific heat (F-D statistics)
Cu   29 Jul 07 (ATP) Find E_F using weights by bisection, not INTNOS
Cu   02 Jan 06 return qval (valence charge and moment)
Cu   17 Jan 05 Use 100s digit norder as flag to treat all cases as metal,
Cu             whether or not a gap is present
Cu    1 May 04 When insulator, write gap to log file
Cu   09 May 04 When insulator, write gap to log file
Cu   01 Jul 03 When insulator, prints highest occ and lowest unocc state
Cu   24 Oct 02 Patch for weird cases when idos doesn't encompass
Cu             zval, where emin, emax found by efrang.
Cu   22 Sep 01 Returns dosef now.  Altered argument list.
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      logical metal,tetra
      integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,
     .idtet(5,ntet)
      double precision zval,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),
     .wtkb(nevx,nsp,nkp),efmax,efermi,sumev,dosef(2),qval(2),ent
C Local variables
      integer:: it , itmax , n , nptdos , nspx , nbmxx , nevxx , ib 
     ., ikp , ipr , job , nbpw , i1mach , nev , stdo , stdl , nglob 
     ., mkdlst , ifi , fopna , i , j , lry , procid , mpipid , master 
     ., nulli , isw
      real(8) ,allocatable :: dos_rv(:)
      real(8) ,allocatable :: bot_rv(:)
      real(8) ,allocatable :: top_rv(:)
      integer ,allocatable :: bmap_iv(:)
      real(8) ,allocatable :: wk_rv(:)

      double precision emin,emax,e1,e2,dum,tol,e,elo,ehi,sumwt,
     .dmin,dmax,dval,egap,amom,dsum,cv,tRy
      character outs*100
      logical cmdopt,efrng2,lfill
      real(8) ,allocatable :: tlst_rv(:)

      parameter (nulli=-99999)
      integer w(1)

C --- Locate band limits: lfill => insulator ---
      procid = mpipid(1)
      master = 0
      call tcn('bzwts')
      call getpr(ipr)
      qval(1) = 0
      qval(2) = 0
      ent = 0
      n = isign(1,norder) * mod(iabs(norder),100)
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
Changenglob      stdl = nglob('stdl')
      stdl = globalvariables%stdl
#ifdef KINODEBUG
      write(*,*) 'kino nevx,nsp',nevx,nsp
#endif
      allocate(bot_rv(nevx*nsp))

      allocate(top_rv(nevx*nsp))
#ifdef KINODEBUG
      write(*,*) 'kino allocate top_rv',nevx*nsp
#endif

      nspx  = 1
      nevxx = nevx*nsp
      nbmxx = nbmx*nsp
C     job = 1 for non spin pol, -1 for spin pol
      job = 3-2*nsp
      dosef(1) = 0
      dosef(nspx) = 0
      egap = nulli
      if (nsp .eq. 2 .and. nspc .eq. 1 .and. cmdopt('--oldbz',7,0,outs))
     .then
        nspx  = nsp
        nevxx = nevx
        nbmxx = nbmx
        job = 1
      endif

C     Force coupled spins: find range
      if (nspx .ne. nsp .and. nspc .eq. 1) then
        nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
        allocate(bmap_iv(nevx*nsp*nkp/nbpw+1))
        bmap_iv(:)=0

        allocate(wk_rv(nevx*nsp))
#ifdef KINODEBUG
        write(*,*) 'kino allocate bmap_iv, wk_rv',nevx*nsp*nkp/nbpw+1,nevx*nsp
#endif

        call ebcpl ( 0 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv 
     .  , wk_rv , eb )

        lfill = efrng2 ( nspx , nkp , nbmxx , nevxx , zval * 2 , eb , 
     .  bot_rv , top_rv , elo , ehi , emin , emax )

        call ebcpl ( 1 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv 
     .  , wk_rv , eb )

C     Spins not coupled: find range
      else
        lfill = efrng2 ( nspx , nkp , nbmxx , nevxx , nspc * zval , eb 
     .  , bot_rv , top_rv , elo , ehi , emin , emax )

      endif
C ... Bands never filled if 100s digit norder set
      if (.not. tetra .and. iabs(norder) .ge. 100) lfill = .false.
      if (allocated(wk_rv)) deallocate(wk_rv)
      if (allocated(bmap_iv)) deallocate(bmap_iv)
      if (allocated(top_rv)) deallocate(top_rv)
      if (allocated(bot_rv)) deallocate(bot_rv)

C ... Case an insulator: put efermi at emin + tiny number
      if (lfill) then
        efermi = emin + 1d-10
C ... Do the best we can should be a metal, but assumption that it isn't
      elseif (.not. metal) then
        efermi = (emin + emax) / 2
      endif

C ... Pretend as though spin-pol bands are coupled to find E_f
      if (nsp .eq. 2 .and. nspc .eq. 1 .and. job .eq. -1) then
        nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
        allocate(bmap_iv(nevx*nsp*nkp/nbpw+1))
        bmap_iv(:)=0

        allocate(wk_rv(nevx*nsp))
#ifdef KINODEBUG
        write(*,*) 'kino allocate bmap_iv, wk_rv=',nevx*nsp*nkp/nbpw+1, nevx*nsp
#endif

        call ebcpl ( 0 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv 
     .  , wk_rv , eb )

      endif

C --- BZ weights, sumev and E_f for an insulator  ---
      if ( .not. metal ) then
        if (.not. lfill .and. ipr .gt. 10) then
          print *, ' BZWTS : partly filled bands encountered; '
          print *, ' expect serious errors from Fermi cut-off. '
          print *, ' **** re-start with METAL=T in control file **** '
        endif
        sumwt = 0d0
        sumev = 0d0
        nev = nint(zval/2)
        call rxx(nev .gt. nbmx,'BZWTS: zval too big')
        do  ikp = 1, nkp
          do  ib = 1, nev*nsp
            e = eb(ib,1,ikp)
            if (e .le. efermi) then
              sumwt = sumwt +   abs(wtkp(ikp))/nsp
              sumev = sumev + e*abs(wtkp(ikp))/nsp
            endif
          enddo
        enddo
        egap = emax-emin
        if (ipr .ge. 20) then
          write(stdo,'(/'' BZWTS : --- Non-metal sampling ---'')')
          call awrit3(' Fermi energy: %,6;6d;  %d electrons;  Sum '//
     .    'occ. bands: %,6;6d',' ',80,stdo,efermi,sumwt,sumev)
          call info5(20,0,0,' VBmax = %,6;6d  CBmin = %,6;6d  gap = '//
     .    '%,6;6d Ry = %,5;5d eV',emin,emax,emax-emin,
     .    (emax-emin)*13.6d0,0)
        endif

C --- BZ weights, sumev and E_f by tetrahedron method (Blochl wts) ---
      else if (tetra) then
        if (ipr .ge. 30) write (stdo,103)
  103   format(/' BZWTS : --- Tetrahedron Integration ---')
        if (lfill) then
          egap = emax-emin
          call info2(10,0,0,' ... only filled or empty bands'//
     .    ' encountered:  ev=%,6;6d  ec=%,6;6d',emin,emax)
          call info5(20,0,0,' VBmax = %,6;6d  CBmin = %,6;6d  gap = '//
     .    '%,6;6d Ry = %,5;5d eV',emin,emax,emax-emin,
     .    (emax-emin)*13.6058d0,0)
          goto 2
        endif
        nptdos = 101
        allocate(dos_rv(nspx*nptdos))
#ifdef KINODEBUG
        write(*,*)'kino allocate dos_rv',nspx*nptdos
#endif

        tol = 1d-6
C       Preliminary check that dos lies within emin,emax
C       Widen emin,emax if not
        if (.not. lfill) then
          call bzints ( n1 , n2 , n3 , eb , dum , nkp , nevxx , nbmxx , 
     .    nspx , emin , emax , dos_rv , nptdos , efermi , job , ntet 
     .    , idtet , sumev , qval )

          dmin = dval ( dos_rv , 1 )

          if ( nspx .eq. 2 ) dmin = dmin + dval ( dos_rv , nptdos + 
     .    1 )

          dmax = dval ( dos_rv , nptdos )

          if ( nspx .eq. 2 ) dmax = dmax + dval ( dos_rv , nptdos + 
     .    nptdos )

          if (dmin .gt. zval) then
            emin = 3*emin-2*emax
            call info5(10,0,0,
     .      ' (warning): initial NOS ( %1,6;6g %1,6;6g ) does'//
     .      ' not encompass Q=%1;6g' ,dmin,dmax,zval,0,0)
          elseif (dmax .lt. zval) then
            emax = 3*emax-2*emin
            call info5(10,0,0,
     .      ' (warning): initial NOS ( %1,6;6g %1,6;6g ) does'//
     .      ' not encompass Q=%1;6g' ,dmin,dmax,zval,0,0)
          endif
        endif
        if (ipr .ge. 35) print 101
  101   format(9x,'Est E_f ',10x,'Window',8x,'Tolerance',2x,'n(E_f)')

        itmax = 5
        do   it = 1, itmax
          call bzints ( n1 , n2 , n3 , eb , dum , nkp , nevxx , nbmxx , 
     .    nspx , emin , emax , dos_rv , nptdos , efermi , job , ntet 
     .    , idtet , sumev , qval )

          call fermi ( zval , dos_rv , nptdos , emin , emax , nspx , 
     .    efermi , emin , emax , dosef )

          if (ipr .ge. 35)
     .    write(stdo,100) efermi,emin,emax,emax-emin,dosef(1)
  100     format(7x,6(f10.6,1x))
          if (emax-emin .lt. tol) goto 1
        enddo
        if (ipr .gt. 10)
     .  call awrit2(' BZWTS (warning): Fermi energy not converged: '//
     .  '%g > tol=%g',' ',80,i1mach(2),emax-emin,tol)
    1   continue
        if (allocated(dos_rv)) deallocate(dos_rv)

    2   continue
        call bzints(n1,n2,n3,eb,wtkb,nkp,nevxx,nbmxx,
     .  nspx,emin,emin,emin,1,efermi,2*job,ntet,idtet,sumev,qval)
      else

C --- BZ weights, sumev and E_f by Methfessel-Paxton sampling ---
C       if (ipr .gt. 30) write (*,102) n,width,nint(rnge)
C       if (ipr .gt. 30) write (fopn('LOG'),102) n,width
C  102   format(/' BZWTS : --- Brillouin Zone sampling; N =',i2,
C     .    ', W =',f6.3:'  range =',i2,'*W ---')
        call info2(31,1,0,' BZWTS : --- Brillouin Zone sampling; '//
     .  'N=%i, W=%10z%;4g%0z ---',n,width)
C   ... Temporarily remove spin degeneracy if spins are coupled
        if (nsp .eq. 2 .and. nspx .eq. 1) call dscal(nkp,.5d0,wtkp,1)
C   ... Find Fermi level, sampling
        if ((.not. lfill) .or. (metal .and. (nkp .eq. 1))) then
          e1 = elo - rnge*width/2
          e2 = ehi + rnge*width/2
          efermi = 0.5d0*(e1 + e2)
          itmax = 1000
          do  it = 1, itmax
            call pshprt(0)
            call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efermi,
     .      .true.,sumev,wtkb,qval,ent,dosef,cv)
            call poppr
            if (dabs(zval - qval(1)) .lt. 1d-12) then
              call info5(41,0,0,
     .        ' Fermi energy, %d found after %i bisections,'//
     .        ' %d electrons, DOS(E_f)=%;4g',efermi,it,qval(1),dosef,0)
              goto 3
            endif
            if (qval(1) .gt. zval) then
              e2 = efermi
            else
              e1 = efermi
            endif
            efermi = 0.5d0*(e1 + e2)
          enddo
          call info0(10,0,0,' BZWTS (warning): '//
     .    'cannot find E_F by bisection, using INTNOS')
          allocate(dos_rv(nspx*npts))
#ifdef KINODEBUG
          write(*,*)'kino allocate dos_rv',nspx*npts
#endif

          emin = elo - rnge*width/2
          emax = emax + rnge*width/2
          call maknos ( nkp , nevxx , nbmxx , nspx , wtkp , eb , n , width 
     .    , - rnge , emin , emax , npts , dos_rv )

          if ( nspx.eq.2 ) call dpsadd ( dos_rv , dos_rv , npts , 
     .    1 , npts + 1 , 1d0 )

          call intnos ( npts , dos_rv , emin , emax , zval , efermi 
     .    , dosef , sumev )

          if (allocated(dos_rv)) deallocate(dos_rv)

    3     continue
        else
          dosef(1) = 0
          dosef(2) = 0
          egap = emax-emin
          call info2(10,0,0,' ... only filled or empty bands'//
     .    ' encountered:  ev=%,6;6d  ec=%,6;6d',emin,emax)
          call info5(20,0,0,' VBmax = %,6;6d  CBmin = %,6;6d  gap = '//
     .    '%,6;6d Ry = %,5;5d eV',emin,emax,emax-emin,
     .    (emax-emin)*13.6d0,0)
        endif

C   ... (optional) Tabulate specific heat in file for list of T's
        if ((cmdopt('--cv:',5,0,outs) .or. cmdopt('--cvK:',6,0,outs))
     .  .and. n .lt. 0 .and. metal) then
          if (procid .eq. master) then
            if (cmdopt('--cvK:',6,0,outs)) then
              lRy = 0
              i = 7
            else
              lRy = 1
              i = 6
            endif
            itmax = mkdlst(outs(i:),1d-8,-1,w)
            if (itmax .gt. 0) then
              allocate(tlst_rv(itmax))

              call word(outs,1,it,j)
              call info2(20,0,0,' Writing CV(T) to file for %i vals of T: '
     .        //outs(i:j)//' %?#(n==1)#(Ry)#(K)#',itmax,lRy)
              it = mkdlst ( outs ( i: ) , 1d - 8 , itmax + 1 , tlst_rv )

              if (it .ne. itmax) call rx('bzwts: bug in mkdlst')
              ifi = fopna('cv',-1,0)
              rewind ifi
              call awrit1('%% rows %i cols 4%N'//
     .        '#    T(K)    T(Ry)   S(k_B)   TdS/dT(k_B)',' ',80,ifi,it)
              do  it = 1, itmax
                try = dval ( tlst_rv , it )

                if (lRy .eq. 0) then
                  tRy = tRy/0.1579d6
                endif
                call pshpr(1)
                call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,tRy,efermi,
     .          metal,sumev,wtkb,qval,ent,dosef,cv)
                call poppr
                call awrit4('%;9,1D%;9,5D  %;9F  %;9F',' ',80,ifi,
     .          0.1579d6*tRy,tRy,ent,cv)
              enddo
              call fclose(ifi)
            endif
          endif
        endif

C   ... Make weights, sampling
        call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efermi,
     .  (.not. lfill) .or. (metal .and. (nkp .eq. 1)),
     .  sumev,wtkb,qval,ent,dosef,cv)

C   ... Put back spin degeneracy if removed
        if (nsp .eq. 2 .and. nspx .eq. 1) call dscal(nkp,2d0,wtkp,1)

      endif

C ... Restore to uncoupled bands; ditto with weights
      if (nsp .eq. 2 .and. nspc .eq. 1 .and. job .eq. -1) then
        call ebcpl ( 1 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv 
     .  , wk_rv , eb )

        if ( metal ) call ebcpl ( 1 , nevx , nevx , nsp , nspc , nkp 
     .  , nbpw , bmap_iv , wk_rv , wtkb )

        if (allocated(tlst_rv)) deallocate(tlst_rv)
        if (allocated(wk_rv)) deallocate(wk_rv)
        if (allocated(bmap_iv)) deallocate(bmap_iv)

      endif

C ... Magnetic moment
      amom = 0
      if (nsp .eq. 2 .and. nspc .ne. 2 .and. metal) then
        do  ikp = 1, nkp
          amom = amom + dsum(nevx,wtkb(1,1,ikp),1) -
     .    dsum(nevx,wtkb(1,2,ikp),1)
        enddo
C        if (ipr .gt. 0) write(stdo,922) amom
C  922   format(9x,'Mag. moment:',f15.6)
C        if (ipr .gt. 0) then
C          call awrit1(' BZWTS: magnetic moment = %d',' ',128,i1mach(2),
C     .                 amom)
C        endif
      endif
      qval(2) = amom

C ... log file
      if (ipr .gt. 0) then
        outs = ' '
        call awrit5(
     .  'bz met %l  tet %l  ef %,6;6d  sev %,6;6d  zval %,1;6d',
     .  outs,len(outs),0,metal,tetra,efermi,sumev,zval)
        call awrit5('%a qval %,1;6d%?#n# amom %,1;6d#%j#'//
     .  '%?#n# gap %,4;4d eV##',outs,len(outs),-stdl,
     .  qval,int(amom*10000),amom,isw(egap.gt.0),egap*13.6d0)
      endif

C ... Check that efmax is above efermi ...
      e = efermi
      if (.not. lfill .and. .not. tetra) e = efermi + rnge*width/2
      if (e .gt. efmax) call fexit2(-1,111,' Exit -1 BZWTS: '//
     .'efmax=%1;6d < %1;6d',efmax,e)

      call tcx('bzwts')

      end

      subroutine ebcpl(mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap,wk,eb)
C- Gather spin-polarized bands into a single group, or redistribute
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0, gather; 1, scatter
Ci   nbmx  :leading dimension of b and dimensions wk
Ci   nevx  :number of eigenvalues
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nq    :number of k-points for which eb are calculated
Ci   nbpw  :a number no larger than the number of bits per integer word
Ci   bmap  :an integer work array with dimension at nbmx*nsp*nq/nbpw
Ci         :see Remarks
Ci   wk    :a work array with dimension nbmx*nsp
Cio Inputs/Outputs
Cio   eb    :energy bands:
Cio         :mode=0 input spin-split, output merged to a single vector
Cio         :mode=1 input merged to a single vector, output spin-split
Cr Remarks
Cio   Call ebcpl with mode=1 to undo call of ebcpl with mode=0.
Cio   bmap must be preserved for mode=1 call.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap(1)
      double precision eb(nbmx,nsp,nq),wk(nbmx*nsp)
      integer ib,iq,ib1,ib2,iqb,iget

      if (nsp .eq. 1 .or. nspc .eq. 2) return

C --- Gather bands at each qp ---
      if (mode .eq. 0) then
        iqb = 0
        do  10  iq = 1, nq

C   ... Gather and order +,- bands at this qp into one column
          ib1 = 1
          ib2 = 1
          do  20  ib = 1, nevx*nsp
            iqb = iqb+1
            if (eb(min(ib1,nevx),1,iq) .lt. eb(min(ib2,nevx),2,iq)
     .      .and. ib1 .le. nevx .or. ib2 .gt. nevx) then
              wk(ib) = eb(ib1,1,iq)
              ib1 = ib1+1
            else
              wk(ib) = eb(ib2,2,iq)
              ib2 = ib2+1
              call mark1(bmap, nbpw, iqb)
            endif
C         call awrit6(' iq=%,2i  ib=%,2i  evl=%d  down=%i  ib1=%i  '//
C    .      'ib2=%i',' ',80,i1mach(2),iq,ib,wk(ib),
C    .      iget(bmap,nbpw,iqb),ib1,ib2)
   20     continue
          call dcopy(nevx*nsp,wk,1,eb(1,1,iq),1)
          if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
   10   continue
      endif

C --- Disperse bands at each qp ---
      if (mode .eq. 1) then
        iqb = 0
        do  110  iq = 1, nq

C   ... Disperse bands into +,- for this qp according to bmap
          ib1 = 1
          ib2 = 1
          call dcopy(nevx*nsp,eb(1,1,iq),1,wk,1)
          do  120  ib = 1, nevx*nsp
            iqb = iqb+1
            if (iget(bmap,nbpw,iqb) .eq. 0) then
              eb(ib1,1,iq) = wk(ib)
C            call awrit4(' iq=%,2i  ib=%,2i  evl=%d  down ib1=%i',
C     .         ' ',80,i1mach(2),iq,ib,wk(ib),ib1)
              ib1 = ib1+1
            else
              eb(ib2,2,iq) = wk(ib)
C            call awrit4(' iq=%,2i  ib=%,2i  evl=%d    up ib2=%i',
C     .         ' ',80,i1mach(2),iq,ib,wk(ib),ib2)
              ib2 = ib2+1
            endif
C         call awrit6(' iq=%,2i  ib=%,2i  evl=%d  down=%i  ib1=%i  '//
C    .      'ib2=%i',' ',80,i1mach(2),iq,ib,wk(ib),iget(bmap,nbpw,iqb),
C    .      ib1,ib2)
  120     continue
          if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
  110   continue
      endif

      end
      subroutine mark1(bitmap, nbpw, n)
C- put a one in the nth bit of bitmap.
C ----------------------------------------------------------------
Ci Inputs
Ci   bitmap, n
Cr Remarks
Cr    nbpw: a number no larger than the number of bits per integer word
C ----------------------------------------------------------------
C     implicit none
      integer bitmap(1), nbpw, n
C Local parameters
      integer nword,nbit,i

      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      i = 2**(nbpw-nbit-1)
      bitmap(nword+1) = bitmap(nword+1) + i*(1-mod(bitmap(nword+1)/i,2))
      end

      integer function iget(bitmap, nbpw, n)
C- Return 0 or 1, depending on the value of the nth bit of bitmap
C ----------------------------------------------------------------
Cr Remarks
Cr   See mark1
C ----------------------------------------------------------------
C     implicit none
      integer bitmap(1), nbpw, n
C Local parameters
      integer nword,nbit
      nword = (n-1)/nbpw
      nbit = mod(n-1,nbpw)
      iget = mod(bitmap(nword+1)/2**(nbpw-nbit-1),2)
      end

