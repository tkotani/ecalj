subroutine bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,&
  metal,tetra,norder,npts,width,rnge,wtkp,eb,efermi,sumev,&
  wtkb,dosef,qval,ent,lfill)
  use m_ftox
  use m_lmfinit,only: stdo,stdl
  use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
  use m_MPItk,only: master_mpi
  !- BZ integration for fermi level, band sum and qp weights
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbmx  :leading dimension of eb
  !i   nevx  :leading dimension of wtkb and max number of evals calculated
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
  !i   n1..n3:number of divisions for the k-point mesh
  !i   nkp   :number of inequivalent k-points (bzmesh.f)
  !i   ntet  :number of inequivalent tetrahedra (tetirr.f)
  !i   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
  !i         :corners of tetrahedron;
  !i         :idtet(0,i) number of tetrahedra of the i'th kind
  !i   zval  :valence charge
  !i   metal :T => metal, F => nonmetal
  !i   tetra :T => tetrahedron integration
  !i   norder:(sampling) polynomial order in Methfessel-Paxton integration
  !i         :100s digit norder flags that metals treatment should apply
  !i         :     regardless of whether a gap is present or not
  !i   width :(sampling) gaussian width in Methfessel-Paxton integration
  !i   npts  :(sampling) number of points in DOS mesh
  !i   rnge  :(sampling) range over which sampling delta function is assumed
  !i         :to vanish, in units of width
  !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
  !i   eb    :energy bands; alias eband
  !i   eb    :energy bands
  !o Outputs
  !o   efermi:Fermi energy
  !o   sumev :sum of eigenvalues
  !o   wtkb  :integration weights (not generated for nonmetal case)
  !o   dosef :DOS at Fermi level
  !o   qval  :qval(1) = total charge; qval(2) = magnetic moment
  !o   ent   :entropy term (actually TS)
  !l Local variables
  !l   lfill :true => insulator
  !u Updates
  !u   12 Jul 08 (ATP) bzwts now returns entropy term (actually kTS)
  !u   04 Jun 08 (ATP) Handles metal case when nkp=1
  !u    4 Aug 07 bzwts can make and tabulate specific heat (F-D statistics)
  !u   29 Jul 07 (ATP) Find E_F using weights by bisection, not INTNOS
  !u   02 Jan 06 return qval (valence charge and moment)
  !u   17 Jan 05 Use 100s digit norder as flag to treat all cases as metal,
  !u             whether or not a gap is present
  !u    1 May 04 When insulator, write gap to log file
  !u   09 May 04 When insulator, write gap to log file
  !u   01 Jul 03 When insulator, prints highest occ and lowest unocc state
  !u   24 Oct 02 Patch for weird cases when idos doesn't encompass
  !u             zval, where emin, emax found by efrang.
  !u   22 Sep 01 Returns dosef now.  Altered argument list.
  ! ----------------------------------------------------------------------
  implicit none
  ! Passed parameters
  logical metal,tetra
  integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,&
       idtet(5,ntet),ifile_handle
  double precision zval,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),&
       wtkb(nevx,nsp,nkp),efermi,sumev,dosef(2),qval(2),ent
  ! Local variables
  integer:: it , itmax , n , nptdos , nspx , nbmxx , nevxx , ib &
       , ikp , ipr , job , nbpw , i1mach , nev&
       , mkdlst , ifi , i , j , lry , procid , mpipid , master &
       , nulli , isw
  real(8) ,allocatable :: dos_rv(:)
  real(8) ,allocatable :: bot_rv(:)
  real(8) ,allocatable :: top_rv(:)
  integer ,allocatable :: bmap_iv(:)
  real(8) ,allocatable :: wk_rv(:)
  double precision emin,emax,e1,e2,dum,tol,e,elo,ehi,sumwt,&
       dmin,dmax,dval,egap,amom,dsum,cv,tRy
  character outs*100,ryy*3
  logical cmdopt,efrng2,lfill
  real(8) ,allocatable :: tlst_rv(:)
  parameter (nulli=-99999)
  integer:: iprint, w(1)

  ! --- Locate band limits: lfill => insulator ---
  procid = mpipid(1)
  master = 0
  call tcn('bzwts')
  ipr=iprint() !call getpr(ipr)
  qval(1) = 0
  qval(2) = 0
  ent = 0
  n = isign(1,norder) * mod(iabs(norder),100)
  !hangenglob      stdo = nglob('stdo')
  !      stdo = globalvariables%stdo
  !hangenglob      stdl = nglob('stdl')
  !      stdl = globalvariables%stdl
  allocate(bot_rv(nevx*nsp))
  allocate(top_rv(nevx*nsp))
  nspx  = 1
  nevxx = nevx*nsp
  nbmxx = nbmx*nsp
  !     job = 1 for non spin pol, -1 for spin pol
  job = 3-2*nsp
  dosef(1) = 0
  dosef(nspx) = 0
  egap = nulli
  if (nsp .eq. 2 .and. nspc .eq. 1 .and. cmdopt('--oldbz',7,0,outs))&
       then
     nspx  = nsp
     nevxx = nevx
     nbmxx = nbmx
     job = 1
  endif

  !     Force coupled spins: find range
  if (nspx .ne. nsp .and. nspc .eq. 1) then
     nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
     allocate(bmap_iv(nevx*nsp*nkp/nbpw+1))
     bmap_iv(:)=0

     allocate(wk_rv(nevx*nsp))
     !#ifdef KINODEBUG
     !        write(*,*) 'kino allocate bmap_iv, wk_rv',nevx*nsp*nkp/nbpw+1,nevx*nsp
     !#endif

     call ebcpl ( 0 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv &
          , wk_rv , eb )

     lfill = efrng2 ( nspx , nkp , nbmxx , nevxx , zval * 2 , eb , &
          bot_rv , top_rv , elo , ehi , emin , emax )

     call ebcpl ( 1 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv &
          , wk_rv , eb )

     !     Spins not coupled: find range
  else
     lfill = efrng2 ( nspx , nkp , nbmxx , nevxx , nspc * zval , eb &
          , bot_rv , top_rv , elo , ehi , emin , emax )

  endif
  ! ... Bands never filled if 100s digit norder set
  if (.not. tetra .and. iabs(norder) .ge. 100) lfill = .false.
  if (allocated(wk_rv)) deallocate(wk_rv)
  if (allocated(bmap_iv)) deallocate(bmap_iv)
  if (allocated(top_rv)) deallocate(top_rv)
  if (allocated(bot_rv)) deallocate(bot_rv)

  ! ... Case an insulator: put efermi at emin + tiny number
  if (lfill) then
     efermi = emin + 1d-10
     ! ... Do the best we can should be a metal, but assumption that it isn't
  elseif (.not. metal) then
     efermi = (emin + emax) / 2
  endif

  ! ... Pretend as though spin-pol bands are coupled to find E_f
  if (nsp .eq. 2 .and. nspc .eq. 1 .and. job .eq. -1) then
     nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
     allocate(bmap_iv(nevx*nsp*nkp/nbpw+1))
     bmap_iv(:)=0
     allocate(wk_rv(nevx*nsp))
     call ebcpl ( 0 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv &
          , wk_rv , eb )
  endif

  ! --- BZ weights, sumev and E_f for an insulator  ---
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
        write(stdo,ftox)' BZWTS : --- Non-metal sampling ---'
        write(stdo,ftox)' Fermi energy:',ftof(efermi),'electrons: num=',ftof(sumwt),&
             'occ. bands =',ftof(sumev)
        write(stdo,ftox)' VBmax=',ftof(emin),'CBmin=',ftof(emax)&
             &  ,'gap=',ftof(emax-emin), 'Ry =',ftof((emax-emin)*13.6058d0),'eV'
     endif

     ! --- BZ weights, sumev and E_f by tetrahedron method (Blochl wts) ---
  else if (tetra) then
     if (ipr .ge. 30) write (stdo,103)
103  format(/' BZWTS : --- Tetrahedron Integration ---')
     if (lfill) then
        egap = emax-emin
        if(ipr>=30)write(stdo,ftox)' ... only filled or empty bands encountered: ev=',&
             ftof(emin),'ec=',ftof(emax)
        if(ipr>=30)write(stdo,ftox)' VBmax=',ftof(emin),'CBmin=',ftof(emax),'gap =',&
             ftof(emax-emin),'Ry = ',ftof((emax-emin)*13.6058d0),'eV'
        goto 2
     endif
     nptdos = 101
     allocate(dos_rv(nspx*nptdos))
     tol = 1d-6
     !  Preliminary check that dos lies within emin,emax.  Widen emin,emax if not
     if (.not. lfill) then
        call bzints ( n1 , n2 , n3 , eb , dum , nkp , nevxx , nbmxx , &
             nspx , emin , emax , dos_rv , nptdos , efermi , job , ntet &
             , idtet , sumev , qval )
        dmin = dval ( dos_rv , 1 )
        if ( nspx .eq. 2 ) dmin = dmin + dval ( dos_rv , nptdos +  1 )
        dmax = dval ( dos_rv , nptdos )
        if ( nspx .eq. 2 ) dmax = dmax + dval ( dos_rv , nptdos +  nptdos )
        if (dmin .gt. zval) then
           emin = 3*emin-2*emax
           write(stdo,"(' (warning): initial NOS ( ',d14.6,x,d14.6,' ) does'//&
                ' not encompass Q=',f14.6)") dmin,dmax,zval
        elseif (dmax .lt. zval) then
           emax = 3*emax-2*emin
           write(stdo,"(' (warning): initial NOS ( ',d14.6,x,d14.6,' ) does'//&
                ' not encompass Q=',f14.6)") dmin,dmax,zval
        endif
     endif
     if (ipr .ge. 35) print 101
101  format(9x,'Est E_f ',10x,'Window',8x,'Tolerance',2x,'n(E_f)')
     itmax = 5
     do   it = 1, itmax
        call bzints ( n1 , n2 , n3 , eb , dum , nkp , nevxx , nbmxx , &
             nspx , emin , emax , dos_rv , nptdos , efermi , job , ntet &
             , idtet , sumev , qval )
        call fermi ( zval , dos_rv , nptdos , emin , emax , nspx , &
             efermi , emin , emax , dosef )
        if (ipr .ge. 35)&
             write(stdo,100) efermi,emin,emax,emax-emin,dosef(1)
100     format(7x,6(f10.6,1x))
        if (emax-emin .lt. tol) goto 1
     enddo
     if (ipr .gt. 10)&
          write(stdo,ftox)' BZWTS (warning): Fermi energy not converged: '&
          ,ftof(emax-emin),' > tol=',ftof(tol)
1    continue
     if (allocated(dos_rv)) deallocate(dos_rv)
2    continue
     call bzints(n1,n2,n3,eb,wtkb,nkp,nevxx,nbmxx,&
          nspx,emin,emin,emin,1,efermi,2*job,ntet,idtet,sumev,qval)
  else
     ! --- BZ weights, sumev and E_f by Methfessel-Paxton sampling ---
     write(stdo,"(a,i0,a,f15.6)")' BZWTS : --- Brillouin Zone sampling; N=',n,' W=',width
     !   ... Temporarily remove spin degeneracy if spins are coupled
     if (nsp .eq. 2 .and. nspx .eq. 1) call dscal(nkp,.5d0,wtkp,1)
     !   ... Find Fermi level, sampling
     if ((.not. lfill) .or. (metal .and. (nkp .eq. 1))) then
        e1 = elo - rnge*width/2
        e2 = ehi + rnge*width/2
        efermi = 0.5d0*(e1 + e2)
        itmax = 1000
        do  it = 1, itmax
           call pshpr(0)
           call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efermi,&
                .true.,sumev,wtkb,qval,ent,dosef,cv)
           call poppr
           if (dabs(zval - qval(1)) .lt. 1d-12) then
              write(stdo,ftox)' Fermi energy, ',ftof(efermi),&
                   &' found after ',it,' bisections,',ftof(qval(1)),&
                   &' electrons, DOS(E_f)=',ftof(dosef(1))
              goto 3
           endif
           if (qval(1) .gt. zval) then
              e2 = efermi
           else
              e1 = efermi
           endif
           efermi = 0.5d0*(e1 + e2)
        enddo
        write(stdo,ftox)' BZWTS (warning): cannot find E_F by bisection, using INTNOS'
        allocate(dos_rv(nspx*npts))
        emin = elo - rnge*width/2
        emax = emax + rnge*width/2
        call maknos ( nkp , nevxx , nbmxx , nspx , wtkp , eb , n , width &
             , - rnge , emin , emax , npts , dos_rv )
        if ( nspx.eq.2 ) call dpsadd ( dos_rv , dos_rv , npts , 1 , npts + 1 , 1d0 )
        call intnos ( npts , dos_rv , emin , emax , zval , efermi , dosef , sumev )
        if (allocated(dos_rv)) deallocate(dos_rv)
3       continue
     else
        dosef(1) = 0
        dosef(2) = 0
        egap = emax-emin
        write(stdo,ftox)' ... only filled or empty bands'//&
             ' encountered:  ev=',ftof(emin),' ec=',ftof(emax)
        write(stdo,ftox)' VBmax = ',ftof(emin),' CBmin = ',ftof(emax),' gap = ',&
             ftof(emax-emin),'Ry = ',ftof((emax-emin)*13.6058d0,3),'eV'
     endif

     !   ... (optional) Tabulate specific heat in file for list of T's
     if ((cmdopt('--cv:',5,0,outs) .or. cmdopt('--cvK:',6,0,outs))&
          .and. n .lt. 0 .and. metal) then
        if (procid .eq. master) then
           if (cmdopt('--cvK:',6,0,outs)) then
              lRy = 0
              i = 7
              ryy='K'
           else
              lRy = 1
              i = 6
              ryy='Ry'
           endif
           itmax = mkdlst(outs(i:),1d-8,-1,w)
           if (itmax .gt. 0) then
              allocate(tlst_rv(itmax))
              call word(outs,1,it,j)
              write(stdo,ftox)' Writing CV(T) to file for ',itmax,' vals of T: '&
                   ,outs(i:j),trim(ryy) !' %?#(n==1)#(Ry)#(K)#',itmax,lRy)
              it = mkdlst ( outs ( i: ) , 1d-8 , itmax + 1 , tlst_rv )
              if (it .ne. itmax) call rx('bzwts: bug in mkdlst')
              ifi = ifile_handle()
              open(ifi,file='cv.'//trim(sname)) 
              rewind ifi
              write(ifi,ftox)'% rows ',it,' cols 4 #   T(K)    T(Ry)   S(k_B)   TdS/dT(k_B)'
              do  it = 1, itmax
                 try = dval ( tlst_rv , it )
                 if (lRy .eq. 0) then
                    tRy = tRy/0.1579d6
                 endif
                 call pshpr(1)
                 call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,tRy,efermi,&
                      metal,sumev,wtkb,qval,ent,dosef,cv)
                 call poppr
                 write(ifi,ftox)ftof(0.1579d6*tRy),ftof(tRy),ftof(ent),ftof(cv)
              enddo
              close(ifi)
           endif
        endif
     endif

     !   ... Make weights, sampling
     call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efermi,&
          (.not. lfill) .or. (metal .and. (nkp .eq. 1)),&
          sumev,wtkb,qval,ent,dosef,cv)

     !   ... Put back spin degeneracy if removed
     if (nsp .eq. 2 .and. nspx .eq. 1) call dscal(nkp,2d0,wtkp,1)

  endif

  ! ... Restore to uncoupled bands; ditto with weights
  if (nsp .eq. 2 .and. nspc .eq. 1 .and. job .eq. -1) then
     call ebcpl ( 1 , nbmx , nevx , nsp , nspc , nkp , nbpw , bmap_iv &
          , wk_rv , eb )

     if ( metal ) call ebcpl ( 1 , nevx , nevx , nsp , nspc , nkp &
          , nbpw , bmap_iv , wk_rv , wtkb )

     if (allocated(tlst_rv)) deallocate(tlst_rv)
     if (allocated(wk_rv)) deallocate(wk_rv)
     if (allocated(bmap_iv)) deallocate(bmap_iv)

  endif

  ! ... Magnetic moment
  amom = 0
  if (nsp .eq. 2 .and. nspc .ne. 2 .and. metal) then
     do  ikp = 1, nkp
        amom = amom + dsum(nevx,wtkb(1,1,ikp),1) -&
             dsum(nevx,wtkb(1,2,ikp),1)
     enddo
     !        if (ipr .gt. 0) write(stdo,922) amom
     !  922   format(9x,'Mag. moment:',f15.6)
  endif
  qval(2) = amom
  if (ipr .gt. 0) then
     write(stdl,ftox)'bzmet',metal,'tet',tetra,'ef',ftof(efermi),&
          'sev',ftof(sumev),'zval',ftof(zval)
     write(stdl,ftox)'qval',ftof(qval(1)),'amom',ftof(amom),'egap(eV)',ftof(egap,3)
!     call awrit5('%a qval %,1;6d%?#n# amom %,1;6d#%j#'//&
!          '%?#n# gap %,4;4d eV##',outs,len(outs),-stdl,&
!          qval,int(amom*10000),amom,isw(egap.gt.0),egap*13.6d0)
  endif
  e = efermi
  if (.not. lfill .and. .not. tetra) e = efermi + rnge*width/2
  call tcx('bzwts')
end subroutine bzwts

subroutine ebcpl(mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap,wk,eb)
  !- Gather spin-polarized bands into a single group, or redistribute
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0, gather; 1, scatter
  !i   nbmx  :leading dimension of b and dimensions wk
  !i   nevx  :number of eigenvalues
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
  !i   nq    :number of k-points for which eb are calculated
  !i   nbpw  :a number no larger than the number of bits per integer word
  !i   bmap  :an integer work array with dimension at nbmx*nsp*nq/nbpw
  !i         :see Remarks
  !i   wk    :a work array with dimension nbmx*nsp
  !io Inputs/Outputs
  !io   eb    :energy bands:
  !io         :mode=0 input spin-split, output merged to a single vector
  !io         :mode=1 input merged to a single vector, output spin-split
  !r Remarks
  !io   Call ebcpl with mode=1 to undo call of ebcpl with mode=0.
  !io   bmap must be preserved for mode=1 call.
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  integer mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap(1)
  double precision eb(nbmx,nsp,nq),wk(nbmx*nsp)
  integer ib,iq,ib1,ib2,iqb,iget

  if (nsp .eq. 1 .or. nspc .eq. 2) return

  ! --- Gather bands at each qp ---
  ifx:if (mode .eq. 0) then
     iqb = 0
     do10:do    iq = 1, nq
        !   ... Gather and order +,- bands at this qp into one column
        ib1 = 1
        ib2 = 1
        do20: do ib = 1, nevx*nsp
           iqb = iqb+1
           if (eb(min(ib1,nevx),1,iq) .lt. eb(min(ib2,nevx),2,iq)&
                .and. ib1 .le. nevx .or. ib2 .gt. nevx) then
              wk(ib) = eb(ib1,1,iq)
              ib1 = ib1+1
           else
              wk(ib) = eb(ib2,2,iq)
              ib2 = ib2+1
              call mark1(bmap, nbpw, iqb)
           endif
           !         call awrit6(' iq=%,2i  ib=%,2i  evl=%d  down=%i  ib1=%i  '//
           !    .      'ib2=%i',' ',80,i1mach(2),iq,ib,wk(ib),
           !    .      iget(bmap,nbpw,iqb),ib1,ib2)
        enddo do20
        call dcopy(nevx*nsp,wk,1,eb(1,1,iq),1)
        if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
     enddo do10
  endif ifx

  ! --- Disperse bands at each qp ---
  if (mode .eq. 1) then
     iqb = 0
     d110  : do  iq = 1, nq
        !   ... Disperse bands into +,- for this qp according to bmap
        ib1 = 1
        ib2 = 1
        call dcopy(nevx*nsp,eb(1,1,iq),1,wk,1)
        do120  :     do  ib = 1, nevx*nsp
           iqb = iqb+1
           if (iget(bmap,nbpw,iqb) .eq. 0) then
              eb(ib1,1,iq) = wk(ib)
              !            call awrit4(' iq=%,2i  ib=%,2i  evl=%d  down ib1=%i',
              !     .         ' ',80,i1mach(2),iq,ib,wk(ib),ib1)
              ib1 = ib1+1
           else
              eb(ib2,2,iq) = wk(ib)
              !            call awrit4(' iq=%,2i  ib=%,2i  evl=%d    up ib2=%i',
              !     .         ' ',80,i1mach(2),iq,ib,wk(ib),ib2)
              ib2 = ib2+1
           endif
           !         call awrit6(' iq=%,2i  ib=%,2i  evl=%d  down=%i  ib1=%i  '//
           !    .      'ib2=%i',' ',80,i1mach(2),iq,ib,wk(ib),iget(bmap,nbpw,iqb),
           !    .      ib1,ib2)
        enddo do120
        if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
     enddo d110
  endif
end subroutine ebcpl
   
subroutine mark1(bitmap, nbpw, n)
  !- put a one in the nth bit of bitmap.
  ! ----------------------------------------------------------------
  !i Inputs
  !i   bitmap, n
  !r Remarks
  !r    nbpw: a number no larger than the number of bits per integer word
  ! ----------------------------------------------------------------
  !     implicit none
  integer bitmap(1), nbpw, n
  ! Local parameters
  integer nword,nbit,i

  nword = (n-1)/nbpw
  nbit = mod(n-1,nbpw)
  i = 2**(nbpw-nbit-1)
  bitmap(nword+1) = bitmap(nword+1) + i*(1-mod(bitmap(nword+1)/i,2))
end subroutine mark1

integer function iget(bitmap, nbpw, n)
  !- Return 0 or 1, depending on the value of the nth bit of bitmap
  ! ----------------------------------------------------------------
  !r Remarks
  !r   See mark1
  ! ----------------------------------------------------------------
  !     implicit none
  integer bitmap(1), nbpw, n
  ! Local parameters
  integer nword,nbit
  nword = (n-1)/nbpw
  nbit = mod(n-1,nbpw)
  iget = mod(bitmap(nword+1)/2**(nbpw-nbit-1),2)
end function iget

