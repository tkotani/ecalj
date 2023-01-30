module  m_bzwts
  public bzwtsf,bzwtsf2
  private
contains
  subroutine bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,&
       metal,tetra,norder,npts,width,rnge,wtkp,eb, efermi,sumev,wtkb,dosef,qval,ent,lfill)
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
    logical metal,tetra
    integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,&
         idtet(5,ntet),ifile_handle
    double precision zval,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),&
         wtkb(nevx,nsp,nkp),efermi,sumev,dosef(2),qval(2),ent
    integer:: it , itmax , n , nptdos , nspx , nbmxx , nevxx , ib &
         , ikp , ipr , job , nbpw , i1mach , nev&
         , mkdlst , ifi , i , j , lry , procid , mpipid , master &
         , nulli , isw
    real(8) ,allocatable :: dos_rv(:,:)
    real(8) ,allocatable :: bot_rv(:)
    real(8) ,allocatable :: top_rv(:)
    integer ,allocatable :: bmap_iv(:)
    real(8) ,allocatable :: wk_rv(:)
    double precision emin,emax,e1,e2,dum,tol,e,elo,ehi,sumwt,&
         dmin,dmax,egap,amom,cv,tRy
    character outs*100,ryy*3
    logical cmdopt0,lfill
    real(8) ,allocatable :: tlst_rv(:)
    parameter (nulli=-99999)
    integer:: iprint, w(1)
    procid = mpipid(1)
    master = 0
    call tcn('bzwts')
    ipr=iprint()
    qval(1) = 0
    qval(2) = 0
    ent = 0
    n = isign(1,norder) * mod(iabs(norder),100)
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
    !  if (nsp .eq. 2 .and. nspc .eq. 1 .and. cmdopt('--oldbz',7,0,outs))&
    !       then
    !     nspx  = nsp
    !     nevxx = nevx
    !     nbmxx = nbmx
    !     job = 1
    !  endif
    !     Force coupled spins: find range
    if (nspx .ne. nsp .and. nspc .eq. 1) then
       nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
       allocate(bmap_iv(nevx*nsp*nkp/nbpw+1))
       bmap_iv(:)=0
       allocate(wk_rv(nevx*nsp))
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
103    format(/' bzwts: --- Tetrahedron Integration ---')
       if (lfill) then
          egap = emax-emin
          if(ipr>=30)write(stdo,ftox)' ... only filled or empty bands encountered: ev=',&
               ftof(emin),'ec=',ftof(emax)
          if(ipr>=30)write(stdo,ftox)' VBmax=',ftof(emin),'CBmin=',ftof(emax),'gap =',&
               ftof(emax-emin),'Ry = ',ftof((emax-emin)*13.6058d0),'eV'
          goto 2
       endif
       nptdos = 101
       allocate(dos_rv(nptdos,nspx))
       tol = 1d-6
       !  Preliminary check that dos lies within emin,emax.  Widen emin,emax if not
       if (.not. lfill) then
          call bzints ( n1 , n2 , n3 , eb , dum , nkp , nevxx , nbmxx , &
               nspx , emin , emax , dos_rv , nptdos , efermi , job , ntet &
               , idtet , sumev , qval )
          dmin = sum(dos_rv(1,1:nspx))/nspx
          dmax = sum(dos_rv(nptdos,1:nspx))/nspx
          !if ( nspx .eq. 2 ) dmin = dmin + dos_rv(1,nspx)
          !if ( nspx .eq. 2 ) dmax = dmax + dos_rv(nptdos,nspx)
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
101    format(9x,'Est E_f ',10x,'Window',8x,'Tolerance',2x,'n(E_f)')
       itmax = 5
       do   it = 1, itmax
          call bzints ( n1 , n2 , n3 , eb , dum , nkp , nevxx , nbmxx , &
               nspx , emin , emax , dos_rv , nptdos , efermi , job , ntet &
               , idtet , sumev , qval )
          call fermi ( zval , dos_rv , nptdos , emin , emax , nspx , &
               efermi , emin , emax , dosef(1) )
          if (ipr .ge. 35)&
               write(stdo,100) efermi,emin,emax,emax-emin,dosef(1)
100       format(7x,6(f10.6,1x))
          if (emax-emin .lt. tol) goto 1
       enddo
       if(ipr>10) write(stdo,ftox)' BZWTS (warning): Fermi energy not converged: '&
            ,ftof(emax-emin),' > tol=',ftof(tol)
1      continue
       if (allocated(dos_rv)) deallocate(dos_rv)
2      continue
       call bzints(n1,n2,n3,eb,wtkb,nkp,nevxx,nbmxx,&
            nspx,emin,emin,emin,1,efermi,2*job,ntet,idtet,sumev,qval)
    else
       ! --- BZ weights, sumev and E_f by Methfessel-Paxton sampling ---
       if(ipr>0) write(stdo,"(a,i0,a,f15.6)")' BZWTS : --- Brillouin Zone sampling; N=',n,' W=',width
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
          allocate(dos_rv(npts,nspx))
          emin = elo - rnge*width/2
          emax = emax + rnge*width/2
          call maknos ( nkp , nevxx , nbmxx , nspx , wtkp , eb , n , width &
               , - rnge , emin , emax , npts , dos_rv )
          if ( nspx.eq.2 ) dos_rv(:,1)=dos_rv(:,2)+dos_rv(:,1)
          call intnos ( npts , dos_rv , emin , emax , zval , efermi , dosef , sumev )
          if (allocated(dos_rv)) deallocate(dos_rv)
3         continue
       else
          dosef(1) = 0
          dosef(2) = 0
          egap = emax-emin
          if(ipr>0) write(stdo,ftox)' ... only filled or empty bands'//&
               ' encountered:  ev=',ftof(emin),' ec=',ftof(emax)
          if(ipr>0) write(stdo,ftox)' VBmax = ',ftof(emin),' CBmin = ',ftof(emax),' gap = ',&
               ftof(emax-emin),'Ry = ',ftof((emax-emin)*13.6058d0,3),'eV'
       endif
       ! ... (optional) Tabulate specific heat in file for list of T's
       if (cmdopt0('--cvK:') .and. n .lt. 0 .and. metal) then
          if (procid .eq. master) then
             lRy = 0
             ryy='K'
             itmax=8
             allocate(tlst_rv(itmax))
             tlst_rv=[10,20,40,80,160,320,640,1280] !fixed now
             write(stdo,ftox)'Writing CV(T) to file for ',itmax,'vals of T:',ftof(tlst_rv,1),trim(ryy) 
             !open(newunit=ifi,file='cv.'//trim(sname)) 
             write(stdo,ftox)'% rows ',it,' cols 4 #   T(K)    T(Ry)   S(k_B)   TdS/dT(k_B)'
             do  it = 1, itmax
                tRy = tlst_rv(it)/0.1579d6
                call pshpr(1)
                call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,tRy,efermi,&
                     metal,sumev,wtkb,qval,ent,dosef,cv)
                call poppr
                write(stdo,ftox)ftof(0.1579d6*tRy),ftof(tRy),ftof(ent),ftof(cv)
             enddo
             !close(ifi)
          endif
       endif
       !   ... Make weights, sampling
       call splwts(nkp,nevxx,nbmxx,nspx,wtkp,eb,n,width,efermi,&
            (.not. lfill) .or. (metal .and. (nkp .eq. 1)), sumev,wtkb,qval,ent,dosef,cv)
       if(nsp .eq. 2 .and. nspx .eq. 1) call dscal(nkp,2d0,wtkp,1)
    endif
    ! ... Restore to uncoupled bands; ditto with weights
    if (nsp .eq. 2 .and. nspc .eq. 1 .and. job .eq. -1) then
       call ebcpl(1 , nbmx , nevx , nsp , nspc , nkp , nbpw, bmap_iv , wk_rv , eb )
       if(metal) call ebcpl(1 , nevx , nevx , nsp , nspc, nkp, nbpw , bmap_iv , wk_rv , wtkb )
       if(allocated(tlst_rv)) deallocate(tlst_rv)
       if(allocated(wk_rv)) deallocate(wk_rv)
       if(allocated(bmap_iv)) deallocate(bmap_iv)
    endif
    amom = 0d0 ! ... Magnetic moment
    if(nsp==2 .and. nspc/=2 .and. metal) amom = sum(wtkb(1:nevx,1,1:nkp)- wtkb(1:nevx,2,1:nkp))
    qval(2) = amom
    if (ipr .gt. 0) then
       write(stdl,ftox)'bzmet',metal,'tet',tetra,'ef',ftof(efermi),'sev',ftof(sumev),'zval',ftof(zval)
       write(stdl,ftox)'qval',ftof(qval(1)),'amom',ftof(amom),'egap(eV)',ftof(egap,3)
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
    integer ib,iq,ib1,ib2,iqb

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


  logical function efrng2(nsp,nkp,nbmax,nband,zval,eband,ebbot,ebtop,elo,ehi,e1,e2)
    !- Find range of Fermi energy.
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nkp   :number of irreducible k-points (bzmesh.f)
    !i   nbmax :leading dimension of eband
    !i   nband :number of bands
    !i   zval  :no. of valence electrons
    !i   eband :energy bands
    !o Outputs
    !o   e1,e2: e1 < ef < e2
    !o   elo, ehi:  lowest and highest band found
    !o   efrng2:flags whether metal or insulator.
    !o         :F metal (highest occ band crossed lowest unoccupied one)
    !o         :T insulator (highest occ band did not cross lowest unocc)
    !r Remarks
    !r    For an even no. of electrons ef is above the bottom of the
    !r    zval/2+1'th band and below the top of the zval/2 'th band. If the
    !r    former is higher that the latter we have an insulator, with
    !r    these two numbers estimates for the conduction band minimum and
    !r    valence band maximum, respectively.
    !r    For an odd no. of electrons ef is between the bottom and top
    !r    of the (zval+1)/2 'th band.
    !r
    !r    For spin pol case:
    !r      bottom of the zval+1'th band < ef < top of the zval'th band.
    !r      If the bottom is higher that the top then we have an insulator,
    !r      with the bottom an estimate for the conduction band minimum and
    !r      to and estimate for the valence band maximum.
    !r      and e1=e2.
    !u Updates
    !u   08 Jul 08 Extend to case where number of bands can be q dependent
    !u             eband=99999 => not calculated: ignore
    !u   01 Apr 03 Set e1=e2 if every band is filled
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer nsp,nkp,nbmax,nband
    double precision zval,e1,e2,eband(nbmax,nsp,nkp),ebbot(nband,nsp),ebtop(nband,nsp),elo,ehi
    ! Local parameters
    double precision xx,d1mach,enull
    integer ikp,isp,iba,nval,nbbot,nbtop,nfound
    parameter (enull=99999d0)
    elo = enull
    ehi = -enull
    nval = zval + 1d-7
    ! --- Find bottom and top of each band ---
    ebbot=elo 
    ebtop=ehi 
    do  ikp = 1, nkp
       do  isp = 1, nsp
          do  iba = 1, nband
             if (eband(iba,isp,ikp) .ne. enull) then
                ebbot(iba,isp) = min(ebbot(iba,isp),eband(iba,isp,ikp))
                ebtop(iba,isp) = max(ebtop(iba,isp),eband(iba,isp,ikp))
             endif
          enddo
       enddo
    enddo
    !     Set all -enull to enull to float to top when sorted
    do  isp = 1, nsp
       do  iba = 1, nband
          if (ebtop(iba,isp) .eq. -enull) ebtop(iba,isp) = enull
       enddo
    enddo
    !     Sort bands irrespective of spin
    call dshell(nband*nsp,ebbot)
    call dshell(nband*nsp,ebtop)
    nfound = nband*nsp
10  continue
    if (ebtop(nfound,1).eq.enull .or. ebbot(nfound,1).eq.enull) then
       nfound = nfound-1
       if (nfound .eq. 0) call rx('efrng2: no bands')
       goto 10
    endif
    ! --- Find limits ---
    nbtop = (nval+2-nsp)/(3-nsp)
    if (zval .gt. nval) nbtop = nbtop+1
    nbbot = nval/(3-nsp) + 1
    if (nbtop .gt. nfound) nbtop = nfound
    if (nbbot .gt. nfound) nbbot = nfound
    elo = ebbot(1,1)
    ehi = ebtop(nfound,1)
    if (elo .eq. enull) call rx('efrng2: no bands')
    e1  = ebbot(nbbot,1)
    e2  = ebtop(nbtop,1)
    efrng2 = .false.
    if (e1-e2 > epsilon(0d0)) then
       xx = e1
       e1 = e2
       e2 = xx
       efrng2 = .true.
    endif
  end function efrng2

  subroutine dshell(n,array)
    implicit none
    integer n
    double precision array(n)
    integer i,j,k,inc
    double precision v
    ! ... Get the largest increment
    if (n .le. 1) return
    inc = 1
10  continue
    inc = 3*inc+1
    if (inc .lt. n) goto 10
    ! ... Loop over partial sorts
12  continue
    inc = inc/3
    !   ... Outer loop of straight insertion
    do  11  i = inc+1, n
       v = array(i)
       j = i
       !     ... Inner loop of straight insertion
20     continue
       if (array(j-inc) .gt. v) then
          array(j) = array(j-inc)
          j = j-inc
          if (j .le. inc) goto 21
          goto 20
       endif
21     continue
       array(j) = v
11  enddo
    if (inc .gt. 1) goto 12
  end subroutine dshell

  subroutine fixef0(zval,nsp,nspc,nevx,ndev,evl,dosw,ef0)
    use m_lgunit,only:stdo
    !- Corrects estimate for Fermi level
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   zval  :valence charge
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
    !i   nevx  :max number of eigenvalues calculated
    !i   ndev  :leading dimension of evl
    !i   evl   :eigenvalues
    !i   dosw  :dos window
    ! o Inputs/Outputs
    ! o  ef0   :on input, estimate for Fermi energy
    ! o        :on output, revised estimate, if ef0 outside bounds
    ! o  dosw  :on input dos window
    ! o        :on output, revised if ebot<dosw(1) or dosw(2)<ef0
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: nsp,nspc,ndev,nevx
    double precision :: zval,ef0,evl(ndev*2),dosw(2)!takao evl(ndev)->evl(ndev*2)
    ! ... Local parameters
    integer:: i , i1 , ipr , nbpw , i1mach
    integer ,allocatable :: bmap_iv(:)
    real(8) ,allocatable :: wk_rv(:)

    double precision :: w2,xx,doso(2)
    ! ... Heap

    call getpr(ipr)
    !      stdo = lgunit(1)
    if (nsp == 2) then
       nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
       allocate(bmap_iv(abs(-(nevx*nsp/nbpw+1))))
       if (-(nevx*nsp/nbpw+1)<0) bmap_iv(:)=0

       allocate(wk_rv(nevx*nsp))

       call ebcpl ( 0 , ndev , nevx , nsp , nspc , 1 , nbpw , bmap_iv &
            , wk_rv , evl )

    endif
    i = max(1,int(zval)/(3-nsp))
    if (ef0 < evl(i)) then
       i1 = (zval + 0.001d0)/(3-nsp)
       w2 = zval/(3-nsp)-i1
       xx = (1-w2)*evl(i1)+w2*evl(i1+1)
       if (ipr >= 10) write(stdo,"(' Est Ef = ',f15.8,' < evl(',f15.8,i0,')=',f15.8, &
            ' ... using qval=',f15.8,' revise to ',f15.8)") ef0,i,evl(i),zval,xx
       ef0 = xx
    endif

    if (nsp == 2) then
       call ebcpl ( 1 , ndev , nevx , nsp , nspc , 1 , nbpw , bmap_iv &
            , wk_rv , evl )

       if (allocated(wk_rv)) deallocate(wk_rv)
       if (allocated(bmap_iv)) deallocate(bmap_iv)

    endif

    if (dosw(1) > evl(1) .OR. dosw(2) < ef0) then
       doso(1) = dosw(1)
       doso(2) = dosw(2)
       dosw(1) = evl(1) - 0.5d0
       dosw(2) = ef0  + 0.5d0
       if (ipr >= 10) write(stdo,"(' DOS window (',f15.7,x,f15.7,')', &
            ' reset to (',f15.7,x,f15.7,')')") doso(1),doso(2),dosw(1),dosw(2)
    endif
  end subroutine fixef0
  !!== FSMOMMETHOD=1 == June2011 takao
  !! This is for molecules.
  !! 1st step : Set initial bias magnetic field under assuming all eigenvalues are discrete.
  !!            Search elumo1,ehomo1 (for up spin) and elumo2,ehomo2 (for down spin) below.
  !! 2nd step : refine the bias field for given temperature.
  !!
  !!  (Takao think this bzwtsf2 may need to be modified for solids).
  subroutine bzwtsf2(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
       fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,lswtk, &
       swtk,efermi,sumev,wtkb,qval,lwtkb,lfill,vnow)
    use m_lmfinit,only: stdo
    use m_ftox

    !- BZ integration for fermi level, band sum and qp weights, fixed-spin
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
    !i   fmom  :fixed spin moment.  If zero, no constraint is applied.
    !i   metal :T => metal, F => nonmetal
    !i   tetra :T => tetrahedron integration
    !i   norder,npts,width,rnge: parameters for sampling integr. (maknos)
    !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
    !i   eb    :energy bands; alias eband
    !i   eb    :energy bands
    !i   lswtk :Flags indicating whether 'spin weights' swtk are available
    !i   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
    !i         :where z are eigenvectors, sigma is the Pauli spin matrix
    !i         :Supplies information about spin moment in noncoll case.
    !i         :Used when lswtk is set
    ! oxxx  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
    ! o        :0 weights are not available; no moment calculation
    ! o        :if 1, weights were generated with no constraint
    ! o        :In this case, print moment, and if fmom ne 0 remake weights
    ! o        :with constraint; set to lwtkb=2 on output.
    ! o        :if 2, weights were generated with constrained global moment
    ! o        :if -1, same as 1
    !o Outputs
    !o   efermi:Fermi energy
    !o   sumev :sum of eigenvalues
    !o   wtkb  :integration weights (not generated for nonmetal case)
    !o   qval  :qval(1) = total charge; qval(2) = magnetic moment
    !u Updates
    !u   12 Jul 08 change arg list in bzwts -- now returns entropy term
    !u   02 Jan 06 return qval (valence charge and moment)
    !u   22 Sep 01 Adapted from bzwts.
    ! ----------------------------------------------------------------------
    implicit none
    ! Passed parameters
    logical :: metal,tetra
    integer :: nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet, &
         idtet(5,ntet),lswtk,lwtkb
    double precision :: zval,fmom,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp), &
         wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efermi,sumev,qval(2)
    ! Local variables
    integer :: ikp,ib,ipr,itmax,iter,iprint
    double precision :: amom,dosef(2),vhold(12),vnow,dvcap,dv,ef0,ent
    parameter (dvcap=.2d0,itmax=50)

    logical:: agreemom
    real(8),parameter::    NULLR =-99999
    integer:: nmom1,nmom2
    real(8):: ehomo1,ehomo2,elumo1,elumo2
    real(8),allocatable:: ebs(:,:,:)
    real(8):: ele1,ele2
    integer:: itermx
    logical:: quitvnow,lfill

    ! --- Fermi level without spin constraint ---
    call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
         metal,tetra,norder,npts,width,rnge,wtkp,eb,efermi, &
         sumev,wtkb,dosef,qval,ent,lfill)
    if (nsp == 1) return

    call getpr(ipr)
    ! angenglob      stdo = nglob('stdo')
    !      stdo = globalvariables%stdo
    !     stdl = nglob('stdl')
    ! --- Make and print out magnetic moment ---
    if ((lswtk == 1 .OR. nspc == 1) .AND. metal) then
       call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
       if (ipr >= 20) write(stdo,922) amom
922    format(9x,'Mag. moment:',f15.6)
       qval(2) = amom
    else
       write(stdo,*) 'spin weights not available ... no spin moment calculated'
       return
    endif
    !! --- Setup for fixed-spin moment method ---
    if (fmom==NULLR .OR. lwtkb == 0) return
    !      if (fmom==NULLR .or. lwtkb .eq. 0) call rx('bzwtsf2: check logic of lwtkb')

    call tcn('bzwtsf2')

    ele1 = (zval+fmom)/2d0
    ele2 = (zval-fmom)/2d0
    nmom1 = ele1+1d-8
    nmom2 = ele2+1d-8
    elumo1= minval(eb(nmom1+1,1,:))
    elumo2= minval(eb(nmom2+1,2,:))
    if(nmom1/=0) then
       ehomo1= maxval(eb(nmom1,1,:))
    else
       ehomo1= elumo1 - 0.5d0
    endif
    if(nmom2/=0) then
       ehomo2= maxval(eb(nmom2,2,:))
    else
       ehomo2= elumo2 - 0.5d0
    endif
    vnow = (ehomo1+elumo1)/2d0 -(ehomo2+elumo2)/2d0
    efermi= ((ehomo1+elumo1)/2d0 +(ehomo2+elumo2)/2d0)/2d0 !/2d0 bug fix Jun26,2014 this only affects to the message.
    write(stdo,"('bzwtsf2: zval fmom nmon1 nmom2=',2f12.8,2x,2i3)")zval,fmom,nmom1,nmom2
    write(stdo,"('bzwtsf2: HOMOup LUMOup Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo1,elumo1,elumo1-ehomo1
    write(stdo,"('bzwtsf2: HOMOdn LUMOdn Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo2,elumo2,elumo2-ehomo2
    write(stdo,"('bzwtsf2: Set Bias initial cond. -Vup+Vdn=',f20.15)")vnow

    !!= takao interted a block taken from original version of bzwtsf.F June-2 2011.=
    vhold= 0d0
    ef0  = efermi
    write(stdo,*)' Seek potential shift for fixed-spin mom ...'

    !!== do loop for new guess at potential shift ==
    ! bisection method takao
    itermx=100
    quitvnow=.false.
    do 10 iter=1,itermx
       !! Potential shift
       allocate(ebs(nevx,2,nkp))
       if (nspc == 2) then
          ebs = eb + vnow/2*swtk
       else
          ebs(:,1,:) = eb(:,1,:) - vnow/2d0
          ebs(:,2,:) = eb(:,2,:) + vnow/2d0
       endif
       !! Fermi level with dv shift
       if( .NOT. quitvnow) call pshpr(ipr-50)
       call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
            metal,tetra,norder,npts,width,rnge,wtkp,ebs,efermi, &
            sumev,wtkb,dosef,qval,ent,lfill)
       if (iprint()>= 20) then
          call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk, &
               amom)
          write(stdo,922) amom
       endif
       deallocate(ebs)
       if ( .NOT. quitvnow) call poppr
       if(quitvnow) exit
       !!=== Magnetic moment ===
       call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
       if(ipr>=41)write(stdo,ftox)' -Vup+Vdn=',ftof(vnow,8),'yields ef=',ftof(efermi), &
            'amom',ftof(amom),'when seeking mom=',ftof(fmom)
       ! takao for molecule Dec1 2010
       agreemom= abs(amom-fmom) < 1d-6 ! 1d-6 on June-2 2011
       if(iprint()>60) print *,'ttttt amom fmom=',amom,fmom,agreemom
       call dvdos(vnow,amom,dosef,vhold,fmom,dvcap,dv)
       if(agreemom) vhold(12)=0
       quitvnow=.false.
       if (abs(dv) < 1d-6 .OR. agreemom) then
          if (vhold(12) == -2 .OR. vhold(12) == -3 .OR. &
               vhold(12) ==  0 .OR. vhold(12) ==  1) then
             if (ipr >= 10) &
                  write(stdo,ftox)' BZWTSF2: potential shift bracketed.', &
                  'Unconstrained efermi=',ftof(ef0), &
                  'constraint fmom=',ftof(fmom),'actual mmom=',amom, &
                  'ef=',efermi,'-Vup+Vdn=',vnow
             quitvnow=.true.
          endif
       else if (iter == itmax) then
          if (ipr >= 10) &
               write(stdo,ftox)' BZWTSF2: failed to converge potential shift', &
               'after',iter,'iterations.', &
               'constraint fmom=',ftof(fmom),'actual amom=',ftof(amom), &
               'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vnow)
          quitvnow=.true.
       endif
10  enddo
    ele1 = (zval+amom)/2d0
    ele2 = (zval-amom)/2d0
    sumev=sumev+ ele1*(vnow/2d0) - ele2*(vnow/2d0) !sumev correction takao
    if(iprint()>20) write(stdo,ftox)' bzwtsf2(METHOD=1): Set Bias field -Vup+Vdn=',ftof(vnow,8)
    if (lswtk == 1 .AND. lwtkb == 1) then
       call rx('bzwtsf2:111 tk think not used here')
    elseif (lswtk == 1 .AND. lwtkb == 2) then
       call rx('bzwtsf2:222 tk think not used here')
    endif
    call tcx('bzwtsf2')
  end subroutine bzwtsf2


  !!== FSMOMMETHOD=0 ogiginal version(modified version. fmom=0 is allowed.)==
  subroutine bzwtsf(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
       fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb,lswtk, &
       swtk,efermi,sumev,wtkb,qval,lwtkb,lfill,vnow)
    use m_lmfinit,only: stdo
    use m_ftox

    !- BZ integration for fermi level, band sum and qp weights, fixed-spin
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
    !i   fmom  :fixed spin moment.  Even zero is allowed.
    !i   metal :T => metal, F => nonmetal
    !i   tetra :T => tetrahedron integration
    !i   norder,npts,width,rnge: parameters for sampling integr. (maknos)
    !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
    !i   eb    :energy bands; alias eband
    !i   eb    :energy bands
    !i   lswtk :Flags indicating whether 'spin weights' swtk are available
    !i   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
    !i         :where z are eigenvectors, sigma is the Pauli spin matrix
    !i         :Supplies information about spin moment in noncoll case.
    !i         :Used when lswtk is set
    ! oxxx  lwtkb :Used in connection w/ fixed spin-moments method.  On input:
    ! o        :0 weights are not available; no moment calculation
    ! o        :if 1, weights were generated with no constraint
    ! o        :In this case, print moment, and if fmom ne 0 remake weights
    ! o        :with constraint; set to lwtkb=2 on output.
    ! o        :if 2, weights were generated with constrained global moment
    ! o        :if -1, same as 1
    !o Outputs
    !o   efermi:Fermi energy
    !o   sumev :sum of eigenvalues
    !o   wtkb  :integration weights (not generated for nonmetal case)
    !o   qval  :qval(1) = total charge; qval(2) = magnetic moment
    !u Updates
    !u   12 Jul 08 change arg list in bzwts -- now returns entropy term
    !u   02 Jan 06 return qval (valence charge and moment)
    !u   22 Sep 01 Adapted from bzwts.
    ! ----------------------------------------------------------------------
    implicit none
    ! Passed parameters
    logical :: metal,tetra
    integer :: nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet, &
         idtet(5,ntet),lswtk,lwtkb
    double precision :: zval,fmom,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp), &
         wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),efermi,sumev,qval(2)
    ! Local variables
    integer :: ikp,ib,ipr,itmax,iter,iprint
    double precision :: amom,dosef(2),vhold(12),vnow,dvcap,dv,ef0,ent
    parameter (dvcap=.2d0,itmax=50)

    real(8):: ele1,ele2
    integer:: itermx
    logical:: agreemom
    real(8),parameter::    NULLR =-99999
    real(8),allocatable:: ebs(:,:,:)
    logical:: quitvnow,lfill

    !!== Fermi level without spin constraint ==
    call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
         metal,tetra,norder,npts,width,rnge,wtkp,eb,efermi, &
         sumev,wtkb,dosef,qval,ent,lfill)
    if (nsp == 1) return

    call getpr(ipr)
    !      stdo = globalvariables%stdo

    !!== Make and print out magnetic moment ==
    if ((lswtk == 1 .OR. nspc == 1) .AND. metal) then
       call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
       if (ipr >= 20) write(stdo,922) amom
922    format(9x,'Mag. moment:',f15.6)
       qval(2) = amom
    else
       write(stdo,*)'spin weights not available ... no spin moment calculated'
       return
    endif

    !!== Setup for fixed-spin moment method ==
    !      if (fmom .eq. 0 .or. lwtkb .eq. 0) return
    if (fmom==NULLR .OR. lwtkb == 0) return
    !      if (fmom==NULLR .or. lwtkb .eq. 0) call rx('bzwtsf:check logic of lwtkb')
    call tcn('bzwtsf')
    call dpzero(vhold,12)
    vnow = 0
    ef0 = efermi
    write(stdo,*)' Seek potential shift for fixed-spin mom ...'

    !!== do loop for new guess at potential shift ==
    ! bisection method takao
    itermx=100
    do 10 iter=1,itermx
       !!=== Magnetic moment ===
       call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk,amom)
       if(ipr>=41) write(stdo,ftox)' -Vup+Vdn=',ftof(vnow,8),'yields ', &
            'ef=',ftof(efermi),'amom=',ftof(amom),'when seeking',ftof(fmom)
       ! takao for molecule Dec1 2010
       agreemom= abs(amom-fmom) < 1d-3
       if(iprint()>60) print *,'ttttt amom fmom=',amom,fmom,agreemom
       !      if(original_dvdos) then
       call dvdos(vnow,amom,dosef,vhold,fmom,dvcap,dv)
       if(agreemom) vhold(12)=0
       !      if (abs(dv) .lt. 1d-6) then
       quitvnow=.false.
       if (abs(dv) < 1d-6 .OR. agreemom) then
          !       A root was found
          if (vhold(12) == -2 .OR. vhold(12) == -3 .OR. &
               vhold(12) ==  0 .OR. vhold(12) ==  1) then
             if (ipr >= 10) write(stdo,ftox)' BZWTSF: potential shift bracketed.', &
                  'Unconstrained efermi=',ftof(ef0), &
                  'constraint fmom=',ftof(fmom),'actual mmom=',ftof(amom), &
                  'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vnow,8)
             quitvnow=.true.
          endif
       else if (iter == itmax) then
          if(ipr>=10)then
             write(stdo,ftox)' BZWTSF: failed to converge potential shift after',iter,'iterations.'
             write(stdo,ftox)' constraint fmom=',ftof(fmom),'actual amom=',ftof(amom), &
                  'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vnow,8)
          endif
          quitvnow=.true.
       endif

       !! Potential shift
       allocate(ebs(nevx,2,nkp))
       if (nspc == 2) then
          ebs = eb + vnow/2*swtk
       else
          ebs(:,1,:) = eb(:,1,:) - vnow/2
          ebs(:,2,:) = eb(:,2,:) + vnow/2
       endif

       !! Fermi level with dv shift
       if( .NOT. quitvnow) call pshpr(ipr-50)
       call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
            metal,tetra,norder,npts,width,rnge,wtkp,ebs,efermi, &
            sumev,wtkb,dosef,qval,ent,lfill)
       if (iprint()>= 20) then
          call bzwtsm(lswtk.eq.1.and.nspc.eq.2,nkp,nsp,nevx,wtkb,swtk, &
               amom)
          write(stdo,922) amom
       endif
       deallocate(ebs)
       if ( .NOT. quitvnow) call poppr
       if(quitvnow) exit
10  enddo

    ele1 = (zval+fmom)/2d0
    ele2 = (zval-fmom)/2d0
    sumev=sumev+ ele1*(vnow/2d0) - ele2*(vnow/2d0) !sumev correction takao
    if(iprint()>20) write(stdo,"(' bzwtsf: Set Bias field -Vup+Vdn=',f20.15)")vnow
    if (lswtk == 1 .AND. lwtkb == 1) then
       !        lwtkb = 2
       call rx('bzwtsf:111 tk think not used here')
    elseif (lswtk == 1 .AND. lwtkb == 2) then
       !        lwtkb = 1
       call rx('bzwtsf:222 tk think not used here')
    endif
    call tcx('bzwtsf')
  end subroutine bzwtsf


  subroutine bzwtsm(lswtk,nkp,nsp,nevx,wtkb,swtk,amom)
    use m_ftox
    !- Determine the magnetic moment, collinear or noncollinear case
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   lswtk :if true, swtk is used.  Otherwise, collinear case assumed:
    !i         :swtk(*,1,*) = 1  and swtk(*,2,*) = -1
    !i   nkp   :number of irreducible k-points (bzmesh.f)
    !i   nevx  :Maximum number of bands
    !i   wtkb  :band weights
    !i   swtk  :'spin weights': diagonal part of  (z)^-1 sigmz z
    !i         :where z are eigenvectors, sigma is the Pauli spin matrix
    !i         :Used when lswtk is set
    !o Outputs
    !o   amom  :magnetic moment
    logical :: lswtk
    integer :: nkp,nevx,nsp
    double precision :: wtkb(nevx,nsp,nkp),swtk(nevx,nsp,nkp),amom
    integer :: ikp,k
    double precision :: dsum
    if (nsp == 1) return
    if (lswtk) then
       amom = sum(wtkb(:,1,:)*swtk(:,1,:) + wtkb(:,2,:)*swtk(:,2,:))
    else
       amom = sum(wtkb(:,1,:) - wtkb(:,2,:))
    endif
  end subroutine bzwtsm

  subroutine fermi(qval,dos,ndos,emin,emax,nsp,eferm,e1,e2,dosef)
    !- Makes fermi energy from integrated density
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   qval, number of electrons to fermi level; dos(i) integrated
    !i   density at bin i; ndos, number of bins + 1; emin, emax, energy
    !i   window.
    !o Outputs
    !o   Eferm, Fermi energy; e1, e2, confidence limits on Fermi energy
    !o   i.e., Fermi energy lies between e1 and e2.
    !o   dosef:  density of states at fermi level
    !r Remarks
    !r   emin and e1 (and emax and e2) may point to the same address.
    !r   This version uses idos decomposed into spin up, down for nsp=2
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer :: ndos,nsp
    double precision :: qval,dos(ndos,1),emin,emax,eferm,e1,e2,dosef
    ! Local parameters
    integer :: i1,ie,i2
    double precision :: de,q1,q2,d1mach,wt
    ! External procedures

    ! --- Check bounds of DOS ---
    wt = 1d0/2
    i2 = 1
    if (nsp == 2) then
       wt = 1
       i2 = 2
    endif
    q1 = (dos(1,1)+dos(1,i2))*wt
    q2 = (dos(ndos,1)+dos(ndos,i2))*wt
    if (q1 > qval .OR. q2 < qval) call fexit3(-1,111, &
         ' Exit -1 FERMI: NOS ( %1,6;6g %1,6;6g ) does not '// &
         'encompass Q = %1;6g',q1,q2,qval)

    ! --- Find bin that boxes E_f ---
    de = (emax-emin)/(ndos-1)
    i1 = 1
    q1 = (qval + d1mach(3))/wt
    do  1  ie = 2, ndos
       if (dos(ie,1)+dos(ie,i2) > q1) goto 2
       i1 = ie
1   enddo
    call rx('bug in FERMI')
2   continue
    ! --- Linear interpolation for the Fermi level ---
    q1 = (dos(i1,1)+dos(i1,i2))*wt
    q2 = (dos(i1+1,1)+dos(i1+1,i2))*wt
    e1 = emin + de*(i1-1)
    e2 = e1 + de
    eferm = e1 + (qval-q1)/(q2-q1)*de
    dosef = (q2-q1)/de
  end subroutine fermi

end module m_bzwts
