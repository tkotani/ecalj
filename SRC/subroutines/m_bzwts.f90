module  m_bzwts ! BZ integration
  public bzwtsf2
  private
contains
  subroutine bzwtsf2(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, & !== FSMOMMETHOD=0 ogiginal version(modified version. fmom=0 is allowed.)==
       fmom,metal,tetra,norder,npts,width,rnge,wtkp,eb, & !- BZ integration for fermi level, band sum and qp weights, fixed-spin
       efermi,sumev,wtkb,qval,lfill,vmag,wtsf2) !,lwtkb lswtk, swtk,
    use m_lgunit,only:stdo
    use m_lmfinit,only:lso
    use m_ftox
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
    logical :: metal,tetra,wtsf2
    integer :: nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet, idtet(5,ntet)!,lswtk!,lwtkb
    double precision :: zval,fmom,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp), &
         wtkb(nevx,nsp,nkp),efermi,sumev,qval(2) !,swtk(nevx,nsp,nkp)
    integer :: ikp,ib,ipr,itmax,iter,iprint
    double precision :: amom,dosef(2),vhold(12),vmag,dvcap,dv,ef0,ent
    parameter (dvcap=.2d0,itmax=50)
    real(8):: ele1,ele2
    integer:: itermx
    logical:: agreemom
    real(8),parameter::    NULLR =-99999
    integer:: nmom1,nmom2
    real(8):: ehomo1,ehomo2,elumo1,elumo2
    real(8),allocatable:: ebs(:,:,:)
    logical:: quitvmag,lfill
    !!== Fermi level without spin constraint ==
    call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
         metal,tetra,norder,npts,width,rnge,wtkp,eb,efermi, &
         sumev,wtkb,dosef,qval,ent,lfill)
    if (nsp == 1) return
    call getpr(ipr)
    if( lso/=1 .AND. metal) then     !only for lso/=1
       amom = sum(wtkb(:,1,:) - wtkb(:,2,:))
       if(ipr >= 20) write(stdo,"(9x,'Mag. moment:',f15.6)") amom !magnetic moment
       qval(2) = amom
    else 
       return
    endif
    vmag = 0d0
    if(fmom==NULLR) return  
    !!== Setup for fixed-spin moment method ==
    call tcn('bzwtsf2')
    if(wtsf2) then
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
       vmag = (ehomo1+elumo1)/2d0 -(ehomo2+elumo2)/2d0
       efermi= ((ehomo1+elumo1)/2d0 +(ehomo2+elumo2)/2d0)/2d0 !/2d0 bug fix Jun26,2014 this only affects to the message.
       write(stdo,"('bzwtsf2: zval fmom nmon1 nmom2=',2f12.8,2x,2i3)")zval,fmom,nmom1,nmom2
       write(stdo,"('bzwtsf2: HOMOup LUMOup Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo1,elumo1,elumo1-ehomo1
       write(stdo,"('bzwtsf2: HOMOdn LUMOdn Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo2,elumo2,elumo2-ehomo2
       write(stdo,"('bzwtsf2: Set Bias initial cond. -Vup+Vdn=',f20.15)")vmag
    endif
    vhold=0d0
    ef0 = efermi
    if(ipr>0) write(stdo,*)' Seek potential shift for fixed-spin mom ...'
    itermx=100
    quitvmag=.false.
    do 10 iter=1,itermx !do loop for new guess at potential shift ==     ! bisection method takao
       if(.not.wtsf2) then
          amom = sum(wtkb(:,1,:) - wtkb(:,2,:)) !!=== Magnetic moment === 
          if(ipr>=41)write(stdo,ftox)&
               ' -Vup+Vdn=',ftof(vmag,8),'yields ef=',ftof(efermi),'amom=',ftof(amom),'when seeking',ftof(fmom)
          agreemom= abs(amom-fmom) < 1d-3
          call dvdos(vmag, amom,dosef(1),vhold,fmom,dvcap, dv)
          if(agreemom) vhold(12)=0        !      if (abs(dv) .lt. 1d-6) then
          quitvmag=.false.
          if (abs(dv) < 1d-6 .OR. agreemom) then           !       A root was found
             if (vhold(12) == -2 .OR. vhold(12) == -3 .OR. &
                  vhold(12) ==  0 .OR. vhold(12) ==  1) then
                if (ipr >= 10) write(stdo,ftox)' BZWTSF: potential shift bracketed.', &
                     'Unconstrained efermi=',ftof(ef0), &
                     'constraint fmom=',ftof(fmom),'actual mmom=',ftof(amom), &
                     'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vmag,8)
                quitvmag=.true.
             endif
          else if (iter == itmax) then
             if(ipr>=10)then
                write(stdo,ftox)' BZWTSF: failed to converge potential shift after',iter,'iterations.'
                write(stdo,ftox)' constraint fmom=',ftof(fmom),'actual amom=',ftof(amom), &
                     'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vmag,8)
             endif
             quitvmag=.true.
          endif
       endif
       !! Potential shift
       allocate(ebs(nevx,2,nkp))
       ebs(:,1,:) = eb(:,1,:) - vmag/2d0
       ebs(:,2,:) = eb(:,2,:) + vmag/2d0
       !! Fermi level with dv shift
       if( .NOT. quitvmag) call pshpr(ipr-50)
       if(iprint()>0) write(stdo,ftox) ' Second call bzwts in bzwtsf for fsmom mode'
       call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
            metal,tetra,norder,npts,width,rnge,wtkp,ebs,efermi, sumev,wtkb,dosef,qval,ent,lfill)
       if (iprint()>= 20) then
          amom = sum(wtkb(:,1,:) - wtkb(:,2,:)) 
          write(stdo,922) amom
       endif
       deallocate(ebs)
       if(.NOT.quitvmag) call poppr
       if(quitvmag) exit
       if(wtsf2) then
          amom = sum(wtkb(:,1,:) - wtkb(:,2,:)) !=== Magnetic moment ===  
          if(ipr>=41)write(stdo,ftox)' -Vup+Vdn=',ftof(vmag,8),'yields ef=',ftof(efermi), &
               'amom',ftof(amom),'when seeking mom=',ftof(fmom)
          agreemom = abs(amom-fmom) < 1d-6 
          call dvdos(vmag,amom,dosef(1),vhold,fmom,dvcap,dv)
          if(agreemom) vhold(12)=0
          quitvmag=.false.
          if(abs(dv) < 1d-6 .OR. agreemom) then
             quitvmag=.true.
             if(vhold(12) == -2 .OR. vhold(12) == -3 .OR. vhold(12) ==  0 .OR. vhold(12) ==  1) then
                if(ipr>=10) write(stdo,ftox)' BZWTSF2: potential shift bracketed. Unconstrained efermi=',ftof(ef0), &
                     'constraint fmom=',ftof(fmom),'actual mmom=',amom,'ef=',efermi,'-Vup+Vdn=',vmag
             endif
          elseif(iter == itmax) then
             quitvmag=.true.
             if(ipr>=10)write(stdo,ftox)' BZWTSF2: failed to converge potential shift after',iter,'iterations.', &
                  'constraint fmom=',ftof(fmom),'actual amom=',ftof(amom),'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vmag)
          endif
       endif
10  enddo
    ele1 = (zval+fmom)/2d0
    ele2 = (zval-fmom)/2d0
    sumev=sumev+ ele1*(vmag/2d0) - ele2*(vmag/2d0) !sumev correction takao
    if(iprint()>20) write(stdo,"(' bzwtsf2: Set Bias field -Vup+Vdn=',f20.15)")vmag
    call tcx('bzwtsf2')
922 format(9x,'Mag. moment:',f15.6)
  end subroutine bzwtsf2
  subroutine bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,& ! BZ integration for fermi level, band sum and qp weights
       metal,tetra,norder,npts,width,rnge,wtkp,eb, efermi,sumev,wtkb,dosef,qval,ent,lfill)
    use m_ftox
    use m_lgunit,only: stdo,stdl
    use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
    use m_MPItk,only: master_mpi
    use m_bzints,only:bzints
    implicit none
    intent(in)::   nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,&
         metal,tetra,norder,npts,width,rnge,wtkp!,eb
    intent(out)::                                  efermi,sumev,wtkb,dosef,qval,ent,lfill
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
    !o   lfill :true => insulator
    logical metal,tetra
    integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet(5,ntet)
    double precision zval,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),&
         wtkb(nevx,nsp,nkp),wtkbx(nevx,nsp,nkp),efermi,sumev,dosef(2),qval(2),ent
    integer:: it,itmax,n,nptdos,nspxx,nbmxx,nevxx,ib &
        ,ikp,ipr,job,nbpw,i1mach,nev, mkdlst,ifi,i,j,lry ,nulli,isw
    real(8) ,allocatable :: dos_rv(:,:), bot_rv(:), top_rv(:)
    integer ,allocatable :: bmap_iv(:)
    real(8) ,allocatable :: wk_rv(:)
    double precision emin,emax,e1,e2,dum(1),tol,e,elo,ehi,sumwt, dmin,dmax,egap,amom,cv,tRy
    character outs*100,ryy*3
    logical cmdopt0,lfill
    real(8) ,allocatable :: tlst_rv(:),eb2(:,:,:)
    
    real(8):: ebx(nevx*nsp,nkp)
    integer:: ibx(nevx*nsp,nkp), isx(nevx*nsp,nkp),ib1,ib2,ib2e,ix

    parameter (nulli=-99999)
    integer:: iprint, w(1)
    call tcn('bzwts')
    ipr=iprint()
    qval = 0
    ent = 0
    n = isign(1,norder) * mod(iabs(norder),100)
    allocate(bot_rv(nevx*nsp),top_rv(nevx*nsp))
    nspxx  = 1
    nevxx = nevx*nsp
    nbmxx = nbmx*nsp
    job = 3-2*nsp !     job = 1 for non spin pol, -1 for spin pol
    dosef(1) = 0
    egap = nulli
    if (nsp==2.and.nspc==1) then ! Force coupled spins: find range
       nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
       allocate(bmap_iv(nevx*nsp*nkp/nbpw+1))
       bmap_iv(:)=0
       allocate(wk_rv(nevx*nsp))
       allocate(eb2,source=eb)
       call ebcpl ( 0,nbmx,nevx,nsp,nspc,nkp,nbpw,bmap_iv, wk_rv,eb2 )
       lfill = efrng2 ( nspxx,nkp,nbmxx,nevxx,zval * 2,eb2, bot_rv,top_rv,elo,ehi,emin,emax )
       deallocate(eb2)
!       call ebcpl ( 1,nbmx,nevx,nsp,nspc,nkp,nbpw,bmap_iv,wk_rv,eb )
    else !     Spins not coupled: find range
       lfill = efrng2 ( nspxx,nkp,nbmxx,nevxx,nspc * zval,eb, bot_rv,top_rv,elo,ehi,emin,emax )
    endif
    ! ... Bands never filled if 100s digit norder set
    if (.not. tetra .and. iabs(norder) .ge. 100) lfill = .false.
    if (allocated(wk_rv)) deallocate(wk_rv)
    if (allocated(bmap_iv)) deallocate(bmap_iv)
    if (allocated(top_rv)) deallocate(top_rv)
    if (allocated(bot_rv)) deallocate(bot_rv)
    if (lfill) then ! ... Case an insulator: put efermi at emin + tiny number
       efermi = emin + 1d-10
    elseif(.not. metal) then ! ... Do the best we can should be a metal, but assumption that it isn't
       efermi = (emin + emax) / 2
    endif
    if (nsp==2 .and. nspc==1 ) then ! ... Pretend as though spin-pol bands are coupled to find E_f
       nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
       allocate(bmap_iv(nevx*nsp*nkp/nbpw+1),source=0)
       allocate(wk_rv(nevx*nsp))
       allocate(eb2,source=eb)
       
       do ikp=1,nkp
          ix = 0
          ib2e=0
          do ib1=1,nevx
             e1 = eb(ib1,1,ikp)
             do ib2= ib2e+1,nevx
                e2 = eb(ib2,2,ikp)
                if(e2>e1) exit
                ix=ix+1
                ebx(ix,ikp)= e2
                ibx(ix,ikp)=ib2
                isx(ix,ikp)=2
                ib2e=ib2
             enddo
             ix=ix+1
             ebx(ix,ikp)= e1
             ibx(ix,ikp)=ib1
             isx(ix,ikp)=1
          enddo
       enddo
       
       call ebcpl ( 0,nbmx,nevx,nsp,nspc,nkp,nbpw,bmap_iv,wk_rv,eb)
    endif
    
    ! --- BZ weights, sumev and E_f for an insulator  ---
    if(.not. metal ) then
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
       if(ipr>=20) then
          write(stdo,ftox)' BZWTS : --- Non-metal sampling ---'
          write(stdo,ftox)' Fermi energy:',ftof(efermi),'electrons: num=',ftof(sumwt),'occ. bands =',ftof(sumev)
          write(stdo,ftox)' VBmax=',ftof(emin),'CBmin=',ftof(emax),'gap=',ftof(emax-emin), 'Ry =',ftof((emax-emin)*13.6058d0),'eV'
       endif
    elseif(tetra) then ! --- BZ weights, sumev and E_f by tetrahedron method (Blochl wts) ---
       if(ipr>=30) write(stdo,ftox)' bzwts: --- Tetrahedron Integration ---'
       if(lfill) then
          egap = emax-emin
          if(ipr>=30)then
             write(stdo,ftox)' ... only filled or empty bands encountered: ev=',ftof(emin),'ec=',ftof(emax)
             write(stdo,ftox)' VBmax=',ftof(emin),'CBmin=',ftof(emax),'gap =',&
                  ftof(emax-emin),'Ry = ',ftof((emax-emin)*13.6058d0),'eV'
          endif
       else
          nptdos = 101
          allocate(dos_rv(nptdos,1))
          tol = 1d-6
          !  Preliminary check that dos lies within emin,emax.  Widen emin,emax if not
          call bzints(n1*n2*n3,eb,dum,nkp,nevxx,nbmxx,nspxx,emin,emax,dos_rv,nptdos,efermi,job,ntet,idtet,sumev,qval(1) )
          dmin = sum(dos_rv(1,1:nspxx))/nspxx
          dmax = sum(dos_rv(nptdos,1:nspxx))/nspxx
          if (dmin .gt. zval) then
             emin = 3*emin-2*emax
             write(stdo,"(' (warning): initial NOS ( ',d14.6,x,d14.6,' ) does not encompass Q=',f14.6)") dmin,dmax,zval
          elseif (dmax .lt. zval) then
             emax = 3*emax-2*emin
             write(stdo,"(' (warning): initial NOS ( ',d14.6,x,d14.6,' ) does not encompass Q=',f14.6)") dmin,dmax,zval
          endif
          if(ipr>=35) write(stdo,"(9x,'Est E_f ',10x,'Window',8x,'Tolerance',2x,'n(E_f)')")
          itmax = 5
          do it = 1, itmax
             call bzints(n1*n2*n3,eb,dum,nkp,nevxx,nbmxx,nspxx,emin,emax,dos_rv,nptdos,efermi,job,ntet,idtet,sumev,qval(1) )
             call fermi ( zval,dos_rv,nptdos,emin,emax,efermi,emin,emax,dosef(1) )
             if(ipr>=35)  write(stdo,"(7x,6(f10.6,1x))") efermi,emin,emax,emax-emin,dosef(1)
             if(emax-emin .lt. tol) goto 111
          enddo
          if(ipr>10) write(stdo,ftox)' BZWTS (warning): Fermi energy not converged: ',ftof(emax-emin),' > tol=',ftof(tol)
111       continue
          deallocate(dos_rv)
       endif
       call bzints(n1*n2*n3,eb,wtkbx,nkp,nevxx,nbmxx, nspxx,emin,emin,dum,1,efermi,2*job,ntet,idtet,sumev,qval(1))
    else ! --- BZ weights, sumev and E_f by Methfessel-Paxton sampling ---
       call rx('this branch is not supported. See before 2024-5-12. Methfessel-Paxton sampling')
    endif
    amom = 0d0 ! ... Magnetic moment
!    eb=eb2
    if (nsp==2.and.nspc==1 ) then ! ... Restore to uncoupled bands; ditto with weights
       eb=eb2
       deallocate(eb2)
       !call           ebcpl(1,nbmx,nevx,nsp,nspc, nkp,nbpw,bmap_iv,wk_rv,eb )
       !if(metal) call ebcpl(1,nevx,nevx,nsp,nspc, nkp,nbpw,bmap_iv,wk_rv,wtkbx )
       if(metal) then
          do ikp=1,nkp
          do ix=1,nevxx
             wtkb(ibx(ix,ikp),isx(ix,ikp),ikp)= wtkbx(ix,1,ikp)
          enddo
          enddo
       endif
       
       if(allocated(tlst_rv)) deallocate(tlst_rv)
       if(allocated(wk_rv)) deallocate(wk_rv)
       if(allocated(bmap_iv)) deallocate(bmap_iv)
       if(metal) amom = sum(wtkb(1:nevx,1,1:nkp)- wtkb(1:nevx,2,1:nkp))
    else
       wtkb=wtkbx
    endif
    
    
    qval(2) = amom
    if(ipr>0)write(stdl,ftox)'bzmet',metal,'tet',tetra,'ef',ftof(efermi),'sev',ftof(sumev),'zval',ftof(zval)
    if(ipr>0)write(stdl,ftox)'qval',ftof(qval(1)),'amom',ftof(amom),'egap(eV)',ftof(egap,3)
    e = efermi
    if(.not. lfill .and. .not. tetra) e = efermi + rnge*width/2
    call tcx('bzwts')
  end subroutine bzwts
  
  subroutine ebcpl(mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap,wk,eb) !- Gather spin-polarized bands into a single group, or redistribute
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
    implicit none
    integer mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap(1)
    double precision eb(nbmx,nsp,nq),wk(nbmx*nsp)
    integer ib,iq,ib1,ib2,iqb
    if (nsp .eq. 1 .or. nspc .eq. 2) return
    ! --- Gather bands at each qp ---
    ifx:if (mode .eq. 0) then
       iqb = 0
       do10:do    iq = 1, nq          !   ... Gather and order +,- bands at this qp into one column
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
          enddo do20
          call dcopy(nevx*nsp,wk,1,eb(1,1,iq),1)
          if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
       enddo do10
    endif ifx
    ! --- Disperse bands at each qp ---
    if (mode .eq. 1) then
       iqb = 0
       d110  : do  iq = 1, nq   !   ... Disperse bands into +,- for this qp according to bmap
          ib1 = 1
          ib2 = 1
          call dcopy(nevx*nsp,eb(1,1,iq),1,wk,1)
          do120  :     do  ib = 1, nevx*nsp
             iqb = iqb+1
             if (iget(bmap,nbpw,iqb) .eq. 0) then
                eb(ib1,1,iq) = wk(ib)
                ib1 = ib1+1
             else
                eb(ib2,2,iq) = wk(ib)
                ib2 = ib2+1
             endif
          enddo do120
          if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
       enddo d110
    endif
  end subroutine ebcpl
  subroutine mark1(bitmap, nbpw, n)    !- put a one in the nth bit of bitmap.
    !i Inputs   !i   bitmap, n
    !r    nbpw: a number no larger than the number of bits per integer word
    !     implicit none
    integer bitmap(*), nbpw, n
    ! Local parameters
    integer nword,nbit,i
    nword = (n-1)/nbpw
    nbit = mod(n-1,nbpw)
    i = 2**(nbpw-nbit-1)
    bitmap(nword+1) = bitmap(nword+1) + i*(1-mod(bitmap(nword+1)/i,2))
  end subroutine mark1
  integer function iget(bitmap, nbpw, n)    !- Return 0 or 1, depending on the value of the nth bit of bitmap
    implicit none
    integer bitmap(*), nbpw, n
    integer nword,nbit
    nword = (n-1)/nbpw
    nbit = mod(n-1,nbpw)
    iget = mod(bitmap(nword+1)/2**(nbpw-nbit-1),2)
  end function iget
  logical function efrng2(nsp,nkp,nbmax,nband,zval,eband,ebbot,ebtop,elo,ehi,e1,e2)  !- Find range of Fermi energy.
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
  
  subroutine fermi(qval,dosi,ndos,emin,emax, eferm,e1,e2,dosef) !Makes fermi energy from integrated density
    use m_ftox
    use m_nvfortran
    !i   qval:    number of electrons to fermi level
    !i   dosi(i): integrated density at bin i;
    !i   ndos: number of bins + 1
    !i   emin, emax: energy window.
    !o Outputs
    !o   Eferm, Fermi energy;
    !o   e1<=Eferm<=e2 : confidence limits on Fermi energy. i.e., Fermi energy lies between e1 and e2.
    !o   dosef:  density of states at fermi level
    !remark: emin and e1 (and emax and e2) may point to the same address.
    implicit none
    integer :: ndos,i1,ie
    double precision :: qval,dosi(ndos),emin,emax,eferm,e1,e2,dosef,de,q1,q2
    if(dosi(1)>qval.OR.dosi(ndos)<qval) call rx('fermi does not encompass qval='//ftof(qval))
    i1 = findloc([(dosi(ie)>qval, ie=1,ndos)],value=.true.,dim=1) - 1
    de = (emax-emin)/(ndos-1)
    e1 = emin + de*(i1-1)
    e2 = emin + de*i1
    q1 = dosi(i1)   
    q2 = dosi(i1+1) 
    eferm = e1 + (qval-q1)/(q2-q1)*de ! Linear interpolation for the Fermi level
    dosef = (q2-q1)/de
  end subroutine fermi
 
  subroutine dvdos(vmag,nosnow,dosnow,vhold,ztarg,dvcap,dv)  !- Estimate shift in potential shift to meet target number-of-states
    !i Inputs
    !i   vmag  :current value of potential shift
    !i   nosnow:current value of charge
    !i   dosnow:current density of states, d nosnow / dv
    !i         :(used only for first iteration to estimate rfalsi step size)
    !i   vhold :vector of 12 numbers, maintained internally by rfalsi
    !i         :Initial call: vhold should be zero.
    !i   ztarg :desired charge
    !i   dvcap :maximum change in potential shift for any step
    !o Outputs
    !o   vmag  :updated value of estimated potential shift
    !o   dv    :change in vmag this step
    !r Remarks
    !r   Routine uses regula falsi to iteratively find target nosnow.
    implicit none
    double precision :: vmag,nosnow,dosnow,vhold(12),ztarg, dvcap,dv, dznow,dxmx
    integer :: ir
    ! ... First order estimate dv = (ztarg-zhave)/slope
    dznow = nosnow-ztarg
    dv = dznow/max(dosnow,1d-5)
    ir = nint(vhold(12)) !  Hang onto vhold(12) because rfalsi destroys it
    if (ir == 0) then
       dxmx = dznow/max(dosnow,1d-5)
       if (abs(dxmx) > dvcap) dxmx = sign(dvcap,dxmx)
       ir = 0
    else
       dxmx = dvcap
    endif
    call pshpr(0)
    call rfalsi(vmag,dznow,5d-8,0d0,5d-8,dxmx,10,vhold(1),ir)
    call poppr
    vhold(12) = ir
    dv = vmag - vhold(1)
  end subroutine dvdos
end module m_bzwts
