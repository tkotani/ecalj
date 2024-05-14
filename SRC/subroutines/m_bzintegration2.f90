!> BZ integration for fermi level, band sum and qp weights, fixed-spin
module  m_bzintegration2 ! BZ integration
  use m_nvfortran
  use m_ftox
  use m_lgunit,only: stdo,stdl
  use m_MPItk,only: master_mpi
  public bzintegration2
  private
contains
  subroutine bzintegration2(eb, efermi,sumev,wtkb,sumqv,vmag)  
    use m_lmfinit,only: lso,bz_fsmommethod,fmom=>bz_fsmom,norder=>bz_n,lmet=>bz_lmet
    use m_lmfinit,only: width=>bz_w,npts=>bz_ndos,nsp,nspc, zbak,NULLR
    use m_mkqp,only: ntet=> bz_ntet, nkabc=> bz_nabc,idtet=>iv_a_oidtet, wtkp=>rv_a_owtkp
    use m_suham,only: nbmx=>ham_ndham 
    use m_qplist,only: nkp
    use m_mkpot,only: qval
    !i   nbmx  : leading dimension of eb
    !i   nsp   :=2 for spin-polarized case, otherwise 1
    !i   nspc  :=2 for lso=1 (spin-up and spin-down channels are coupled), otherwise 1
    !i   n1..n3:number of divisions for the k-point mesh
    !i   nkp   :number of inequivalent k-points (bzmesh.f)
    !i   ntet  :number of inequivalent tetrahedra (tetirr.f)
    !i   idtet :idtet(1..4,i) points to the 4 irreducible k-points defining
    !i         :corners of tetrahedron;
    !i         :idtet(0,i) number of tetrahedra of the i'th kind
    !i   zval = qval-zbak  :valence charge
    !i   fmom  :fixed spin moment.  Even zero is allowed.
    !i   metal :T => metal, F => nonmetal
    !i   tetra :T => tetrahedron integration
    !
    !i   norder,npts,width,rnge: parameters for sampling integr. (maknos)
    !i   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
    !i   eb    :energy bands; alias eband
    !o Outputs
    !o   efermi:Fermi energy
    !o   sumev :sum of eigenvalues
    !o   wtkb  :integration weights (not generated for nonmetal case)
    !o   sumqv  :sumqv(1) = total charge; sumqv(2) = magnetic moment
    !o   vmag: magnetic field to keep give magnetic moment for fsmom (fixed moment) method.
    implicit none
    logical :: metal,tetra,wtsf2, agreemom,quitvmag,lfill
    integer :: nevx,n1,n2,n3 ,ikp,ib,ipr,iter,iprint,itermx,nmom1,nmom2
    real(8) :: zval,eb(nbmx,nsp,nkp),wtkb(nbmx,nsp,nkp),efermi,sumev,sumqv(2) 
    real(8) :: amom,dosef(2),vhold(12),vmag,dv,ef0,ent,ele1,ele2, ehomo1,ehomo2,elumo1,elumo2,rnge
    real(8),parameter:: dvcap=.2d0
    integer,parameter:: itmax=50
    real(8),allocatable:: ebs(:,:,:)
    metal=lmet.ne.0
    wtsf2= bz_fsmommethod==1
    nevx = nbmx
    sumev=0d0
    tetra = ntet>0
    rnge = 8
    if(norder<0) rnge = 16
    n1=nkabc(1)
    n2=nkabc(2)
    n3=nkabc(3)
    zval=qval-zbak
    !!== Fermi level without spin constraint ==
    call bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval, &
         metal,tetra,norder,npts,width,rnge,wtkp,eb,efermi, &
         sumev,wtkb,dosef,sumqv,ent,lfill)
    if (nsp == 1) return
    call getpr(ipr)
    if( lso/=1 .AND. metal) then     !only for lso/=1
       amom = sum(wtkb(:,1,:) - wtkb(:,2,:))
       if(ipr >= 20) write(stdo,"(9x,'Mag. moment:',f15.6)") amom !magnetic moment
       sumqv(2) = amom
    else 
       return
    endif
    vmag = 0d0
    if(fmom==NULLR) return  
    !!== Setup for fixed-spin moment method ==
    call tcn('bzintegration2')
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
       write(stdo,"('bzintegration2: zval fmom nmon1 nmom2=',2f12.8,2x,2i3)")zval,fmom,nmom1,nmom2
       write(stdo,"('bzintegration2: HOMOup LUMOup Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo1,elumo1,elumo1-ehomo1
       write(stdo,"('bzintegration2: HOMOdn LUMOdn Diff=',3f20.15,' (Diff=0.5forNoOccupied)')")ehomo2,elumo2,elumo2-ehomo2
       write(stdo,"('bzintegration2: Set Bias initial cond. -Vup+Vdn=',f20.15)")vmag
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
            metal,tetra,norder,npts,width,rnge,wtkp,ebs,efermi, sumev,wtkb,dosef,sumqv,ent,lfill)
       if (iprint()>= 20) then
          amom = sum(wtkb(:,1,:) - wtkb(:,2,:)) 
          write(stdo,"(9x,'Mag. moment:',f15.6)") amom
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
                if(ipr>=10) write(stdo,ftox)' BZINTEGRATION2: potential shift bracketed. Unconstrained efermi=',ftof(ef0), &
                     'constraint fmom=',ftof(fmom),'actual mmom=',amom,'ef=',efermi,'-Vup+Vdn=',vmag
             endif
          elseif(iter == itmax) then
             quitvmag=.true.
             if(ipr>=10)write(stdo,ftox)' BZINTEGRATION2: failed to converge potential shift after',iter,'iterations.', &
                  'constraint fmom=',ftof(fmom),'actual amom=',ftof(amom),'ef=',ftof(efermi),'-Vup+Vdn=',ftof(vmag)
          endif
       endif
10  enddo
    ele1 = (zval+fmom)/2d0
    ele2 = (zval-fmom)/2d0
    sumev=sumev+ ele1*(vmag/2d0) - ele2*(vmag/2d0) !sumev correction takao
    if(iprint()>20) write(stdo,"(' bzintegration2: Set Bias field -Vup+Vdn=',f20.15)")vmag
    call tcx('bzintegration2')
  end subroutine bzintegration2
  subroutine bzwts(nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,& ! BZ integration for fermi level, band sum and qp weights
       metal,tetra,norder,npts,width,rnge,wtkp,eb, efermi,sumev,wtkb,dosef,sumqv,ent,lfill)
    use m_bzints,only:bzints
    implicit none
    intent(in)::   nbmx,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet,zval,&
         metal,tetra,norder,npts,width,rnge,wtkp!,eb
    intent(out)::                                  efermi,sumev,wtkb,dosef,sumqv,ent,lfill
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
    !o   sumqv  :sumqv(1) = total charge; sumqv(2) = magnetic moment
    !o   ent   :entropy term (actually TS)
    !o   lfill :true => insulator
    logical metal,tetra
    integer nbmx,norder,npts,nevx,nsp,nspc,n1,n2,n3,nkp,ntet,idtet(5,ntet)
    real(8)::zval,eb(nbmx,nsp,nkp),width,rnge,wtkp(nkp),&
         wtkb(nevx,nsp,nkp),efermi,sumev,dosef(2),sumqv(2),ent ,wtkbx(nevx,nsp,nkp),wtkb2(nevx,nsp,nkp)
    integer:: it,itmax,n,nptdos,nspxx,nbmxx,nevxx,ib &
         ,ikp,ipr,job,i1mach,nev, mkdlst,ifi,i,j,lry ,nulli,isw !,nbpw
    real(8) ,allocatable :: dos_rv(:)
    real(8) emin,emax,e1,e2,dum(1),tol,e,elo,ehi,sumwt, dmin,dmax,egap,amom,cv,tRy
    character outs*100,ryy*3
    logical cmdopt0,lfill
    real(8) ,allocatable :: tlst_rv(:),eb2(:,:,:)
    real(8):: ebx(nevx*nsp,nkp),de,q1,q2,bot_rv(nevx*nsp),top_rv(nevx*nsp)
    integer:: ibx(nevx*nsp,nkp), isx(nevx*nsp,nkp),ib1,ib2,ib2e,ix,isp,i1,ie
    parameter (nulli=-99999)
    integer:: iprint, w(1)
    call tcn('bzwts')
    ipr=iprint()
    sumqv = 0
    ent = 0
    n = isign(1,norder) * mod(iabs(norder),100)
    nspxx  = 1
    nevxx = nevx*nsp
    nbmxx = nbmx*nsp
    job = 3-2*nsp !     job = 1 for non spin pol, -1 for spin pol
    dosef(1) = 0
    egap = nulli
    if (nsp==2 .and. nspc==1 ) then ! ... merge up-dn bands
       do ikp=1,nkp
          ix = 0
          ib2e=0
          do ib1=1,nevx+1
             e1 = merge(eb(ib1,1,ikp),9999d99,ib1<=nevx)
             do ib2= ib2e+1,nevx
                e2 = eb(ib2,2,ikp)
                if(e2>e1) exit
                ix=ix+1
                ebx(ix,ikp)= e2
                ibx(ix,ikp)=ib2
                isx(ix,ikp)=2
                ib2e=ib2
             enddo
             if(ib1==nevx+1) exit
             ix=ix+1
             ebx(ix,ikp)= e1
             ibx(ix,ikp)=ib1
             isx(ix,ikp)=1
          enddo
       enddo
    else
       ebx=reshape(eb,shape(ebx))
    endif
    lfill = efrng2 ( nspxx,nkp,nbmxx,nevxx, nsp*zval, ebx, bot_rv,top_rv,elo,ehi,emin,emax )
    if (.not. tetra .and. iabs(norder) .ge. 100) lfill = .false.
    if (lfill) then ! ... Case an insulator: put efermi at emin + tiny number
       efermi = emin + 1d-10
    elseif(.not. metal) then ! ... Do the best we can should be a metal, but assumption that it isn't
       efermi = (emin + emax) / 2
    endif
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
             e = ebx(ib,ikp)
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
          allocate(dos_rv(nptdos))
          tol = 1d-6
          !  Preliminary check that dos lies within emin,emax.  Widen emin,emax if not
          call bzints(n1*n2*n3,ebx,dum,nkp,nevxx,nbmxx,nspxx,emin,emax,dos_rv,nptdos,efermi,job,ntet,idtet,sumev,sumqv(1) )
          dmin = dos_rv(1)
          dmax = dos_rv(nptdos)
          if (dmin .gt. zval) then
             emin = 3*emin-2*emax
             write(stdo,"(' (warning): initial NOS ( ',d14.6,x,d14.6,' ) does not encompass Q=',f14.6)") dmin,dmax,zval
          elseif (dmax .lt. zval) then
             emax = 3*emax-2*emin
             write(stdo,"(' (warning): initial NOS ( ',d14.6,x,d14.6,' ) does not encompass Q=',f14.6)") dmin,dmax,zval
          endif
          if(ipr>=35) write(stdo,"(9x,'Est E_f ',10x,'Window',8x,'Tolerance',2x,'n(E_f)')")
          itmax = 5
          GetFermienergy:do it = 1, itmax
             call bzints(n1*n2*n3,ebx,dum,nkp,nevxx,nbmxx,nspxx,emin,emax,dos_rv,nptdos,efermi,job,ntet,idtet,sumev,sumqv(1) )
             !   !i   qvalx:    number of electrons to fermi level
             !   !i   dosi(i): integrated density at bin i;
             !   !i   ndos: number of bins + 1
             !   !i   emin, emax: energy window.
             !   !o   Eferm, Fermi energy;
             !   !o   e1<=Eferm<=e2 : confidence limits on Fermi energy. i.e., Fermi energy lies between e1 and e2.
             !   !o   dosef:  density of states at fermi level
             associate(dosi=>dos_rv(:),qvalx=>zval,ndos=>nptdos) !Makes fermi energy from integrated density
               if(dosi(1)>qvalx.OR.dosi(ndos)<qvalx) call rx('fermi does not encompass qval='//ftof(qvalx))
               i1 = findloc([(dosi(ie)>qvalx, ie=1,ndos)],value=.true.,dim=1) - 1
               de = (emax-emin)/(ndos-1)
               e1 = emin + de*(i1-1)
               e2 = emin + de*i1
               q1 = dosi(i1)   
               q2 = dosi(i1+1) 
               efermi = e1 + (qvalx-q1)/(q2-q1)*de ! Linear interpolation for the Fermi level
               dosef(1) = (q2-q1)/de
               emin=e1
               emax=e2
             endassociate
             if(ipr>=35)  write(stdo,"(7x,6(f10.6,1x))") efermi,emin,emax,emax-emin,dosef(1)
             if(emax-emin .lt. tol) goto 111
          enddo GetFermienergy
          write(stdo,ftox)' BZWTS (warning): Fermi energy not converged: ',ftof(emax-emin),' > tol=',ftof(tol)
111       continue
          deallocate(dos_rv)
       endif
       call bzints(n1*n2*n3,ebx,wtkbx,nkp,nevxx,nbmxx, nspxx,emin,emin,dum,1,efermi,2*job,ntet,idtet,sumev,sumqv(1))
    else ! --- BZ weights, sumev and E_f by Methfessel-Paxton sampling --- not maintained well...
       if(ipr>0) write(stdo,"(a,i0,a,f15.6)")' BZWTS : --- Brillouin Zone sampling; N=',n,' W=',width
       if(nsp==2) call dscal(nkp,.5d0,wtkp,1) !   ... Temporarily remove spin degeneracy if spins are coupled
       if((.not. lfill) .or. (metal .and. (nkp .eq. 1))) then !   ... Find Fermi level, sampling
          e1 = elo - rnge*width/2
          e2 = ehi + rnge*width/2
          efermi = 0.5d0*(e1 + e2)
          itmax = 1000
          do it = 1, itmax
             call pshpr(0)
             call splwts(nkp,nevxx,nbmxx,nspxx,wtkp,ebx,n,width,efermi, .true.,sumev,wtkbx,sumqv(1),ent,dosef(1),cv)
             call poppr
             if (dabs(zval-sumqv(1))<1d-12) then
                if(ipr>0)write(stdo,ftox)' Fermi energy, ',ftof(efermi),' found after ',it,' bisections,',ftof(sumqv(1)),&
                     ' electrons, DOS(E_f)=',ftof(dosef(1))
                goto 333
             endif
             if (sumqv(1) > zval) then
                e2 = efermi
             else
                e1 = efermi
             endif
             efermi = 0.5d0*(e1 + e2)
          enddo
          write(stdo,ftox)' BZWTS (warning): cannot find E_F by bisection, using INTNOS'
          allocate(dos_rv(npts))
          emin = elo - rnge*width/2
          emax = emax + rnge*width/2
          call maknos ( nkp,nevxx,nbmxx,nspxx,wtkp,ebx,n,width,- rnge,emin,emax,npts,dos_rv )
          call intnos ( npts,dos_rv,emin,emax,zval,efermi,dosef(1),sumev )
          deallocate(dos_rv)
333       continue
       else
          dosef = 0d0
          egap = emax-emin
          if(ipr>0) write(stdo,ftox)' ... only filled or empty bands encountered:  ev=',ftof(emin),' ec=',ftof(emax)
          if(ipr>0) write(stdo,ftox)' VBmax = ',ftof(emin),' CBmin = ',ftof(emax),' gap = ',&
               ftof(emax-emin),'Ry = ',ftof((emax-emin)*13.6058d0,3),'eV'
       endif
       if(master_mpi.and.(cmdopt0('--cvK:') .and. n<0 .and. metal)) then ! ... (optional) Tabulate specific heat in file for list of T's
          lRy = 0
          ryy='K'
          itmax=8
          allocate(tlst_rv(itmax))
          tlst_rv=[10,20,40,80,160,320,640,1280] !fixed now
          write(stdo,ftox)'Writing CV(T) to file for ',itmax,'vals of T:',ftof(tlst_rv,1),trim(ryy) 
          write(stdo,ftox)'% rows ',it,' cols 4 #   T(K)    T(Ry)   S(k_B)   TdS/dT(k_B)'
          do  it = 1, itmax
             tRy = tlst_rv(it)/0.1579d6
             call pshpr(1)
             call splwts(nkp,nevxx,nbmxx,nspxx,wtkp,ebx,n,tRy,efermi, metal,sumev,wtkbx,sumqv(1),ent,dosef(1),cv)
             call poppr
             write(stdo,ftox)ftof(0.1579d6*tRy),ftof(tRy),ftof(ent),ftof(cv)
          enddo
       endif
       call splwts(nkp,nevxx,nbmxx,nspxx,wtkp,ebx,n,width,efermi,& !   ... Make weights, sampling
            (.not. lfill) .or. (metal .and. (nkp .eq. 1)), sumev,wtkbx,sumqv(1),ent,dosef(1),cv)
       if(nsp==2) call dscal(nkp,2d0,wtkp,1)
    endif
    amom = 0d0 ! ... Magnetic moment
    if(metal) then
       if (nsp==2.and.nspc==1 ) then ! ... Restore to uncoupled bands; ditto with weights
          forall(ikp=1:nkp,ix=1:nevxx) wtkb(ibx(ix,ikp),isx(ix,ikp),ikp)= wtkbx(ix,1,ikp)
          amom = sum(wtkb(1:nevx,1,1:nkp)- wtkb(1:nevx,2,1:nkp))
       else
          wtkb=wtkbx !call dcopy(nevx*nsp*nkp,wtkbx,1,wtkb,1) 
       endif
    endif
    sumqv(2) = amom
    if(ipr>0) write(stdl,ftox)'bzmet',metal,'tet',tetra,'ef',ftof(efermi),'sev',ftof(sumev),'zval',ftof(zval)
    if(ipr>0) write(stdl,ftox)'sumqv',ftof(sumqv(1)),'amom',ftof(amom),'egap(eV)',ftof(egap,3)
    e = efermi
    if(.not. lfill .and. .not. tetra) e = efermi + rnge*width/2
    call tcx('bzwts')
  end subroutine bzwts
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
    real(8) zval,e1,e2,eband(nbmax,nsp,nkp),ebbot(nband,nsp),ebtop(nband,nsp),elo,ehi
    ! Local parameters
    real(8) xx,d1mach,enull
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
    real(8) array(n)
    integer i,j,k,inc
    real(8) v
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
    real(8) :: vmag,nosnow,dosnow,vhold(12),ztarg, dvcap,dv, dznow,dxmx
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

!!!!!! Folloings are for Methfessel-Paxton branch !!!!
  subroutine intnos(ndos,dos,emin,emax,qval,efermi,dosef,eband)!Finds E_F from a tabulated number of states function
    !i  Input
    !i    ndos : number of tabulated points; dos : integrated DOS
    !i    emin, emax : energy range of tabulation;
    !i    qval : number of valence electrons
    !o  Output
    !o    efermi : Fermi energy, dosef : DOS at E_f
    !-----------------------------------------------------------------------
    !     implicit none
    integer :: ndos
    real(8) :: dos(0:ndos-1),emin,emax,qval,efermi,eband,dosef
    integer :: i,meshpt,iprint
    real(8) :: step,sum,q,q1,q2,e1,eps,d1mach
    eps = d1mach(3)
    ! --- make Fermi energy ---
    step = (emax - emin) / (ndos - 1)
    meshpt = 0
    q = qval + eps
    do  1  i = 1, ndos-1
       if ( dos(i) >= q ) goto 2
       meshpt = i
1   enddo
2   continue
    if (meshpt == ndos-1) &
         call rx('INTNOS : Fermi energy lies above emax')
    ! E_F lies between mesh points meshpt and meshpt+1 -- interpolate :
    q1 = dos(meshpt)
    q2 = dos(meshpt+1)
    e1 = emin + step * meshpt
    efermi = e1 + ( qval-q1 ) / ( q2-q1 ) * step
    dosef = (q2 - dos(meshpt-1)) / (2*step)
    ! --- make band energy by partial integration ---
    sum = .5d0 * q1
    do  3  i = 1, meshpt-1
       sum = sum + dos(i)
3   enddo
    sum = sum * step
    sum = sum + .5d0 * (efermi - e1) * (qval + q1)
    eband = efermi * qval - sum
    if (iprint() >= 30) then
       !        do  12  i = 1, 1
       write(stdo,"(' INTNOS: Fermi energy=,f15.8', &
            '  band energy=',f15.8,' DOS(E_f)=',f15.8)") efermi, eband, dosef
       !   12   continue
       !        write(*,10) efermi, eband, dosef
       !        write(fopn('LOG'),10) efermi, eband, dosef
       !   10 format(' INTNOS: Fermi energy =',f10.6,'; band energy =',f11.6/
       !     .       '          DOS(E_f) =',f11.6)
    endif
  end subroutine intnos
  subroutine maknos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax, ndos,dos)!- Make density of states from bands
    !i  Input
    !i    nqp : number of q-points; nband : number of bands;
    !i    nsp : 2 for spin polarised bands, 1 otherwise;
    !i    wgts, evl : weights and bands (eigenvalues);
    !i    nbmx : first dimension of evl ;
    !i    n   : n>0 Methfessel-Paxton polynomial order
    !i        : n<0 sampling done with Fermi-Dirac statistics
    !i    w   : n>0 gaussian width in Methfessel-Paxton integration (Ry)
    !i        : n<0 Temperature for Fermi distribution (Ry)
    !i    tol : allowed error in DOS due to truncating the gaussian,
    !i          if negative on entry, range is set to -tol*W
    !i    emin, emax, ndos; energy range and number of energy mesh points
    !o  Ouput
    !o    dos: integrated DOS
    !u Updates
    !u   2 Nov 1995 (JEK) returns spin-polarized integrated dos
    !     implicit none
    integer :: nqp,nband,nbmx,nsp,n,ndos
    real(8) :: wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp), &
         w,emin,emax,tol,wt,emesh
    integer :: i,isp,iband,iq,meshpt,mesh1,mesh2,mrange,iprint,i1mach
    real(8) :: e,x,range,test,step,d,s,xx
    !  external delstp
    mrange=-9999999
    call dpzero(dos,nsp*ndos)
    step = (emax - emin) / (ndos - 1)
    if ( tol > 0d0 ) then
       do  2  i = 0, ndos-1
          x = i * step / w
          call delstp(n,x,test,s,xx)
          if ( test < tol ) then
             mrange = i + 1
             goto 3
          endif
2      enddo
       if (iprint() > 30) print *,'maknos (warning) : tol too small'
3      continue
       range = 2 * mrange * step
       test = tol
    else
       range = -tol * w
       mrange = range / ( 2 * step )
       call delstp(n,-tol/2,test,s,xx)
    endif
    if (iprint() >= 40) write(stdo,ftox) ' MAKNOS: range=',ftof(range/w), &
         ' (',2*mrange,'bins) DOS error estimate=',ftof(test),'per state'
    do  7  iq = 1, nqp
       wt = abs(wgts(iq)) / nsp
       do  61  iband = 1, nband
          do  6  isp = 1, nsp
             e = evl(iband,isp,iq)
             meshpt = (e - emin) / step
             mesh1 = meshpt - mrange
             mesh2 = meshpt + mrange
             if (mesh2 >= ndos) mesh2 = ndos-1
             call rxx(mesh1 .lt. 0,'MAKNOS: emin too large')
             do  4  meshpt = mesh1, mesh2
                emesh = emin + meshpt * step
                x = (emesh - e) / w
                call delstp(n,x,d,s,xx)
                dos(meshpt,isp) = dos(meshpt,isp) + wt * (1d0 - s)
4            enddo
             do  5  meshpt = mesh2+1, ndos-1
                dos(meshpt,isp) = dos(meshpt,isp) + wt
5            enddo
6         enddo
61     enddo
7   enddo
  end subroutine maknos
  subroutine splwts(nqp,nband,nbmx,nsp,wgts,evl,n,w,efermi, & !make sampling weights for integrals under the Fermi surface
       metal,sumev,bndwts,wtot,entrpy,dosef,cv)
    !i  Input
    !i    nqp : number of q-points; nband : number of bands
    !i    wgts: band weights
    !i    evl : energy eigenvalues
    !i    n   : n>0 Methfessel-Paxton polynomial order
    !i        : n<0 sampling done with Fermi-Dirac statistics
    !i    w   : n>0 gaussian width in Methfessel-Paxton integration (Ry)
    !i        : n<0 Temperature for Fermi distribution (Ry)
    !i    nbmx : first dimension of evl ;
    !i    metal : if F, weights unity below E_f and zero above.
    !i    efermi : Fermi energy
    !o  Output
    !o    bndwts : band and E_F - dependent k-point weights for integration
    !o    wtot   : sum of all weights (charge) qval
    !o    entrpy : electron entropy
    !o    dosef  : DOS at Fermi energy
    !o    cv     : electronic specific heat
    !o             (only evaluated with Fermi-Dirac statistics)
    !r  Remarks
    !r    sum of occupied eigenvalues = sum_n,k  w_nk E_nk
    !r    w_nk are generalised occupation numbers;
    !r    see Needs et al. Phys Rev B 33 (1986) 3778, eqs 1 & 2.
    !u Updates
    !u   16 Jul 08 returns entropy as TS for all n
    !u   04 Aug 07 Generates dos(efermi), cv(T=w) for F-D statistics
    !u   02 May 07 (MvS) prints entropy to stdout
    !u   21 Jun 06 (ATP) generates entrpy as output
    !u   17 Jan 05 Output wtot
    !-----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: nqp,nband,nbmx,nsp,n,ix,isplwts,i_copy_size
    logical :: metal
    real(8) :: wgts(nqp),evl(nbmx,nsp,nqp),w,efermi,sumev, &
         bndwts(nband,nsp,nqp),wtot,entrpy,dosef,cv
    ! ... Local parameters
    integer :: iqp,iband,isp,iprint,i1mach
    real(8) :: e,s,d,wt,x,xx,dsdt,tdsdt
    logical :: fractional
    sumev = 0d0
    wtot = 0d0
    entrpy = 0d0
    dosef = 0
    tdsdt = 0
    fractional=.false.
    do  3  iqp = 1, nqp
       do  2  iband = 1, nband
          do  1  isp = 1, nsp
             e = evl(iband,isp,iqp)
             if (metal) then

                !             debugging: check derivative numerically
                !              x = (efermi - e) / (w+1d-7)
                !              call delstp(n,x,d,s,sp)
                !              x = (efermi - e) / (w-1d-7)
                !              call delstp(n,x,d,s,sm)
                !              dsdt1 = (sp-sm)/2d-7

                x = (efermi - e) / w
                !             call delstp(n,x,d,s,xx)
                call delstd(n,x,d,s,xx,dsdt)
                if (abs(x) < 36) then
                   dsdt = -(efermi-e)/w**2 * dsdt
                else
                   dsdt = 0
                endif

                !C             debugging: compare analytical, numerical derivative
                !              if (abs(x) .lt. 30) then
                !                print 222, x,dsdt1,dsdt,dsdt-dsdt1
                !  222           format(3f14.8,1pe12.3)
                !              endif
             else
                s = 1d0
                if (e <= efermi) s = 0d0
                xx = 0
                d = 0
             endif
             wt = abs(wgts(iqp)) * (1d0 - s) / nsp
             bndwts(iband,isp,iqp) = wt
             if(0.1d0<wt .AND. wt<0.9d0) then
                fractional=.true.
             endif
             dosef = dosef + d*abs(wgts(iqp))/w/nsp
             wtot = wtot + wt
             sumev = sumev + e * wt
             entrpy = entrpy + xx  * abs(wgts(iqp)) / nsp
             tdsdt  = tdsdt + dsdt * abs(wgts(iqp)) / nsp
1         enddo
2      enddo
3   enddo
    tdsdt = tdsdt*w
    entrpy = entrpy*w
    if (n < 0) cv = tdsdt
    ! ... Print out band weights, if only 1 kp
    if (iprint() > 30 .AND. nqp == 1) then
       isplwts=1093
       open(isplwts,file="BandWeight.dat")
       do  isp = 1, nsp
          write(stdo,ftox)'SPLWTS: band weights .. Spin=',isp,'       eval      weight'
          if(fractional) then
             write(isplwts,*)"! Fractional occupation (criterion 0.1<wgt<0.9)"
             write(stdo,*)   "! Fractional occupation (criterion 0.1<wgt<0.9)"
          endif
          ix=0
          do  iband = 1, nband
             if(bndwts(iband,isp,1)<1d-7) ix=ix+1
             if(ix==10) then
                write (stdo,"('     ... ')")
                exit
             endif
             write (stdo,20) iband,evl(iband,isp,1),bndwts(iband,isp,1)
             write (isplwts,20) iband,evl(iband,isp,1),bndwts(iband,isp,1)
20           format (4x,i5,2f10.6)
          enddo
       enddo
       close(isplwts)
    endif
    ! ... Print out various k-point integrations
    if (iprint() >= 10) then
       if (n >= 0) then
          !          call awrit6(' N=%i, W=%d, E_F=%d, sumev=%d, entropy term:'
          !     .    //' %d, %d electrons',' ',256,i1mach(2), n,w,efermi,sumev,entrpy,wtot)
          write(stdo,"(a,i5,5d13.5)")'N W E_F sumev TS Nele=', &
               n,w,efermi,sumev,entrpy,wtot
       else
          !          call awrit5(' T=%dK, E_F=%d, sumev=%d, TS=%;3g,'//' %d electrons',' ',256,i1mach(2),
          write(stdo,"(a,5d13.5)")'T(K) E_F sumev TS Nele=', 0.1579d6*w,efermi,sumev,entrpy,wtot
          write(stdo,"(a,2d13.5)")'Entropy S, specific heat TdS/dT=',entrpy/w,tdsdt
          write(stdo,"(a)")'Fermi-Dirac;sampling'
       endif
       !      call info5(10,0,0,' SPLWTS: Fermi energy:%;6d;'//
       !     .  '  band energy=%;6d;  %;6d electrons  DOS(E_f)=%;4g',
       !     .    efermi,sumev,wtot,dosef,0)
    endif
  end subroutine splwts
  subroutine delstd(n,x,d,s,e,ep) !- Returns generalised delta and step functions (Methfessel & Paxton)
    !i  Inputs
    !i    n  : order of approximant; see Remarks
    !i       : n>=0 returns Methfessel-Paxton broadening
    !i       : n<0  returns Fermi-Dirac broadening
    !i    x  : (efermi - e) / width
    !i       : width should be gaussian width (n>=0)
    !i       : or temperature (n<0)
    !o  Outputs
    !o    D_n (x): smeared delta-function
    !o    S_n (x): smeared heaviside function
    !o    e_n (x): entropy
    !o    ep     : de/dx (Fermi-Dirac case only)
    !r  Remarks
    !r    For Methfessel-Paxton (generalized gaussian) broadening
    !r    (see Phys Rev B40, 3616 (1989))
    !r      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
    !r      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
    !r      where H is a Hermite polynomial and
    !r      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
    !r    For Fermi-Dirac broadening
    !r      s = 1/(exp(x)+1)   (fermi function)
    !r      d = ds/dx = exp(x)*s*s
    !r      e = -( s*log(s) + (1-s)*log(1-s) )
    !r     ep = log(s) - log(1-s)
    !u Updates
    !u   04 Aug 07 extended delstp to returns ep in Fermi-dirac case
    !u   23 May 00 extended to handle Fermi-Dirac broadening
    !-----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: n
    real(8) :: x,d,s,e,ep
    ! ... Local parameters
    integer :: i,k
    real(8) :: a,h1,h2,h3,s0,ex2,derfc,srpi
    !      intrinsic dsqrt,datan,dexp
    srpi = dsqrt(4d0*datan(1d0))
    ! ... Fermi-Dirac broadening
    if (n < 0) then
       if (x < -36d0) goto 91
       if (x >  36d0) goto 92
       s = 1d0/(dexp(x)+1d0)
       d = dexp(x)*s*s
       e = -( s*dlog(s) + (1-s)*dlog(1-s) )
       ep = (dlog(s) - dlog(1-s)) * s**2 * exp(x)
       return
    endif

    ! ... Methfessel-Paxton broadening
    if (x < -6d0) goto 91
    if (x >  6d0) goto 92
    ex2 = dexp(-x*x)
    s0 = .5d0 * erfc(x)
    a = 1d0/srpi
    k = 0
    h1 = 1d0
    h2 = 2d0 * x
    s = 0d0
    d = a
    do  1  i = 1, n
       a = -a / ( 4d0*i )
       k = k+1
       h3 = h1
       h1 = h2
       h2 = 2*x*h2 - 2*k*h3
       s = s + a*h1
       k = k+1
       h3 = h1
       h1 = h2
       h2 = 2*x*h2 - 2*k*h3
       d = d + a*h1
1   enddo
    d = d * ex2
    s = s0 + s*ex2
    e = 0.5d0*a*h1*ex2
    ep = 0
    return
    ! ... Branch for very small or very large x
91  s = 1d0
    e = 0d0
    d = 0d0
    ep = 0d0
    return
92  s = 0d0
    e = 0d0
    d = 0d0
    ep = 0d0
    return
  end subroutine delstd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module m_bzintegration2

! subroutine ebcpl(mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap,wk,eb) !- Gather spin-polarized bands into a single group, or redistribute
!   !i Inputs
!   !i   mode  :0, gather; 1, scatter
!   !i   nbmx  :leading dimension of b and dimensions wk
!   !i   nevx  :number of eigenvalues
!   !i   nsp   :2 for spin-polarized case, otherwise 1
!   !i   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
!   !i   nq    :number of k-points for which eb are calculated
!   !i   nbpw  :a number no larger than the number of bits per integer word
!   !i   bmap  :an integer work array with dimension at nbmx*nsp*nq/nbpw
!   !i         :see Remarks
!   !i   wk    :a work array with dimension nbmx*nsp
!   !io Inputs/Outputs
!   !io   eb    :energy bands:
!   !io         :mode=0 input spin-split, output merged to a single vector
!   !io         :mode=1 input merged to a single vector, output spin-split
!   !r Remarks
!   !io   Call ebcpl with mode=1 to undo call of ebcpl with mode=0.
!   !io   bmap must be preserved for mode=1 call.
!   !u Updates
!   ! ----------------------------------------------------------------------
!   implicit none
!   integer mode,nbmx,nevx,nsp,nspc,nq,nbpw,bmap(1)
!   real(8) eb(nbmx,nsp,nq),wk(nbmx*nsp)
!   integer ib,iq,ib1,ib2,iqb
!   if (nsp .eq. 1 .or. nspc .eq. 2) return
!   ! --- Gather bands at each qp ---
!   ifx:if (mode .eq. 0) then
!      iqb = 0
!      do10:do    iq = 1, nq          !   ... Gather and order +,- bands at this qp into one column
!         ib1 = 1
!         ib2 = 1
!         do20: do ib = 1, nevx*nsp
!            iqb = iqb+1
!            if (eb(min(ib1,nevx),1,iq) .lt. eb(min(ib2,nevx),2,iq)&
!                 .and. ib1 .le. nevx .or. ib2 .gt. nevx) then
!               wk(ib) = eb(ib1,1,iq)
!               ib1 = ib1+1
!            else
!               wk(ib) = eb(ib2,2,iq)
!               ib2 = ib2+1
!               call mark1(bmap, nbpw, iqb)
!            endif
!         enddo do20
!         call dcopy(nevx*nsp,wk,1,eb(1,1,iq),1)
!         if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
!      enddo do10
!   endif ifx
!   ! --- Disperse bands at each qp ---
!   if (mode .eq. 1) then
!      iqb = 0
!      d110  : do  iq = 1, nq   !   ... Disperse bands into +,- for this qp according to bmap
!         ib1 = 1
!         ib2 = 1
!         call dcopy(nevx*nsp,eb(1,1,iq),1,wk,1)
!         do120  :     do  ib = 1, nevx*nsp
!            iqb = iqb+1
!            if (iget(bmap,nbpw,iqb) .eq. 0) then
!               eb(ib1,1,iq) = wk(ib)
!               ib1 = ib1+1
!            else
!               eb(ib2,2,iq) = wk(ib)
!               ib2 = ib2+1
!            endif
!         enddo do120
!         if (ib1-1 .ne. nevx .and. ib2-1 .ne. nevx) call rx('bug')
!      enddo d110
!   endif
! end subroutine ebcpl
! subroutine mark1(bitmap, nbpw, n)    !- put a one in the nth bit of bitmap.
!   !i Inputs   !i   bitmap, n
!   !r    nbpw: a number no larger than the number of bits per integer word
!   !     implicit none
!   integer bitmap(*), nbpw, n
!   ! Local parameters
!   integer nword,nbit,i
!   nword = (n-1)/nbpw
!   nbit = mod(n-1,nbpw)
!   i = 2**(nbpw-nbit-1)
!   bitmap(nword+1) = bitmap(nword+1) + i*(1-mod(bitmap(nword+1)/i,2))
! end subroutine mark1
! integer function iget(bitmap, nbpw, n)    !- Return 0 or 1, depending on the value of the nth bit of bitmap
!   implicit none
!   integer bitmap(*), nbpw, n
!   integer nword,nbit
!   nword = (n-1)/nbpw
!   nbit = mod(n-1,nbpw)
!   iget = mod(bitmap(nword+1)/2**(nbpw-nbit-1),2)
! end function iget

! subroutine fermi(qval,dosi,ndos,emin,emax, eferm,e1,e2,dosef) !Makes fermi energy from integrated density
!   use m_ftox
!   use m_nvfortran
!   !i   qval:    number of electrons to fermi level
!   !i   dosi(i): integrated density at bin i;
!   !i   ndos: number of bins + 1
!   !i   emin, emax: energy window.
!   !o Outputs
!   !o   Eferm, Fermi energy;
!   !o   e1<=Eferm<=e2 : confidence limits on Fermi energy. i.e., Fermi energy lies between e1 and e2.
!   !o   dosef:  density of states at fermi level
!   !remark: emin and e1 (and emax and e2) may point to the same address.
!   implicit none
!   integer :: ndos,i1,ie
!   real(8) :: qval,dosi(ndos),emin,emax,eferm,e1,e2,dosef,de,q1,q2
!   if(dosi(1)>qval.OR.dosi(ndos)<qval) call rx('fermi does not encompass qval='//ftof(qval))
!   i1 = findloc([(dosi(ie)>qval, ie=1,ndos)],value=.true.,dim=1) - 1
!   de = (emax-emin)/(ndos-1)
!   e1 = emin + de*(i1-1)
!   e2 = emin + de*i1
!   q1 = dosi(i1)   
!   q2 = dosi(i1+1) 
!   eferm = e1 + (qval-q1)/(q2-q1)*de ! Linear interpolation for the Fermi level
!   dosef = (q2-q1)/de
! end subroutine fermi

