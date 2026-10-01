!> Quantities derived from the space group operations: site tables (symtbl, mptauof), classes (splcls), the closure of a
!> group (grpgen), rotation matrices of real harmonics (rotdlmm). Split from m_mksym_util (2026-10-02). Pure routines, no state.
module m_symderive
  use m_lgunit,only:stdo
  use m_ftox
  use m_nvfortran,only: findloc
  use m_mpi,only: master_mpi
  use m_symop_util,only: asymop,grpeql,latvec
  public symtbl,splcls,grpgen,mptauof,rotdlmm,iclbsjx
  private
  real(8),parameter:: toll=1d-4,tiny=1d-4,epsr=1d-12
contains
  subroutine symtbl(mode,nbas,pos,g,ag,ng,qlat,istab) ! Make symmetry transformation table for posis atoms; check classes
    !i Inputs
    !i   mode  :1st digit
    !i         :0  site ib is transformed into istab(ib,ig) by grp op ig
    !i         :1  site istab(i,ig) is transformed into site i by grp op ig
    !i   nbas  :size of basis
    !i   pos   :pos(i,j) are Cartesian coordinates of jth atom in basis
    !i   g     :point group operations
    !i   ag    :translation part of space group
    !i   ng    :number of group operations
    !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
    !o Outputs  istab :table of site permutations for each group op; see mode
    implicit none
    integer :: nbas,ng,mode, istab(nbas,*),ib,ic,ig,jb,jc,ka
    real(8) :: pos(3,nbas),g(3,3,ng),ag(3,ng),qlat(9)
    character(8):: xt,xn
    do    ig = 1, ng !Make atom transformation table ---
       do ib = 1, nbas          
          jb=findloc( [(latvec(1,toll,qlat, matmul(g(:,:,ig),pos(:,ib))+ag(:,ig)-pos(:,ka)), ka=1,nbas)],dim=1,value=.true.)!ib is mapped to jb by g,ag 
          if(jb == 0) call rx("SYMTBL: no map for atom ib="//trim(xn(ib))//" ig="//trim(xn(ig)))
          if (mode == 0) then
             istab(ib,ig) = jb
          else
             istab(jb,ig) = ib
          endif
       enddo
    enddo
  end subroutine symtbl
  subroutine splcls(bas,nbas,ng,istab,nspec,slabl,nclass,ipc, ics,nrclas) !- Splits species into classes
    !i   bas,nbas: dimensionless basis vectors, and number
    !i   nspec:    number of species
    !io   ipc:      on input, site j belongs to species ipc(j)
    !i   slabl:    on input, slabl is species label
    !i   ng:       number of group operations
    !i   istab:    site ib is transformed into istab(ib,ig) by grp op ig
    !o  Outputs:
    !o   ipc:      site j belongs to class ipc(j)
    !o   ics:      class i belongs to species ics(i)
    !o   nclass:   number of classes
    !o   nrclas:   number of classes per each species
    implicit none
    logical :: nosplt
    integer :: nbas,nspec,nclass,ng,istab(nbas,ng),ipc(nbas), ics(nbas),nrclas(nspec)
    real(8) :: bas(3,*)
    character(8) :: slabl(*)
    integer :: ib,ic,icn,ig,jb,m,i,is,ipr,idx,ispec,j
    logical :: lyetno
    character(80) :: outs,clabl=''
    call getpr(ipr)
    nclass = nspec
    ics = [(i,i=1,nspec)]
    ic = 1
    do while(ic <= nclass) 
       is = ics(ic)
       ib = iclbsjx(ipc,nbas, ic,1)
       if (ib == 0) goto 11 !   ... No sites of this class ... skip
       lyetno = .true.
       do 20  jb = 1, nbas !For each basis atom in this class, do
          if (ipc(jb) == ic) then !class of jb
             if(  any(istab(ib,1:ng) == jb).or.&                    !If there is a g mapping ib->jb, sites are equivalent
                  any(istab(jb,1:ng) == ib).or.&                    !If there is a g mapping jb->ib, sites are equivalent
                  any([(istab(istab(ib,ig),ig)== jb,ig=1,ng)]).or.&   !If there is a g mapping ib->kb,jb, sites are equivalent
                  any([(istab(istab(jb,ig),ig)== ib,ig=1,ng)])) cycle !If there is a g mapping jb->kb,ib, sites are equivalent
             if (lyetno) then !If the classes haven't been split yet, do so
                nclass = nclass+1
                icn  =  nclass
                ics(icn) = is
                nrclas(is) = nrclas(is)+1
                lyetno = .false.
             endif
             if(nclass > nbas) call rx('splcls:  problem with istab')
             icn  =  nclass
             ipc(jb)=  icn !class index
          endif
20     enddo
11     continue
       ic = ic + 1
    enddo
    if(ipr>=30) then
       write(stdo,"(a)")' splcls:  ibas iclass ispec label(ispec)'
       do j=1,nbas
          ic   = ipc(j) !class
          ispec= ics(ic)!spec
          write(stdo,"(a,3i6,a)")"       ",j,ic,ispec,'     '//trim(slabl(ispec))
       enddo
    endif
  end subroutine splcls
  subroutine grpgen(gen,ngen,symops,ng,ngmx) !Generate all point symmetry operations from the generation group
    !i   gen,ngen,ngmx
    !i   if ng>0, add symops to the ng already in list.
    !o Outputs: symops,ng
    !r Remarks   This works for point groups only and is set up for integer  generators.
    implicit none
    integer :: ngen,ng,ngmx
    real(8) :: gen(3,3,ngen),symops(3,3,ngmx), h(3,3),hh(3,3),sig(3,3)
    integer :: igen,ig,itry,iord,nnow,j,ip,i,k,n2,m1,m2,n,m, ipr
    character(80) :: sout
    character(8):: xn
    real(8),parameter:: ee(9)=[1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],e(3,3)= reshape(ee,[3,3]),ae(3)=0d0
    call getpr(ipr)
    sout = ' '
    symops(:,:,1)= e
    if (ng < 1) ng = 1
    igenloop: do  80  igen = 1, ngen
       sig = gen(:,:,igen) !  Extend the group by all products with sig ---
       do ig = 1, ng 
          if (grpeql(symops(:,:,ig),sig) .AND. ipr > 30)  write(stdo,ftox)' Generator ',igen,' already in group as element',ig
          if (grpeql(symops(:,:,ig),sig)) goto 80
       enddo
       h=sig
       do  itry = 1, 100
          iord = itry
          if (grpeql(h,e)) exit
          h=matmul(sig,h) 
       enddo
       nnow = ng
       if(ipr >= 40) write(stdo,ftox) trim(sout),' ',igen,' is',iord
       do j = 1, ng !Products of type  g1 sig**p g2 ---
          h = symops(:,:,j) 
          do   ip = 1, iord-1  
             h = matmul(sig,h) ! h = sig**ip
             do i = 1, ng    
                hh = matmul(symops(:,:,i),h) ! hh = symops_i sig**ip
                if(any([(grpeql(symops(:,:,k),hh),k=1,nnow)])) cycle
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                symops(:,:,nnow)=hh 
             enddo
          enddo
          if (j == 1) n2 = nnow
       enddo
       m1 = ng+1
       m2 = nnow
       do       i = 2, 50 ! --- Products with more than one sandwiched sigma-factor ---
          do    n = ng+1, n2
             do m = m1, m2
                h= matmul(symops(:,:,n),symops(:,:,m)) 
                if(any([(grpeql(symops(:,:,k),h),k=1,nnow)])) cycle
                nnow = nnow+1
                if (nnow > ngmx) goto 99
                symops(:,:,nnow)=h 
             enddo
          enddo
          if (m2 == nnow) exit
          m1 = m2 + 1
          m2 = nnow
       enddo
       ng = nnow
80  enddo igenloop
    if( ipr >= 30) then
       if(sout /= ' ' .AND. ipr >= 60) write(stdo,ftox)' Order of generator '//trim(sout)
       write(stdo,ftox)'GRPGEN:',ng,'symmetry operations from',ngen,'generator(s)'
       if(ipr >= 80 .AND. ng > 1) then
          write(stdo,'('' ig  group op'')')
          do  ig = 1, ng
             call asymop(symops(1,1,ig),ae,' ',sout)
             write(stdo,'(i4,2x,a)') ig,trim(sout)
          enddo
       endif
    endif
    return
99  continue
    call rx('GRPGEN: too many elements nnow ngmx='//trim(xn(nnow))//' '//trim(xn(ngmx)))
  end subroutine grpgen
  subroutine mptauof(symops,ng,plat,nbas,bas, iclass,miat,tiat,invg,delta,afmode) !- Mapping of atomic sites by points group operations.
    use  m_lmfinit,only: iantiferro
    !i  Input
    !i     symops(1,ng),ng,plat,nbas,bas(3,nbas)
    !i     iclass(nbas); denote class for each atom
    !o  Output
    !o    miat(ibas  ,ig); ibas-th atom is mapped to miat-th atom, by the ig-th
    !o    points group operation.  Origin is (0,0,0).
    !o    tiat(k,ibas,ig);
    !o    delta : shifting vector for non-symmorphic group.
    !o            r' = matmul (am, r) + delta
    !r  Remarks
    !r
    !r (1) The ibas-th atom (position at bas(k,ibas) ) is mapped to
    !r
    !r    bas( k,miat(ibas,ig) )+ tiat(k,ibas,ig), k=1~3.
    !r
    !r (2) tiat= unit translation
    implicit none
    integer :: ng,nbas, miat(nbas,ng),iclass(nbas),invg(ng), &
         nbmx, nsymx, ig,igd,i,j,ibas,mi,i1,i2,i3
    double precision :: SYMOPS(9,ng),plat(3,3), &
         tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas), &
         tr1,tr2,tr3,ep, dd1,dd2,dd3,t1,t2,t3
    integer::  iprintx=0
    integer :: ires(3, nbas, ng)
    integer:: ib1,ib2
    real(8) ::tran(3),delta(3,ng)
    logical,optional:: afmode
    ep=1d-3
    if(iprintx>=46) write(6,*)'MPTAUOf: search miat tiat for wave function rotation'
    do 10 ig=1,ng
       do igd=1,ng
          ! seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ) <= ep .AND. &
               abs( symops(2,ig)-symops(4,igd) ) <= ep .AND. &
               abs( symops(3,ig)-symops(7,igd) ) <= ep .AND. &
               abs( symops(4,ig)-symops(2,igd) ) <= ep .AND. &
               abs( symops(5,ig)-symops(5,igd) ) <= ep .AND. &
               abs( symops(6,ig)-symops(8,igd) ) <= ep .AND. &
               abs( symops(7,ig)-symops(3,igd) ) <= ep .AND. &
               abs( symops(8,ig)-symops(6,igd) ) <= ep .AND. &
               abs( symops(9,ig)-symops(9,igd) ) <= ep  ) then
             invg(ig)=igd
             goto 16
          endif
       enddo
16     continue
       do i=1,3
          do j=1,3
             am(i,j)=symops(i+3*(j-1),ig)
          enddo
       enddo
       do 120 ib1=1,nbas ! trial shift vector tran
          do 121 ib2=1,nbas
             tran =  bas(:,ib2)  - matmul(am,bas(:,ib1))
             if(present(afmode)) then
                if(iantiferro(ib1)==0) cycle
                if(iantiferro(ib2)==0) cycle
                if(iantiferro(ib1)+iantiferro(ib2)/=0) cycle
             endif
             do 30 ibas=1,nbas
                !bb1=matmul(am,bas(:,ibas))+trans
                b1=am(1,1)*bas(1,ibas)+am(1,2)*bas(2,ibas)+am(1,3)*bas(3,ibas) +tran(1)
                b2=am(2,1)*bas(1,ibas)+am(2,2)*bas(2,ibas)+am(2,3)*bas(3,ibas) +tran(2)
                b3=am(3,1)*bas(1,ibas)+am(3,2)*bas(2,ibas)+am(3,3)*bas(3,ibas) +tran(3)
                do 40 mi=1,nbas
                   if( iclass(mi) /= iclass(ibas) ) cycle
                   do  i1=-3,3
                      do  i2=-3,3
                         do  i3=-3,3
                            dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                            dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                            dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )
                            t1 = b1 - (bas(1,mi)+dd1)
                            t2 = b2 - (bas(2,mi)+dd2)
                            t3 = b3 - (bas(3,mi)+dd3)
                            if(abs(t1) <= ep .AND. abs(t2) <= ep .AND. abs(t3) <= ep) go to 60
                         enddo
                      enddo
                   enddo
40              enddo
                goto 121 ! seach failed, Not found mi and dd1. Try next (tr).
60              continue
                miat(ibas,ig)  = mi
                tiat(1,ibas,ig)= dd1
                tiat(2,ibas,ig)= dd2
                tiat(3,ibas,ig)= dd3
                ires(1,ibas,ig)= i1
                ires(2,ibas,ig)= i2
                ires(3,ibas,ig)= i3
30           enddo
             goto 21 ! When the do-30 loop has been completed, we get out of do-20 loop
121       enddo
120    enddo
       call rx('mptauof: Can not find miat and tiat')
21     continue
       delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000
       !- have gotten the translation-> check write --------------------
       if(iprintx >= 46) then
          write(6,4658)tran
4658      format('  Obtained translation operation=',3d12.4)
          do 123  ibas=1,nbas
             write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig), &
                  tiat(2,ibas,ig), tiat(3,ibas,ig), &
                  ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
150          format(' iiiiibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
123       enddo
       endif
10  enddo
  end subroutine mptauof

  subroutine rotdlmm(symops,ng,nl ,dlmm) ! Generate rotation matrix D^l_{m,m'} for L-representaiton,
    !  corresponding to points group operations.
    !i symops(9,ng),ng; point ops.
    !i nl; num.of l =lmax+1
    !o dlmm(2*nl-1,2*nl-1,0:nl-1,ng,2); D^l_{m,m'}. Indexes are for Real harmonics.
    !r  dlmmc is used as work area about 200kbyte used for  s,p,d,f -> nl=4
    !-----------------------------------------------------------------
    implicit double precision (a-h,o-z)
    integer:: is,i,ig,ikap,j,l,m,m1,m2,m3,md,mx,ix, ng,nl
    double precision :: SYMOPS(9,ng), am(3,3) ,fac1,fac2
    double precision :: dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
    double precision :: det,osq2
    complex(8):: msc(0:1,2,2), mcs(0:1,2,2),dum(2),&
         dlmmc(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
    complex(8),parameter:: Img=(0d0,1d0)
    integer:: debugmode
    real(8):: ep=1d-3 !ep was 1d-8 before feb2013
    do 10 ig =1,ng
       do  i=1,3
          do  j=1,3
             am(i,j) = symops(i+3*(j-1),ig)
          enddo
       enddo
       ! calculate determinant(signature)
       det= am(1,1)*am(2,2)*am(3,3) &
            -am(1,1)*am(3,2)*am(2,3) &
            -am(2,1)*am(1,2)*am(3,3) &
            +am(2,1)*am(3,2)*am(1,3) &
            +am(3,1)*am(1,2)*am(2,3) &
            -am(3,1)*am(2,2)*am(1,3)
       if(abs(abs(det)-1d0) >= 1d-10) then
          print *,' rotdlmm: det/=1 ig and det=',ig,det
          stop
       endif
       ! seek Euler angle   print *,' goto cbeta',ig,det
       cbeta = am(3,3)/det
       ! added region correction so as to go beyond domain error for functions, dsqrt and acos.
       if(abs(cbeta-1d0) <= 1d-6) cbeta= 1d0
       if(abs(cbeta+1d0) <= 1d-6) cbeta=-1d0
       beta = dacos(cbeta) ! beta= 0~pi
       sbeta= sin(beta)
       if(sbeta <= 1.0d-6) then
          calpha= 1d0
          salpha= 0d0
          alpha = 0d0
          cgamma= am(2,2)/det
          sgamma= am(2,1)/det
       else
          salpha =  am(2,3)/sbeta/det
          calpha =  am(1,3)/sbeta/det
          sgamma =  am(3,2)/sbeta/det
          cgamma = -am(3,1)/sbeta/det
       endif
       co2 = dcos(beta/2d0)
       so2 = dsin(beta/2d0)
       if(abs(calpha-1.0d0) <= 1.0d-6) calpha= 1.0d0
       if(abs(calpha+1.0d0) <= 1.0d-6) calpha=-1.0d0
       if(abs(cgamma-1.0d0) <= 1.0d-6) cgamma= 1.0d0
       if(abs(cgamma+1.0d0) <= 1.0d-6) cgamma=-1.0d0
       alpha=dacos(calpha)
       if(salpha < 0d0) alpha=-alpha
       gamma=dacos(cgamma)
       if(sgamma < 0d0) gamma=-gamma  !print *,'alpha beta gamma det=',alpha,beta,gamma,det
       do l =  0, nl-1
          do md= -l, l
             do m = -l, l
                !  from 'Ele theo. ang. mom. by M. E. Rose 5th 1967 Wisley and Sons.  p.52 (4.13)
                fac1 = dsqrt( igann(l+m)*igann(l-m)*igann(l+md)*igann(l-md) )
                fac2 = 0d0
                do ikap=0,2*l
                   if(l-md-ikap >= 0 .AND. l+m-ikap >= 0 &
                        .AND. ikap+md-m >= 0) then
                      add= dble((-1)**ikap)/( igann(l-md-ikap)*igann(l+m-ikap) &
                           *igann(ikap+md-m)*igann(ikap) )
                      if(2*l+m-md-2*ikap /= 0) add=add*co2**(2*l+m-md-2*ikap)
                      if(md-m+2*ikap /= 0)     add=add*(-so2)**(md-m+2*ikap)
                      fac2 = fac2+add
                   endif
                enddo
                ! l-th rep. is odd or even according to (det)**l
                dlmmc(md,m,l,ig) = fac1*fac2*det**l* cdexp( -Img*(alpha*md+gamma*m) )
             enddo
          enddo
       enddo
       am(1,1)= cos(beta)*cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)
       am(1,2)=-cos(beta)*cos(alpha)*sin(gamma)-sin(alpha)*cos(gamma)
       am(1,3)= sin(beta)*cos(alpha)
       am(2,1)= cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
       am(2,2)=-cos(beta)*sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)
       am(2,3)= sin(beta)*sin(alpha)
       am(3,1)=-sin(beta)*cos(gamma)
       am(3,2)= sin(beta)*sin(gamma)
       am(3,3)= cos(beta)
       if(abs(am(1,1)*det-symops(1,ig))>ep .OR. &
            abs(am(2,1)*det-symops(2,ig))>ep .OR. &
            abs(am(3,1)*det-symops(3,ig))>ep .OR. &
            abs(am(1,2)*det-symops(4,ig))>ep .OR. &
            abs(am(2,2)*det-symops(5,ig))>ep .OR. &
            abs(am(3,2)*det-symops(6,ig))>ep .OR. &
            abs(am(1,3)*det-symops(7,ig))>ep .OR. &
            abs(am(2,3)*det-symops(8,ig))>ep .OR. &
            abs(am(3,3)*det-symops(9,ig))>ep) then
          print *,' rotdlmm: not agree. symgrp and one by eular angle'
          stop
       endif
       if(debugmode()>9) then
          print *;print *;print *,' **** group ops no. ig=', ig
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
          print *,' by Eular angle '
          write(6,1731)am(1,1)*det,am(1,2)*det,am(1,3)*det
          write(6,1731)am(2,1)*det,am(2,2)*det,am(2,3)*det
          write(6,1731)am(3,1)*det,am(3,2)*det,am(3,3)*det
       endif
1731   format (' ',3f9.4)
10  enddo
    ! conversion to real rep. Belows are from csconvs
    !  msc mcs conversion matrix generation 2->m 1->-m for m>0
    osq2 = 1d0/sqrt(2d0)
    do m = 0,1
       Msc(m,1,:)= osq2*[complex(8):: (-1d0)**m, -Img*(-1d0)**m] !spherical to real(cubic)
       Msc(m,2,:)= osq2*[complex(8)::       1d0,            Img]
       Mcs(m,1,:)= osq2*[complex(8):: (-1d0)**m,      1d0]     !inverse
       Mcs(m,2,:)= osq2*[complex(8):: Img*(-1d0)**m, -Img]
    enddo
    converttoreal:do 123 is=1,ng ! convert to real rep.
       llooop:do 23   l =0,nl-1
          do  m2=-l,l
             do  m1= 1,l
                mx    = mod(m1,2)
                dum= [dlmmc(m2, m1,l,is), dlmmc(m2,-m1,l,is)]
                dlmmc(m2,  m1,l,is)= sum(dum(:)*msc(mx,:,1))
                dlmmc(m2, -m1,l,is)= sum(dum(:)*msc(mx,:,2))
             enddo
          enddo
          do m2=  1,l
             do m1= -l,l
                mx=mod(m2,2)
                dum= [dlmmc( m2, m1,l,is),dlmmc(-m2, m1,l,is)]
                dlmmc( m2, m1,l,is)= sum(mcs(mx,1,:)*dum(:))
                dlmmc(-m2, m1,l,is)= sum(mcs(mx,2,:)*dum(:))
             enddo
          enddo
          do m2=-l,l
             do m1=-l,l
                dlmm(m2,m1,l,is)=dreal( dlmmc(m2,m1,l,is) )
                if( abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12 ) &
                     call rx(' rotdlmm: abs(dimag(dlmmc(m2,m1,l,is))) >= 1.0d-12')
             enddo
          enddo
          if( .FALSE. ) then
             print *; print *,'  points ops  ig, l=', is,l,' cubic   '
             do m2=-l,l
                write(6,"(28f10.5)")( dreal(dlmmc (m2, m1,l,is) ), m1=-l,l)
                !    &    , ( dimag(dlmmc (m2, m1,l,is) ), m1=-l,l),( dlmm(m2, m1,l,is), m1=-l,l)
             enddo
          endif
23     enddo llooop
123 enddo converttoreal
    if(debugmode()>1) print *,' end of rotdlmm'
  end subroutine rotdlmm
  !--------------------------------------------
  double precision function igann(i)
    integer:: i,ix
    igann  = 1d0
    do ix =1,i
       igann=igann*dble(ix)
    enddo
  end function igann
  integer function iclbsjx(ipc,nbas, ic,nrbas) !the nrbas-th atom belonging to class ic (ipc(ibas)==ic)
    use m_nvfortran,only: findloc,count
    implicit none
    integer :: ic,nbas,ipc(nbas),nrbas,ib,ibas
    integer:: rrr(nbas),ccc(nbas)
    do ibas=1,nbas
       ccc(ibas)=count([(ipc(ib)==ic,ib=1,ibas)])
    enddo
    iclbsjx = findloc( ccc(1:nbas), value=nrbas,dim=1)
  !!!!not working in nvfortran24.1
  !  ccc(:)=[(countl([(ipc(ib)==ic,ib=1,ibas)]),ibas=1,nbas)]
  !  iclbsjx = findloci( ccc(1:nbas), value=nrbas,dim=1)
  
  ! iclbsjx = findloci([(countl([(ipc(ib)==ic,ib=1,ibas)]), ibas=1,nbas)], value=nrbas,dim=1)

  !  associate(ccc => [(countl([(ipc(ib)==ic,ib=1,ibas)]),ibas=1,nbas)] )
  !    iclbsjx = findloci( ccc(1:nbas), value=nrbas,dim=1)
  !  endassociate
  end function iclbsjx
end module m_symderive
