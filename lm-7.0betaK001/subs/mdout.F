      subroutine mdout(nbas,amass,tstep,f,qmpol,iclas,temp,relax,
     x                 bas,vel,zeta,zacc,zsqacc,itrlx,ekin,tkin,
     x                 etot)
C --- Print out some infprmation about the MD.
C ---------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,iclas(nbas),itrlx
      double precision amass(nbas),f(3,nbas),bas(3,nbas),
     .  vel(3,nbas),temp,ekin,tkin,etot,tstep,zeta,qmpol(9,nbas),
     .  relax,zacc,zsqacc
C Local Variables
      integer id,np,fopn,i1mach
      double precision tzeta,tzetsq,totmom(3),totfrc(3),cons,
     x     autime,autmp
      data autime/0.048377d0/
      data autmp/0.6333328d-5/

c --- Check the health of the simulation 
      zacc= zacc + zeta
      zsqacc= zsqacc + zeta*zeta
c --- 1) zeta should avarage to zero:
      tzeta= zacc / itrlx
c --- 2) mean fake kinetic energy:
      tzetsq= zsqacc / itrlx
c --- 3) conserved quantity:
      cons= etot + ekin + 3*(nbas-1)*temp*zacc*tstep  
     x     + 0.5d0*3*(nbas-1)*temp*relax**2*zeta**2
c --- 4) center of mass (at the half step!) and total force:
      do 150 id=1,3
        totmom(id)=0.d0
        totfrc(id)=0.d0
        do 140 np=1,nbas
          totmom(id)= totmom(id) + amass(iclas(np))*vel(id,np)
          totfrc(id)= totfrc(id) + f(id,np)
 140    continue
 150  continue

      write(fopn('MD'),100)itrlx*tstep*autime,
     x   tkin/autmp,(totmom(id),id=1,3),zeta,tzeta,tzetsq
      write(fopn('ENG'),200)itrlx*tstep*autime,cons,etot+ekin,
     x     etot,ekin,(totfrc(id),id=1,3)
      write(fopn('MPOL1'),300)itrlx*tstep*autime,
     x     (qmpol(id,1) , id=1,9)
      write(fopn('MPOL2'),300)itrlx*tstep*autime,
     x     (qmpol(id,nbas) , id=1,9)
        
      call awrit4('%N ITER= %;5i    ETB= %;8d     EKIN= %;8d'//
     x     '     SUM= %;9d Ry %N',' ',120,i1mach(2),
     x     itrlx,etot,ekin,etot+ekin)

      call awrit3(' MDOUT: conserved quantity = %;8d Ry %N'//
     x            '        temperature        = %;3d K '//
     x            ' target temperature= %;3d K ',' '
     x       ,240,i1mach(2),cons,tkin/autmp
     x       , temp/autmp)
      call awrit6('        total momentum     = %;8d %;8d %;8d %N'//
     x            '        total force        = %;8d %;8d %;8d ',' '
     x       ,240,i1mach(2),totmom(1),totmom(2),totmom(3)
     x       ,totfrc(1),totfrc(2),totfrc(3))

      if (relax.lt.1d10)
     x       call awrit3('        zeta = %;8d '//
     x                   '       <zeta>= %;8d '//
     x                   '       <zeta>^2= %;8d ',' '
     x       ,240,i1mach(2),zeta,tzeta,tzetsq)
 
 100  format(f9.2,f10.2,6e12.3)
 200  format(f9.2,4f14.8,3f12.8)
 300  format(f9.2,9f11.6)
      end
        
        

