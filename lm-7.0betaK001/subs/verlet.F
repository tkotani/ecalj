      subroutine verlet(ensbl,it,nbas,nf,amass,tstep,f,ipc,dclabl,temp,
     .                  etot,tpv,vol,pext,lq,qmpol,taup,taub,mdtime,
     .                  time,alat,pos,vel,zeta,eps,veps,logs,zacc,
     .                  zsqacc,ekin,tkin,cons)
C- Velocity Verlet
C ----------------------------------------------------------------------
Ci Inputs: 
Ci   ensbl: ensemble ('NVE', 'NVT' or 'NPT')
Ci   it: iteration number of current MD run
Ci   nbas, nf: number of atoms, number of degrees of freedom 
Ci                              set in initv
Ci   amass: atomic mass for each species (in a.u.)
Ci   tstep: time step, dt (in a.u.)
Ci   alat: lattice constant (scale for pos into a.u.)
Ci   f: force at time t; accelerations at time t are f/amass
Ci   ipc: pointers from atom to species 
Ci   dclabl: d.p. translations of atom labels, passed through to mdwrit
Ci   temp: target temperature (actually kT in Ry)
Ci   etot: potential energy at time t
Ci   tpv: 3pV, p is the "internal virial" at time t (see tbtote)
Ci             that is, p = (1/3) \sum r_{ij} . f_{ij}
Ci             (Allen and Tildesley, p. 46)
Ci   vol: volume of simulation cell (a.u.)
Ci   pext: external (applied) pressure (in Ry/bohr^3)
Ci   lq: logical switch, if T qmpol is passed
Ci   qmpol: multipole moments (for writing to xyz file)
Ci   taup: thermostat relaxation time
Ci   taub: barostat relaxation time
Ci   mdtime: total mdtime to be completed, including from previous runs
Ci   time: time t (in a.u.), including time accumulated from previous
Ci         runs
Ci   pos: positions at time t
Ci   vel: velocities at time t - dt/2
Ci   zeta: Hoover thermostat friction at time t - dt/2 (this is the
Ci         quantity p_\xsi/Q in Martyna et al.)
Ci   eps: (1/3) ln V at time t (used if NPT)
Ci   veps: barostat friction at time t - dt (used if NPT) (this is 
Ci         the quantity p_\epsilon/W_c in Martyna et al.)
Ci   logs: ln s, "thermostat position;" s is Nose's scaling variable, s,
Ci         at time t - dt/2
Ci   zacc: zeta accumulated over time, t - dt
Ci   zsqacc: zeta^2 accumulated over time, t - dt
Co Outputs:
Ci   time: time t + dt (in a.u.)
Ci   zacc: zeta accumulated over time, t
Ci   zsqacc: zeta^2 accumulated over time, t
Co   ekin: Kinetic energy at time t
Co   tkin: Temperature at time t (actually kT in Ry) 
Co   zeta: Hoover friction at time t + dt/2
Co   logs: ln s at time t + dt/2
Co   pos: positions at time t + dt
Co   vel: velocities at time t + dt/2
Co   cons: conserved quantity at time t
Co   eps: (1/3) ln V/V_0 at time t + dt
Cr Remarks
Cr   qmpol (tight-binding multipole moments) are passed in so they
Cr   can be written to the xyz file if switch lq is set. In the same
Cr   way this segment can be modified to pass through other quantities
Cr   that may be written out for latter processing of files to get
Cr   ergodic averages.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer it,nbas,nf,ipc(nbas)
      character*3 ensbl
      logical lq
      double precision amass(1),tstep,alat,f(3,nbas),temp,taup,taub,
     .                 pos(3,nbas),vel(3,nbas),zeta,logs,ekin,tkin,etot,
     .                 qmpol(9,nbas),dclabl(1),zacc,zsqacc,mdtime,time,
     .                 tpv,eps,veps,pext,vol,cons
C Local Variables
      integer ib,id,ic
      double precision pint,acc,t,t2,sinhtt,et,et2

C --- operate with iL_1/2 (advance velocities to time t) ---
      do  ib = 1, nbas
        ic = ipc(ib)
        do  id = 1, 3
          acc = f(id,ib) / amass(ic)
          vel(id,ib) = vel(id,ib) + 0.5d0 * acc * tstep
        enddo
      enddo

C --- now have velocity at time t. Calculate kinetic energy ---
      call tekin(nbas,ipc,amass,vel,ekin,tkin)

C --- operate with iL_C/2 (advance zeta, logs to t; scale v) ---
      if (ensbl .eq. 'NVT') then
        call vvlc2(nbas,nf,temp,tstep,taup,vel,ekin,zeta,logs)
        call tekin(nbas,ipc,amass,vel,ekin,tkin)
      endif

C --- operate with iL_CP/2 (advance veps, zeta, logs to t; scale v) ---
      pint = 0d0
      if (ensbl .eq. 'NPT') then
        call vvlcp2(nbas,nf,temp,pext,tstep,taup,taub,tpv,vol,vel,
     .              ekin,zeta,veps,logs,pint)
        call tekin(nbas,ipc,amass,vel,ekin,tkin)
      endif

C --- Set to zero the velocity of the centre of mass ---
      call zercmv(nbas,vel,amass,ipc)

C --- now all quantities known at time t ---
      call mdinfo(ensbl,nbas,nf,amass,tstep,f,ipc,temp,taup,taub,vel,
     .            zeta,veps,logs,zacc,zsqacc,time,mdtime,it,ekin,tkin,
     .            etot,pext,pint,vol,cons)

C --- write data files ---
      call mdwrit(nbas,nf,ipc,dclabl,lq,qmpol,it,time,ekin,tkin,etot,
     .            pint,cons,pos,alat,vel)

C --- start a new timestep ---
      time = time + tstep

C --- operate with iL_C/2 (advance zeta, logs to t + dt/2; scale v) ---
      if (ensbl .eq. 'NVT') then
        call vvlc2(nbas,nf,temp,tstep,taup,vel,ekin,zeta,logs)
      endif

C --- operate with iL_CP/2 (advance veps, zeta, logs to t + dt/2; 
C                           scale v) ---
      if (ensbl .eq. 'NPT') then
        call vvlcp2(nbas,nf,temp,pext,tstep,taup,taub,tpv,vol,vel,
     .              ekin,zeta,veps,logs,pint)
      endif

C --- operate with iL_1/2 (advance velocities to time t + dt/2) ---
      do  ib = 1, nbas
        ic = ipc(ib)
        do  id = 1, 3
          acc = f(id,ib) / amass(ic)
          vel(id,ib) = vel(id,ib) + 0.5d0 * acc * tstep
        enddo
      enddo

C --- operate with iL_2 (advance positions to t + dt) ---
      if (ensbl .eq. 'NPT') then
C --- modified Verlet for NPT ---
        eps = eps + veps * tstep
        t = 0.5d0 * veps * tstep
        t2 = t * t
C --- Maclaurin expansion of sinh t / t to 8th order ---
        sinhtt = 1d0 + (1d0/6d0)*t2 * (1d0 + (1d0/20d0)*t2
     .                              * (1d0 + (1d0/49d0)*t2
     .                              * (1d0 + (1d0/72d0)*t2 )))
        et  = exp(t)
        et2 = et * et
        do  ib = 1, nbas
          do  id = 1, 3
            pos(id,ib) = pos(id,ib) * et2
     .                 + et * sinhtt * (vel(id,ib) * tstep / alat)
     .                                   
          enddo
        enddo
      else
        do  ib = 1, nbas
          do  id = 1, 3
            pos(id,ib) = pos(id,ib) + vel(id,ib) * tstep / alat
          enddo
        enddo
      endif

C --- get force and come back ---
      end

      subroutine tekin(nbas,ipc,amass,vel,ekin,tkin)
C- For MD, get instantaneous kinetic energy and temperature
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas,ipc,amass,v (velocities)
Co Outputs:
Co   ekin, tkin (actually kT in Ry)
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,ipc(nbas)
      double precision amass(1),vel(3,nbas),ekin,tkin
C Local Variables
      integer id,ib,ic

      ekin = 0d0
      do  ib = 1, nbas
        ic = ipc(ib)
        do  id = 1, 3
          ekin = ekin + amass(ic)*vel(id,ib)*vel(id,ib)
        enddo
      enddo
      tkin = ekin / (3*nbas - 3)
      ekin = 0.5d0 * ekin
      end

      subroutine vvlc2(nbas,nf,temp,tstep,taup,vel,ekin,zeta,logs)
C- "(1/2) L_C" advance zeta, ln s and scale velocities by dt/2
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas, taup (thermostat relaxation time)
Ci   nf: number of degrees of freedom
Ci   tstep: time step (dt)
Ci   temp: target temperature, T_ext, (actually kT in Ry)
Ci   vel:  velocity at time t
Ci   ekin: instantaneous kinetic energy at time t
Ci   zeta: Hoover "viscosity" at time t - dt/2
Ci         (given the symbol \dot\xsi_1 by Martyna 
Ci         (\chi in some of my notes and DL_POLY manual)
Ci   logs: ln s, or \xsi_1 in Martyna, at time t - dt/2
Co Outputs:
Co   vel;  velocity, still at time t, but scaled by the thermostat
Co   ekin; kinetic energy, still at time t, but scaled by the thermostat
Co   tkin; temperature, still at time t, but scaled by the thermostat
Co   zeta: advanced to t
Co   logs: advanced to t
Cr Remarks
Cr    This applies the Liouville operator (1/2) L_NHC (Martyna et al., 
Cr    Mol Phys, 87 (1996) 1117) 
Cr    The thermostat position, ln s, serves as a check on the 
Cr    conservation of the Nose "hamiltonian": we should have
Cr    ln s = \int^t zeta(t) dt (see mdinfo).
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,nf
      double precision temp,taup,tstep,vel(3,nbas),ekin,zeta,logs
C Local Variables
      double precision g,efac,zdt,Q

      Q = taup**2*nf*temp
C --- G_1, equation (25) in Martyna et al. ---
      g = (2d0*ekin - nf*temp) / Q
      zeta = zeta + 0.25d0 * g * tstep
      zdt = zeta * tstep
      efac = exp(-0.5d0 * zdt)
      call dscal(3*nbas,efac,vel,1)
      logs = logs + 0.5d0 * zdt
      efac = exp(-zdt)
      ekin = ekin * efac
      g = (2d0*ekin - nf*temp) / Q
      zeta = zeta + 0.25d0 * g * tstep
      end

      subroutine vvlcp2(nbas,nf,temp,pext,tstep,taup,taub,tpv,vol,vel,
     .                  ekin,zeta,veps,logs,pint)
C- "(1/2) L_CP" advance veps, zeta, ln s and scale velocities by dt/2
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas, taup, taub (thermostat and barostat relaxation times)
Ci   nf: number of degrees of freedom
Ci   tstep: time step (dt)
Ci   temp: target temperature, T_ext, (actually kT in Ry)
Ci   tpv: 3pV, p is the "internal virial" at time t (see tbtote)
Ci   pext: target pressure (in a.u.)
Ci   vel:  velocity at time t
Ci   ekin: instantaneous kinetic energy at time t
Ci   zeta: Hoover "viscosity" at time t - dt/2
Ci         (given the symbol \dot\xsi_1 by Martyna 
Ci         (\chi in some of my notes and DL_POLY manual)
Ci   logs: ln s, or \xsi_1 in Martyna, at time t - dt/2
Ci   veps: barostat viscosity at time t - dt
Co Outputs:
Co   vel;  velocity, still at time t, but scaled by the thermostat
Co   ekin; kinetic energy, still at time t, but scaled by the thermostat
Co   tkin; temperature, still at time t, but scaled by the thermostat
Co   zeta: advanced to t
Co   logs: advanced to t
Co   veps: advanced to t - dt/2
Co   pint: internal pressure at time t (returned for mdinfo)
Cr Remarks
Cr    This applies the Liouville operator (1/2) L_NHCP (Martyna et al., 
Cr    Mol Phys, 87 (1996) 1117) to the Martyna-Tobias-Klein algorithm
Cr    (J Chem Phys 115, 1678 (2001), ibid 101, 4177 (1994))
Cr    The thermostat position, ln s, serves as a check on the 
Cr    conservation of the Nose "hamiltonian": we should have
Cr    ln s = \int^t zeta(t) dt (see mdinfo).
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,nf
      double precision temp,pext,taup,taub,tstep,vel(3,nbas),ekin,zeta,
     .                 veps,logs,tpv,vol,pint
C Local Variables
      double precision g,ge,Q,Wc,a,gkf

      gkf = 1d0 + 1d0/nbas
      Q = taup**2*nf*temp
      Wc = taub**2*(nf + 3)*temp
      g = (2d0*ekin + Wc*veps**2 - 3*nbas*temp) / Q
      ge = (2d0*ekin*gkf + tpv - 3d0*pext*vol) / Wc
      zeta = zeta + 0.25d0 * g * tstep
      veps = veps*exp(-0.25d0*zeta*tstep)
     .     + 0.25d0*tstep*ge*exp(-0.125d0*zeta*tstep)
      logs = logs + 0.5d0 * zeta * tstep
      a = exp(-0.5d0*tstep*(zeta + gkf*veps))
      call dscal(3*nbas,a,vel,1)
      ekin = ekin * a * a
      ge = (2d0*ekin*gkf + tpv - 3d0*pext*vol) / Wc
      veps = veps*exp(-0.25d0*zeta*tstep)
     .     + 0.25d0*tstep*ge*exp(-0.125d0*zeta*tstep)
      g = (2d0*ekin + Wc*veps**2 - 3*nbas*temp) / Q
      zeta = zeta + 0.25d0 * g * tstep
      pint = (2d0*ekin + tpv) / (3d0 * vol)
      end

      subroutine mdinfo(ensbl,nbas,nf,amass,tstep,f,iclas,temp,taup,
     .                  taub,vel,zeta,veps,logs,zacc,zsqacc,time,mdtime,
     .                  it,ekin,tkin,etot,pext,pint,vol,cons)
C- Print out info and monitor conservation in MD
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ensbl: ensemble ('NVE', 'NVT' or 'NPT')
Ci   nbas, amass, tstep, f, iclas, 
Ci   time: accumulated MD time (a.u.)
Ci   mdtime: time requested in ctrl file plus any from previous runs
Ci   nf: degrees of freedom (3*nbas in absence of constraints)
Ci   temp: target temperature (actually kT in Ry)
Ci   taup: Nose-Hoover relaxation time = \sqrt{Q/LkT}, where Q is the
Ci           "mass" or inertia of the thermostat, kT is target temp,
Ci           L is the degrees of freedom (nf)
Ci   taub: barostat relaxation time
Ci   vel: velocities
Ci   zeta: Hoover thermostat friction at time t
Ci   veps: barostat friction at time t
Ci   logs: ln s, "thermostat position;" s is Nose's scaling variable, s.
Ci   zacc: zeta accumulated over time, t - dt
Ci   zsqacc: zeta^2 accumulated over time, t - dt
Ci   ipc: pointers from atom to species
Ci   it: this iteration number
Ci   ekin: kinetic energy
Ci   tkin: current temperature (actually kT in Ry)
Ci   etot: potential energy (total TB or DFT energy) 
Ci   pext: external (applied) pressure (Ry/bohr^3)
Ci   pint: internal pressure
Ci   vol: cell volume (bohr^3)
Co Outputs:
Ci   zacc: zeta accumulated over time, t
Ci   zsqacc: zeta^2 accumulated over time, t
Co   cons: conserved quantity
Cr Remarks
Cr   For NVT, conserved quantity is equation (2) in (Martyna et al., 
Cr   Mol Phys, 87 (1996) 1117) 
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,nf,iclas(nbas),it
      character*3 ensbl
      double precision amass(1),f(3,nbas),vel(3,nbas),temp,ekin,tkin,
     .                 etot,tstep,zeta,veps,logs,time,mdtime,
     .                 pext,pint,vol,taup,taub,zacc,zsqacc
C Local Variables
      integer id,ib,i1mach,iprint
      double precision tzeta,tzetsq,totmom(3),totfrc(3),cons,
     .                 autime,autmp,Q,Wc
      data autime/0.048377d0/
      data autmp/0.6333328d-5/

C --- thermostat inertial "mass" ---
      Q = taup**2*nf*temp

C --- barostat inertial "mass" ---
      Wc = taub**2*(nf + 3)*temp

C --- accumulate zeta and zeta^2 ---
      zacc = zacc + zeta
      zsqacc = zsqacc + zeta*zeta

C --- conserved quantity ---
      cons = etot + ekin
      if (ensbl .eq. 'NVT') then
        cons = cons  + nf*temp*logs + 0.5d0*Q*zeta**2
      endif
      if (ensbl .eq. 'NPT') then
        cons = cons  + (3*nbas + 1)*temp*logs + 0.5d0*Q*zeta**2
     .               + 0.5d0*Wc*veps**2 + pext*vol
      endif
      if (iprint() .lt. 10) return

C --- zeta should avarage to zero ---
      tzeta = zacc / it

C --- mean kinetic energy of the Nose variable, s ---
      tzetsq = 0.5d0 * Q * zsqacc / it

C --- centre of mass and total force ---
      do  id = 1, 3
        totmom(id) = 0d0
        totfrc(id) = 0d0
        do  ib = 1, nbas
          totmom(id) = totmom(id) + amass(iclas(ib))*vel(id,ib)
          totfrc(id) = totfrc(id) + f(id,ib)
        enddo
      enddo
      call awrit6(' MDINFO: at iteration %i accumulated %;2dfs'//
     .            ' MD time of a total %;2dfs%N'//
     .            '         potential energy    = %;8d Ry %N'//
     .            '         kinetic energy      =  %;8d Ry %N'//
     .            '         T+V                 = %;8d Ry',
     .            ' ',1028,i1mach(2),it,time*autime,
     .            mdtime*autime,etot,ekin,etot+ekin)
      if (ensbl .eq. 'NVE') then
        call awrit2('         conserved quantity  = %;8d Ry %N'//
     .              '         temperature         = %;3d K ',' ',
     .              240,i1mach(2),cons,tkin/autmp)
      else
        call awrit3('         conserved quantity  = %;8d Ry %N'//
     .              '         temperature         = %;3d K %N'//
     .              '         target temperature  = %;3d K ',' ',
     .              240,i1mach(2),cons,tkin/autmp,temp/autmp)
      endif
      if (ensbl .eq. 'NPT') then
      call awrit3('         volume              = %;8d a.u. %N'//
     .            '         pressure            = %d a.u. %N'//
     .            '         target pressure     = %;3d a.u. ',' ',
     .            240,i1mach(2),vol,pint,pext)
      endif        
      call awrit6('         total momentum      = %;8d %;8d %;8d %N'//
     .            '         total force         = %;8d %;8d %;8d ',' ',
     .            240,i1mach(2),totmom(1),totmom(2),totmom(3),
     .            totfrc(1),totfrc(2),totfrc(3))
      if (ensbl .eq. 'NVT') then
        call awrit6('   Hoover viscosity, zeta    = %;8d %N'//
     .              '                    <zeta>   = %;8d %N'//
     .              '            (1/2) Q <zeta>^2 = %;8d %N'//
     .              '            Sum zeta(t) dt   = %;8d %N'//
     .              '                     log s   = %;8d %N'//
     .              '  Nose t-scale variable, s   = %;8d %N',
     .              ' ',1028,i1mach(2),zeta,tzeta,tzetsq,
     .              zacc*tstep,logs,exp(logs))
      endif
      if (ensbl .eq. 'NPT') then
        call awrit7('   barostat viscosity, veps  = %;8d %N'//
     .              ' thermostat viscosity, zeta  = %;8d %N'//
     .              '                    <zeta>   = %;8d %N'//
     .              '            (1/2) Q <zeta>^2 = %;8d %N'//
     .              '            Sum zeta(t) dt   = %;8d %N'//
     .              '                     log s   = %;8d %N'//
     .              '  Nose t-scale variable, s   = %;8d %N',
     .              ' ',1028,i1mach(2),veps,zeta,tzeta,tzetsq,
     .              zacc*tstep,logs,exp(logs))
      endif
      end

      subroutine mdwrit(nbas,nf,ipc,dclabl,lq,qmpol,it,time,ekin,tkin,
     .                  etot,pint,cons,pos,alat,vel)
C- Write statistics to MD, MV and XYZ files
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas, nf (number of atoms, degrees of freedom)
Ci   dclabl: d.p. translations of atom labels, passed through to ioxyz
Ci   lq: switch to write out multipole moments
Ci   it: iteration number
Ci   time: accumulated MD time (a.u.)
Co   ekin: Kinetic energy at time t
Co   tkin: Temperature at time t (actually kT in Ry) 
Co   pint: internal pressure at time t (a.u.)
Co   zeta: Hoover friction at time t
Ci   pos: positions at time t
Ci   vel: velocities at time t
Co Outputs:
Co   None
Cr Remarks
Cr   command-line arguments --md=# --mv=# --xyz=# cause output to be
Cr   written for each file at every # iterations (# need not be the 
Cr   same for each file). The md.ext file lists energy, temperature
Cr   and pressure as a function of time for simple plotting; the mv.ext
Cr   has the format of an XBS mv file and the xyz.ext file is for
Cr   programs like Aten or jmol (this file outputs pos in Angstroms).
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,nf,it,ipc(1)
      logical lq
      double precision time,ekin,tkin,etot,cons,pos(3,nbas),vel(3,nbas),
     .                 pint,alat,dclabl(1),qmpol(9,nbas)
C Local Variables
      integer i,ii,id,ib,ifi1,ifi2,ifi3,ip,skip,ipr
      integer fopn,iprint,parg,i1mach
      logical cmdopt
      character*72 outs
      double precision autime,autmp,aupres

C --- first pass through mdwrit ---
      integer first
      data first /1/
      save first

      data autime/0.048377d0/ autmp/0.6333328d-5/ aupres/147.116d0/

      ipr = iprint()

C --- write data to md.ext file ---
      if (cmdopt('--md=',5,0,outs)) then
        ip = 5
        call skipbl(outs,len(outs),ip)
        i = parg(' ',2,outs,ip,len(outs),' ',1,1,ii,skip)
        call rxx(i .ne. 1,' MDWRIT: error parsing --md=')
        ifi1 = fopn('MD')
        if (cmdopt('--st',4,0,outs) .and. it .eq. 1) then
          if (ipr .gt. 10) then
            call awrit1
     .        (' Starting new md file. Write every %i iterations',' ',
     .         128,i1mach(2),skip)
          endif
          if (ipr .gt. 0) then
            call awrit0('#  t (fs)  T (K)  P (Mbar)  V (Ry)  K (Ry)  '//
     .                  'H (Ry)',' ',128,ifi1)
            call awrit0('%% cols 6',' ',128,ifi1)
          endif
        else
          if (ipr .gt. 10) then
            call awrit1
     .        (' Appending to md file. Write every %i iterations',' ',
     .         128,i1mach(2),skip)
          endif
          if (first .eq. 1) then
            call poseof(ifi1)
          endif
        endif
        if (mod(it-1,skip) .eq. 0 .and. ipr .gt. 0) then
          write (ifi1,1) time*autime,tkin/autmp,pint*aupres,etot,ekin,
     .                  cons
        endif
      endif
C --- write to mv.ext file for xbs ---
      if (cmdopt('--mv=',5,0,outs)) then
        ip = 5
        call skipbl(outs,len(outs),ip)
        i = parg(' ',2,outs,ip,len(outs),' ',1,1,ii,skip)
        call rxx(i .ne. 1,' MDWRIT: error parsing --mv=')
        ifi2 = fopn('MV')
        if (cmdopt('--st',4,0,outs) .and. it .eq. 1) then
          if (ipr .gt. 10) then
            call awrit1
     .        (' Starting new mv file. Write every %i iterations',' ',
     .         128,i1mach(2),skip)
          endif
        else
          if (ipr .gt. 10) then
            call awrit1
     .        (' Appending to mv file. Write every %i iterations',' ',
     .         128,i1mach(2),skip)
          endif
          if (first .eq. 1) then
            call poseof(ifi2)
          endif
        endif
        if (mod(it-1,skip) .eq. 0 .and. ipr .gt. 0) then
          call awrit5('frame t= %;2d fs T=%;0dK V= %;4d Ry'//
     .        ' K.E.+V= %;4d Ry C.Q.= %;6d Ry',' ',128,ifi2,time*autime,
     .        tkin/autmp,etot,ekin+etot,cons)
          write (ifi2,2) ((alat*pos(id,ib),id=1,3),ib=1,nbas)
        endif
      endif
C --- write to xyz.ext file for jmol ---
      if (cmdopt('--xyz=',6,0,outs)) then
        ip = 6
        call skipbl(outs,len(outs),ip)
        i = parg(' ',2,outs,ip,len(outs),' ',1,1,ii,skip)
        call rxx(i .ne. 1,' MDWRIT: error parsing --xyz=')
        ifi3 = fopn('XYZ')
        if (cmdopt('--st',4,0,outs) .and. it .eq. 1) then
          if (ipr .gt. 10) then
            call awrit1
     .        (' Starting new xyz file. Write every %i iterations',' ',
     .         128,i1mach(2),skip)
          endif
        else
          if (ipr .gt. 10) then
            call awrit1
     .        (' Appending to xyz file. Write every %i iterations',' ',
     .         128,i1mach(2),skip)
          endif
          if (first .eq. 1) then
            call poseof(ifi3)
          endif
        endif
        if (mod(it-1,skip) .eq. 0 .and. ipr .gt. 0) then
          call awrit1('%i',' ',20,ifi3,nbas)
          call awrit4('t=%;2dfs T=%;0dK V=%;4dRy'//
     .                ' K.E.+V=%;4dRy',' ',72,ifi3,time*autime,
     .                tkin/autmp,etot,ekin+etot)
          call ioxyz(nbas,ipc,pos,alat,vel,qmpol,dclabl,lq,.true.,ifi3)
        endif
      endif
      first = 0
C --- flush buffers ---
      if (cmdopt('--flush',7,0,outs)) then
        if (ipr .gt. 10) then
          call awrit0
     .        (' MDWRIT: flush buffers ..',' ',128,i1mach(2))
        endif
        call flushs(-1)
      endif
    1 format(f12.3,5f12.6)
    2 format(8f10.5)
      end

      subroutine ioxyz(nbas,ipc,pos,alat,vel,qmpol,dclabl,sc,md,ifi)
C- write to XYZ file for Jmol etc.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas,ipc,pos,alat,vel,qmpol,dclabl,sc,md,ifi
Co Outputs:
Co   none
Cr Remarks
Cr   Writes to a file xyz.ext for use by visualisation programs like
Cr   Aten and Jmol; and for use by programs that analyse trajectories
Cr   (correlation functions etc.) In the case of MD, velocities are
Cr   written. For self consistent TB, the multipole moments are also
Cr   written. Positions are converted to Angstroms, all other
Cr   quantities in a.u.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,ifi,ipc(1)
      logical sc,md
      double precision pos(3,nbas),alat,vel(3,nbas),qmpol(9,nbas),
     .                 dclabl(1)
C Local Variables
      integer ib,is
      double precision apos(3),avel(3)
      character*8 clabl

      do  ib = 1, nbas
        is = ipc(ib)
        call r8tos8(dclabl(is),clabl)
        call dpcopy(pos(1,ib),apos,1,3,alat*0.529177249d0)
        if (md) then
          call dpcopy(vel(1,ib),avel,1,3,1d4)
          if (sc) then
            call awrit3(clabl//'%3:-2,4;4#4d %3:-2,4;4#4d %9:-2,6;6#4d',
     .                  ' ',256,ifi,apos,avel,qmpol(1,ib))
          else
            call awrit2(clabl//'%3:-2,4;4#4d %3:-2,4;4#4d',
     .                  ' ',256,ifi,apos,avel)
          endif
        else
          if (sc) then
            call awrit2(clabl//'%3:-2,4;4#4d %9:-2,6;6#4d',
     .                  ' ',256,ifi,apos,qmpol(1,ib))
          else
            call awrit1(clabl//'%3:-2,4;4#4d',' ',256,ifi,apos)
          endif
        endif
      enddo
      end

      subroutine initv(ensbl,nbas,nf,tstep,ipc,amass,temp,pext,taup,
     .                 taub,nspec,lsp8,spid,dclabl,zeta,zacc,zsqacc,
     .                 veps,vel)
C- Set up the velocities on a Gaussian distribution.
C-----------------------------------------------------------------------
Ci Inputs
Ci   nbas
Ci   tstep: time step
Ci   ipc: class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   amass: atomic mass (in a.u.)
Ci   temp: target temperature, (actually kT in a.u.)
Ci   pext: target pressure (Ry/bohr^3)
Ci   taup: relaxation time of the thermostat (a.u.)
Ci   taub: relaxation time of the barostat (a.u)
Ci   spid / dclabl : species label in char or real*8
Ci   lsp8  : 1 if to use dclabl, 0 if to use spid (merge mol with tbe)
Co Outputs
Co   nf: number of degrees of freedom (3*nbas in absence of constraints)
Co   zeta: "hoover viscosity"
Co   zacc,zsqacc: accumulators which monitor the conservation of the
Co                Nose "hamiltonian"
Co   veps: barostat viscosity (v_\epsilon in Martyna et al.)
Co   vel: velocities
Cr Remarks
C-----------------------------------------------------------------------
C     Passed parameters
C     implicit none
      integer nbas,nf,ipc(nbas),nspec,lsp8
      character*3 ensbl
      double precision vel(3,nbas),temp,pext,amass(1),dclabl(1),
     .                 zeta,tstep,taup,taub,zacc,zsqacc,veps
      character*8 spid(1), clabl
C     Local Variables
      integer n,id,i1mach,iprint,ic
      real r1,r2,ran1
      double precision vscale,scale,tkin,tmp0,cs,speed,twopi,ekin,
     .     autime,autmp,Q,Wc
      data autime/0.048377d0/
      data autmp/0.6333328d-5/

      twopi = 8d0*datan(1d0)

C --- Initialize the Hoover viscosity and the accumulators ---
      zeta = 0d0
      zacc = zeta
      zsqacc = zeta*zeta

C --- Initialise barostat viscosity ---
      veps = 0

C --- Give the particles a Maxwellian velocity distribution
      call ran1in(1)
      do  id = 1, 3
        do  n = 1, nbas
          vscale = dsqrt(temp/amass(ipc(n)))
          r1 = ran1()
          r2 = ran1()
          cs = dcos(twopi*dble(r1))
          speed = dsqrt(-2d0*dble(log(r2)))
          vel(id,n) = vscale*speed*cs
        enddo
      enddo

C --- Compute the kinetic energy ---
      call tekin(nbas,ipc,amass,vel,ekin,tkin)
      nf = 3 * nbas - 3
      
C --- Scale the velocities ---
      if (ekin .gt. 0d0) then
        tmp0 = 2d0*ekin/(3d0*(nbas-1))
        scale = dsqrt(temp/tmp0)
      else
        scale = 0d0
      endif
      do  id = 1, 3
        do  n = 1, nbas
          vel(id,n) = scale*vel(id,n)
        enddo
      enddo
      call tekin(nbas,ipc,amass,vel,ekin,tkin)

C --- Set to zero the velocity of the centre of mass ---
      call zercmv(nbas,vel,amass,ipc)

C --- print out ---
      if (iprint() .ge. 10) then
        call awrit2(' INITV: time-step = %;3dfs (%d a.u.)',' ',120,
     .              i1mach(2),tstep*autime,tstep)
        write(*,200)
C --- 'atomic weight' is weight in grams of Avogadro's no. of atoms ---
        do  ic = 1, nspec
          if (lsp8 .eq. 0) then
            clabl = spid(ic)
          else
            call r8tos8(dclabl(ic),clabl)
          endif
          write(*,100) clabl,amass(ic),amass(ic)*1.09716d-3
        enddo
        if (ensbl .eq. 'NVT') then
          Q = taup**2*nf*temp
          call awrit6('          Thermostat ON       Barostat OFF %N'//
     .    '       initial Maxwell temperature = %;3dK,%N'//
     .    '       target temperature = %;3dK %N'//
     .    '       initial Hoover viscosity zeta = %;8d%N'//
     .    '       Nose thermostat inertia, Q=%g a.u., relaxation '//
     .    'time tau=%;0dfs (%d a.u.)%N ',' ',1028,i1mach(2),tkin/autmp,
     .      temp/autmp,zeta,Q,taup*autime,taup)
        endif
        if (ensbl .eq. 'NPT') then
          Q = taup**2*nf*temp
          Wc = taub**2*(nf + 3)*temp
          call awrit6('          Thermostat ON       Barostat ON %N'//
     .    '       initial Maxwell temperature = %;3dK,%N'//
     .    '       target temperature = %;3dK %N'//
     .    '       initial Hoover viscosity zeta = %;8d%N'//
     .    '       Nose thermostat inertia, Q=%g a.u., relaxation '//
     .    'time tau=%;0dfs (%d a.u.)%N ',' ',1028,i1mach(2),tkin/autmp,
     .      temp/autmp,zeta,Q,taup*autime,taup)
          call awrit5('       target pressure = %;3d a.u. %N'//
     .    '       initial barostat viscosity veps = %;8d%N'//
     .    '       barostat inertia, W_c=%g a.u., relaxation '//
     .    'time tau=%;0dfs (%d a.u.)%N ',' ',1028,i1mach(2),pext,veps,
     .      Wc,taub*autime,taub)
        endif
        if (ensbl .eq. 'NVE') then 
          call awrit0('          Thermostat OFF       Barostat OFF%N',
     .                ' ',120,i1mach(2))
        endif
      endif
  100 format (11x,a4,2x,f15.6,6x,f10.5/)
  200 format (7x,' species  atomic mass (a.u.)  ''atomic weight''')
      end

