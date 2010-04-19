      subroutine zercmv(nbas,vel,amass,iclas)
C- calculate centre of mass velocity and subtract from all v(i)'s
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas, vel, amass, iclas
Co Outputs:
Co   vel: velocities shifted to zero the total momentum
Cr Remarks
Cr  zero the CM velocity (i.e. set total momentum to zero)
C ----------------------------------------------------------------------
C     implicit none

C Passed Parameters
      integer nbas,iclas(nbas)
      double precision vel(3,nbas),amass(nbas)
C Local Variables
      integer n
      double precision v1,v2,v3,wtsum

      v1 = 0d0
      v2 = 0d0
      v3 = 0d0
      wtsum = 0d0
      do  n = 1, nbas
        wtsum = wtsum + amass(iclas(n))
        v1 = v1 + amass(iclas(n))*vel(1,n)
        v2 = v2 + amass(iclas(n))*vel(2,n)
        v3 = v3 + amass(iclas(n))*vel(3,n)
      enddo

      v1 = v1/wtsum
      v2 = v2/wtsum
      v3 = v3/wtsum

      do  n = 1, nbas
        vel(1,n) = vel(1,n) - v1
        vel(2,n) = vel(2,n) - v2
        vel(3,n) = vel(3,n) - v3
      enddo
      end

