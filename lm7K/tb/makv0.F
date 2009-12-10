      subroutine makv0(d,v0,erep,derep,p)
C- Makes pair repulsive interaction
C ----------------------------------------------------------------------
Ci Inputs
Ci   v0,d
Co Outputs
Co   erep: repulsive interaction and its derivative and pressure
Cr Remarks
Cr   V0 is a vector length 9 allowing the pairwise energy to be a sum of
Cr   three terms of the form A r^B exp(-Cr). If the first power exponent
Cr   is positive, then the pair potential is of Chadi's form, namely
Cr   A1 eps + a2 eps^2: the third number in each set is the equilibrium
Cr   bond length. If the first power exponent is positive and the third
Cr   number of the set (C) is negative then the Goodwin-Skinner-Pettifor
Cr   form is used: A (r0/d)^m exp[m (-{d/rc}^mc + {r0/rc}^mc)].
Cr   In the GSP there is a remaining two parameters. These are used
Cr   to add an additional power or exponential depending whether
Cr   V0(3,3) is positive or negative
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer i
      double precision d,v0(3,3),erep,derep,p
C Local paramters
      double precision e,de,eps,m,mc,r0,rc

C      print 333, v0
C  333 format(' makv0: V0=',9f8.2)

      erep = 0
      derep = 0d0

      if (v0(2,1) .gt. 0) then
        if (v0(3,1) .lt. 0) then
          m  = v0(1,2)
          mc = v0(2,2)
          r0 = v0(3,2)
          rc = v0(1,3)
          e = v0(1,1)*r0**m*exp(m*(r0/rc)**mc)
          erep = e*exp(-m*(d/rc)**mc)/d**m
          derep = -m*erep*(1d0 + mc*(d/rc)**mc)/d
          if (v0(2,3) .gt. 0) then
            if (v0(3,3) .gt. 0) then
              e = v0(2,3)*dexp(-v0(3,3)*d)
              de = -e*v0(3,3)
            else
              e = v0(2,3)**d**v0(3,3)
              de = e*v0(3,3)/d
            endif
            erep = erep + e
            derep = derep + de
          endif
        else
          eps   = (d - v0(3,1))/v0(3,1)
          erep  = v0(1,1)*eps + v0(1,2)*eps**2
          derep = (v0(1,1) + 2d0*v0(1,2)*eps)/v0(3,1)
        endif
      else
        do  2  i = 1, 3
          e = v0(1,i)*d**v0(2,i)*dexp(-v0(3,i)*d)
          erep = erep + e
          derep = derep + e*(v0(2,i)/d - v0(3,i))
    2   continue
      endif
      p = derep*d
      end
