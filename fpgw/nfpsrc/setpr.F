C   $Id: setpr.f,v 1.2 1998/04/28 18:45:51 coa Exp $
c  defines:   fct iprint(), setpr(i), getpr(i), pushpr(i), poppr
      integer function iprint()
      implicit real*8 (a-h,p-z), integer(o)
      integer setpr, getpr, pushpr, poppr
      parameter( nx=25 )
      dimension kpr(nx)
      save
      data ip /1/, kpr /nx*30/
      iprint=kpr(ip)
      return

      entry setpr(i3)
      ip=1
      kpr(1)=i3
      return

      entry getpr(i1)
      i1=kpr(ip)
      return

      entry pushpr(i2)
      ip=ip+1
      if(ip.gt.nx) call rx('pushpr: max depth exceeded')
      kpr(ip)=i2
      return

      entry poppr
      ip=ip-1
      if(ip.le.0) call rx('poppr: stack empty')
      return

      end
