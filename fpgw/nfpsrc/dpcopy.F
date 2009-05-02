      subroutine dpcopy(afrom,ato,n1,n2,fac)
      real*8 afrom(1),ato(1),fac
      do 1 i=n1,n2
        ato(i)=fac*afrom(i)
    1 continue
      return
      end
