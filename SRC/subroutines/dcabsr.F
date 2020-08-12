      double precision function dcabsr(zr,zi)
C- real-arithmetic analog of dcabs1
      double precision zr,zi
      dcabsr = dabs(zr) + dabs(zi)
      end
      subroutine cpy(tr,ti,dr,di,t1,t2)
C- complex multiply (t1,t2) = (tr,ti) * (dr,di) 
      double precision tr,ti,dr,di,t1,t2      
      double precision tmp
      tmp = tr*dr - ti*di
      t2  = tr*di + ti*dr
      t1 = tmp
      end

