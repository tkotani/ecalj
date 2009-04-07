      subroutine packv(lpack,scall,struc,listi,is1,is2,x1)
C- Pack/unpack vector elements in a struc
Ci ----------------------------------------------------------------
Ci Inputs
Ci   lpack  1s digit 
Ci           0 unpacks, 1 pack
Ci         10s digit 
Ci           0 x1 is same for all species (packing only)
Ci           1 x1 is different for each species
Ci           2 x1 and its size is different for each species
Ci             (not implemented)
Ci   scall   subroutine name that returns offset, cast and no. elts.
Ci   listi   which elements to unpack
Ci   is1,is2 range of species over which to copy
Cio Inputs/Outputs 
Cio  struc   structure
Cio  x1      contains elements to pack/unpack
Cr Remarks
Cr   packv packs/upacks a entries over a group of species in
Cr   structure 'struc'.  The structure-specific information it requires, 
Cr   (offset, cast and number of entries for a given element), by
Cr   calling routine 'scall', whose name is passed to packv.  scall
Cr   must have a standard argument list, is
Cr     subroutine scall(struc,listi,nlisti,is,offset,cast,nelt)
C ----------------------------------------------------------------
C     implicit none
      integer lpack,is1,is2,listi,x1(1)
      double precision struc(1)
      external scall
C Local variables
      logical pack
      integer offi,casti,nelti,off,i,ssize,nsmx,is,incx,nel,offx
      parameter (nsmx=2048)
      double precision xx(nsmx)

C ... Get the list of offsets,casts,number of entries
      call scall(struc,listi,1,is1,offi,casti,nelti)
      ssize = nint(struc(1))
      pack = mod(lpack,10) .eq. 1
      off = offi + ssize*(is1-1)

C ... Case nelti is 1
      if (nelti .eq. 1) then
        incx = 0
        if (mod(lpack/10,10) .eq. 1) incx=1
        nel = is2-is1+1
        if (pack) then
          if (casti .eq. 4) then
            call dcopy(nel,x1,incx,struc(off),ssize)
          else
            i = nel
            if (incx .eq. 0) i = 1
            call idscop(x1,xx,i,1,1)
            call dcopy(nel,xx,incx,struc(off),ssize)
          endif
        else
          if (incx .eq. 0) nel = 1
          if (casti .eq. 4) then
            call dcopy(nel,struc(off),ssize,x1,incx)
          else
            if (nel .gt. nsmx) call rxi('packv: increase nsmx to',nel)
            call dcopy(nel,struc(off),ssize,xx,incx)
            call discop(xx,x1,nel,1,1,0)
          endif
        endif
C ... If nelti is not one, copy elements one spec at a time
      else
        if (pack) then
          offx = 1
          do  10  is = is1, is2
            off = offi + ssize*(is-1)
            if (mod(lpack/10,10) .eq. 1) offx = 1 + (is-is1)*nelti
            if (casti .eq. 4) then
              call dpscop(x1,struc,nelti,offx,off,1d0)
            else
              call idscop(x1,struc,nelti,offx,off)
            endif
   10     continue
        else
          do  20  is = is1, is2
            off = offi + ssize*(is-1)
            offx = 1 + (is-is1)*nelti
            if (casti .eq. 4) then
              call dpscop(struc,x1,nelti,off,offx,1d0)
            else
              call discop(struc,x1,nelti,off,offx,0)
            endif
   20     continue
        endif
      endif
      end
      

