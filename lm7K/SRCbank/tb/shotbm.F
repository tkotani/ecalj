      subroutine shotbm(tabme,decay,v0,nset,npair,nterm,nl,nsp,nlmesp,
     .  memode,nclass,dclabl,iam,ipair,strn)
C- Print table of tight-binding matrix elements and pair potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   tabme,decay,v0,nterm,memode,iam: see subroutine rdtbh
Ci   nset: number of matrix element sets
Ci   npair: number of matrix element pairs
Ci   nl,nsp,nlmesp,nclass,dclabl
Ci   ipair: 1, print pair potential; 0, do not
Ci   strn: identifier string for MEs (printed)
Cu Updates
Cu    8 Jun 07 (MvS) Merged Klepeis's additions to TB package
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nset,npair,nterm,nl,nsp,nlmesp,memode,nclass,ipair
      integer iam(3,1)
      double precision tabme(nterm,nlmesp,nset),v0(3,3,nset),
     .  decay(nlmesp,nset),dclabl(nclass)
      character*(*) strn
C Local variables
      integer ip,np,ic,jc,k,l,m,nme,nlme,lgunit
      double precision wk(10)
      character clabli*8,clablj*8,strng*40
C      integer nlme(3)
C      data nlme /1,4,10/

      if (nset .eq. 0) return

C --- Print pair potential ---
      if (ipair .eq. 1) then
        call awrit1('%N SHOTBM: pair potential, %i rules:'
     .    ,' ',80,lgunit(1),nset)
        write (*,100)
        do  20  k = 1, nset
          np = 0
          do  10  ip = 1, npair
C       ... Find the rule connecting this pair
            if (iam(3,ip) .eq. k) then
              np = np + 1
              if (np .eq. 1) then
                ic = iam(1,ip)
                jc = iam(2,ip)
              endif
            endif
   10     continue
          if (np .eq. 0) then
            call awrit1(' ... no pairs matching rule %i',' ',80,
     .        lgunit(1),k)
            goto 20
          endif
          call r8tos8(dclabl(ic),clabli)
          call r8tos8(dclabl(jc),clablj)
          call awrit2(' '//clabli//'%a,'//clablj//
     .      '%a:%?#n>1# (%i more)',' ',80,lgunit(1),np,np-1)
          if (v0(2,1,k) .gt. 0) then
            if (v0(3,1,k) .lt. 0) then
              write(*,115) v0(1,1,k),(v0(l,2,k), l = 1, 3),v0(1,3,k)
            else
              write(*,120) v0(1,1,k),v0(1,2,k),v0(3,1,k)
            endif
          else
            write(*,130) ((v0(l,m,k), l = 1, 3), m = 1, 3)
          endif
   20   continue
      endif

  100 format(
     .  24x,'a1',4x,'b1',4x,'c1',6x,'a2',4x,'b2',4x,'c2',6x,'a3',
     .   4x,'b3',4x,'c3')
C 110 format(1x,a,',',a,':',:,'  (',i4,' more)')
  115 format('   A(r0/d)^m e^...:',1x,'A=',f7.3,2x,'m=',f6.2,
     .  2x,'mc=',f6.2,2x,'r0=',f7.3,2x,'rc=',f7.3,2x,'GSP')
  120 format('   a1 ep + a2 ep^2:',f8.3,12x,f8.3,16x,'d0=',f7.3)
  130 format('   a d^b exp(-c d):',3(f8.3,f6.2,f6.3))

      if (memode .ge. 10) return

C --- Print MEs ---
      strng = strn
      call awrit1('%N SHOTBM: %i rules for '//strng//'%a',' ',80,
     .  lgunit(1),nset)
      if (memode .eq. 4) then
        print 150
      elseif (memode .eq. 5) then
        print 155
      elseif (memode .eq. 6) then
        print 156
      elseif (memode .eq. 7) then
        print 157
      else
        if (nl .ge. 3) then
          print 160
        elseif (nl .ge. 2) then
          print 161
        else
          print 162
        endif
      endif
      nme = nlme(nl)
      if (memode .eq. 1) call ume(nme,1d0,wk)
      do  40  k = 1, nset
        np = 0
        do  30  ip = 1, npair
          if (iam(3,ip) .eq. k) then
            np = np + 1
            if (np .eq. 1) then
              ic = iam(1,ip)
              jc = iam(2,ip)
            endif
          endif
   30   continue
        if (np .eq. 0) then
          call awrit1(' ... no pairs matching rule %i',' ',80,
     .      lgunit(1),k)
          goto 40
        endif
        call r8tos8(dclabl(ic),clabli)
        call r8tos8(dclabl(jc),clablj)
        call awrit2(' '//clabli//'%a,'//clablj//
     .    '%a:%?#n>1# (%i more)',' ',80,lgunit(1),np,np-1)
        if (memode .eq. 4) then
          write (*,170) (tabme(m,1,k),m=1,9)
          if (nsp .eq. 2) then
            write (*,180) (tabme(m,1+1*nme,k),m=1,9),
     .                    (tabme(m,1+2*nme,k),m=1,9)
          endif
          if (nl .ge. 2) then
            if (nsp .eq. 2) then
              write (*,190) ((tabme(m,l+0*nme,k),m=1,9),
     .                       (tabme(m,l+1*nme,k),m=1,9),
     .                       (tabme(m,l+2*nme,k),m=1,9),l=2,4)
            else
              write (*,200) ((tabme(m,l,k),m=1,9),l=2,4)
            endif
          endif
          if (nl .ge. 3) then
            if (nsp .eq. 2) then
              write (*,210) ((tabme(m,l+0*nme,k),m=1,9),
     .                       (tabme(m,l+1*nme,k),m=1,9),
     .                       (tabme(m,l+2*nme,k),m=1,9),l=5,10)
            else
              write (*,220) ((tabme(m,l,k),m=1,9),l=5,10)
            endif
          endif
        elseif (memode .eq. 5) then
          write (*,175) (tabme(m,1,k),m=1,5)
          if (nsp .eq. 2) then
            write (*,185) (tabme(m,1+1*nme,k),m=1,5),
     .                    (tabme(m,1+2*nme,k),m=1,5)
          endif
          if (nl .ge. 2) then
            if (nsp .eq. 2) then
              write (*,195) ((tabme(m,l+0*nme,k),m=1,5),
     .                       (tabme(m,l+1*nme,k),m=1,5),
     .                       (tabme(m,l+2*nme,k),m=1,5),l=2,4)
            else
              write (*,205) ((tabme(m,l,k),m=1,5),l=2,4)
            endif
          endif
          if (nl .ge. 3) then
            if (nsp .eq. 2) then
              write (*,215) ((tabme(m,l+0*nme,k),m=1,5),
     .                       (tabme(m,l+1*nme,k),m=1,5),
     .                       (tabme(m,l+2*nme,k),m=1,5),l=5,10)
            else
              write (*,225) ((tabme(m,l,k),m=1,5),l=5,10)
            endif
          endif
        elseif (memode .eq. 7) then
          write (*,177) (tabme(m,1,k),m=1,4)
          if (nsp .eq. 2) then
            write (*,187) (tabme(m,1+1*nme,k),m=1,4),
     .                    (tabme(m,1+2*nme,k),m=1,4)
          endif
          if (nl .ge. 2) then
            if (nsp .eq. 2) then
              write (*,197) ((tabme(m,l+0*nme,k),m=1,4),
     .                       (tabme(m,l+1*nme,k),m=1,4),
     .                       (tabme(m,l+2*nme,k),m=1,4),l=2,4)
            else
              write (*,207) ((tabme(m,l,k),m=1,4),l=2,4)
            endif
          endif
          if (nl .ge. 3) then
            if (nsp .eq. 2) then
              write (*,217) ((tabme(m,l+0*nme,k),m=1,4),
     .                       (tabme(m,l+1*nme,k),m=1,4),
     .                       (tabme(m,l+2*nme,k),m=1,4),l=5,10)
            else
              write (*,227) ((tabme(m,l,k),m=1,4),l=5,10)
            endif
          endif
        elseif (memode .eq. 1) then
          write (*,230) (wk(l),l=1,nme)
          if (nsp .eq. 2) then
            write (*,240) (wk(l),l=1,nme)
            write (*,250) (wk(l),l=1,nme)
          endif
        elseif (memode .eq. 6) then
          write (*,230) (tabme(1,l,k),l=1,8)
          if (nsp .eq. 2) then
            write (*,240) (tabme(1,l+1*8,k),l=1,8)
            write (*,250) (tabme(1,l+2*8,k),l=1,8)
          endif
        else
          write (*,230) (tabme(1,l,k),l=1,nme)
          if (nsp .eq. 2) then
            write (*,240) (tabme(1,l+1*nme,k),l=1,nme)
            write (*,250) (tabme(1,l+2*nme,k),l=1,nme)
          endif
          if (memode .eq. 2 .or. memode .eq. 3 .or. memode .eq. 7) then
            write (*,260) (decay(l,k),l=1,nme)
            if (nsp .eq. 2) then
              write (*,240) (decay(l+1*nme,k),l=1,nme)
              write (*,250) (decay(l+2*nme,k),l=1,nme)
            endif
          endif
        endif
   40 continue

* 140 format(/' SHOTBM: found ',i4,' sets of ',a,' matrix elements:')
  150 format(13x,'a1',5x,'b1',5x,'c1',7x,'a2',5x,'b2',5x,'c2',7x,'a3',
     .   5x,'b3',5x,'c3')
  155 format(13x,'V',8x,'n',7x,'nc',7x,'r0',7x,'rc')
  156 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp',4x,'usp',4x,'uxy',3x,
     .       'upps',3x,'uppp')
  157 format(13x,'a',8x,'b',8x,'c',7x,'d0')
  160 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp',4x,'sds',4x,'pds',4x,
     .       'pdp',4x,'dds',4x,'ddp',4x,'ddd')
  161 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp')
  162 format(13x,'sss')
C  163 format(13x,'sss',4x,'sps',4x,'pps',4x,'ppp',4x,'usp',4x,'uxy',
C     .       4x,'upps',4x,'uppp')
  170 format('   sss:',3(f9.3,f7.2,f7.3))
  175 format('   sss:',5f9.3)
  177 format('   sss:',4f9.3)
  180 format('    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3))
  185 format('    --:',5f9.3/'    +-:',5f9.3)
  187 format('    --:',4f9.3/'    +-:',4f9.3)
  190 format('   sps:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   pps:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3)/
     .       '   ppp:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3))
  195 format('   sps:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   pps:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3/
     .       '   ppp:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3)
  197 format('   sps:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   pps:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3/
     .       '   ppp:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3)
  200 format('   sps:',3(f9.3,f7.2,f7.3)/'   pps:',3(f9.3,f7.2,f7.3)/
     .       '   ppp:',3(f9.3,f7.2,f7.3))
  205 format('   sps:',5f9.3/'   pps:',5f9.3/
     .       '   ppp:',5f9.3)
  207 format('   sps:',4f9.3/'   pps:',4f9.3/
     .       '   ppp:',4f9.3)
  210 format('   sds:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   pds:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3)/
     .       '   pdp:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   dds:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3)/
     .       '   ddp:',3(f9.3,f7.2,f7.3)/'    --:',3(f9.3,f7.2,f7.3)/
     .       '    +-:',3(f9.3,f7.2,f7.3)/'   ddd:',3(f9.3,f7.2,f7.3)/
     .       '    --:',3(f9.3,f7.2,f7.3)/'    +-:',3(f9.3,f7.2,f7.3))
  215 format('   sds:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   pds:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3/
     .       '   pdp:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   dds:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3/
     .       '   ddp:',5f9.3/'    --:',5f9.3/
     .       '    +-:',5f9.3/'   ddd:',5f9.3/
     .       '    --:',5f9.3/'    +-:',5f9.3)
  217 format('   sds:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   pds:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3/
     .       '   pdp:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   dds:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3/
     .       '   ddp:',4f9.3/'    --:',4f9.3/
     .       '    +-:',4f9.3/'   ddd:',4f9.3/
     .       '    --:',4f9.3/'    +-:',4f9.3)
  220 format('   sds:',3(f9.3,f7.2,f7.3)/'   pds:',3(f9.3,f7.2,f7.3)/
     .       '   pdp:',3(f9.3,f7.2,f7.3)/'   dds:',3(f9.3,f7.2,f7.3)/
     .       '   ddp:',3(f9.3,f7.2,f7.3)/'   ddd:',3(f9.3,f7.2,f7.3))
  225 format('   sds:',5f9.3/'   pds:',5f9.3/
     .       '   pdp:',5f9.3/'   dds:',5f9.3/
     .       '   ddp:',5f9.3/'   ddd:',5f9.3)
  227 format('   sds:',4f9.3/'   pds:',4f9.3/
     .       '   pdp:',4f9.3/'   dds:',4f9.3/
     .       '   ddp:',4f9.3/'   ddd:',4f9.3)
  230 format('   coeff:',10f7.2)
  240 format('      --:',10f7.2)
  250 format('      +-:',10f7.2)
  260 format('   decay:',10f7.2)

      end
