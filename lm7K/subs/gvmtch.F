Cdelw1 var ['ogg', '0', 'ogg', 'ngs', 'real(8)', 'rv_w_', 'gvmtch']
Cdelw1 var ['okk', '0', 'okk', 'ngs', 'real(8)', 'rv_w_', 'gvmtch']
Cdelw1 var ['oiwk', '0', 'oiwk', 'ngs', 'real(8)', 'rv_w_', 'gvmtch']
Cdelw1 var ['ogvs', '-1', 'undef', '0', 'undef', 'x', 'gvmtch']
Cdelw1 var ['ogvb', '-1', 'undef', '0', 'undef', 'x', 'gvmtch']
Cdelw1 var ['okvb', '-1', 'undef', '0', 'undef', 'x', 'gvmtch']
Cdelw1 undel [['call', 'undef', 'wref'], 'ogvs', 'gvmtch']
Cdelw1 undel [['call', 'undef', 'wref'], 'ogvb', 'gvmtch']
Cdelw1 undel [['call', 'undef', 'wref'], 'okvb', 'gvmtch']
      subroutine gvmtch(ngs,ogvs,okvs,ngb,ogvb,okvb)
C-  Ensure lower block of big list matches lower block of small list
C     implicit none
Cdelw1 do not change  [['call', 'undef', 'wref'], 'ogvs', 'gvmtch']
Cdelw1 do not change  [['call', 'undef', 'wref'], 'ogvb', 'gvmtch']
Cdelw1 do not change  [['call', 'undef', 'wref'], 'okvb', 'gvmtch']
      integer ngs,ngb,ogvs,okvs,ogvb,okvb
Cdelw1       integer ogg,okk,oiwk
       real(8) ,allocatable :: rv_w_ogg(:)
       real(8) ,allocatable :: rv_w_okk(:)
       real(8) ,allocatable :: rv_w_oiwk(:)

      integer w(1)
      common /w/ w

Cdelw1       call defrr(ogg,   ngs)
       allocate(rv_w_ogg(ngs))
       if (ngs<0) rv_w_ogg(:)=0.0d0

Cdelw1       call defrr(okk,   ngs)
       allocate(rv_w_okk(ngs))
       if (ngs<0) rv_w_okk(:)=0.0d0

Cdelw1       call defrr(oiwk,  ngs)
       allocate(rv_w_oiwk(ngs))
       if (ngs<0) rv_w_oiwk(:)=0.0d0

Cdelw1 warning(1) , probably  ogvs  is not defined yet at linenumber= 12
Cdelw1 do not change ogvs because of ['call', 'undef', 'wref']
Cdelw1 warning(1) , probably  ogvb  is not defined yet at linenumber= 12
Cdelw1 do not change ogvb because of ['call', 'undef', 'wref']
Cdelw1 warning(1) , probably  okvb  is not defined yet at linenumber= 12
Cdelw1 do not change okvb because of ['call', 'undef', 'wref']
Cdelw1 warning(1) , probably  ogvs  is not defined yet at linenumber= 12
Cdelw1 do not change ogvs because of ['call', 'undef', 'wref']
Cdelw1 warning(1) , probably  ogvb  is not defined yet at linenumber= 12
Cdelw1 do not change ogvb because of ['call', 'undef', 'wref']
Cdelw1 warning(1) , probably  okvb  is not defined yet at linenumber= 12
Cdelw1 do not change okvb because of ['call', 'undef', 'wref']
Cdelw1       call pgvmat(w(ogg),w(okk),w(oiwk),ngs,w(ogvs),ngb,w(ogvb),w(okvb))
       call pgvmat ( rv_w_ogg , rv_w_okk , rv_w_oiwk , ngs , w ( ogvs 
     . ) , ngb , w ( ogvb ) , w ( okvb ) ) 

Cdelw1 rlse name= ogg old_list= ogg okk oiwk 
Cdelw1 rlse name= ogg new_list= (None)
Cdelw1       call rlse(ogg)
       if (allocated(rv_w_oiwk)) deallocate(rv_w_oiwk)
       if (allocated(rv_w_okk)) deallocate(rv_w_okk)
       if (allocated(rv_w_ogg)) deallocate(rv_w_ogg)


Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ogvs ogvb okvb]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end
      subroutine pgvmat(gg,kk,iwk,ngs,gvs,ngb,gvb,kvb)
C     implicit none
      integer ngs,ngb
      integer kvb(ngb,3),kk(ngs),iwk(ngs)
      double precision gvs(ngs,3),gvb(ngb,3),gg(ngs)
      integer ig,low,high,jg,mm
      double precision xx,tol,tol2
      parameter (tol=1d-6,tol2=1d-9)

C      call prm('gvs',gvs,ngs,ngs,3)
C      call prm('gvb',gvb,ngb,ngb,3)

C ... Generate length of g and ensure length matches for both
      do  10  ig = 1, ngs
   10 gg(ig) = gvs(ig,1)**2 + gvs(ig,2)**2 + gvs(ig,3)**2

C     call prm('gg',gg,ngs,ngs,1)


C ... Sanity check.  Make sure |G| is same for each vector in both lists
C      do  20  ig = 1, ngs
C        xx = gvb(ig,1)**2 + gvb(ig,2)**2 + gvb(ig,3)**2
C        if (abs(xx-gg(ig)) .gt. tol) then
C          print *, xx,gg(ig)
C          call rxi('bug in gvmatch, ig=',ig)
C        endif
C   20 continue

C --- For each G in small, find matching G in big list ---
      xx = -1d0
      low = 0
      high = 0
      do  30  ig = 1, ngs
        iwk(ig) = -1
C   ... Find first and last g-vector list with same length
        if (abs(xx-gg(ig)) .gt. tol) then
          call huntx(gg,ngs,gg(ig)+tol,0,high)
          low = ig
          high = min(high,ngs)
          xx = gg(ig)
        endif

        do  32  jg = low, high
          do  34  mm = 1, 3
   34     if (abs(gvb(jg,mm)-gvs(ig,mm)) .gt. tol2) goto 32
C     ... Found a match
          iwk(ig) = jg
   32   continue
C   ... Sanity check
        if (iwk(ig) .eq. -1) call rxi('bug in gvmatch, ig=',ig)
   30 continue

C ... Rearrange gvb, kvb
      do  40  mm = 1, 3
        do  42  ig = 1, ngs
          jg = iwk(ig)
          gg(ig) = gvb(jg,mm)
          kk(ig) = kvb(jg,mm)
   42   continue
        do  44  ig = 1, ngs
          gvb(ig,mm) = gg(ig)
          kvb(ig,mm) = kk(ig)
   44   continue

   40 continue

c      call prm('gvb',gvb,ngb,ngb,3)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end


