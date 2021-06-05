      subroutine zpmblk(nblk,ipm,ioffo,ioffn,lds,SRC,mode,alfa,nlma,ldd,
     .DST,nmax)
C- Clips and optionally permutes subblocks of a double complex array
C ----------------------------------------------------------------------
Ci Inputs
Ci   nblk  :number of subblocks to permute
Ci   ipm   :permutation of subblocks; see Remarks
Ci   ioffo :ioffo(1:nblk+1) markers to the start of subblocks in SRC:
Ci          Subblock i starts at ioffo(i)+1 and ends at ioffo(i+1)
Ci          Thus SRC subblock has has size ioffo(i+1)-ioffo(i)
Ci   ioffn :ioffn(1:nblk+1) markers to the start of subblocks in DST:
Ci          Subblock i starts at ioffn(i)+1 and ends at ioffn(i+1)
Ci          Thus DST subblock has size ioffn(i+1)-ioffn(i)
Ci          NB: ioffn is the same as ioffo, but for DST.
Ci   lds   :leading dimension of SRC
Ci   SRC   :source matrix
Ci   mode:  one's digit: 0 SRC offsets to block i computed by ioffo(i)
Ci                     : 1 SRC offsets are relative:  offset to block i
Ci                     :       computed by ioffo(i)-ioffo(1)
Ci                     : 2 SRC is not permuted at all: columns read
Ci                     :       sequentially starting from col 1
Ci                     :       In this case, ioffo is not used.
Ci                     : 3     same as 2, but cols are read starting
Ci                     :       at col ioffo(1).
Ci                     :       In this case, only ioffo(1) is used.
Ci                     : NB: only 1 or 10s digit may be set to 2 or 3
Ci          ten's digit: 0 DST offsets to block i computed by ioffn(i)
Ci                     : 1 DST offsets are relative:  offset to block i
Ci                     :       computed by ioffn(i)-ioffn(1)
Ci                     : 2 DST is not permuted at all: columns written
Ci                     :       sequentially starting from col ioffn(1).
Ci                     :       In this case, ioffn is not used.
Ci                     : 3     same as 2, but cols are written starting
Ci                     :       at col ioffn(1)
Ci                     :       In this case, only ioffn(1) is used.
Ci                     : NB: only 1 or 10s digit may be set to 2 or 3
Ci          100's digit: 0 no permutation of subblock order
Ci                         In this case, ipm is not used
Ci                       1 permute sequence of subblocks in SRC by ipm
Ci                         i.e. SRC subblock i endpoints are computed
Ci                         from ioffo(j) and ioffo(j+1), j=ipm(i)
Ci                         NB: not compatible with 1s digit >=2
Ci                       2 permute sequence of subblocks in DST by ipm
Ci                         i.e. DST subblock i endpoints are computed
Ci                         from ioffo(j) and ioffn(j+1), j=ipm(i)
Ci                         NB: not compatible with 10s digit >=2
Ci         1000's digit: 1 add into new array, not copy
Ci                       2 copy (add) alfa(j)*SRC(*,j) into DST(*,j')
Ci                       3 both 1 and 2
Ci                       4 add 4 to exchange roles of rows and columns,
Ci                         i.e. permute rows for each nlma columns
Ci                       8 permute subblocks in both rows and cols
Ci                        (see Remarks).
Ci                       9 like 8, but add into DST.
Ci   alfa  :copy or add alfa(j)*SRC(*,j) into DST(*,j'), if mode set
Ci   nlma  :number of rows for which to permute columns
Ci          (not used if rows are also permuted; see 1000s digit mode)
Ci   ldd   :leading dimension of DST
Co Outputs
Co   DST   :destination matrix
Co   nmax  :largest column index in DST affected by call to pmblk.
Cr Remarks
Cr  Columns are grouped into a sequence of subblocks.  pmblk copies
Cr  entire subblocks of SRC: columns within a subblock retain their
Cr  order when copied.  The subblock order, however, may be permuted,
Cr  as specified by the 10s digit of mode.
Cr
Cr  * The permutations may occur by permuting the order in which the
Cr    subblocks are copied from SRC (if 1000s digit mode is 1) or by
Cr    permuting the order in which subblocks are copied to DST (if 1000s
Cr    digit mode is 2).  By default subblocks are copied in their
Cr    natural order; thus there is no permutation of subblocks if 1000s
Cr    digit mode is 0.
Cr
Cr  * If both SRC and DST are being permuted, the size of the ith
Cr    subblock copied is the smaller of sizes of the subblock i in SRC
Cr    and subblock i in DST.  Thus, if the size of DST subblock i is
Cr    smaller than SRC, only the first portion of the SRC subblock is
Cr    copied; if the size DST is larger than SRC, only the first columns
Cr    in the DST subblock will be affected.
Cr
Cr  * By default pmblk permutes columns, doing the permutation for each of
Cr    nlma rows.  If mode 1000's digit >=8, pmblk permutes both columns
Cr    and rows.
Cu Updates
Cu   16 May 03 Added options 2 and 3 to 1 and 10s digit mode
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nblk,lds,ldd,nlma,mode,ioffn(nblk),ioffo(nblk),
     .ipm(nblk),nmax
      double complex dst(ldd,1),src(lds,1),alfa(1)
C ... Local parameters
      integer ilma,ib,iofo,ilm,iofn,nlmi,ipn,ipo,offo0,offn0,
     .nn,no,mode1,mode2,mode3,mode4,jb,jofo,jlm,jofn,nlmj,jpn,jpo
      logical:: isanrg, l_dummy_isanrg

      mode1 = mod(mode,10)
      mode2 = mod(mode/10,10)
      mode3 = mod(mode/100,10)
      mode4 = mod(mode/1000,10)

      offn0 = 0
      offo0 = 0
      if (mod(mode1,2) .eq. 1) offo0 = ioffo(1)
      if (mod(mode2,2) .eq. 1) offn0 = ioffn(1)
      if (mode1.ge.2 .and. mode2.ge.2)
     .call rx('zpmblk: must permute one of SRC or DST')
      if (mode1.ge.2 .and. mode3.eq.1)
     .call rxi('zpmblk: imcompatible switches, mode=',mode)
      if (mode2.ge.2 .and. mode3.eq.2)
     .call rxi('zpmblk: imcompatible switches, mode=',mode)
      l_dummy_isanrg=isanrg(mode3,0,2,'zpmblk:','100s digit mode', .true.)

      nmax = -1
      nlmi = 0
      iofo = offo0
      iofn = offn0
      if (mode4 .le. 7) then
        do  10  ib = 1, nblk
          ipo = ib
          if (mode3 .eq. 1) ipo = ipm(ib)
          ipn = ib
          if (mode3 .eq. 2) ipn = ipm(ib)
C         sizes and offsets to SRC and DST blocks
          if (mode1 .ge. 2) then
            iofo = iofo + nlmi
            iofn = ioffn(ipn) - offn0
            nlmi = ioffn(ipn+1)-ioffn(ipn)
          elseif (mode2 .ge. 2) then
            iofn = iofn + nlmi
            iofo = ioffo(ipo) - offo0
            nlmi = ioffo(ipo+1)-ioffo(ipo)
          else
            iofn = ioffn(ipn) - offn0
            iofo = ioffo(ipo) - offo0
            nn = ioffn(ipn+1)-ioffn(ipn)
            no = ioffo(ipo+1)-ioffo(ipo)
            nlmi = min(nn,no)
          endif

          if (mode4 .eq. 0) then
            do  ilm = 1, nlmi
              do  ilma = 1, nlma
                 dst(ilma,ilm+iofn) = src(ilma,ilm+iofo)
              enddo
           enddo

          elseif (mode4 .eq. 1) then
            do   ilm = 1, nlmi
              do  ilma = 1, nlma
                 dst(ilma,ilm+iofn) = dst(ilma,ilm+iofn)
     .                + src(ilma,ilm+iofo)
              enddo
             enddo
             
          elseif (mode4 .eq. 2) then
            do    ilm = 1, nlmi
              do  ilma = 1, nlma
                 dst(ilma,ilm+iofn) = src(ilma,ilm+iofo)*alfa(ilm+iofn)
              enddo
           enddo
                 
          elseif (mode4 .eq. 3) then
            do    ilm = 1, nlmi
              do  ilma = 1, nlma
                 dst(ilma,ilm+iofn) = dst(ilma,ilm+iofn)
     .                + src(ilma,ilm+iofo)*alfa(ilm+iofn)
              enddo
             enddo
          elseif (mode4 .eq. 4) then
            do   ilma = 1, nlma
              do  ilm = 1, nlmi
                 dst(ilm+iofn,ilma) = src(ilm+iofo,ilma)
              enddo
             enddo
          elseif (mode4 .eq. 5) then
            do    ilma = 1, nlma
              do   ilm = 1, nlmi
                 dst(ilm+iofn,ilma) = dst(ilm+iofn,ilma)
     .      + src(ilm+iofo,ilma)
              enddo
             enddo
          elseif (mode4 .eq. 6) then
            do   ilma = 1, nlma
            do  ilm = 1, nlmi
               dst(ilm+iofn,ilma) = src(ilm+iofo,ilma)*alfa(ilm+iofn)
            enddo
            enddo

          else
            do   ilma = 1, nlma
              do  ilm = 1, nlmi
                 dst(ilm+iofn,ilma) = dst(ilm+iofn,ilma)
     .                + src(ilm+iofo,ilma)*alfa(ilm+iofn)
              enddo
           enddo
          endif
          nmax = max(nn+iofn,nmax)
   10   continue
      else
        do  100  ib = 1, nblk
          ipo = ib
          if (mode3 .eq. 1) ipo = ipm(ib)
          ipn = ib
          if (mode3 .eq. 2) ipn = ipm(ib)
C         sizes and offsets to SRC and DST blocks, col index
          if (mode1 .ge. 2) then
            iofo = iofo + nlmi
            iofn = ioffn(ipn) - offn0
            nlmi = ioffn(ipn+1)-ioffn(ipn)
          elseif (mode2 .ge. 2) then
            iofn = iofn + nlmi
            iofo = ioffo(ipo) - offo0
            nlmi = ioffo(ipo+1)-ioffo(ipo)
          else
            iofn = ioffn(ipn) - offn0
            iofo = ioffo(ipo) - offo0
            nn = ioffn(ipn+1)-ioffn(ipn)
            no = ioffo(ipo+1)-ioffo(ipo)
            nlmi = min(nn,no)
          endif
          nlmj = 0
          jofo = offo0
          jofn = offn0
          do  110  jb = 1, nblk
            jpo = jb
            if (mode3 .eq. 1) jpo = ipm(jb)
            jpn = jb
            if (mode3 .eq. 2) jpn = ipm(jb)
C           sizes and offsets to SRC and DST blocks, row index
            if (mode1 .ge. 2) then
              jofo = jofo + nlmj
              jofn = ioffn(jpn) - offn0
              nlmj = ioffn(jpn+1)-ioffn(jpn)
            elseif (mode2 .ge. 2) then
              jofn = jofn + nlmj
              jofo = ioffo(jpo) - offo0
              nlmj = ioffo(jpo+1)-ioffo(jpo)
            else
              jofn = ioffn(jpn) - offn0
              jofo = ioffo(jpo) - offo0
              nn = ioffn(jpn+1)-ioffn(jpn)
              no = ioffo(jpo+1)-ioffo(jpo)
              nlmj = min(nn,no)
            endif
            if (mode4 .eq. 8) then
              do   jlm = 1, nlmj
                 do   ilm = 1, nlmi
                    dst(ilm+iofn,jlm+jofn) = src(ilm+iofo,jlm+jofo)
                 enddo
              enddo
            else
               do    jlm = 1, nlmj
                  do  ilm = 1, nlmi
                     dst(ilm+iofn,jlm+jofn) = dst(ilm+iofn,jlm+jofn) +
     .                    src(ilm+iofo,jlm+jofo)
                  enddo
               enddo
            endif
  110     continue
          nmax = max(nn+iofn,nmax)
  100   continue
      endif
      end
#if TEST
      subroutine fmain
      implicit none
      integer nblk,ioffs(5),ioffd(5),iprm(4)
      double precision s(10),d(10)
      integer nmax,i,mode

      ioffs(1) = 0
      ioffs(2) = 1
      ioffs(3) = 3
      ioffs(4) = 6
      ioffs(5) = 10

      ioffd(1) = 0
      ioffd(2) = 2
      ioffd(3) = 6
      ioffd(4) = 9
      ioffd(5) = 10

      nblk = 4

      iprm(1) = 2
      iprm(2) = 4
      iprm(3) = 3
      iprm(4) = 1

      do  i = 1, 10
        s(i) = 1.1d0*i
      enddo

      print 332, ioffs
  332 format(10i3)
      print 332, ioffd
      print 332, iprm

C     Do the permutation
C     NB: 4000 should have no effect; just tests different branch.
      call dpzero(d,10)
      mode = 4100
      call pmblk(nblk,iprm,ioffs,ioffd,1,s,mode,0d0,1,1,d,nmax)

      print 333, (s(i), i=1,10)
      print 333, (d(i), i=1,10)
  333 format(10f8.1)

C     Undo the permutation
      mode = 0200
      call dpzero(s,10)
      call pmblk(nblk,iprm,ioffd,ioffs,1,d,mode,0d0,1,1,s,nmax)
      print 333, (s(i), i=1,10)

C     Same permutation, but use mode for sequential DST order
      call dpzero(d,10)
      mode = 4120
      call pmblk(nblk,iprm,ioffs,ioffd,1,s,mode,0d0,1,1,d,nmax)
      print 333, (d(i), i=1,10)

      end
#endif

