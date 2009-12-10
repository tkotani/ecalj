      subroutine shosg(ivl,nds,nsite,nbas,plat,bas,
     .                  iax,ntab,ntabg,sg,io,scal)
C- Print out value-Laplacian structure constants 
C ----------------------------------------------------------------------
Ci Inputs
Ci   ivl   :identifies functions used to build the value-Laplacian set
Ci   nds   :leading dimension of sg
Ci   nsite :number of pairs in neighbor table
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   ntabg :ntabg(ib)= # of neighbors for cluster ib (pairg.f)
Ci   sg    :real-space structure constants, in order of iax table
Ci   io:    0, print out all structure constants
Ci          1, print out only sc's whose ss element differs from last
Ci   scal  : Multiply structure constants by scal in printout
Co Outputs
Co   strx sg for the auxiliary ("Gaussian") basis are printed to stdo
Cv Verbosity
Cv  >=10: prints out L=0 sg for each site
Cv  >=30: prints out and diagonal strux sg_LL for each site
Cv  >=40: prints out and full strux sg_LL for each site
Cv  
Cv  >30: prints out iax for each site
Cv   40: prints out structure constant matrix for any site whose
Cv       ss element differs from previous
Cv  >40: prints out full structure constant matrix
Cb Bugs
Cr Remarks
Cu Updates
Cu   18 Jan 08 (S. Lozovoi) Adapted from shostr.f
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ivl,nds,nsite,nbas,io,ntab(nbas+1),ntabg(nbas)
      integer niax,nkaps
      parameter (niax=10, nkaps=2)
      integer iax(niax,nsite)
      double precision sg(nds,nds,nkaps,nkaps,*),scal,
     .                 plat(3,3),bas(3,*)
C Local parameters
      double precision dx
      integer i,j,k,il1,il2,ib,jb,iat2,ix,ipr,it,iclus,ig,
     .  lgunit,iic,jjc,nlma,nlmb,ip,i0,stdo,isw,
     .  ikap,jkap
      character outs*120,out2*120
      character*7 cgh(0:2)
      data cgh/'G1     ','Hsm    ','Hsm-dot'/

      dx(ix,iat2,i,j,k) = bas(ix,iat2) - bas(ix,ib) +
     .                    plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k

      call getpr(ipr)
      stdo = lgunit(1)

      call info2(10,1,0,
     .  ' SHOSG:  '//
     .  'Strux for value-Laplacian basis'//
     .  ' prepared from G0 and '//trim(cgh(ivl))//
     .  '%?#n#, scaled by %d##',
     .  isw(scal.ne.1d0),scal)
      call info2(10,0,0,
     .  '%9f%i pairs out of %i total',
     .  nsite,ntab(nbas+1)) 


      ig = 0
      do  it = 1, nbas
        do  iclus = 1, ntabg(it)
          ig = ig+1
          i = ntab(it) + iclus
          outs = ' '
          ip = 3

          if (iclus .eq. 1) then
            i0 = i-1
            out2 = ' '
            ip = 3
            ib  = iax(1,i)
            nlma = iax(9,i)
            jb = 0
            call info5(0,1,0,' Site  %i  has  %i  neighbors',
     .        ib, ntabg(ib),ipr,0,0)
          endif

        jb = jb+1
C       ii = i
        if (ipr .ge. 30) then
          write(stdo,300) i, i-i0,iax(2,i), iax(1,i), (iax(j,i),j=3,5),
     .      (dx(ix,iax(2,i),iax(3,i),iax(4,i),iax(5,i)), ix=1,3)
  300     format(' pair',i5,'(',i2,')',' R''=',i3,
     .           ' R=',i3,' plat=',3i3,' tau=',3f9.4)
          call awrit0('%a',outs,80,-stdo)
        endif

C       jj = max(i-1,1)
        iic = ig
        jjc = max(ig-1,1)
        if (io .eq. 0 .or. i .eq. 1 .or. ipr .gt. 40 .or.
     .      dabs(sg(1,1,1,1,iic) - sg(1,1,1,1,jjc)) .gt. 1.d-8) then
c 110     continue
          nlmb = iax(9,iic)
          if (ipr .ge. 40) then
            do  ikap = 1, nkaps
              do  jkap = 1, nkaps
                do  il1 = 1, nlmb
                  if (ipr .gt. 50) then
                    print 302, (scal*sg(il1,il2,ikap,jkap,iic),
     .                il2=1,nlma)
                  else
                    print 301, (scal*sg(il1,il2,ikap,jkap,iic),
     .                il2=1,nlma)
                  endif
                enddo
                if (nkaps .gt. 1) write(stdo,'('' ------'')')
              enddo
            enddo
          elseif (ipr .ge. 30) then
            do  ikap = 1, nkaps
              do  jkap = 1, nkaps
                outs = 'SG_LL'
                ip = 5
                do  il2 = 1, min(nlmb,nlma)
                  call bin2a('F6',1,0,scal*
     .              sg(il2,il2,ikap,jkap,iic),4,0,len(outs),outs,ip)
                enddo
                print '(a)', outs(1:ip)
              enddo
            enddo
          elseif (ipr .ge. 10) then
            do  ikap = 1, nkaps
              do  jkap = 1, nkaps
                call bin2a('F8',1,0,scal*sg(1,1,ikap,jkap,iic),4,0,
     .            len(outs),out2,ip)
                if (ip .ge. 72) then
                  print '(a)', out2(1:ip)
                  out2 = ' '
                  ip = 3
                endif
              enddo
            enddo
          endif
        endif

        enddo
      enddo
      if (i .ne. 1 .and. ipr .le. 30 .and. ipr .ge. 10
     .    .and. out2 .ne. ' ') print '(a)', out2(1:ip)

  301 format(1x,9f8.4)
  302 format(1x,9f12.8)

      end

