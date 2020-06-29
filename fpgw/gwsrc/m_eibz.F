      module m_eibz
      use m_read_bzdata,only: Read_bzdata, !<--- 'call read_bzdata' sets up following data.
     &   ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,qbas,ginv,
     &   dq_,qbz,wbz,qibz,wibz,
     &     ntetf,idtetf,ib1bz, qbzw,nqbzw !for tetrahedron
c     &     idteti, nstar,irk,nstbz
      use m_genallcf_v3,only: Genallcf_v3,
     &     nclass,natom,nspin,nl,nn,ngrp,
     &     nlmto,nlnmx, nctot,!niw, !nw_input=>nw,
     &     alat, deltaw,symgrp,clabl,iclass, !diw,dw,delta,
     &     invg, il, in, im, nlnm, 
     &     plat, pos, ecore, symgg 
      use m_qbze,only: Setqbze,
     &     nqbze,nqibze,qbze,qibze
      use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      use m_rdpp,only: Rdpp,      !"call rdpp" generate following data.
     &     nxx,lx,nx,mdimx,nbloch,cgr,ppbrd 
      use m_pbindex,only: PBindex !,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
      use m_readqgcou,only: Readqgcou
      implicit none

      integer,allocatable,protected:: nwgt(:,:)
      integer,allocatable,protected:: neibz(:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
      
      contains

      subroutine Seteibz(iqxini,iqxend,iprintx)
      logical:: eibzmode,eibz4x0,tiii,iprintx
      integer:: iqxini,iqxend,iqxendx,timereversal
      eibzmode = eibz4x0()
!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
!!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
      if(eibzmode) then
!! commentout block inversion Use iqxendx=iqxend because of full inversion
         iqxendx=iqxend
         allocate( nwgt(nqbz,iqxini:iqxendx), !qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
     &        igx(ngrp*2,nqbz,iqxini:iqxendx),igxt(ngrp*2,nqbz,iqxini:iqxendx),
     &        eibzsym(ngrp,-1:1,iqxini:iqxendx))
         iprintx=.false.
         write(6,*)
         write(6,"('=== Goto eibzgen === TimeRevesal switch =',l1)")timereversal() 
         call eibzgen(nqibz,symgg,ngrp,qibze(:,iqxini:iqxend),iqxini,iqxendx,qbz,nqbz,
     i        timereversal(),ginv,iprintx,
     o        nwgt,igx,igxt,eibzsym,tiii)
         write(6,"('Used timeRevesal for EIBZ = ',l1)") tiii
         call cputid(0)
         call PBindex(natom,lx,2*(nl-1),nx) !this returns required index stored in arrays in m_pbindex.
           ! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
         call readqgcou() ! no input. Read QGcou and store date into variables.
      else 
         iqxendx=iqxend
         allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &    ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      endif
      end subroutine
!!
      end module
