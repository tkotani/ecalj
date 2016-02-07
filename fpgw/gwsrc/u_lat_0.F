c --- u_lat_0
      subroutine u_lat_0(s_lat, alat,vol,awald0,awald,tol,nkdmx,nkqmx,
     .  nkd,nkq,odlv,oqlv)
C- Extracts from s_lat scalar parameters relating to lattice
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(60)
      alat   = s_lat(1)
      vol   = s_lat(9)
      awald  = s_lat(10)
      awald0 = s_lat(14)
      tol    = s_lat(15)
      nkdmx  = nint(s_lat(16))
      nkqmx  = nint(s_lat(17))
      odlv   = nint(s_lat(18))
      oqlv   = nint(s_lat(19))
      nkd    = nint(s_lat(20))
      nkq    = nint(s_lat(21))
      end

c --- p_lat_0
      subroutine p_lat_0(s_lat, alat,vol,awald0,awald,tol,nkdmx,nkqmx,
     .  nkd,nkq,odlv,oqlv)
C- puts into s_lat scalar parameters relating to lattice
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(60)
      s_lat(1)  = alat
      s_lat(9)  = vol
      s_lat(10) = awald
      s_lat(14) = awald0
      s_lat(15) = tol
      s_lat(16) = nkdmx
      s_lat(17) = nkqmx
      s_lat(18) = odlv
      s_lat(19) = oqlv
      s_lat(20) = nkd
      s_lat(21) = nkq
      end

c --- u_lat_vecs
      subroutine u_lat_vecs(s_lat, alat,plat,qlat)
      implicit real*8 (a-h,p-z), integer (o)
      dimension plat(9),qlat(9),s_lat(50)
      alat=s_lat(1)
      do 10 i=1,9
        plat(i)=s_lat(i+30)
        qlat(i)=s_lat(i+40)
  10  continue
      end

      subroutine p_lat_vecs(s_lat, alat,plat,qlat)
      implicit real*8 (a-h,p-z), integer (o)
      dimension plat(9),qlat(9),s_lat(50)
      s_lat(1) = alat
      do 10 i=1,9
        s_lat(i+30) = plat(i)
        s_lat(i+40) = qlat(i)
  10  continue
      end

c --- u_lat_v0: get plat0 (lattice basis before distortion)
      subroutine u_lat_v0(s_lat, plat0)
      implicit real*8 (a-h,p-z), integer (o)
      dimension plat0(9),s_lat(50)
      do 10 i=1,9
        plat0(i)=s_lat(i+50)
  10  continue
      end

c --- u_lat_dist
      subroutine u_lat_dist(s_lat, gx,gy,gz,gam)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(50)
      gx     = s_lat(22)
      gy     = s_lat(23)
      gz     = s_lat(24)
      gam    = s_lat(25)
      end

c --- u_lat_FT
      subroutine u_lat_FT(s_lat,n1,n2,n3)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(50)

      n1 = nint(s_lat(11))
      n2 = nint(s_lat(12))
      n3 = nint(s_lat(13))

      end

c --- u_lat_qmesh
      subroutine u_lat_qmesh(s_lat,job,nkx,nky,nkz)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      job = nint(s_lat(71))
      nkx = nint(s_lat(72))
      nky = nint(s_lat(73))
      nkz = nint(s_lat(74))

      end

c --- u_lat_grp, p_lat_grp
      subroutine u_lat_grp(s_lat,ngrp,og,oag)
C- unpacks from s_lat the symmetry group operations
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      ngrp = nint(s_lat(66))
      og   = nint(s_lat(67))
      oag  = nint(s_lat(68))
      end

      subroutine p_lat_grp(s_lat,ngrp,og,oag)
C- packs into s_lat parameters the symmetry group operations
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      s_lat(66) = ngrp
      s_lat(67) = og
      s_lat(68) = oag
      end

c --- u_lat_gv, p_lat_gv
      subroutine u_lat_gv(s_lat,gmax,ng,ogv,okv)
C- unpacks from s_lat parameters for list of recips in cutoff
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      gmax =      s_lat(60)
      ng   = nint(s_lat(61))
      ogv  = nint(s_lat(62))
      okv  = nint(s_lat(63))
      end

      subroutine p_lat_gv(s_lat,gmax,ng,ogv,okv)
C- packs into s_lat parameters for list of recips in cutoff
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      s_lat(60) = gmax
      s_lat(61) = ng
      s_lat(62) = ogv
      s_lat(63) = okv
      end

      subroutine u_lat_rgv(s_lat,gmax,ng,ogv,okv)
C- unpacks from s_lat parameters for list of recips in cutoff
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      gmax =      s_lat(96)
      ng   = nint(s_lat(97))
      ogv  = nint(s_lat(98))
      okv  = nint(s_lat(99))
      end

      subroutine p_lat_rgv(s_lat,gmax,ng,ogv,okv)
C- packs into s_lat parameters for list of recips in big cutoff
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      s_lat(96) = gmax
      s_lat(97) = ng
      s_lat(98) = ogv
      s_lat(99) = okv
      end

      subroutine u_lat_FTb(s_lat,k1,k2,k3,n1,n2,n3)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(50)

      n1 = nint(s_lat(11)) + 4
      n2 = nint(s_lat(12)) + 4
      n3 = nint(s_lat(13)) + 4
      k1 = n1
      k2 = n2
      k3 = n3

      end

c --- u_lat_sym
      subroutine u_lat_sym(s_lat,oips0,obgv)
C- unpacks from s_lat parameters for symmetrization of mesh density
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      oips0 = nint(s_lat(64))
      obgv  = nint(s_lat(65))
      end

      subroutine p_lat_sym(s_lat,oips0,obgv)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      s_lat(64) = oips0
      s_lat(65) = obgv
      end

c --- u_lat_asym
      subroutine u_lat_asym(s_lat,nclas,oipcl,ndcl,osym)
C- unpacks from s_lat parameters for symmetrization of atomic density
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      nclas = nint(s_lat(85))
      oipcl = nint(s_lat(86))
      ndcl  = nint(s_lat(87))
      osym  = nint(s_lat(88))
      end

      subroutine p_lat_asym(s_lat,nclas,oipcl,ndcl,osym)
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      s_lat(85) = nclas
      s_lat(86) = oipcl
      s_lat(87) = ndcl
      s_lat(88) = osym
      end

c --- u_lat_cg  p_lat_cg
      subroutine u_lat_cg(s_lat,ocg,ojcg,oindxcg,ocy)
C- unpacks from s_lat pointers to clebsch-gordan coeffs and cy
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      ocg       = nint(s_lat(81))
      ojcg      = nint(s_lat(82))
      oindxcg   = nint(s_lat(83))
      ocy       = nint(s_lat(84))
      end

      subroutine p_lat_cg(s_lat,ocg,ojcg,oindxcg,ocy)
C- packs into s_lat pointers to clebsch-gordan coeffs and cy
      implicit real*8 (a-h,p-z), integer (o)
      dimension s_lat(200)
      s_lat(81) = ocg
      s_lat(82) = ojcg
      s_lat(83) = oindxcg
      s_lat(84) = ocy
      end

c ----- v_lat: unpack float element by name ----
      doubleprecision function v_lat(s_lat,name)
      implicit real*8 (a-h,p-z), integer (o)
      parameter ( n_lat=200 )
      dimension s_lat(n_lat)
      character*(*) name
      if (name .eq. 'gmax') then
        v_lat=s_lat(60)
      else if (name .eq. 'alat') then
        v_lat=s_lat(1)
      else if (name .eq. 'vol') then
        v_lat=s_lat(9)
      else if (name .eq. 'tolgv') then
        v_lat=s_lat(69)
      else
        call rxs('v_lat: unknown name: ',name)
      endif
      end

c ----- i_lat: unpack integer element by name ----
      integer function i_lat(s_lat,name)
      implicit real*8 (a-h,p-z), integer (o)
      parameter ( n_lat=200 )
      dimension s_lat(n_lat)
      character*(*) name
      if (name .eq. 'nkq') then
        i_lat=nint(s_lat(20))
      else if (name .eq. 'nkd') then
        i_lat=nint(s_lat(21))
      else
        call rxs('v_lat: unknown name: ',name)
      endif
      end

c ----- p_lat: pack float element by name ----
      subroutine p_lat (s_lat,name,value)
      implicit real*8 (a-h,p-z), integer (o)
      parameter ( n_lat=200 )
      dimension s_lat(n_lat)
      character*(*) name
      if (name .eq. 'gmax') then
        s_lat(60)=value
      else
        call rxs('p_lat_n: unknown name: ',name)
      endif
      end
