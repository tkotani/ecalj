C tests complex matrix inversion
      program test
      implicit double precision (a-h,p-z)
      parameter (nra=3,nca=nra,nrb=nra,ncb=nca,nrc=5,ncc=9)
      complex*16 a,b,c,wk
      dimension a(nra,nca),b(nrb,ncb),c(nrc,ncc)
      dimension aa(0:nra*nca),wk(nra),aaa(nra,nca)
      integer kpvt(nra)
      
      do 100 i=1,nra
      do 100 j=1,i
        a(i,j) = i+nra*(j-1) + i*j*dcmplx(0d0, 1d0)
        a(j,i) = dconjg(a(i,j))
        a(i,i) = dble(a(i,i))
  100 continue     

      do 200 i=1,nrb
      do 200 j=1,ncb
        b(i,j) = 1.d0/(i+nrb*(j-1)) + (i+j)*dcmplx(0d0, 1d0)
        b(j,i) = dconjg(b(i,j))
        b(i,i) = dble(b(i,i))

        b(i,j) = a(i,j)
  200 continue     

      call zmprnt('a:',a,nra,nra,nca)
      call zmprnt('b:',b,nrb,nrb,ncb)

      call zsifa(a,nra,nra,kpvt,info)      
      print *, 'info=', info
      call zmprnt('a:',a,nra,nra,nca)
      call zsidi(a,nra,nra,kpvt,det,wk,1)
      print *, 'det=', det

      do 100 i=1,nra
      do 100 j=
        a(j,i) = dconjg(a(i,j))
  100 continue     

      call zmprnt('a:',a,nra,nra,nca)

      call zmpy(a,2*nra,2,1,b,2*nrb,2,1,c,2*nrc,2,1,nra,nra,nra)
      call zmprnt('c:',c,nrc,nrb,ncb)

      end
      subroutine zmprnt(string,a,na,n,m)
      integer na,n,m
      double precision a(0:1,na,na)
      character*(*) string
      integer i,j

800   format(6(1PG13.6))

      write(*,*) string
820   format(1x,80a)
      do 100 i = 1,n
        write(*,800) (a(0,i,j), j = 1,m)
        write(*,800) (a(1,i,j), j = 1,m)
100   continue
      write(*,*) ' '
      return
      end
