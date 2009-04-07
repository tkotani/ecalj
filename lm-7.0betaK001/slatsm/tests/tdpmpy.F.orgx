C tests dpmpy:
      program test
      implicit double precision (a-h,p-z)
c      parameter (nra=4,nca=5,nrb=nra,ncb=nca,nrc=8,ncc=9)
      parameter (nra=6,nca=nra,nrb=nra,ncb=nca,nrc=8,ncc=9)
      dimension a(nra,nca),b(nrb,ncb),c(nrc,ncc)
      dimension aa(0:nra*nca),wk(nra),aaa(nra,nca)
      
      do 100 i=1,nra
      do 100 j=1,nca
        a(i,j) = i+nra*(j-1)
  100 continue     

      do 200 i=1,nrb
      do 200 j=1,ncb
        b(i,j) = 1.d0/(i+nrb*(j-1))
  200 continue     

      call dmprnt('a:',a,nra,nra,nca)
      call dmprnt('b:',b,nrb,nrb,ncb)

      call dpack(a,min(nra,nca))
      call dupack(a,min(nra,nca))

      call dmprnt('a after packing and unpacking',a,nra,nra,nca)

      n = 2
      m = 3
      l = 4
      print *, 'n,m,l=',n,m,l

      call dmpy(a,nra,1,b,nrb,1,c,nrc,1,n,m,l)
      call dmprnt('c: a*b',c,nrc,n,m)

      call dmpy(a,nra,1,b,1,nrb,c,nrc,1,n,m,l)
      call dmprnt('c: a*b-transpose',c,nrc,n,m)

      call dmpy(a,nra,1,b,1,nrb,c,1,nrc,n,m,l)
      call dmprnt('c: (a*b-transpose)-transpose',c,nrc,n,m)

      call dpack(a,min(nra,nca))      

      call dpmpy(a,b,nrb,1,c,nrc,1,n,m,l)
      call dmprnt('c: a(packed)*b',c,nrc,n,m)

      call dpmpy(a,b,1,nrb,c,nrc,1,n,m,l)
      call dmprnt('c: a(packed)*b-transpose',c,nrc,n,m)

      call dpmpy(a,b,1,nrb,c,1,nrc,n,m,l)
      call dmprnt('c: (a(packed)*b-transpose)-transpose',c,nrc,n,m)

      end
      subroutine dmprnt(string,a,na,n,m)
      integer na,n,m
      double precision a(na,na)
      character*(*) string
      integer i,j

800   format(6(1PG13.6))

      write(*,*) string
820   format(1x,80a)
      do 100 i = 1,n
        write(*,800) (a(i,j), j = 1,m)
100   continue
      write(*,*) ' '
      return
      end
