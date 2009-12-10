--- lmfav7.F.ORG	2009-04-07 10:41:32.000000000 +0900
+++ lmfav7.F	2009-04-30 14:35:57.000000000 +0900
@@ -67,6 +67,10 @@
       character*1,allocatable:: recrd(:)
 
       real(8),parameter::    NULLR =-99999
+#if MPI|MPIK
+      include 'mpif.h'
+      integer ierr
+#endif
 
 C -------------- First executable statement ---------------
 C     Link in routines that cause potential library conflicts
@@ -92,6 +96,9 @@
       vn(2) = vn(2) + vn2(2)/1d5
       stdo = lgunit(1)
       i = dglob('stdo',dble(stdo),1)
+#if MPI|MPIK
+      call mpi_init(ierr)
+#endif
       master = 0
       procid = mpipid(1)
       nproc  = mpipid(0)
