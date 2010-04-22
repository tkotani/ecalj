BEGIN {
  iline=0
  stderr="/dev/stderr"
}
/^#/{
   iline++
   line[iline,0]=$0
   line[iline,-1]=0
   next;
}
/^\t/ {
   line[iline,-1]= line[iline,-1]+1
   iadd= line[iline,-1]
   line[iline,iadd]  = $0
   if (is_dep) {
      #print "exec", $0
   }
   next; 
}
{  
   iline++; 
   line[iline,0]=$0 
   line[iline,-1]=0
   is_dep=0
   if (match($0,":=")==0 && match($0,":")) {
      is_dep=1
   }
   if (is_dep==0) { next; }
   gsub(":"," : ");
   n=split($0,a)
   
   for (i=1;i<=n;i++) {
      if (a[i]==":") {
          nsep=i
          break;
      }
   }

   for (i=1;i<nsep;i++) {
      target[iline,i]=a[i]
   }
   target[iline,0]=nsep-1
   for (i=nsep+1;i<=n;i++) {
      source[iline,i-nsep]=a[i]
   }
   source[iline,0]=n-nsep

#   n=target[iline,0]
#   for (i=1;i<=n;i++) {
#      print iline,i,"target",target[iline,i]
#   }
#   n=source[iline,0]
#   for (i=1;i<=n;i++) {
#      print iline,i,"source", source[iline,i]
#   }
}
END{
  nline=iline

# find duplicated target
  for (iline1=1;iline1<=nline;iline1++) {
  for (iline2=iline1+1;iline2<=nline;iline2++) {
     ntarget1=target[iline1,0]
     ntarget2=target[iline2,0]
     
     idup=0
     for (itarget1=1;itarget1<=ntarget1;itarget1++) {
     for (itarget2=1;itarget2<=ntarget2;itarget2++) {
         name1=target[iline1,itarget1]
         name2=target[iline2,itarget2]
	#print "test dup>",iline1,iline2,name1,name2
         if (name1==name2) {
             idup=1
	     print "#",name1," is the same" > stderr
         }

     }
     }

     if (idup) {
         print "#duplicated target" > stderr
	 print "#line1>" > stderr
         print "#",line[iline1,0] > stderr
         n=line[iline1,-1]
         for (i=1;i<=n;i++) {
		print "#",line[iline1,i] > stderr
	 }
	 print "#line2>" > stderr
         print "#",line[iline2,0] > stderr
         n=line[iline2,-1]
         for (i=1;i<=n;i++) {
		print "#",line[iline2,i] > stderr
	 }
	 print "# use line2=",line[iline2,0] > stderr
	 print "# and execute line1 " > stderr
         n=line[iline1,-1]
         for (i=1;i<=n;i++) {
                print "#",line[iline1,i] > stderr
         }
         line[iline1,0] =  line[iline2,0]
         line[iline2,0]=""
         line[iline2,-1]=0
	 print "#" > stderr

     }
  } 
  }


  for (i=1;i<=nline;i++) {
      nadd=line[i,-1]
      print line[i,0]
      for (j=1;j<=nadd;j++) {
         print line[i,j]
      }
  }

}
