BEGIN{

 for (iarg=1;iarg<ARGC;iarg++) {
     file= ARGV[iarg]
     ifiles=0
     files[iarg,-1]=file
    while (1) {
     ret=getline < file
     if (ret<=0) { break; }
     ifiles++
     files[iarg,ifiles] = $1
     }
     files[iarg,0]=ifiles
 }

# find duplicated files

 for (iarg1=1;iarg1<ARGC; iarg1++) {
 for (iarg2=iarg1+1;iarg2<ARGC; iarg2++) {
    nfiles1=files[iarg1,0]
    nfiles2=files[iarg2,0]

    for (ifile1=1; ifile1<=nfiles1; ifile1++) {

   name1=files[iarg1,ifile1]
   for (ifile2=1;ifile2<=nfiles2;ifile2++) {
       name2=files[iarg2,ifile2]
       if ( name1==name2 ) {
          print "duplicated file", files[iarg1,-1], files[iarg1,ifile1], files[iarg2,-1], files[iarg2,ifile2]
       }
   }
}

}
}

}


