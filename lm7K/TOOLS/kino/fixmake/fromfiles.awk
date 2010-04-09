BEGIN {
  debug=0
  moddir="$(moddir)/"
  modext=".mod"
  IGNORECASE=1
  narg=ARGC-1
  for (i=1;i<=narg;i++) {
     g_files[i]=ARGV[i]
  }

  dep_fromfile(narg)
  fix_use(narg)
  fix_moduse(narg)
  make_dep(narg,"__makefile__")

}

func make_dep(narg,fileout0, fileout,iarg,source,nmodule,nuse,targetname,i,sourcename) {
  fileout=fileout0
  if  (fileout=="") {
     fileout="/dev/stdout"
     print "#output to ",fileout
  }
  for (iarg=1;iarg<=narg;iarg++) {
     source=g_sourcefile[iarg]
     nmodule=g_module[iarg,0]
     nuse=g_use[iarg,0]
     if (nmodule>0 || nuse>0 ) { 
        targetname=fix_dir(source,"t",".o")
        printf("%s ",targetname) >fileout
        for (i=1;i<=nmodule;i++) {
        printf("%s ", moddir g_module[iarg,i] modext) >fileout
        }
        printf(" : ") >fileout
        sourcename=fix_dir(source,"s","")
        printf("%s ",sourcename) >fileout
        for (i=1;i<=nuse;i++) {
          if (length(g_use[iarg,i])>1) {
             printf("%s ",moddir g_use[iarg,i] modext) >fileout
          }
        }
        printf("\n") >fileout
	#command
	#$(FC) $(FFLAGS) -c  -o $(subs_obj_path)/susite.o subs/susite.F
	printf("\t$(FC) $(FFLAGS) -c  -o %s %s\n",targetname,sourcename) >fileout
     }
  }
}

func fix_dir(name,ts,ext) {
   if (ts=="t") {
      gsub("^.*subs/","$(subs_obj_path)/",name)
      gsub("^.*fp/","$(fp_obj_path)/",name)
      gsub("^.*gwd/","$(gwd_obj_path)/",name)
      gsub("^.*slatsm/","$(sla_obj_path)/",name)
   }
   else {
      gsub("^.*subs/","subs/",name)
      gsub("^.*fp/","fp/",name)
      gsub("^.*gwd/","gwd/",name)
      gsub("^.*slatsm/","slatsm/",name)
   }
   if (ext!="") {
      gsub("\\.F$",ext,name)
   }
   return name
}

func fix_moduse(narg, nuse,nmodule,imodule,namemodule,iuse ) {
  for (iarg=1;iarg<=narg;iarg++) {
    nuse=g_use[iarg,0]
    nmodule=g_module[iarg,0]
    if (nuse>0 && nmodule>0) {
       for (imodule=1;imodule<=nmodule;imodule++) {
          namemodule=g_module[iarg,imodule]
          for (iuse=1;iuse<=nuse;iuse++) {
             if (namemodule==g_use[iarg,iuse]) {
                print "#cyclic dependency:", namemodule, " in ",g_sourcefile[iarg]
                g_use[iarg,iuse]=""
             } 
          }
       }
        
    }
  }
}

func fix_use(narg, iarg,nuse,i,file,file2,cmd,iuse,ret ) {
  for (iarg=1;iarg<=narg;iarg++) {
    nuse=g_use[iarg,0]
    if (nuse>0) {
       file="__FIXMAKE__" 
       for (i=1;i<=nuse;i++) {
         print g_use[iarg,i] >file
       }
       close(file)
       file2="__FIXMAKE2__"
       cmd="sort <"file" | uniq > "file2
       system(cmd)
       iuse=0
       while (1) {
         ret= getline < file2
         if (ret<=0) { break; }
         iuse++
         g_use[iarg,iuse]=$1
       }
       close(file2)
       g_use[iarg,0]=iuse
       if (debug) {
          print "use:",g_use[iarg,0]
          for (i=1;i<=g_use[iarg,0];i++) {
              print g_use[iarg,i]
          }
       }
       cmd="rm "file" "file2
	system(cmd)
    }
    
  }
}

func dep_fromfile(narg, iarg,file,ret,n,a,imodule,iuse) {

  for (iarg=1;iarg<=narg;iarg++) {
     file=g_files[iarg]
     g_sourcefile[iarg]=file
     if (debug) { print "read ",file }
     imodule=0
     iuse=0
     while (1) {
        ret=getline < file
        if (ret<=0) { break; }
        if (match($0,"^ +module")){
            gsub(","," ");
	    n=split($0,a)
            if (a[2]!="procedure") {
                imodule++
                g_module[iarg,imodule]=a[2]
                g_module[iarg,0]=imodule
		if (debug) { print "module ",a[2] }
            }
        }
        else if (match($0,"^ +use ")){
	    gsub(","," , ")
            n=split($0,a)
            if (a[2]!="=") {
            	iuse++
	    	g_use[iarg,iuse] = a[2]
            	g_use[iarg,0]=iuse
		if (debug) { print "use ",a[2] }
            }
        }
     }
     close(file)
  }
}

