import sys
from types import *
import re
import datetime

spc="      "
spccont="     ,"

thisprogram="Cstruc_makeinit.py"

def is_comment(aline):
      thisfunc="is_comment"
      if len(aline)>0 and aline[0].isdigit():
          ret=0
          return ret
      if len(aline)>0 and aline[0]!=' ':
          ret= 1
          return ret
      return 0

############################################

def del_null(tok):
	tok2=[]
	for nn in tok:
		if len(nn)>0:
			tok2.append(nn)
	return tok2

############################################

class  DoIt:
    content=[]
    struc_list=[]

    def __init__(self):
	self.content=[]
	self.struc_list=[]

    def read(self):
	for aline in sys.stdin:
		aline=aline[:-1]
		aline=aline.replace('\t','        ')
		self.content.append(aline)

    def strucdef(self):
	n_line=len(self.content)
	i_line=0
	while i_line < n_line:
		aline=self.content[i_line]
		if is_comment(aline)==1:
			idummy=1
		else:
			tok=aline.split(" ")
			tok=del_null(tok)
			if len(tok)==2 and tok[0]=="type":
				# start type 
				#print "start>",tok
				struc_name=tok[1]
				i_line==i_line+1
				struc_def=[]
				while 1:
					aline=self.content[i_line]
					if is_comment(aline)==1:
						idummy=1
					else:
					    pline=re.sub("::"," :: ",aline)
					    tok=pline.split(" ")
					    tok=del_null(tok)
					    if len(tok)>=2 and tok[0]=="end" and tok[1]=="type":
						#print "end>", struc_def
						self.struc_list.append([struc_name,struc_def])
						break
					    elif "::" in tok:
						# var definition
						pline=re.sub(" ","",aline)
						varline=pline.split("::")
						vartype=varline[0].split(",")
						var=varline[1].split("!")
						struc_def.append([vartype,var])
					i_line=i_line+1
		i_line=i_line+1

    def show(self):
	del_routine=[]
        for struc_name_struc_def in self.struc_list:
                if len(struc_name_struc_def)==0:
                        continue
                struc_name,struc_def=struc_name_struc_def
                struc_prefix=struc_name
                struc_prefix=re.sub("s_","u",struc_prefix)
                routine_name=struc_prefix+"_init"
		del_routine.append(routine_name)
		routine_name=struc_prefix+"_show"
                del_routine.append(routine_name)	

	iline=0
	while iline<len(self.content):
		aline=self.content[iline]
		if is_comment(aline):
			print aline
			iline=iline+1
			continue
		else:
			bline=aline
			bline=re.sub("\("," ( ",bline)
			tok=bline.split(' ')	
			tok=del_null(tok)
			if len(tok)==3 and tok[0]=="end" and tok[1]=="module" and tok[2]=="m_struc_func":
				self.make_bcast()
			elif len(tok)>=2 and tok[0]=="subroutine":
				if tok[1] in del_routine:
					# search "end"
					while iline<len(self.content):
						aline=self.content[iline]
						#print iline,aline
						if is_comment(aline):
							iline=iline+1
							continue
						bline=aline
						bline=re.sub("\("," ( ",bline)
						tok=bline.split(' ')
						tok=del_null(tok)
						flag1=len(tok)==1 and tok[0]=="end"
						flag2=len(tok)==3 and tok[0]=="end" and tok[1]=="subroutine"
						if flag1 or flag2:
							iline=iline+1	
							break
						iline=iline+1	
					continue
			print aline
			iline=iline+1


    def make_init(self):
	for struc_name_struc_def in self.struc_list:
		if len(struc_name_struc_def)==0:
                        continue
		struc_name,struc_def=struc_name_struc_def
                struc_prefix=struc_name
		struc_prefix=re.sub("s_","u",struc_prefix)
                routine_name=struc_prefix+"_init"
		print spc,"subroutine",routine_name,"(struc)"
		print spc,"use m_struc_def"
      		print spc,"use m_struc_func"
      		print spc,"implicit none"
		print spc,"type("+struc_name+"):: struc"
		for vartype_var in struc_def:
			vartype,var=vartype_var
			if len(vartype)==1:
				# not a pointer
			    if var[0]=="size":
				print spc,"struc%size=",struc_prefix+"_size()"	
			    else:
				var_simple=re.sub("\(.*\)","",var[0])
				# remove array def
				if vartype[0].find("integer")!=-1:
					print spc,"struc%"+ var_simple+ "=0"
				elif vartype[0].find("real")!=-1:
					print spc,"struc%"+ var_simple+ "=0.0d0"
		print spc,"end subroutine",routine_name
		print

    def make_bcast(self):
	for struc_name_struc_def in self.struc_list:
		if len(struc_name_struc_def)==0:
			continue
		struc_name,struc_def=struc_name_struc_def
		struc_prefix=struc_name
                routine_name="mpibc1_"+ struc_prefix
		print spc,"subroutine "+routine_name+"(struc,mlog,funnam,label)"
		print spc,"use m_struc_def, only:", struc_name
		print spc,"implicit none"
		print "#if MPI|MPIK"
		print "      include 'mpif.h'"
		print "      integer numprocs, ierr"
		print "      integer MAX_PROCS"
		print "      parameter (MAX_PROCS = 100)"
		print "      integer resultlen"
		print "      character*(MPI_MAX_PROCESSOR_NAME) name"
		print "      character*10 shortname(0:MAX_PROCS-1)"
		print "      character*26 datim"
		print "      integer namelen(0:MAX_PROCS-1)"
		print "      character*256 strn"
		print "      logical lgunit"
		print "      integer procid,master,i_data_size"
                print "      integer:: n=0, cast=4"
		print "C use i_data_size in order to cast an address to integer"
		print "#endif"
		print "C ... Passed parameters"
		print "      logical mlog"
		print spc,"type("+struc_name+"):: struc"
		print "      character funnam*(*), label*(*)"

		print "#if MPI|MPIK"
		print spc,"master=0"
		print spc,"call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )"
		print spc,"call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )"
		for vartype_var in struc_def:
			vartype,var=vartype_var
			if len(vartype)!=1:
				print thisprogram, "skip",var[0],"because it is a pointer."
			else:
				# not a pointer
				var_simple=re.sub("\(.*\)","",var[0])
				is_array=0
				if var[0].find("(")!=-1:
					is_array=1
				# remove array def
				if vartype[0]=="integer(8)":
					transfertype="MPI_INTEGER8"
				elif vartype[0]=="real(8)":
					transfertype="MPI_REAL8"
				elif vartype[0]=="integer":
					transfertype="MPI_INTEGER"
				elif vartype[0]=="real":
					transfertype="MPI_REAL"
				else:
					print "unknown type",vartype,var
					sys.exit(10)
				if is_array==1:
					print spc,"i_data_size=size(struc%"+var_simple+")",
					print "!",var[0]
					transfersize="i_data_size"
				else:
					transfersize="1"
				print spc,"call mpi_bcast(struc%"+ var_simple+",",
				print transfersize+","+transfertype
				print spccont,", master, MPI_COMM_WORLD,ierr)"
		print 
		print "      if (mlog) then"
		print "        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)"
		print "        call strcop(shortname(procid),name,10,'.',ierr)"
		print "        namelen(procid) = ierr-1"
		print "        call gettime(datim)"
		print "        strn = ' '//funnam//' '//datim//' Process %i of %i on '"
		print "     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//"
		print "     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'"
		print "        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,"
		print "     .  cast)"
		print "      endif"
		print "#endif"
		print "      end subroutine",routine_name
		print 

#################################################

doit=DoIt()

doit.read()
doit.strucdef()

print thisprogram,", make *_init, *_mpibc1. ",datetime.datetime.today()
doit.show()
print 
#doit.make_bcast()
doit.make_init()

