class ReadFile
  @defline 
  @target
  @source
  def initialize(filename)
     @defline=[]
     fr=File.open(filename,"r")
     iline=0
     fr.each { | line | 
        line.chop!()
        if /^\t/ =~ line then
           @defline[iline][0] +=1
           ith=@defline[iline][0]
           @defline[iline][ith+1]=line
	else
           iline+=1
	   if line.size>0
	   @defline[iline]=[]
           @defline[iline][0]=0
	   @defline[iline][1]=line
           end
	end
     }
     @defline.compact!
     #@defline.each { |line|
     #   print line.class," ",line[1],"\n"
     #}
     nline=@defline.size
     #nline=iline
     iline=0
     while iline<nline
        iline+=1
        if @defline[iline].class==NilClass then
           next
        end
        ncmd=@defline[iline][0]
	#print @defline[iline][1],"\n"
        #i=1
	#while i<=ncmd
	#   i=i+1
        #   print @defline[iline][i],"\n"
        #end 
     end 
     fr.close
  end

  def analyze
     @target=[]
     @source=[]
     nline=@defline.size
     iline=0
     while iline <nline
        iline+=1
        line=@defline[iline]
        if line.class==NilClass  
           next
        end
        defline=line[1]
        if /^#/ =~ defline 
           next 
        end
        if /:=/ =~ defline
           next
        end
        if  /:/ !~ defline
           next
        end

        #print defline,"\n"
           @target[iline]=[]
           @source[iline]=[]
           sp=defline.split(':')
           @target[iline]=sp[0].split(' ')
           @source[iline]=sp[1].split(' ')
           #print @target[iline],"\n"
        #ntarget=@target[iline].size
        #for itarget in Range.new(0,ntarget-1)
        #  print iline," ",ntarget," target ",@target[iline][itarget],"\n"
        #end
        #nsource=@source[iline].size
        #isource=0
        #for isource in Range.new(0,nsource-1)
        #  print iline," ",nsource, " source ",@source[iline][isource] , "\n"
        #end
     end
  end

  def del_dup
     nline=@target.size
     for iline1 in Range.new(0,nline-1)
        if @target[iline1].class==NilClass
           next
        end
     for iline2 in Range.new(iline1+1,nline-1)
        if @target[iline2].class==NilClass
           next
        end
        idup=0
        @target[iline1].each { | target1 |
        @target[iline2].each { | target2 |
           if target1 == target2 
              idup=1
              break
           end
        }
           if (idup==1) 
              break
           end
        }
        if idup==1
           message="#dup line(s): "
           print message,"\n"
           print "#def ",iline1," ", @defline[iline1][1],"\n"
           print "#def ",iline2," ", @defline[iline2][1],"\n"

	   # use the first cmd, use the last def
	   @target[iline1]=@target[iline2]
           @defline[iline1][0]=@defline[iline2][0]
           @defline[iline1][1]=@defline[iline2][1]
           @target[iline2]=nil
           @source[iline2]=nil
           @defline[iline2]=nil
           nsize=@defline[iline1].size
           for i in Range.new(1,nsize-1)
              print "#use ",@defline[iline1][i],"\n"
           end
        end
     end
     end
  end 

  def print_modified
     nline=@defline.size
     for iline in Range.new(0,nline-1)
        line=@defline[iline]
        if line.class==NilClass
           next
        end
        target=@target[iline]
        if target.class==NilClass
	   if /^##/ =~ line[1]
              print "\n"
           end
           print line[1],"\n"
           next
        end
        target.each  { |t|
           print t," "
        }
        print ": "
        source = @source[iline]
        if source.class==NilClass
           next
        end
        source.each  { |s|
           print s," "
        }
        print "\n"
	nsize=@defline[iline].size
	for  i in Range.new(2,nsize-1)
	   print @defline[iline][i],"\n"
        end
        print "\n"
     end
  end

end 

file=ReadFile.new("Make.inc.ifort")
file.analyze
file.del_dup
file.print_modified

