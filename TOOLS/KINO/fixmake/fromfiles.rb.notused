class FromFile
  @filelist
  @uselist
  @moduleilst

      @@modpath="$(moddir)/"
      @@fppath= "$(fp_obj_path)/"
      @@subpath="$(subs_obj_path)/"
      @@gwdpath="$(gwd_obj_path)/"
      @@slapath="$(sla_obj_path)/"
#-----------------------
  def initialize(cmd)
      @filelist=[]
      @modulelist=[]
      @uselist=[]
      cmd.each { |file|
         do_file(file)
      }
  end
#-----------------------
  def do_file(file)
      fr=File.open(file)
      modulelist=[]
      uselist=[]
      fr.each { |line|
         if /^ +module / =~ line
             a=line.split(' ')
             modulelist.push(a[1])
         end
         if /^ +use / =~line
            line2=line.split(",")
            line2=line2[0].split(":")
            a=line2[0].split(' ')
            uselist.push(a[1])
         end
      }
      #modulelist.each { |a| 
      #  print a,"\n"
      #}
      #uselist.each { |a| 
      #  print a,"\n"
      #}
      @filelist.push(file)
      ith=@filelist.size
      modulelist.uniq!
      @modulelist[ith-1]=modulelist
      uselist.uniq!
      # delete items in list that are in modulelist
      uselist = uselist-modulelist
      @uselist[ith-1]=uselist
  end

  def show
      nsize=@filelist.size
      for i in Range.new(0,nsize-1)
          if @modulelist[i].size>0 or @uselist[i].size>0 
          print @filelist[i]," , "
          @modulelist[i].each { |m|
             print m," "
          }
          print " : "
          @uselist[i].each { |u|
             print u," "
          }
          print "\n"
          end 
      end
  end

  def show_makedep
      print "show_makedep"
      nsize=@filelist.size
      for i in Range.new(0,nsize-1)
          if @modulelist[i].size>0 or @uselist[i].size>0 
          srcfile=@filelist[i].dup
          srcfile.sub!(/^.*\/fp\//,"fp/")
          srcfile.sub!(/^.*\/subs\//,"subs/")
          srcfile.sub!(/^.*\/gwd\//,"gwd/")
          srcfile.sub!(/^.*\/slatsm\//,"slatsm/")
          objfile=srcfile.dup
          objfile.sub!(/^fp\//,@@fppath)
          objfile.sub!(/^subs\//,@@subpath)
          objfile.sub!(/^gwd\//,@@gwdpath)
          objfile.sub!(/^slatsm\//,@@slapath)
          objfile.sub!(/\.F$/,".o")
          # dependency
          print objfile," "
          @modulelist[i].each { |m|
              print @@modpath,m,".mod "
          }
          print ": ",srcfile, " "
          @uselist[i].each { |u|
              print @@modpath,u,".mod "
          }
          print "\n"
	  # command
          print "\t$(FC) $(FFLAGS) -o ",objfile," -c ",srcfile,"\n"
          end 
      end
  end 
#-----------------------
end

fromfile=FromFile.new(ARGV)
fromfile.show_makedep

