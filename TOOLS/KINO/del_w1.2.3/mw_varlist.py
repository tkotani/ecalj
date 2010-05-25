class Sw_varlist:
    v=[]
    thisprogram=''
    def __init__(self,thisprogram0='*'):
	v=[]
	self.thisprogram=thisprogram0

    def len(self):
	return len(self.v)

    def append(self,w):
	self.v.append(w)

    def append_uniq(self,w_var,linen,tok,printit=0):
        thisfunc="w_varlist_append_uniq"
        name=w_var.name
        for w in self.v:
                if name==w.name:
                        if printit==1:
                                print thisprogram,thisfunc,"warning, duplicated symbol:",name
                        #print thisprogram,"line=",linen,tok
                        return  self.v
        self.v.append(w_var)

    def show(self):
        for nn in self.v:
                print nn.name,
        print ''

    def addall(self):
        nn=''
        for n in self.v:
                nn=nn+n.name+' '
        return nn

    def release(self,name,level,linenumber,printit=1):
        thisfunc="w_varlist_release"
        idisable_change=0
        list_rel=[]
        if printit==1:
                print self.thisprogram,"rlse name=",name,"old_list=",
                self.show()

        if len(self.v)==0:

                if printit==1:
                        print ''
                        print 'error, try to release name=',name, ',but list=null at linenumber=',linenumber
                        print 'list=',
                        self.show()
                        print >>sys.stderr,'error, try to release name=',name, ',but list=null, at linenumber=',linenumber

                idisable_change=1
                return [list_rel,idisable_change]

        idel=-1
        for i in range(len(self.v)):
                a=self.v[i].name
                if a == name:
                        idel=i
                        break

        if idel>=0:
                list_rel=self.v[idel:]
                tok2= []
                for i in range(idel):
                        tok2.append(self.v[i])
                self.v=tok2
        else:
                if printit==1:
                        print ''
                        print 'error, try to release name=',name, ',but list does not have ',name
                        print "list=",
                        self.show()
                        print ''

        if printit==1:
                print self.thisprogram,"rlse name=",name,"new_list=",
                self.show()
        return [list_rel,idisable_change]


    def contains(self,name,funcname='*'):
        thisfunc='w_var_toundel_contains'
        #print thisfunc,'search', name,funcname 
        for i in range(len(self.v)):
                nn=self.v[i]
                #print thisfunc,'?',name,funcname,nn
                if nn.name==name and (nn.execfuncname==funcname or funcname=='*'):
                        #print thisfunc,name,funcname,'ret=',i
                        return i
        return -1

