FC=mpif90
#OPTS=-132 -O3 -parallel
OPTS=
#FFLAGS=-I${HOME}/include
FFLAGS=  -ffixed-line-length-132  -O3
#LDFLAGS=-L${HOME}/lib
LDFLAGS=  -fopenmp /usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf 
LIBS=/usr/lib/libfftw3.so.3 -fftw3 
#PROGS= calc_wfn calc_wannier
PROGS= wanplot

#OBJS=m_DATA4GW.o m_QG.o m_LMTO.o m_FFT3D.o lin3.o mymath.o myYlm.o wfn2dx.o m_MLWF.o
#OBJS2=m_DATA4GW.o m_QG.o m_LMTO.o m_FFT3D.o lin3.o mymath.o myYlm.o wfn2dx_2.o m_MLWF.o cubeformat.o xsfformat2.o  #expand_mesh.o  wfn2dx_abc.o wfnrho_abc.o write_mesh2d.o 
OBJS2=m_DATA4GW.o m_QG.o m_FFT3D.o mymath.o myYlm.o wfn2dx_2.o m_MLWF.o cubeformat.o xsfformat2.o  expand_mesh.o  wfn2dx_abc.o wfnrho_abc.o write_mesh2d.o 

all : cpmod ${PROGS}

#calc_wfn : ${OBJS} calc_wfn.o
#	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS} calc_wfn.o -o $@ ${LIBS}

#calc_wannier : ${OBJS} calc_wannier.o
#	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS} calc_wannier.o -o $@ ${LIBS}

#calc_wfn2 : ${OBJS2} calc_wfn2.o
#	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wfn2.o	-o $@ ${LIBS}

# calc_wannier2 : ${OBJS2} calc_wannier2.o
# 	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier2.o	-o $@ ${LIBS}
# calc_wannier3 : ${OBJS2} calc_wannier3.o
# 	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier3.o	-o $@ ${LIBS}
# calc_wannier4 : ${OBJS2} calc_wannier4.o
# 	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier4.o	-o $@ ${LIBS}
# calc_wannier5 : ${OBJS2} calc_wannier5.o
# 	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} calc_wannier5.o	-o $@ ${LIBS}
wanplot : ${OBJS2} wanplot.o ${GW0}
	${FC}  ${OPTS} ${FFLAGS} ${LDFLAGS} ${OBJS2} $(GW0) $(COMM) ${LDFLAGS} wanplot.o	-o $@ ${LIBS}


#calc_wannier3.o: cubeformat.o xsfformat.o 
#calc_wannier4.o: cubeformat.o xsfformat.o 
#calc_wannier5.o: cubeformat.o xsfformat.o  keyvalue.o 
wanplot.o: cubeformat.o xsfformat2.o  

cpmod:
	cp ${top}/exec/*.mod .
clean :
	rm -f *.o *~ ${PROGS} *.mod

.SUFFIXES: .F
.F.o   :
	${FC} ${OPTS} ${FFLAGS} $< -c

top=../../fpgw
nfpsrc  = ${top}/nfpsrc/
gwsrc   = ${top}/gwsrc/
para_g=.o

COMM=$(nfpsrc)rxx.o \
$(gwsrc)mopen.o 

GW0 = \
 $(gwsrc)m_hamindex.o\
 $(gwsrc)readpomat.o \
 $(gwsrc)keyvalue.o \
 $(gwsrc)rppovl.o \
 $(gwsrc)nocctotg.o \
 $(gwsrc)ppbafp.fal$(para_g) \
 $(gwsrc)psi2b_v2$(para_g) \
 $(gwsrc)psi2b_v3$(para_g) \
 $(gwsrc)wfacx.o \
 $(gwsrc)sortea.o \
 $(gwsrc)rydberg.o \
 $(gwsrc)polinta.o \
 $(gwsrc)efsimplef.o \
 $(gwsrc)extension.o \
 $(gwsrc)rangedq.o \
 $(gwsrc)nword.o \
 $(gwsrc)scg.o \
 $(gwsrc)matm.o \
 $(gwsrc)rdpp.o \
 $(gwsrc)mptauof.o \
 $(gwsrc)genallcf_mod.o \
 $(gwsrc)rgwinf_mod.o \
 $(gwsrc)rotdlmm.o \
 $(gwsrc)iopen.o \
 $(gwsrc)cputid.o \
 $(gwsrc)rw.o \
 $(gwsrc)ext.o \
 $(gwsrc)ext2.o \
 $(gwsrc)cross.o \
 $(gwsrc)mate.o \
 $(gwsrc)mate1.o \
 $(gwsrc)icopy.o \
 $(gwsrc)bib1.o \
 $(gwsrc)index.o \
 $(gwsrc)idxk.o \
 $(gwsrc)maxnn.o \
 $(gwsrc)reindx.o \
 $(gwsrc)iprint.o \
 $(gwsrc)bz.o \
 $(gwsrc)bzmesh.o \
 $(gwsrc)genqbz.o \
 $(gwsrc)linpackdummy.o \
 $(gwsrc)switches.o \
 $(gwsrc)rwbzdata.o \
 $(gwsrc)llnew.o  \
 $(gwsrc)readeigen.o \
 $(gwsrc)readqg.o \
 $(gwsrc)iqindx.o  \
 $(gwsrc)alloclist.o \
 $(gwsrc)MPI_fpgw.o
