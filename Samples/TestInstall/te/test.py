from comp import test1_check,test2_check,runprogs,rmfiles
def test(args,bindir,testdir,workdir):
    lmfa= f'mpirun -np 1 {bindir}/lmfa '
    lmf = f'mpirun -np {args.np} {bindir}/lmf '
    message1='''
    # Case te: molecular statics in an open structure
    # --- Test 1.  Basic check of programs lmfa,lmf ---
    Four pre-generated ctrlG variants are switched in via cp before each call:
       ctrlG.te.lmfa.toml   (default,       used by lmfa)
       ctrlG.te.float.toml  (nbas=12,       w/ floating orbitals)
       ctrlG.te.pw.toml     (KMXA=5,PWMODE=11,  with PWs)
       ctrlG.te.mto.toml    (MTO only)
'''
    print(message1)
    outfile='out.lmf.te'
    rmfiles(workdir,[outfile])
    cp = lambda tag: f"cp {testdir}/ctrlG.te.{tag}.toml ctrlG.te.toml"
    runprogs([
        cp('lmfa'),  lmfa+ "te > "+outfile,
        cp('float'), lmf + "te >> "+outfile,
        "rm -f *mixm.te",
        "cp rst.te rst.te.bk",
        cp('pw'),    lmf + "te >> "+outfile,
        "rm -f *mixm.te",
        "cp rst.te.bk rst.te",
        cp('mto'),   lmf + "te >> "+outfile,
    ])
    tall=test1_check(testdir+'/'+outfile, workdir+'/'+outfile)
    return tall
