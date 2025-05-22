import os,sys,time,shutil,re
import subprocess
from pathlib import Path
from datetime import datetime
import argparse

def replace_in_file(filepath, replacements):
    """Replace keys in file with values from replacements dict."""
    with open(filepath, 'r') as f:
        content = f.read()
    for key, value in replacements.items():
        content = content.replace(key, value)
    with open(filepath, 'w') as f:
        f.write(content)

def rsync(local, remote, user, remotehost, to_remote=True, includes=None, verbose=False, retries=200, timeout=300):
    """Rsync between local and remote with retry and timeout."""
    base_cmd = ["rsync", "-e", "ssh -x", "-avz", "--timeout", str(timeout)]
    if includes:
        for inc in includes:
            base_cmd += ["--include", inc]
        base_cmd += ["--exclude", "*"]
    if not to_remote:
        base_cmd.append("--delete") # delete files not in remote
    print()
    if to_remote:
        print("rsync to remote")
        src = str(local) + "/"
        dst = f"{user}@{remotehost}:{remote}/"
    else:
        print("rsync from remote")
        src = f"{user}@{remotehost}:{remote}/"
        dst = str(local) + "/"
    cmd = base_cmd + [src, dst]
    for attempt in range(retries):
        try:
            if verbose:
                subprocess.run(cmd, check=True, timeout=timeout)
            else:
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)
            print("rsync done:", cmd)
            return
        except subprocess.TimeoutExpired:
            print(f"rsync timeout (attempt {attempt+1}/{retries}), retrying...")
            time.sleep(1)
        except subprocess.CalledProcessError as e:
            print(f"rsync failed (attempt {attempt+1}/{retries}): {e}, retrying...")
            time.sleep(1)
    print("rsync failed after retries.")

def ssh_cmd(remote_dir, cmd, user, remotehost, retries=1000, timeout=60):
    """Run a command on the remote server with retry and timeout. Returns first number found in output."""
    ssh_command = f"ssh -x {user}@{remotehost} 'cd {remote_dir} && {cmd}'"
    #print("ssh_cmd:", ssh_command)
    for attempt in range(retries):
        try:
            result = subprocess.run(
                ssh_command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            parts = result.stdout.strip().split()
            # 数字部分を抽出
            #print("ssh_cmd output:", parts)
            #print("ssh_cmd output:", result.stdout.strip(),result.stderr.strip(),len(result.stderr.strip()))
            for p in parts:
                m = re.search(r'\d+', p)
                if m:
                    return m.group()
            continue #submission failed
        except subprocess.TimeoutExpired:
            print(f"ssh_cmd timeout (attempt {attempt+1}/{retries}), retrying...")
            time.sleep(5)
    print("ssh_cmd failed after retries.")
    return None

def get_now_str():
    return datetime.now().strftime("%Y%m%d_%H%M%S")

def read_quelist(quelist_path):
    if not quelist_path.exists():
        return []
    with open(quelist_path) as f:
        return [line.strip() for line in f if line.strip()]

def write_quelist(quelist_path, lines):
    with open(quelist_path, "w") as f:
        for line in lines:
            f.write(line + "\n")
            
def ssh_mkdir(remote_dir, user, remotehost):
    """Create remote directory if it does not exist."""
    ssh_command = f"ssh -x {user}@{remotehost} 'mkdir -p {remote_dir}'"
    subprocess.run(ssh_command, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(
        description=
        """Automatic remote qsub control from local machine. Run after LOCAL/foobar is initialized.
           For example >python initposcar.py --ldir testSGA0
        """,
        epilog="example: python auto/jobmon.py --pythonpath=/home/takao/.pyenv/versions/3.9.13/bin/python --binpath=/home/takao/bin --ldir=testSGA0 --rdir=/home/takao/testSGA0 --maxqsub=3 --maxcore=2"
                             )
    parser.add_argument("--ldir", required=True, help="local dir   LOCAL/foobar, in which we do initial set up.")
    parser.add_argument("--rdir", required=True, help="remode dir REMOTE/foobar, whih is generated automatically.")
    parser.add_argument("--binpath",    required=True,help="all binaries including auto things")
    parser.add_argument("--pythonpath", required=True,help="python path which we use")
    parser.add_argument("--maxcore",    type=int, default=2,help="max core (mpi process) per qsub job")
    parser.add_argument("--initonly",   action="store_true", help="quit after initial sync")
    parser.add_argument("--maxqsub", type=int, default=2, help="max qsub count (default: 2)")
    parser.add_argument("--user",   default="takao", help="remote user name (default: takao)")
    parser.add_argument("--remote", default="ucgw", help="remote machine name (default: ucgw)")
    parser.add_argument("--quiet", action="store_true", help="verbose for rsync")
    parser.add_argument("--qsubheader", required=True,help="machine-dependent qsub header")
    parser.add_argument("--datedir", required=False, help="specify date directory to continue previous run")
    parser.add_argument("--interval", required=False, default=5,help="specify date directory to continue previous run")
    args = parser.parse_args()
    local_dir = Path(args.ldir).resolve()
    remote_dir = args.rdir
    user = args.user
    remotehost = args.remote
    quiet= args.quiet
    
    print()
    print(f"=================  initialization ==============")
    print(args)

    # make remote directory
    ssh_mkdir(remote_dir, user, remotehost)
        # date
    if(args.datedir): 
        date_str = args.datedir
    else:
        date_str = get_now_str()
        
    local_date_dir = local_dir / date_str
    remote_date_dir = f"{remote_dir}/{date_str}"
    print(f"LOCAL : {local_date_dir}")
    print(f"REMOTE: {remote_date_dir}")
    
    # initial sync LOCAL and REMOTE
    print(f"Replacing keywords in qsub files in {local_date_dir} and syncing to {remote_date_dir}")
    if not local_date_dir.exists():
        shutil.copytree(local_dir / "init", local_date_dir)
        # qsub.{numbers} are substituted recursively
        qsubheader_content = Path(args.qsubheader).read_text()
        for qsub_file in local_date_dir.rglob("qsub.*"):
            name = qsub_file.name
            if qsub_file.is_file() and name.startswith("qsub.") and name[5:].isdigit() and name.count(".") == 1:
                remote_qsub_file_dir = Path(str(qsub_file.parent).replace(str(local_date_dir), remote_date_dir, 1))
                #print(f"RRReplacing keywords in {qsub_file}")
                #print(f"RRReplacing keywords in {remote_qsub_file_dir}")
                replace_in_file(qsub_file, {
                    "__HEADER__": qsubheader_content,
                    "__jobname__": remote_qsub_file_dir.name,
                    "__binpath__": args.binpath,
                    "__pythonpath__": args.pythonpath,
                    "__cwd__": str(remote_qsub_file_dir),
                    "__maxcore__": str(args.maxcore),
                })
        ssh_mkdir(remote_date_dir, user, remotehost)
        rsync(local_date_dir, remote_date_dir, user, remotehost, to_remote=True)
    else: # skip when restarted
        print(f"date directory {local_date_dir} already exists. Skipping initial sync.")

# get files already
    rsync(local_date_dir, remote_date_dir, user, remotehost, to_remote=False, includes=
    ["*/","out*", "*finished*", "llmf*","ctrl*","save*","atmpnu*","rst*","QPU*","QPD*"],verbose=not quiet)

    if args.initonly:
        print("Initialization of  qsub.0 (--initonly). Exiting.")
        return
    
    quelist_path  = local_date_dir / "quelist"
    if not quelist_path.exists():
        with open(quelist_path, "w") as f: pass # create empty quelist
        quelist = []
    else:
        quelist = read_quelist(quelist_path) 
        quelist = [l.replace("started@", "spre@") for l in quelist]
       
    iter=0
    interval=int(args.interval)     #qstatcommand = args.qstatcommand
    while True:
        print()
        print(f"=================  watch ============== timing= {iter*interval} ")
        print(f"user@machine   = {user}@{remotehost} ")      
        print(f"quelist        = {quelist_path}")            
        print(f"remote_date_dir= {remote_date_dir}")         
        print(f"local_date_dir = {local_date_dir}")
        # rsync from remote to local          
        rsync(local_date_dir, remote_date_dir, user, remotehost, to_remote=False, includes=["*/","save*","*finished*","QPU*","QPD*"])#,verbose=not quiet)
         
        # Add submittable qsub to the end of quelist. Check dependency.
        for dep_file in local_date_dir.rglob("qsub.dependency.*"):
            if dep_file.name.endswith("~"): continue
            id = dep_file.name.split(".")[-1]
            lwork=Path(dep_file).parent
            if lwork == local_date_dir: continue
            qsub_file = lwork / f"qsub.{id}"
            #print()
            #print(f"Checking dependency file: {dep_file}")
            #print(f"Checking       qsub file: {qsub_file}")
            depok = True
            if dep_file.stat().st_size > 0: # dependency file is not empty
                with open(dep_file) as f:
                    for line in f: #Check if all files in dependency file exist
                        fname = local_date_dir/lwork/line.strip()
                        if(line.startswith("not " )): 
                            fnameN=local_date_dir/lwork/line.strip()[4:]
                            if fnameN.exists():
                                depok = False
                        elif not fname.exists():
                            #print(f"  Checking content: {fname} ---> not exists. ")
                            depok = False
                            break
            if depok :
                if not any(str(qsub_file) in line for line in quelist):
                    quelist.append(str(qsub_file))

        # Add finished to quelist
        finished = list(local_date_dir.rglob("qsub.finished.*"))
        #print("finished files:", len(list(finished)))
        #print(f"finished files:, {finished}")
        nfinished=0
        for i,fdir in enumerate(finished):
            nfinished+=1
            fqdir = str(fdir).replace('.finished','')
            #print("finished file:",fdir,fqdir)
            quelist = [l + f" finished@{get_now_str()}" if l.startswith(fqdir) and not 'finished' in l
                       else l  
                       for l in quelist]
        #for line in quelist:
        #    print(line)
        # number of submitted already
        submitted=0
        for line in quelist: #count numfer of lines in quelist which contains "finished" or "started"
            if "finished" not in line and "started" in line: submitted+=1
        
        # new submission
        print()
        isub=0
        if submitted < args.maxqsub:
            nsub=args.maxqsub - submitted
            isub=0
            for i, line in enumerate(quelist): #print("quelist line:",i,line)
                if "finished" not in line and "started" not in line: # and "finished" not in line and "failed" not in line:
                    qsub_file  = line.split(" ")[0]
                    remote_qsub_file = str(qsub_file).replace(str(local_date_dir), remote_date_dir, 1)
                    submit_cmd = f"qsub {remote_qsub_file}"
                    jobid = ssh_cmd(remote_date_dir, submit_cmd, user, remotehost) #jobid : qsub jobid
                    if(jobid==None): 
                        print(f'submission failed. skip {submit_cmd}')
                        continue
                    print("qsub submission: ",jobid, remote_qsub_file, '\n                  LOCAL=',qsub_file)
                    quelist[i] = line + f" started@{get_now_str()}@{jobid}"
                    isub=isub + 1
                    if(isub >= nsub): break
        if(isub>0): print()
        print("quelist          :", len(quelist))
        print("quelist done     :", nfinished)
        print("qsubmax            :", args.maxqsub)
        print("  submitted        :", submitted)
        print("  newly submitted  :", isub)         #print("running actually:", running)
        write_quelist(quelist_path, quelist) 
        
        if all("finished" in line for line in quelist): 
            rsync(local_date_dir, remote_date_dir, user, remotehost, to_remote=False, includes=
              ["*/","*out*", "*finished*", "llmf*","ctrl*","save*","atmpnu*","rst*"])
            print('REMOTE:', remote_date_dir)
            print('LOCAL :', str(local_date_dir)+"/quelist:")
            print(f'--- OK! all jobs finished. rsync done. Check {quelist_path} !')
            break #exit while loop
        print('waiting...')
        time.sleep(interval)
        iter=iter+1 #        if(iter==10): break

if __name__ == "__main__":
    main()
