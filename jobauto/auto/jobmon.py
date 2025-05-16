import os
import sys
import time
import shutil
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


def rsync(local, remote, user, remotehost, to_remote=True, includes=None, verbose=False):
    """Rsync between local and remote."""
    base_cmd = ["rsync", "-avz"]
    if includes:
        for inc in includes:
            base_cmd += ["--include", inc]
        base_cmd += ["--exclude", "*"]
    if to_remote:
        print("rsync to remote")
        src = str(local) + "/"
        dst = f"{user}@{remotehost}:{remote}/"
    else:
        print("rsync from remote")
        src = f"{user}@{remotehost}:{remote}/"
        dst = str(local) + "/"
    cmd = base_cmd + [src, dst]
    if verbose:
        subprocess.run(cmd, check=True)
    else:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print()
    #print("rsync done:", cmd)

# def rsync(local, remote, user, remotehost, to_remote=True, includes=None):
#     """Rsync between local and remote."""
#     #base_cmd = ["rsync", "-avz", "--quiet"]
#     base_cmd = ["rsync", "-avz"]
#     if includes:
#         for inc in includes:
#             base_cmd += ["--include", inc]
#         base_cmd += ["--exclude", "*"]
#     if to_remote:
#         print("rsync to remote")
#         src = str(local) + "/"
#         dst = f"{user}@{remotehost}:{remote}/"
#     else:
#         print("rsync from remote")
#         src = f"{user}@{remotehost}:{remote}/"
#         dst = str(local) + "/"
#     cmd = base_cmd + [src, dst]
#     subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# def rsync(local, remote, user, remotehost, to_remote=True, includes=None):
#     """Rsync between local and remote."""
#     base_cmd = ["rsync", "-avz"]
#     if includes:
#         for inc in includes:
#             base_cmd += ["--include", inc]
#         base_cmd += ["--exclude", "*"]
#     if to_remote:
#         print("rsync to remote")
#         src = str(local) + "/"
#         dst = f"{user}@{remotehost}:{remote}/"
#     else:
#         print("rsync from remote")
#         src = f"{user}@{remotehost}:{remote}/"
#         dst = str(local) + "/"
#     cmd = base_cmd + [src, dst]
#     subprocess.run(cmd, check=True)

import time
import re

def ssh_cmd(remote_dir, cmd, user, remotehost, retries=3, timeout=5):
    """Run a command on the remote server with retry and timeout. Returns first number found in output."""
    ssh_command = f"ssh {user}@{remotehost} 'cd {remote_dir} && {cmd}'"
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
            for p in parts:
                m = re.search(r'\d+', p)
                if m:
                    return m.group()
            return None
        except subprocess.TimeoutExpired:
            print(f"ssh_cmd timeout (attempt {attempt+1}/{retries}), retrying...")
            time.sleep(1)
    print("ssh_cmd failed after retries.")
    return None

def ssh_cmd0(remote_dir, cmd, user, remotehost, retries=3, timeout=5):
    """Run a command on the remote server with retry and timeout. Returns output as list of words."""
    ssh_command = f"ssh {user}@{remotehost} 'cd {remote_dir} && {cmd}'"
    for attempt in range(retries):
        try:
            result = subprocess.run(
                ssh_command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=timeout
            )
            return result.stdout.strip().split()
        except subprocess.TimeoutExpired:
            print(f"ssh_cmd0 timeout (attempt {attempt+1}/{retries}), retrying...")
            time.sleep(1)
    print("ssh_cmd0 failed after retries.")
    return []


# def ssh_cmd(remote_dir, cmd, user, remotehost):
#     """Run a command on the remote server."""
#     ssh_command = f"ssh {user}@{remotehost} 'cd {remote_dir} && {cmd}'"
#     #print("ssh_cmd remote_dir:", ssh_command)
#     result = subprocess.run(ssh_command, shell=True, capture_output=True, text=True)
#     #print("ssh_cmd result:", result)
#     parts= result.stdout.strip().split()
#     parts= " ".join(parts).split('.')
#     #print('ppp111',parts)
#     first_number = next((p for p in parts if p.isdigit()), None)
#     #print('ppp222',first_number)
#     #print('ssh_cmd: ',ssh_command)
#     return first_number

# def ssh_cmd0(remote_dir, cmd, user, remotehost):
#     """Run a command on the remote server."""
#     ssh_command = f"ssh {user}@{remotehost} 'cd {remote_dir} && {cmd}'"
#     #print("ssh_cmd0 remote_dir:", ssh_command)
#     result = subprocess.run(ssh_command, shell=True, capture_output=True, text=True)
#     #print("ssh_cmd0 result:", result)
#     return result.stdout.strip().split()

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
    ssh_command = f"ssh {user}@{remotehost} 'mkdir -p {remote_dir}'"
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
    parser.add_argument("--qstatcommand", required=True,help="qstat command")
    args = parser.parse_args()
    print(args)
    
    local_dir = Path(args.ldir).resolve()
    remote_dir = args.rdir
    user = args.user
    remotehost = args.remote
    quiet= args.quiet
# --- make remote directory
    ssh_mkdir(remote_dir, user, remotehost)
    
    # date
    date_str = get_now_str()
    local_date_dir = local_dir / date_str
    remote_date_dir = f"{remote_dir}/{date_str}"
    print(f"LOCAL : {local_date_dir}")
    print(f"REMOTE: {remote_date_dir}")
    
    def local2remote(local_path):
        """Convert local path to remote path."""
        return str(local_path).replace(str(local_dir), remote_dir, 1)
    def remote2local(remote_path):
        """Convert remote path to local path."""
        return str(remote_path).replace(remote_dir, str(local_dir), 1)

    # initial sync LOCAL and REMOTE
    print(f"Replacing keywords in qsub files in {local_date_dir} and syncing to {remote_date_dir}")
    #local_date_dir = local_dir / date_str
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

    if args.initonly:
        print("Initial sync done (--initonly). Exiting.")
        return
    
    quelist_path  = local_date_dir / "quelist"
    qsub_now_path = local_date_dir / "qsub.now"
    qsub_log_path = local_date_dir / "qsub.log"

    
    iter=0
    interval=5
    while True:
        # 2. リモートfileをローカルにrsync
        print()
        print(f"=================  watch ============== timing= {iter*interval} sec")
        print(f"remote_date_dir= {user}@{remotehost} {remote_date_dir}")
        rsync(local_date_dir, remote_date_dir, user, remotehost, to_remote=False, includes=
              ["*/","*out*", "*finished*", "qsub.failed.*","llmf*","ctrl*","save*","atmpnu*","rst*"],
              verbose=not quiet)

        # Add submittable qsub to the end of quelist. Check dependency.
        for dep_file in local_date_dir.rglob("qsub.dependency.*"):
            id = dep_file.name.split(".")[-1]
            lwork=Path(dep_file).parent
            if lwork == local_date_dir: continue
            qsub_file = lwork / f"qsub.{id}"
            #print(f"Checking dependency file: {dep_file}")
            #print(f"Checking       qsub file: {qsub_file}")
            depok = True
            if dep_file.stat().st_size > 0:
                with open(dep_file) as f:
                    for line in f:
                        fname = line.strip()
                        if fname and not (local_date_dir / fname).exists():
                            depok = False
                            break
            if depok :
                quelist = read_quelist(quelist_path)
                if not any(str(qsub_file) in line for line in quelist):
                    with open(quelist_path, "a") as f:
                        f.write(f"{qsub_file}\n")
        # print("quelist:")
        # quelist = read_quelist(quelist_path)
        # for line in quelist:
        #     print('quelist:',line)
        
        # qsub.nowのIDをqstatで監視
        qstatcommand = args.qstatcommand
        running = 0
        if qsub_now_path.exists():
            with open(qsub_now_path,"r") as f:
                now_lines = [line.strip() for line in f if line.strip()]
            running = 0    
            qnow=[]        
            for line in now_lines:
                jobid, qsub_file = line.split()[:2]
                qstat_out = ssh_cmd0('~/.', qstatcommand+f" {jobid}", user, remotehost)
                qstat_out = " ".join(qstat_out)
                nojob = "not exist" in qstat_out or "has finished" in qstat_out or len(qstat_out)<2
                #print('jjjj qstat ',jobid, qsub_file,'runnning=',not nojob,'output='+qstat_out[0:15])
                finished_file = qsub_file.replace("qsub.", "qsub.finished.")
                #print('jjjj finished_file:',finished_file)
                quelist = read_quelist(quelist_path)
                if nojob or Path(finished_file).exists():  # Note finished_file can apperar before qstat_out 
                    if Path(finished_file).exists():       # finished
                        quelist = [l if not l.startswith(str(qsub_file)) else l + f" finished {get_now_str()}" for l in quelist]
                    else: # failed
                        quelist = [l if not l.startswith(str(qsub_file)) else l + f" failed {get_now_str()}" for l in quelist]
                    continue
                running += 1
                print("qsub running:",line)
                qnow.append(line)
            with open(qsub_now_path, "w") as f:
                for line in qnow:
                    f.write(line + "\n")
        isub=0            
        if running < args.maxqsub:
            quelist = read_quelist(quelist_path)
            nsub=args.maxqsub - running
            isub=0
            for i, line in enumerate(quelist):
                #print("quelist line:",i,line)
                if "started" not in line and "finished" not in line and "failed" not in line:
                    qsub_file  = line.split(" ")[0]
                    remote_qsub_file = qsub_file.replace(str(local_date_dir), remote_date_dir, 1)
                    submit_cmd = f"qsub {remote_qsub_file}"
                    jobid = ssh_cmd(remote_date_dir, submit_cmd, user, remotehost) #jobid : qsubの戻り値
                    print("qsub submission: ",jobid, remote_qsub_file, ' LOCAL=',qsub_file)
                    quelist[i] = line + f" started@{get_now_str()}"
                    with open(qsub_now_path, "a") as f:
                        f.write(f"{jobid} {qsub_file} ! jobid qsub_file\n")
                    with open(qsub_log_path, "a") as f:
                        f.write(f"{jobid} {qsub_file} started@{get_now_str()}\n")
                    isub=isub + 1
                    if(isub >= nsub): break
            if(isub>0): write_quelist(quelist_path, quelist)
            print("quelist:", len(quelist))
        print("running  :", running)
        print("submitted:", isub)

        #with open(qsub_now_path) as f:
        #    lines=f.readlines()
        #    print("qqqqqqq222 qsub_now_path:",qsub_now_path,len(lines),lines)

        # 6. rsync from remote to local
        print()
        if(len(qsub_now_path.read_text())==0): # 5. qsub.nowが空なら終了
            rsync(local_date_dir, remote_date_dir, user, remotehost, to_remote=False, includes=["*/","*out*", "*finished*"])
            print('REMOTE:',str(remote_date_dir))
            print('LOCAL:',str(local_date_dir)+"/quelist:")
            quelist = read_quelist(quelist_path)
            for line in quelist:
                print('quelist:',line)
            break         
        # 7. rsync from local to remote
        time.sleep(interval)
        iter=iter+1 #        if(iter==10): break

if __name__ == "__main__":
    main()
