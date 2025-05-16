## Usage:

### preparation at kugui
We use kugui. (or ucgw in osaka-u)
At kugui set following lines at .bashrc.
```
if [ -z "$PS1" ]; then
    source /etc/profile.d/modules.sh
fi
export PATH=/opt/pbs/bin/:$PATH
```

### usage
>cp -r testSGA0/ testSGA   #note testSGA is empty
>./jobtestSGA   #look into jobtestSGA. 


### Mechanism memo. How it works.
See jobautoplan.md.

#### 0. starting.
We have initial testSGA0/ as
```
testSGA/
└── init
    ├── POSCARALL
    │   ├── POSCAR.mp-1330
    │   ├── POSCAR.mp-149
               ...
    ├── qsub.0            :qsub template
    └── qsub.dependency.0    :dependency filempty, always)
```

#### 1. initial set up stage
Take a look at testSGA0/. We have testSGW0/init/qsub.0,
which is a template of qsub. testSGW0/init/qsub.0 is machine-independent.
Note variables `__foobar__` in it.
`__foobar__` are replaced by the augments given to `jobmon.py`.

At first, we do `cp -r testSGA0/ at testSGA/` 
(note: we need testSGA0/==testSGA/. Careful when testSGA/ already).

Then we appy `python initposcar.py --ldir=testSGA` at the beginnig of jobtestSGA.
This just initialize the directory of testSGA/
That is, to distribute qsub.0, qsub.dependency.0, and POSCARALL/POSCAR* under testSGA/*

After `python initposcar.py --ldir=testSGA`', we have directories as
```
testSGA/
└   ...
    ├── mp-1330
    │   ├── POSCAR.mp-1330
    │   ├── qsub.0
    │   └── qsub.dependency.0
...
```

#### 2. The second step of ./jobtestSGA
This is the main part of job control. See the bottom of jobtestSGA, where we run
```bash
python auto/jobmon.py --user k041303 --remote kugui --pythonpath=/home/k0413/k041303/.local/share/mise/installs/python/3.13.0/bin/python --binpath=/home/k0413/k041303/bin --ldir=testSGA --rdir=/home/k0413/k041303/testSGA111 --maxqsub=3 --maxcore=4 --quiet --qsubheader=qsubheader.kugui2.temp --qstatcommand='qstat '
```
(note: python auto/jobmon.py --help gives a help).

* qsubheader.kugui2.temp is given by Here Document in jobtestSGA. We need another here documents for other machine,
  and set `--qsubheader=qsubheader.kugui2.temp`).
