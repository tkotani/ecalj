## Usage:

### preparation at kugui
We use kugui. (or ucgw in osaka-u)
At kugui, set following lines at `.bashrc`. This is usually configured by default.
```
if [ -f /etc/profile ]; then
  . /etc/profile
fi
```

### usage
At ecalj/jobauto/, run `ls -R testSGA`
Then you see init0/ directory.
qsub.0 is default qsub script.
All POSCAR is under POSCARALL/

Look into the bash script `jobtestSGA`.
You have to custumize this for your purpose.
After that, run the following command:
```
> ./jobtestSGA 
```


### Mechanism memo. How it works.

#### 0. starting.
We have initial testSGA0/ as
```
testSGA/
└── init
    ├── POSCARALL
    │   ├── POSCAR.mp-1330
    │   ├── POSCAR.mp-149
               ...
    ├── qsub.0               :qsub template
    └── qsub.dependency.0    :dependency filempty, always)
```

#### 1. initial set up stage 
* We have testSGW/init0/qsub.0,
which is a template of qsub. testSGW/init/qsub.0 is not machine-dependent.
`__foobar__` are replaced with the augments given to `jobmon.py`.
See `jobtestSGA`.

* We appy `python initposcar.py --ldir=testSGA` at the beginnig of jobtestSGA.
This create the directory `testSGA/init` from `init0/`.
`initposcar.py` distribute `qsub.0, qsub.dependency.0, and POSCARALL/POSCAR*` under `testSGA/mp*`

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

#### 2. The main step of ./jobtestSGA
See the bottom of jobtestSGA, where we run
```bash
python auto/jobmon.py --user k041303 --remote kugui --pythonpath=...'
```
(python auto/jobmon.py --help gives a help).

* qsubheader.kugui2.temp is given by the Here Document in jobtestSGA. We include machine-dependent here documents in it.
  Set `--qsubheader=qsubheader.kugui2.temp` to specify which one.

* See [the algorism of jobmon.py](./jobautoplan.md)
