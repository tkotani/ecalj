# 2024-3-12 new test system in python ===

## Usage
For tests, we have testecalj.py, which call comp.py and difnum0.py.

>./testecalj [-np mpi_size] [list of tests]

To know what test we have, run
>./testecalj show

To run only copt and si_gwsc test with mpi_size=8, run
>./testecalj -np 8 copt si_gwsc

To run all tests, 
>./testecalj -np 8

-----------
## How the testecalj.py work?
For each test, we had performed calculations and keep files such as
{out.foobar, EPS*,QPU,QPD,log.*}.
In the manner described in testecalj.py, we reproduce such calculations,
and compare these filed.

## How to add your test
If you like to add new test, you do run calculations
and keep some files. Then you add your own test as
```python
    elif(tname=='yourtest'):
        ...
		runlist
		...
		how to compare
		...
```


