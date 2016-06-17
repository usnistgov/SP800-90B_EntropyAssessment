all:
	(time python3 iid_main.py urand.bin 8 -p 4 -v) 1>> results.log 2>> times.log 
	(time python3 iid_main.py urand.bin 8 -v) 1>> results.log 2>> times.log
