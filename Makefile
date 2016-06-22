all:
	(time python3 iid_main.py urand.bin 8 -p 4 -v) 1>> logs/results.log 2>> logs/times.log 
	(time python3 iid_main.py urand.bin 8 -v) 1>> logs/results.log 2>> logs/times.log
