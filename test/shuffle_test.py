import time
import random
import numpy

if __name__ == '__main__':
	
	with open('../bin/urand.bin', 'rb') as file:

		start_time = time.time()

		# Read in and cast data
		bytes_in = bytearray(file.read())
		data = [b & 255 for b in bytes_in]

		read_time = time.time() - start_time

		print("Read time: %fs" % read_time)
		print("Shuffling %d bytes" % len(data))

		shuffles = 10000

		print("Shuffling %d times" % shuffles)

		for i in range(shuffles):
			numpy.random.shuffle(data)

		shuffle_time = time.time() - start_time - read_time

		print("Shuffle time: %fs" % shuffle_time)
		print(data[1:20])

