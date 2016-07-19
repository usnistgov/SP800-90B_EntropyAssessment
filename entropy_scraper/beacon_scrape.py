# pip3 install requests
import requests

# pip3 install xmltodict
import xmltodict

import time
import sys
import binascii

# Work slightly in the past in case clocks aren't sync'ed
cur_time = int(time.time()) - 120

# Round down to a multiple of 60 and make sure it has been published
cur_time = cur_time - (cur_time % 60) - 60

# 64 bytes are produced per minute
# We need 10 days, 20 hours, 25 minutes for exactly 1mil bytes.
# This is 937500 seconds (or 15625 API hits)
start_time = cur_time - 937500

# Open file for writing
f = open('beacon_rand2.bin', 'ab')

# Loop backwards from the current time
for t in range(cur_time, start_time, -60):

	url = 'https://beacon.nist.gov/rest/record/' + str(t)
	
	# This loop avoids the API cap on the beacon.
	# We are restricted to 100 hits per minute, so if a 
	# request isn't satisfied then keep trying until it is
	successful = False
	while not successful:

		# Try to get a random value
		try:

			# Hit the API and decode the XML
			response = requests.get(url)
			xml = xmltodict.parse(response.text)

			# Decode ascii-hex characters to binary and write
			rand = xml['record']['outputValue']
			f.write(binascii.unhexlify(rand))
			
			# If no error is raised in the previous line
			# move on to the next time to read
			successful = True

		# If no value gets pulled, wait a bit and continue
		except KeyError:
			
			# Lets not DDOS them...
			time.sleep(60)

# Done
f.close()