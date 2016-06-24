import binascii

f = open('beacon_rand.txt', 'rb')

data = binascii.unhexlify(f.read())

f.close()

g = open('beacon_rand.bin', 'wb')

g.write(data)

g.close()