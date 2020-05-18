import numpy as np

filename = 'test_1.xyz'

ts = []
labels = []
trajs = []

with open(filename, 'r') as input_file:

	N = int( input_file.readline() )
	for i in range(N): trajs.append([])

	for line in input_file:
		if len(line.split()) == 1:
			continue
		if line.split()[0].split('.')[-1] == 'xyz':
			t = float( line.split()[3] )
			ts.append( t )
			counter = 0
		else:
			if len(ts) == 1:
				label = line.split()[0]
				labels.append(label)
			position = np.array( [ float(line.split()[i]) for i in range(1, 4) ] )
			trajs[counter].append(position)
			counter += 1

# print(ts)
# print(labels)
print(np.array(trajs))
