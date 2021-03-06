# learn imprecise hmm model using local IDM
import sys
from array import array

if len(sys.argv) < 6:
	print "Usage: python learn_idm.py N M observations_filename tags_filename model_filename"

N=int(sys.argv[1]) # nr. of states
M=int(sys.argv[2]) # nr. of observations/symbols

# uncomment to choose different prior
#prior = "laplace"
#sN = 1.0*N # laplace priors
#sM = 1.0*M
prior = "perks"
#sN = 0.5 # Perks priors
#sM = 0.5
sN = 1.0 # Perks priors
sM = 1.0
#sN = 2.0 # Perks priors
#sM = 2.0
#prior = "jeffreys"
#sN = N/2.0 # Jeffreys priors
#sM = M/2.0

print 'N:',N
print 'M:',M

a = [] # transition probabilities
b = [] # emission probabilities
p = [] # prior probabilities
O = [] # to count unobserved symbols

step = int(0.2*N)+1 # show amount of initialized array at % steps

print 'initializing count arrays... '
# assign zeros to all possible combinations
for i in xrange(N):
	if i % step == 0:
		print '%d%% ..' % (int(i*100.0/N)),
		sys.stdout.flush()
	a.append( array('I',[0 for j in xrange(N)]) )
	b.append( array('I',[0 for j in xrange(M)]) )
print '100%'
p = [0 for j in xrange(N)]
O = [0 for j in xrange(M)]

# read training  data
o = file(sys.argv[3]).read().splitlines() # observations file
q = file(sys.argv[4]).read().splitlines() # states file

# compute counts from data
if len(o) != len(q):
	print 'Dataset dimension mismatch!'
	exit(-1)
print 'computing counts...'
step = int(0.15*len(o))+1 # show amount of sequences read at % steps
for i in xrange(len(o)): # for each sequence in dataset
	if i % step == 0: 
		print '%d%% ..' % (int(i*100.0/len(o))),
		sys.stdout.flush()
	
	tokens = o[i].split()
	tags = q[i].split()
	if len(tokens) != len(tags):
		print 'Sequence length mismatch!', 'Sequence:', i
		exit(-1)
	# increase prior probability count
	# ignore states grater or equal than N
	if int(tags[0]) < N:
		p[int(tags[0])]+=1
	# increase transition prob count for each token in seq
	for t in xrange(0,len(tokens)-1):
		# ignore states grater or equal than N
		if int(tags[t]) < N and int(tags[t+1]) < N:
			a[int(tags[t])][int(tags[t+1])] += 1

	# increase emission prob count for each token in seq
	for t in xrange(len(tokens)):
		# ignore states and observations not in range
		if int(tags[t]) < N and int(tokens[t]) < M:
			b[int(tags[t])][int(tokens[t])] += 1
			O[int(tokens[t])] += 1
print '100%'

# compute nr. of unseen symbols
unseen = 0
for j in xrange(M):
	if O[j] == 0:
		unseen += 1
print 'there were %d unseen words. MAR will be assumed.' % unseen

imprecision = [1.0,0.0] # min and max imprecision in transition sets
eimprecision = [1.0,0.0] # min and max imprecision in emission sets
print 'saving model to file...'
# save model to file
out = open(sys.argv[5],'w')
# header
out.write(str(N) + '\n')
out.write(str(M) + '\n')
# lower transition probabilities
for i in xrange(N):
	Ni = 0.0
	for j in xrange(N):
		Ni += a[i][j]

	if Ni == 0:
		print '*** found non occurring state', i, 'in transition. Using uniform distribution here.'
		for j in xrange(N):
			prob = 1.0*sN/(Ni+sN)
			out.write(str(prob) + ' ')
	else:
		for j in xrange(N):
			prob = 1.0*a[i][j]/(Ni+sN)
			out.write(str(prob) + ' ')

		imprec = sN/(Ni+sN)

		if imprec < imprecision[0]:
			imprecision[0] = imprec
		if imprec > imprecision[1]:
			imprecision[1] = imprec

out.write('\n')
# upper transition probabilities
for i in xrange(N):
	Ni = 0.0
	for j in xrange(N):
		Ni += a[i][j]

	for j in xrange(N):
		prob = 1.0*(a[i][j]+sN)/(Ni+sN)
		out.write(str(prob) + ' ')
out.write('\n')
# lower emission probabilities
for i in xrange(N):
	Ni = 0.0
	for j in xrange(M):
		Ni += b[i][j]

	if Ni == 0:
		print '*** found non occurring state', i, 'in emission. Using uniform here.'
		for j in xrange(M):
			if O[j] == 0:
				# assume MAR for unseen words
				out.write(str(1.0) + ' ')
			else:
				prob = 1.0/(Ni+sM)
				out.write(str(prob) + ' ')
		
	else:
		imprec = sM/(Ni+sM)

		if imprec < eimprecision[0]:
			eimprecision[0] = imprec
		if imprec > eimprecision[1]:
			eimprecision[1] = imprec		
		for j in xrange(M):
			if O[j] == 0:
				# assume MAR for unseen words
				out.write(str(1.0) + ' ')
			else:
				prob = 1.0*b[i][j]/(Ni+sM)
				out.write(str(prob) + ' ')
out.write('\n')
# upper emission probabilities
for i in xrange(N):
	Ni = 0.0
	for j in xrange(M):
		Ni += b[i][j]

	for j in xrange(M):
		if O[j] == 0:
			# assume MAR for unseen words
			out.write(str(1.0) + ' ')
		else:
			prob = 1.0*(b[i][j]+sM)/(Ni+sM)
			out.write(str(prob) + ' ')
out.write('\n')
# lower prior probabilities

Ni = 0
for i in xrange(N):
	if p[i] == 0:
		print '* found non occurring state', i, 'in prior.'
	Ni += p[i]

pimprecision = sN/(Ni+sM)

for i in xrange(N):
	prob = 1.0*p[i]/(Ni+sN)
	out.write(str(prob) + ' ')
out.write('\n')
# upper prior probabilities
for i in xrange(N):
	prob = 1.0*(p[i]+sN)/(Ni+sN)
	out.write(str(prob) + ' ')
out.write('\n')
out.close()
# imprecision stats
print "minimum imprecision in transition sets:", imprecision[0]
print "maximum imprecision in transition sets:", imprecision[1]
print "  minimum imprecision in emission sets:", eimprecision[0]
print "  maximum imprecision in emission sets:", eimprecision[1]
print "       imprecision in prior credal set:", pimprecision

