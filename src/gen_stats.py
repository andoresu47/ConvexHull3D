import subprocess
import csv
import sys

MAX_PROC_POW = 2
MIN_INPUT_POW = 12
MAX_INPUT_POW = 22
MAX_NUM_TESTS = 10

input_family = sys.argv[1]
output_file = '../output/stats/' + input_family + '.csv'

with open(output_file,'w') as f:
	writer=csv.writer(f, delimiter=',',lineterminator='\n',)
	writer.writerow(['p/size', '2^12', '2^13', '2^14', '2^15', '2^16', '2^17', '2^18', '2^19', '2^20', '2^21', '2^22'])
	
	for i in range(0, MAX_PROC_POW + 1):
		n = 1 << i
		time = [n]
		for j in range(MIN_INPUT_POW, MAX_INPUT_POW + 1):
			input_file = input_family + '_' + str(1 << j) + '.in'
			command = 'mpirun -np {0} ./divide_and_conquer_MPI {1} -sort'.format(n, input_file)
		
			tokenized = command.split()
			print "Running test {0}/{1} with {2} processors".format(j - MIN_INPUT_POW, MAX_INPUT_POW - MIN_INPUT_POW, n)
			cumsum = 0.0
			for k in range(0, MAX_NUM_TESTS):
				try:
					cumsum += float(subprocess.check_output(tokenized))
				except subprocess.CalledProcessError:
					continue
			time.append(cumsum / MAX_NUM_TESTS)
		
		print "Writing to file..."
		writer.writerow(time)
