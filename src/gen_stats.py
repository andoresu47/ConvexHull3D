import subprocess
import csv
import sys

input_file = sys.argv[1]
output_file = '../output/stats/' + input_file.split('.')[0] + '.csv'

with open(output_file,'w') as f:
	writer=csv.writer(f, delimiter=',',lineterminator='\n',)
	
	for i in range(0, 4):
		n = 1 << i
		command = 'mpirun -np {0} ./divide_and_conquer_MPI {1} -sort'.format(n, input_file)
		tokenized = command.split()
		time = [n]
		print "Running with {0} processors".format(n)
		for j in range(1,11):
			time.append(float(subprocess.check_output(tokenized)))
		
		print "Writing to file..."
		writer.writerow(time)
