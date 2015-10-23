

#---------------------------------------
from MG_design_module import RBS_Design
import csv, random, string, sys

RunNow = RBS_Design()

# input csv with batch to run and number of iterations on each
input_file = "test_batch5.csv"
iter_per_seq = 10

# print list of output files
name = input_file.strip('.csv')
output_filename = "Batch_Result_" + "".join([random.choice(string.digits) for x in range(3)]) + "_" + name + ".txt"

print "Output file will be %s" % output_filename

with open(input_file, 'rb') as csvfile:
    Input_info = csv.reader(csvfile, delimiter=' ', quotechar='|')

    # create dictionary for seq name and output file string
    output_dict = {}
    
    # create dictionary for each sequence and condition
    for row in Input_info:
        lame_string = ''.join(row)
        new_list = lame_string.split(",")
        new_string = ' '.join(new_list)
        identifier = new_string[:2]
        arg_string = new_string[2:]
        
        for i in range(0,iter_per_seq):
            print "working on version " + str(i) + " " + str(arg_string)
            RunNow.init_var(arg_string)
            result_filename = RunNow.run()
            spec_id = identifier + "_" + str(i)
            output_dict[spec_id] = result_filename
            print "result: " + str(i) + " " + str(result_filename)
            i += 1

output_file = open(output_filename, "w")
output_file.write("Sequence and file names: \n")

for k, v in output_dict.iteritems():
    output_file.write(str(k) + " " + str(v) + "\n")

output_file.close()
    

print "Batch %s is done" % output_filename
