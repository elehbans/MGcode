# Created by Evan Henrich for Matrix Genetics, 16 Jul 2015


# This script allows the user to input files in FASTA format with RNA/DNA sequences less
# than 500nt and then run the RBS calculator code by Salis Lab on that dna.
# It is recommended to use only plus / minus 50bp of the start codon, if the 
# start codon is known.

# Finally it returns as output a csv file with the results of the RBS Calculator
# in the same directory. The results csv includes the sequence name from the 
# FASTA file, the number of bp between the middle A in the GGAGG shine-dalgarno
# sequence from E. coli and the start codon, and the predicted expression level.

#--------------------------------------------------------------------------------
# NOTE: I modified RBS_Calculator to enable us to pull aligned spacing info
from RBS_CalculatorEH import RBS_Calculator
import sys, math, random, string, csv

#-----------define-class-and-functions------------------------------------------

class RunSeq:
    
    # Variables used by multiple functions within class
    filename = []
    seq_dict = {}
    spacing_ids_dict = {}
    results_file = []
    input_file = []
    expr_list = []
    start_pos_list = []
    
    def GetInputData(self):
      
        input_file = open(filename)
        input_lines = open(filename).readlines()
        
        line_number = 0
        
        #Take in data from format of odd rows with '>seqname and even rows sequence, put this in a dictionary with name as key, seq as value
        while line_number < len(input_lines):
            self.seq_name = input_lines[line_number].strip('\r\n')
            self.seq_name = self.seq_name.replace('>','')
            seq_itself = input_lines[line_number+1].strip('\r\n')
            seq_itself = seq_itself.replace(' ','')
            self.seq_dict[self.seq_name] = seq_itself
            line_number += 2
        
        #print self.seq_dict
    
    def AssignSpacingValues(self):
        
        # The spacing_ids csv was developed by testing the RBS calculator with
        # sequences with known spacing and then deriving what calculated dG 
        # had been given to these spacings.
        reader = csv.reader(open('spacing_ids.csv', 'r'))
        
        # Write the spacing, energy pairings to a dictionary
        for row in reader:
            k, v = row
            self.spacing_ids_dict[k] = v
   
   
    def RunRBSCALC(self,seq_name,seq_itself):
       
        # Run the RBS_Calculator on the contents of the input data from seq_dict
        start_range = [0, len(seq_itself)]
        calcObj = RBS_Calculator(seq_itself, start_range, self.input_file)
        calcObj.calc_dG()
        
        # Pull list of energies and start position of translation from RBS_CalculatorEH.py
        self.dG_total_list = calcObj.dG_total_list[:]
        self.dG_spacing_list = calcObj.dG_spacing_list[:]
        self.dG_mRNA_list = calcObj.dG_mRNA_list[:]
        self.dG_mRNA_rRNA_list = calcObj.dG_mRNA_rRNA_list[:]
        self.dG_start_energy_list = calcObj.dG_start_energy_list[:]
        self.dG_standby_site_list = calcObj.dG_standby_site_list[:]
        self.start_pos_list = calcObj.start_pos_list[:]
       
        # Create new lists to hold the number of spaces between start site / RBS and expression
        self.dG_num_spaces_list = []
        self.expr_list = []
        
        # Calculate expression from dG total
        for dG in self.dG_total_list:
            self.expr_list.append(calcObj.K * math.exp(-dG/calcObj.RT_eff))
        
        # Pulling aligned spacing was a pain since it's stuck inside a method
        # So instead, I pull the dG from the spacing list and convert that
        # back to the bp difference in aligned spacing later.
        for dGS in self.dG_spacing_list:
            dGS = float(dGS)
            dGS = "%.4f" % round(dGS, 4)
            
            if dGS in self.spacing_ids_dict:
                self.dG_num_spaces_list.append(self.spacing_ids_dict[dGS])
            else:
                self.dG_num_spaces_list.append(dGS)
                
            
    def WriteOutputData(self,seq_name):
        with open(self.results_file, 'a') as csvfile:
    
            self.target = csv.writer(csvfile, 
                            delimiter=',',
                            quotechar='|', 
                            quoting=csv.QUOTE_MINIMAL,
                            )
        
            # Write rows of sequence name, start position, and expression level from RBSCalc lists
            for (spacing,start_pos,expr,dG_total,dG_start,dG_mRNA_rRNA,dG_Spacing,dG_Standby,dG_mRNA) in zip(self.dG_num_spaces_list,
                            self.start_pos_list,
                            self.expr_list,
                            self.dG_total_list,
                            self.dG_start_energy_list,
                            self.dG_mRNA_rRNA_list,
                            self.dG_spacing_list,
                            self.dG_standby_site_list,
                            self.dG_mRNA_list
                            ):
                if start_pos == 35 or start_pos == 34:
                    self.target.writerow([seq_name,
                    str(spacing),
                    str(start_pos + 1),
                    str(round(expr, 3)),
                    str(round(dG_total, 3)),
                    str(round(dG_start, 3)),
                    str(round(dG_mRNA_rRNA, 3)),
                    str(round(dG_Spacing, 3)),
                    str(round(dG_Standby, 3)),
                    str(round(dG_mRNA, 3))
                    ])
            
            

#-------------Execute-code------------------------------------------------------

print "Before running this script, be sure each nucleotide sequence is less than 500nt \n"

Do_work = RunSeq()

filename = raw_input("what is the name of your file? ")
        
Do_work.GetInputData()
Do_work.AssignSpacingValues()

filenameforcsv = filename[:-4]
         
# Create a new temporary results filename
Do_work.results_file = "RBS_results_" + "".join([random.choice(string.digits) for x in range(6)]) + "_" + filenameforcsv + ".csv"

with open(Do_work.results_file, 'a') as csvfile:
    
    Do_work.target = csv.writer(csvfile, 
                    delimiter=',',
                    quotechar='|', 
                    quoting=csv.QUOTE_MINIMAL,
                    )
    
    # Write Headers
    Do_work.target.writerow(['Sequence Name',
                'Spacing nt',
                'Translation Start',
                'Expression',
                'dG_Total',
                'dG_Start_Site',
                'dG_mRNA_rRNA',
                'dG_Spacing',
                'dG_Standby_Site',
                'dG_mRNA'
                ])

for key, value in Do_work.seq_dict.iteritems():
    Do_work.RunRBSCALC(key,value)
    Do_work.WriteOutputData(key)

print "Your file is named %s" % Do_work.results_file
    
