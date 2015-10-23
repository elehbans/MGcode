# Created by Evan Henrich on August 10, 2015

#------------------------------------------------------------------------------
from __future__ import division
from RBS_CalculatorEH import RBS_Calculator
from collections import defaultdict
import sys, math, random, string, csv

#------------------------------------------------------------------------------

class RBS_Design:
    
    # setup data structures
    Bad_RBS_list_hyb = {}
    Bad_RBS_list_mRNA = {}
    Bad_CDS_list = {}
    spacing_ids_dict = {}
    
    
    # Instanting index call, (the index number of the CDS Translation start site
    # in the list of start sites). It will be changed later as the RBSCalc
    # starts to run
    index_call = 0
        
        
    def mRNA_Setup(self, preseq, CDS):
        self.start_mRNA = preseq + CDS
        self.CDS_loc = len(preseq)
    
    
    def Make_best_dict(self):

        self.best_mRNA = defaultdict()
        
        for k, v in self.curr_mRNA.iteritems():
            self.best_mRNA[k] = v
        
        
    def first_run_setup(self, Nuc_mode):
        self.curr_mRNA = {}
        
        self.curr_mRNA['seq_name'] = self.seqname_input
        self.curr_mRNA['iteration'] = -1
        self.curr_mRNA['made_best'] = 0
        
        # I added the '-3' because the spacing list is to the middle A in the ggagg consensus sequence
        # And we want to have the generator go all the way to the 3' "G" in terms of iterations
        if Nuc_mode == "All":
            self.curr_mRNA['spacing'] = 0
            self.curr_mRNA['RBS_length'] = 35
        elif Nuc_mode == "Only":
            self.curr_mRNA['spacing'] = int(self.dG_num_spaces_list[self.index_call]) - 3
            self.curr_mRNA['RBS_length'] = 9
        else:
            Nuc_list = self.Nuc_change_input.split(":")
            self.curr_mRNA['RBS_length'] = Nuc_list[1]
            self.curr_mRNA['spacing'] = 35 - (int(Nuc_list[0]) + int(Nuc_list[1]))
        
        self.curr_mRNA['CDS_loc'] = self.CDS_loc
        
        # The 'RBS_start' key holds the value of nt in RBS_length plus nt of space between
        # 3' of RBS and 5' of start codon
        self.curr_mRNA['RBS_start'] = self.curr_mRNA['spacing'] + self.curr_mRNA['RBS_length']
        
        # 'x' is the nt where the RBS actually starts (e.g. if ATG is at 36, spacing is 7, and RBS 
        # length is 9, then x = 36 - 7 - 9 = 20
        self.curr_mRNA['RBS_5p'] = self.curr_mRNA['CDS_loc'] - self.curr_mRNA['RBS_start']
        
        # Pre RBS is anything before 'x'
        self.curr_mRNA['pre_RBS'] = self.Preseq_input[:self.curr_mRNA['RBS_5p']]
                    
        # Post RBS is anything between 3' of RBS and ATG
        self.curr_mRNA['post_RBS'] = self.Preseq_input[(self.curr_mRNA['RBS_5p'] + self.curr_mRNA['RBS_length']):]
        
        self.curr_mRNA['CDS_seq'] = self.CDS_input
            
                        
    def Update_working_dict(self, rbsseq, cdsseq):
        
        # If the start position is the same as the CDS loc, then update working dictionary
        for b in self.start_pos_list:
            
            if b == self.CDS_loc:
                # Set the index_call variable to the index of the desired start codon
                self.index_call = self.start_pos_list.index(b)
                
                # There is a funky edge case where the dG_num_spaces can be 10000000 if non-existent, so this 
                # bypasses that scenario
                if float(self.dG_num_spaces_list[self.index_call]) < 25:
                    
                    self.curr_mRNA['curr_expr'] = self.expr_list[self.index_call]
                    self.curr_mRNA['dG_Hyb'] = self.dG_mRNA_rRNA_list[self.index_call]
                    self.curr_mRNA['dG_mRNA'] = self.dG_mRNA_list[self.index_call]
                    self.curr_mRNA['RBS_seq'] = rbsseq
                    self.curr_mRNA['CDS_seq'] = cdsseq
                    self.curr_mRNA['Num_Spaces'] = self.dG_num_spaces_list[self.index_call]
                    self.curr_mRNA['iteration'] += 1
                    
                    
                    # Full_preseq is the preRBS, RBS, and post-RBS joined together
                    self.curr_mRNA['full_preseq'] = self.curr_mRNA['pre_RBS'] + self.curr_mRNA['RBS_seq'] + self.curr_mRNA['post_RBS']
                
                    # Measure distance of current mRNA dG Hyb to target Hyb
                    self.curr_mRNA['dG_Hyb_dist'] = self.calc_dist(self.Target_Hyb,self.curr_mRNA['dG_Hyb'])
        
                    # Measure the distance of the current best RBS seq dG_mRNA from the target dG_mRNA before running calc on rebuilt RBS
                    self.curr_mRNA['dG_mRNA_dist'] = self.calc_dist(self.curr_mRNA['dG_mRNA'],self.Target_dG_mRNA)
   
                    self.curr_mRNA['tot_vector'] = self.calc_vector(self.curr_mRNA['dG_Hyb_dist'],self.curr_mRNA['dG_mRNA_dist'])
                    
                else:
                    pass
        
    
    def AssignSpacingValues(self):
        # The spacing_ids csv was developed by testing the RBS calculator with
        # sequences with known spacing and then deriving what calculated dG 
        # had been given to these spacings.
        reader = csv.reader(open('spacing_ids.csv', 'r'))
        
        # Write the spacing, energy pairings to a dictionary
        for row in reader:
            k, v = row
            self.spacing_ids_dict[k] = v
    
            
    def RunCalc(self,seq_name,seq_itself):
        # Run the RBS_Calculator on the contents of the input data from seq_dict
        start_range = [0, len(seq_itself)]
        calcObj = RBS_Calculator(seq_itself, start_range, seq_name)
        calcObj.calc_dG()
        
        # Pull list of energies and start position of translation from RBS_CalculatorEH.py
        self.dG_total_list = calcObj.dG_total_list[:]
        self.dG_spacing_list = calcObj.dG_spacing_list[:]
        self.start_pos_list = calcObj.start_pos_list[:]
        self.dG_mRNA_list = calcObj.dG_mRNA_list[:]
        self.dG_mRNA_rRNA_list = calcObj.dG_mRNA_rRNA_list[:]
        
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
                
    
    def dG_Hyb_stat_determine(self, curr_Hyb):
        if self.Hyb_Mode == "P":
            dG_Hyb_status = (curr_Hyb < (1.5 * self.Target_Hyb)) or (curr_Hyb > (0.5 * self.Target_Hyb))
            return dG_Hyb_status
        else:
            dG_Hyb_status = (self.Target_Hyb < curr_Hyb)
            return dG_Hyb_status
            
    
    def dG_mRNA_stat_determine(self, curr_dG_mRNA):
        if self.dG_mRNA_Mode == "P":
            dG_mRNA_status = (curr_dG_mRNA < (1.2 * self.Target_dG_mRNA)) or (curr_dG_mRNA > (0.8 * self.Target_dG_mRNA))
            return dG_mRNA_status
        else:
            dG_mRNA_status = (self.Target_dG_mRNA > curr_dG_mRNA)
            return dG_mRNA_status
   
    
    def setup_codon_dict(self):
        translation_list = [('Isoleucine','ATT'),
                        ('Isoleucine','ATC'),
                        ('Isoleucine','ATA'),
                        ('Leucine','CTT'),
                        ('Leucine','CTC'),
                        ('Leucine','CTA'),
                        ('Leucine','CTG'),
                        ('Leucine','TTA'),
                        ('Leucine','TTG'),
                        ('Valine','GTT'),
                        ('Valine','GTC'),
                        ('Valine','GTA'),
                        ('Valine','GTG'),
                        ('Phenylalanine','TTT'),
                        ('Phenylalanine','TTC'),
                        ('Methionine','ATG'),
                        ('Cysteine','TGT'),
                        ('Cysteine','TGC'),
                        ('Alanine','GCT'),
                        ('Alanine','GCC'),
                        ('Alanine','GCG'),
                        ('Alanine','GCA'),
                        ('Proline','CCT'),
                        ('Proline','CCG'),
                        ('Proline','CCC'),
                        ('Proline','CCA'),
                        ('Threonine','ACT'),
                        ('Threonine','ACA'),
                        ('Threonine','ACG'),
                        ('Threonine','ACC'),
                        ('Serine','TCT'),
                        ('Serine','TCC'),
                        ('Serine','TCA'),
                        ('Serine','TCG'),
                        ('Serine','AGT'),
                        ('Serine','AGC'),
                        ('Tyrosine','TAT'),
                        ('Tyrosine','TAC'),
                        ('Tryptophan','TGG'),
                        ('Glutamine','CAA'),
                        ('Glutamine','CAG'),
                        ('Asparagine','AAT'),
                        ('Asparagine','AAC'),
                        ('Histidine','CAC'),
                        ('Histidine','CAT'),
                        ('Glutamic Acid','GAA'),
                        ('Glutamic Acid','GAG'),
                        ('Aspartic Acid','GAC'),
                        ('Aspartic Acid','GAT'),
                        ('Lysine','AAA'),
                        ('Lysine','AAG'),
                        ('Arginine','CGT'),
                        ('Arginine','CGC'),
                        ('Arginine','CGA'),
                        ('Arginine','CGG'),
                        ('Arginine','AGA'),
                        ('Arginine','AGG')]
        
        self.codon_dict = defaultdict(list)
        
        for k,v in translation_list:
            self.codon_dict[k].append(v)
     
     
    def make_cod_freq_dict(self,codon_filename):
        self.Cod_freq_dict = {}
        self.dict_cod_list = defaultdict(list)
        
        reader = csv.reader(open(codon_filename, 'r'))
        
        # Write the spacing, energy pairings to a dictionary
        for row in reader:
            k, v = row
            self.Cod_freq_dict[k] = v
        
        for k, v in self.codon_dict.iteritems():
            for y in v:
                if y in self.Cod_freq_dict.iterkeys():
                    val = int(self.Cod_freq_dict[y])
                    for b in range(0,val):
                        self.dict_cod_list[k].append(y)
                        b = b + 1
            
    
    def calc_dist(self,x,y):
        a = max(x,y)
        b = min(x,y)
        dist = a - b
        return dist
    
    
    def calc_vector(self,x,y):
        vector = math.sqrt((x*x) + (y*y))
        return vector
    
    
    def quick_check(self,Hyb_to_check,mRNA_to_check,RBSseq):
        
        if self.dG_Hyb_stat_determine(Hyb_to_check) == False and self.dG_mRNA_stat_determine(mRNA_to_check) == False:
            self.curr_mRNA['made_best'] = 4
            self.Make_best_dict()
            return self.results_file
        else:
            pass

    
    def Create_output_csv(self,name):
        # Create a results file name
        self.results_file = name + "_RBS_design_" + "".join([random.choice(string.digits) for x in range(6)]) + ".csv"
        
        with open(self.results_file, 'a') as csvfile:
            
            self.target = csv.writer(csvfile, 
                            delimiter=',',
                            quotechar='|', 
                            quoting=csv.QUOTE_MINIMAL,
                            )
            
            # Write Headers
            self.target.writerow(['Sequence Name',
                        'Pre_sequence',
                        'CDS_sequence',
                        'Iteration',
                        'Made_best',
                        'Expression',
                        'dG_Hybridization',
                        'dG_mRNA',
                        'Spacing',
                        'Total_Vector'
                        ])
                        
    
    def WriteOutputData(self,seq_name,dict_choice):
        if dict_choice == "Curr":
            dict_name = self.curr_mRNA
        else:
            dict_name = self.best_mRNA
        
        with open(self.results_file, 'a') as csvfile:
        
            self.target = csv.writer(csvfile, 
                            delimiter=',',
                            quotechar='|', 
                            quoting=csv.QUOTE_MINIMAL,
                            )
            
            self.target.writerow([dict_name['seq_name'],
                        dict_name['full_preseq'],
                        dict_name['CDS_seq'],
                        dict_name['iteration'],
                        dict_name['made_best'],
                        dict_name['curr_expr'],
                        dict_name['dG_Hyb'],
                        dict_name['dG_mRNA'],
                        dict_name['Num_Spaces'],
                        dict_name['tot_vector']
                        ])
                
    def init_var(self,arg_string):
        arg_list = arg_string.split()
        self.Nuc_change_input = arg_list[0]
        self.dG_Hyb_input = arg_list[1]
        self.dG_mRNA_input = arg_list[2]
        self.MaxIter_input = int(arg_list[3])
        self.Preseq_input = arg_list[4]
        
        if len(self.Preseq_input) < 35:
            print "Presequence is less than 35nt"
            sys.exit()
        
        self.CDS_input = arg_list[5]
        self.seqname_input = arg_list[6]
        self.codon_file = arg_list[7]

        Hyb_list = self.dG_Hyb_input.split(":")
        self.Target_Hyb = float(Hyb_list[1])
        self.Hyb_Mode = Hyb_list[0]
        
        input_list = self.dG_mRNA_input.split(":")
        self.Target_dG_mRNA = float(input_list[1])
        self.dG_mRNA_Mode = input_list[0]
    
    def run(self):
        
        # Run the spacing values method to generate dictionary to translate dG spacing to number of spaces b/t CDS and 3' end of RBS
        self.AssignSpacingValues()
        
        # Setup a general codon translation dictionary for AA name to three-nucleotide code
        self.setup_codon_dict()
        
        #Can enter the specific organism's codon frequency table here, c1 = three letter codon, c2 = val (0 to 100)
        self.make_cod_freq_dict(self.codon_file)
        
        # Setup the basic variables of the preseq and CDS
        self.mRNA_Setup(self.Preseq_input,self.CDS_input)
        
        # Run the RBS calc on the original sequence first to get starting point
        self.RunCalc(self.seqname_input,self.start_mRNA)
        
        self.first_run_setup(self.Nuc_change_input)
        
        init_RBS_seq = self.Preseq_input[self.curr_mRNA['RBS_5p']:(self.curr_mRNA['RBS_5p'] + self.curr_mRNA['RBS_length'])]
        
        # Update the working dictionary so that iteration can begin
        self.Update_working_dict(init_RBS_seq, self.CDS_input)
        
        self.Create_output_csv(self.seqname_input)
        
        self.WriteOutputData(self.seqname_input,"Curr")
        
        # See if the original sequence satisfies the target energies first
        self.quick_check(self.dG_mRNA_rRNA_list[self.index_call],self.dG_mRNA_list[self.index_call],self.curr_mRNA['RBS_seq'])
        
        # Make the current iteration the best one to start with, necessary because BEST is what gets iterated on
        self.Make_best_dict()
        
        # Setup a counter
        iterator = 0
        
        # Create list of potential nucleotides to change that favor RBS and nearby nucleotides
        Nucleotide_List = []
        
        if self.Nuc_change_input == "All":
            for i in range(1,(int(self.MaxIter_input/35) + 1)):
                for f in range(0,34):
                    Nucleotide_List.append(f)

                    
        else:
            for i in range(1,(int(self.MaxIter_input/(self.best_mRNA['RBS_length']) + 1))):
                for z in range(0,self.best_mRNA['RBS_length']):
                    Nucleotide_List.append(z)
            
        
        # Iterate toward the target dG hybrid / dG mRNA combo and stop either at a certain number of iterations or when within +/- 5% of target
        while ((self.dG_Hyb_stat_determine(self.best_mRNA['dG_Hyb'])) or (self.dG_mRNA_stat_determine(self.best_mRNA['dG_mRNA']))) and (iterator < self.MaxIter_input):
            
            # select random number from nucleotide list to mutate
            nt_to_mutate = random.choice(Nucleotide_List)
            
            # mutate with other nucleotide and save new mRNA
            NucGroup = ['a','t','c','g']
            r = random.choice(NucGroup)
            
            z = self.best_mRNA['RBS_seq']
            
            # break the RBS seq into a list of single characters
            f = list(z)
            
            # pull out a nucleotide from the RBS Seq       
            t = f[nt_to_mutate].lower()
            
           # If the replacement nucleotide is the same as the selected one, try another replacement
            while t == r:
                r = random.choice(NucGroup)
            
            # replace the selected nucleotide with the replacement one and rebuild the string of RBS seq
            f[nt_to_mutate] = r
            m = ''.join(f)
            
            # If the rebuilt RBS seq is NOT in the list bad RBS already tried, then test it
            if m not in self.Bad_RBS_list_hyb.keys():
                
                # Generate a new full sequence for running by joining the rebuilt RBS seq to other pieces
                New_Full_Seq = self.best_mRNA['pre_RBS'] + m + self.best_mRNA['post_RBS'] + self.CDS_input
               
                self.RunCalc(self.seqname_input, New_Full_Seq)
                
                self.Update_working_dict(m, self.CDS_input)
                self.quick_check(self.curr_mRNA['dG_Hyb'],self.curr_mRNA['dG_mRNA'],self.curr_mRNA['RBS_seq'])
                
                self.Bad_RBS_list_hyb[(self.curr_mRNA['RBS_seq'])] = self.curr_mRNA['dG_Hyb_dist']
                self.Bad_RBS_list_mRNA[(self.curr_mRNA['RBS_seq'])] = self.curr_mRNA['dG_mRNA_dist']
                
                iterator = iterator + 1
                
                if self.curr_mRNA['tot_vector'] <= self.best_mRNA['tot_vector']:
                    self.curr_mRNA['made_best'] = 1
                    self.Make_best_dict()
                    
                else:
                    var_R = random.random()
                    var_W = math.exp(-self.curr_mRNA['tot_vector']/(self.MaxIter_input/10))
                    var_I = (1 - ((iterator + 1) / self.MaxIter_input))
                    
                    if (var_I * var_W) > var_R:
                        self.curr_mRNA['made_best'] = 2
                        self.Make_best_dict()
                    
                    else:
                        self.curr_mRNA['made_best'] = 0
            
                # write output row now
                self.WriteOutputData(self.seqname_input,"Curr")
            
            else:
                pass
          
        # reset counter
        iterator = 0
        
        Acceptable_RBS = {}
        
        # create dictionary of acceptable RBS as key, dg mRNA vector as value
        for k,v in self.Bad_RBS_list_hyb.iteritems():
            if v <= 1.5:
                Acceptable_RBS[k] = self.Bad_RBS_list_mRNA[k]
        
        if self.curr_mRNA['dG_Hyb_dist'] <= 1.5:
            Acceptable_RBS[self.curr_mRNA['RBS_seq']] = self.curr_mRNA['dG_mRNA_dist']
        
        if not Acceptable_RBS:
            return self.results_file
            
        Best_RBS_val = min(Acceptable_RBS.itervalues())
        
        for k,v in Acceptable_RBS.iteritems():
            if v == Best_RBS_val:
                Best_RBS = k
        
        temp_seq = self.curr_mRNA['pre_RBS'] + Best_RBS + self.curr_mRNA['post_RBS'] + self.CDS_input
        self.RunCalc(self.seqname_input, temp_seq)
        self.Update_working_dict(Best_RBS, self.CDS_input)
        self.curr_mRNA['made_best'] = 3
        self.WriteOutputData(self.seqname_input,"Curr")
        self.Make_best_dict()

        if self.dG_mRNA_stat_determine(self.best_mRNA['dG_mRNA']) == False:
            return self.results_file
            
        # while the mRNA secondary structure is too great and the iterator is less than input
        while (self.dG_mRNA_stat_determine(self.best_mRNA['dG_mRNA'])) and (iterator < self.MaxIter_input):
            
            CDS_codon_list = []
            
            last_two = self.best_mRNA['CDS_seq'][-2:]
            
            for i in range(1,12):
                codon = self.best_mRNA['CDS_seq'][(3*i-3):(3*i)]
                CDS_codon_list.append(codon)
            
            # select random number 1 to 10 for using later to find index for changing codon
            codon_number = random.randint(1,10)
            
            # Determine codon's corresponding amino acid
            codon_to_mutate = CDS_codon_list[codon_number]
            amino_acid = ""
            
            for key, value in self.codon_dict.iteritems():
                for v in value:
                    if v == codon_to_mutate:
                        amino_acid = key
                        
            # make the set of options equal to the variations of codons for the chosen amino acid
            codon_lottery = self.dict_cod_list[amino_acid]
            
            # select random codon to try
            Rand_codon = random.choice(codon_lottery)
            
            while Rand_codon == codon_to_mutate:
                Rand_codon = random.choice(codon_lottery)
                if codon_to_mutate == 'ATG':
                    break
            
            CDS_codon_list[codon_number] = Rand_codon
                
            New_CDS = ''.join(CDS_codon_list) + last_two
            
            # Generate a new full sequence for running by joining the rebuilt CDS seq to other pieces
            temp2 = self.best_mRNA['full_preseq'] + New_CDS 
            
            self.RunCalc(self.seqname_input, temp2)
            self.Update_working_dict(Best_RBS, New_CDS)
            self.Bad_CDS_list[(self.curr_mRNA['CDS_seq'])] = self.curr_mRNA['dG_mRNA_dist']
            
            if self.dG_mRNA_stat_determine(self.curr_mRNA['dG_mRNA']) == False:
                self.curr_mRNA['made_best'] = 5
                self.WriteOutputData(self.seqname_input,"Curr")
                self.Make_best_dict()
                return self.results_file
            
            else:
                
                iterator = iterator + 1
                
                if self.curr_mRNA['dG_mRNA_dist'] <= self.best_mRNA['dG_mRNA_dist']:
                    self.curr_mRNA['made_best'] = 1
                    self.Make_best_dict()
                    
                else:
                    var_R = random.random()
                    var_W = math.exp(-self.curr_mRNA['dG_mRNA_dist']/(self.MaxIter_input/10))
                    var_I = (1 - ((iterator + 1) / self.MaxIter_input))
                    
                    if (var_I * var_W) > var_R:
                        self.curr_mRNA['made_best'] = 2
                        self.Make_best_dict()
                        
                    else:
                        pass
              
                self.WriteOutputData(self.seqname_input,"Curr")
                
        else:
            pass
        
        # Select dG mRNA from all attempted CDS with smallest vector to target value
        Best_CDS_val = min(self.Bad_CDS_list.itervalues())
        
        # If Best CDS is less than the final iteration, replace it as presented sequence
        if Best_CDS_val < self.curr_mRNA['dG_mRNA_dist']:
            for k,v in self.Bad_CDS_list.iteritems():
                if v == Best_CDS_val:
                    self.curr_mRNA['CDS_seq'] = k
            
            temp3 = self.best_mRNA['full_preseq'] + self.curr_mRNA['CDS_seq']
            self.RunCalc(self.seqname_input, temp3)
            self.Update_working_dict(self.best_mRNA['RBS_seq'], self.curr_mRNA['CDS_seq'])
            self.curr_mRNA['made_best'] = 3
            self.WriteOutputData(self.seqname_input,"Curr")
            self.Make_best_dict()
        
        return self.results_file

#------------------------------------------------------------------------------

#
#  def print_output(self, dict_choice, file_to_output):
#        
#        if dict_choice == "Curr":
#            dict_name = self.curr_mRNA
#        else:
#            dict_name = self.best_mRNA
#            
#        print "RBS Sequence: " + str(dict_name['pre_RBS'] + dict_name['RBS_seq'] + dict_name['post_RBS'])
#        print "Full Sequence: " + str(dict_name['pre_RBS'] + dict_name['RBS_seq'] + dict_name['post_RBS'] + dict_name['CDS_seq'])
#        print "Current Expression: " + str(dict_name['curr_expr'])
#        print "Curren dG Hybridization: " + str(dict_name['dG_Hyb'])
#        print "Current dG mRNA: " + str(dict_name['dG_mRNA'])
#        print "\n"
 #       
  #      if file_to_output == True:
   #         print "Report File Name: " + self.results_file

         
        # This code is not currently functional for weighted nucleotide list creation            
        #elif self.Nuc_change_input == "W":
        #    for i in range(1,(int(self.MaxIter_input/35) + 1)):
        #        for f in range(0,34):
        #            Nucleotide_List.append(f)
        #        for g in range((self.best_mRNA['RBS_5p'] - 7),(self.best_mRNA['RBS_5p'] + 17)):
        #            Nucleotide_List.append(g)
        #        for h in range((self.best_mRNA['RBS_5p']),(self.best_mRNA['RBS_5p'] + 10)):
        #            Nucleotide_List.append(h)