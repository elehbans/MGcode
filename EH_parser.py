# To Do:
# read through gb file
# when hit CDS, create line in FASTA file with ">" name
# go to full genome and pull -35 of TSS to 35+ TSS based on CDS name
# 


#----------------------------------------------------------------------------

import sys, csv, random, re


#----------------------------------------------------------------------------
class parser():

    def io_files(self,gb,fasta):
        
        self.genbank_filename = gb
        fasta_filename = fasta
        self.opened_fasta = open(fasta_filename, "rU")
        fasta_temp = self.opened_fasta.readlines()
        fasta_temp = fasta_temp[1:]
        self.fasta_data = ''.join(line.rstrip() for line in fasta_temp)
        #print self.fasta_data[:1000]
        self.results_filename = self.genbank_filename[:-3] + str(random.randint(0,9)) + ".txt"
        self.output_data = open(self.results_filename, 'a')
    
    def parse(self):
        
        with open(self.genbank_filename, "rU") as self.gb_data:
            print "gb file opened"
            
            for line in self.gb_data:
                if " CDS " in line and "," not in line and "::" not in line and "of" not in line and "in" not in line and ">" not in line:
                    if "complement" in line:
                        working_string = ''.join(str(x) for x in line)
                        new_string_object = re.split('complement', working_string)
                        working_string = new_string_object[1]
                        length = len(working_string)
                        working_string = working_string[((length/2) + 1):-2]
                        start_site = int(working_string)
                        sense_seq = self.fasta_data[(start_site - 35):(start_site + 35)]
                        sense_seq_list = list(sense_seq)
                        reversed_sense_seq = sense_seq_list[::-1]
                        reverse_comp_list = []
                        for i in reversed_sense_seq:
                            if i == 'A':
                                reverse_comp_list.append('t')
                            if i == 'T':
                                reverse_comp_list.append('a')
                            if i == 'G':
                                reverse_comp_list.append('c')
                            if i == 'C':
                                reverse_comp_list.append('g')
                        reverse_comp = ''.join(reverse_comp_list)
                        print reverse_comp    
                        self.output_data.write(">" + str(start_site) + "_complement" + "\n")
                        self.output_data.write(reverse_comp + "\n")
                    else:
                        working_string = ''.join(str(x) for x in line)
                        new_string_object = re.search('(?<=..)\d+', working_string)
                        start_site = int(new_string_object.group(0)) - 1
                        print "this is " + str(start_site)
                        
                        sense_seq = self.fasta_data[(start_site - 35):(start_site + 35)]
                        print sense_seq
                        self.output_data.write(">" + str(start_site) + "\n")
                        self.output_data.write(sense_seq + "\n") 
                else:
                    pass
#----------------------------------------------------------------------------------
RunParser = parser()

RunParser.io_files("Anabaena.gb","Anabaena.fasta")

print "files opened"

RunParser.parse()

RunParser.gb_data.close()
RunParser.opened_fasta.close()
RunParser.output_data.close()

print "File " + RunParser.results_filename + " is done"