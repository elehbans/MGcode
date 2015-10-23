from __future__ import division
from collections import defaultdict
import sys, math, random, string, csv


class CodOpt:
    
    def GetInputData(self):
        
        #filename = raw_input("what is the name of your file")
        filename = "GFP_only.txt"
        with open (filename, "r") as myfile:
            data = myfile.read()
            self.cdsseq = data.upper()
        
        
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
                        ('Arginine','AGG'),
                        ('Glycine','GGA'),
                        ('Glycine','GGC'),
                        ('Glycine','GGG'),
                        ('Glycine','GGT')]
        
        self.codon_dict = defaultdict(list)
        
        for k,v in translation_list:
            self.codon_dict[k].append(v)
     
     
    def make_cod_freq_dict(self,codon_filename):
        self.Cod_freq_dict = {}
        self.dict_cod_list = defaultdict(list)
        
        reader = csv.reader(open(codon_filename, 'r'))
        
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

#-----------------------------------------------------------------------------

RunNow = CodOpt()

RunNow.GetInputData()

RunNow.setup_codon_dict()

cod_file = "Chlamy_codons.csv"

RunNow.make_cod_freq_dict(cod_file)

num_codons = len(RunNow.cdsseq)/3
if num_codons - int(num_codons) != 0:
    print "Excess nucleotides. fix it"
    sys.exit()
num_codons = int(num_codons)

CDS_codon_list = []   
for i in range(1,(num_codons+1)):
    codon = RunNow.cdsseq[(3*i-3):(3*i)]
    CDS_codon_list.append(codon)


for i in range(0,(num_codons-1)):
    curr_codon = CDS_codon_list[i]
    
    amino_acid = ""
    
    for key, value in RunNow.codon_dict.iteritems():
        for v in value:
            if v == curr_codon:
                amino_acid = key
    #print i
    #print curr_codon
    #print amino_acid

    # make the set of options equal to the variations of codons for the chosen amino acid
    codon_lottery = RunNow.dict_cod_list[amino_acid]
    
    # select random codon to try
    Rand_codon = random.choice(codon_lottery)

    while Rand_codon == curr_codon:
        Rand_codon = random.choice(codon_lottery)
        if curr_codon == 'ATG' or 'TGG':
            break
    
    CDS_codon_list[i] = Rand_codon
        
    RunNow.New_CDS = ''.join(CDS_codon_list)
    
print "This is your new CDS \n"
print RunNow.New_CDS