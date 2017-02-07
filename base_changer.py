
import math as math
import pandas as pd

#mutations = {'E585K', 'D44A'}

# Sequence goes here (pol6 wild-type)
my_seq = "ATGCATCACCATCATCATCACCACCACAGCGGCGGTTCCGACAAACACACGCAGTACGTCAAAGAGCATAGCTTCAATTATGACGAGTATAAGAAAGCGAATTTCGACAAGATCGAGTGCCTGATCTTTGACACCGAGAGCTGCACGAATTATGAGAACGATAATACCGGTGCACGTGTTTACGGTTGGGGTCTTGGCGTCACCCGCAACCACAATATGATCTACGGCCAAAATCTGAATCAGTTTTGGGAAGTATGCCAGAACATTTTCAATGATTGGTATCACGACAACAAACATACCATTAAGATTACCAAGACCAAGAAAGGCTTCCCGAAACGTAAGTACATTAAGTTTCCGATTGCAGTTCACAATTTGGGCTGGGATGTTGAATTCCTGAAGTATAGCCTGGTGGAGAATGGTTTCAATTACGACAAGGGTCTGCTGAAAACTGTTTTTAGCAAGGGTGCGCCGTACCAAACCGTGACCGATGTTGAGGAACCGAAAACGTTCCATATCGTCCAGAATAACAACATCGTTTATGGTTGTAACGTGTATATGGACAAATTCTTTGAGGTCGAGAACAAAGACGGCTCTACCACCGAGATTGGCCTGTGCTTGGATTTCTTCGATAGCTATAAGATCATCACGTGTGCTGAGAGCCAGTTCCACAATTACGTTCATGATGTGGATCCAATGTTCTACAAAATGGGTGAAGAGTATGATTACGATACTTGGCGTAGCCCGACGCACAAGCAGACCACCCTGGAGCTGCGCTACCAATACAATGATATCTATATGCTGCGTGAAGTCATCGAACAGTTTTACATTGACGGTTTATGTGGCGGCGAGCTGCCGCTGACCGGCATGCGCACCGCTTCCAGCATTGCGTTCAACGTGCTGAAAAAGATGACCTTTGGTGAGGAAAAGACGGAAGAGGGCTACATCAACTATTTTGAATTGGACAAGAAAACCAAATTCGAGTTTCTGCGTAAGCGCATTGAAATGGAATCGTACACCGGTGGCTATACGCACGCAAATCACAAAGCCGTTGGTAAGACTATTAACAAGATCGGTTGCTCTTTGGACATTAACAGCAGCTACCCTTCGCAGATGGCGTACAAGGTCTTTCCGTATGGCAAACCGGTTCGTAAGACCTGGGGTCGTAAACCAAAGACCGAGAAGAACGAAGTTTATCTGATTGAAGTTGGCTTTGACTTCGTGGAGCCGAAACACGAAGAATACGCGCTGGATATCTTTAAGATTGGTGCGGTGAACTCTAAAGCGCTGAGCCCGATCACCGGCGCTGTCAGCGGTCAAGAGTATTTCTGTACGAACATTAAAGACGGCAAAGCAATCCCGGTTTACAAAGAACTGAAGGACACCAAATTGACCACTAACTACAATGTCGTGCTGACCAGCGTGGAGTACGAGTTCTGGATCAAACACTTCAATTTTGGTGTGTTTAAGAAAGACGAGTACGACTGTTTCGAAGTTGACAATCTGGAGTTTACGGGTCTGAAGATTGGTTCCATTCTGTACTACAAGGCAGAGAAAGGCAAGTTTAAACCTTACGTGGATCACTTCACGAAAATGAAAGTGGAGAACAAGAAACTGGGTAATAAGCCGCTGACGAATCAGGCGAAGCTGATTCTGAACGGTGCGTACGGCAAATTCGGCACCAAACAAAACAAAGAAGAGAAAGATTTGATCATGGATAAGAACGGTTTGCTGACCTTCACGGGTAGCGTCACGGAATACGAGGGTAAAGAATTCTATCGTCCGTATGCGAGCTTCGTTACTGCCTATGGTCGCCTGCAACTGTGGAACGCGATTATCTACGCGGTTGGTGTGGAGAATTTTCTGTACTGCGACACCGACAGCATCTATTGTAACCGTGAAGTTAACAGCCTCATTGAGGATATGAACGCCATTGGTGAAACCATCGATAAAACGATTCTGGGTAAATGGGACGTGGAGCATGTCTTTGATAAGTTTAAGGTCCTGGGCCAGAAGAAGTACATGTATCATGATTGCAAAGAAGATAAAACGGACCTGAAGTGTTGCGGTCTGCCGAGCGATGCCCGTAAGATTATCATTGGTCAAGGTTTCGACGAGTTTTATCTGGGCAAAAATGTCGAAGGTAAGAAGCAACGCAAAAAAGTGATCGGCGGTTGCCTGCTGCTGGACACCCTGTTTACGATCAAGAAAATCATGTTCTAA"

# G278D (GGT --> GAT)
primer_seq = 'TTACATTGACgatTTATGTGGCG'
temp_seq = 'AATGTAACTGCCAAATACACCGC'


# P594A
#primer_seq = 'ATTCTATCGTgcgTATGCGAGCT'
#temp_seq = 'TAAGATAGCAGGCATACGCTGCA'

# G278E
#primer_seq = 'TTACATTGACgaaTTATGTGGCG'
#temp_seq = 'AATGTAACTGCCAAATACACCGC'

#G278E_R
#primer_seq = 'CAACACTTCAGGTCCGTTTTATCC'
#temp_seq = 'GTTGTGAAGTCCAGGCAAAATAGG'

#Y692F
#primer_seq = 'CGGTCTGCCGgaaGATGCCCGTA'
#temp_seq = 'GCCAGACGGCATACTACGGGCAT'

# Dictionary of preferred codons for E. coli
gencode_preferred = {
    'A':'gcg', 'R':'cgt', 'N':'aac', 'D':'gat', 'C':'tgc', 'Q':'cag', 'E':'gaa',
    'G':'ggc', 'H':'cac', 'I':'atc', 'L':'ctg', 'K':'aaa', 'M':'atg', 'F':'ttt',
    'P':'ccg', 'S':'agc', 'T':'acc', 'W':'tgg', 'Y':'tat', 'V':'gtg'
}

codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
# Dictionary of nearest neighbor enthalpy (delta_h) and entropy (delta_s) values for Watson-Crick pairs
dna_NN = {
    'init': (0.2, -5.7), 'init_A/T': (2.2, 6.9), 'init_G/C': (0, 0),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),
    'AA/TT': (-7.6, -21.3), 'AT/TA': (-7.2, -20.4), 'TA/AT': (-7.2, -20.4),
    'CA/GT': (-8.5, -22.7), 'AC/TG': (-8.5, -22.7), 'GT/CA': (-8.4, -22.4), 'TG/AC': (-8.4, -22.4), 
    'CT/GA': (-7.8, -21.0), 'TC/AG': (-7.8, -21.0), 'GA/CT': (-8.2, -22.2), 'AG/TC': (-8.2, -22.2),
    'CG/GC': (-10.6, -27.2), 'GC/CG': (-9.8, -24.4), 'GG/CC': (-8.0, -19.0)
}

# Dictionary of mismatch thermodynamics;  'Nearest neighbor mismatch': (delta_h, delta_s)
dna_mismatch = {
    'AG/TT': (1.0, 0.9), 'AT/TG': (-2.5, -8.3), 'CG/GT': (-4.1, -11.7),
    'CT/GG': (-2.8, -8.0), 'GG/CT': (3.3, 10.4), 'GG/TT': (5.8, 16.3),
    'GT/CG': (-4.4, -12.3), 'GT/TG': (4.1, 9.5), 'TG/AT': (-0.1, -1.7),
    'TG/GT': (-1.4, -6.2), 'TT/AG': (-1.3, -5.3), 'AA/TG': (-0.6, -2.3),
    'AG/TA': (-0.7, -2.3), 'CA/GG': (-0.7, -2.3), 'CG/GA': (-4.0, -13.2),
    'GA/CG': (-0.6, -1.0), 'GG/CA': (0.5, 3.2), 'TA/AG': (0.7, 0.7),
    'TG/AA': (3.0, 7.4),
    'AC/TT': (0.7, 0.2), 'AT/TC': (-1.2, -6.2), 'CC/GT': (-0.8, -4.5),
    'CT/GC': (-1.5, -6.1), 'GC/CT': (2.3, 5.4), 'GT/CC': (5.2, 13.5),
    'TC/AT': (1.2, 0.7), 'TT/AC': (1.0, 0.7),
    'AA/TC': (2.3, 4.6), 'AC/TA': (5.3, 14.6), 'CA/GC': (1.9, 3.7),
    'CC/GA': (0.6, -0.6), 'GA/CC': (5.2, 14.2), 'GC/CA': (-0.7, -3.8),
    'TA/AC': (3.4, 8.0), 'TC/AA': (7.6, 20.2),
    'AA/TA': (1.2, 1.7), 'CA/GA': (-0.9, -4.2), 'GA/CA': (-2.9, -9.8),
    'TA/AA': (4.7, 12.9), 'AC/TC': (0.0, -4.4), 'CC/GC': (-1.5, -7.2),
    'GC/CC': (3.6, 8.9), 'TC/AC': (6.1, 16.4), 'AG/TG': (-3.1, -9.5),
    'CG/GG': (-4.9, -15.3), 'GG/CG': (-6.0, -15.8), 'TG/AG': (1.6, 3.6),
    'AT/TT': (-2.7, -10.8), 'CT/GT': (-5.0, -15.8), 'GT/CT': (-2.2, -8.4),
    'TT/AT': (0.2, -1.5),
}

dna_NN_2 = { 
   'init': (0, 0), 'init_A/T': (0, 0), 'init_G/C': (0, 0), 
   'init_oneG/C': (0, -16.8), 'init_allA/T': (0, -20.1), 'init_5T/A': (0, 0), 
   'sym': (0, -1.3), 
   'AA/TT': (-9.1, -24.0), 'AT/TA': (-8.6, -23.9), 'TA/AT': (-6.0, -16.9), 
   'CA/GT': (-5.8, -12.9), 'GT/CA': (-6.5, -17.3), 'CT/GA': (-7.8, -20.8), 
   'GA/CT': (-5.6, -13.5), 'CG/GC': (-11.9, -27.8), 'GC/CG': (-11.1, -26.7), 
   'GG/CC': (-11.0, -26.6)} 
    
# Sugimoto et al. (1996), Nuc Acids Res 24 : 4501-4505 
dna_NN_3 = { 
   'init': (0.6, -9.0), 'init_A/T': (0, 0), 'init_G/C': (0, 0), 
   'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0), 
   'sym': (0, -1.4), 
   'AA/TT': (-8.0, -21.9), 'AT/TA': (-5.6, -15.2), 'TA/AT': (-6.6, -18.4), 
   'CA/GT': (-8.2, -21.0), 'GT/CA': (-9.4, -25.5), 'CT/GA': (-6.6, -16.4), 
   'GA/CT': (-8.8, -23.5), 'CG/GC': (-11.8, -29.0), 'GC/CG': (-10.5, -26.4), 
   'GG/CC': (-10.9, -28.4)} 
    
# Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594 
dna_NN_4 = { 
   'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8), 
   'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0), 
   'sym': (0, -1.4), 
   'AA/TT': (-7.9, -22.2), 'AT/TA': (-7.2, -20.4), 'TA/AT': (-7.2, -21.3), 
   'CA/GT': (-8.5, -22.7), 'GT/CA': (-8.4, -22.4), 'CT/GA': (-7.8, -21.0), 
   'GA/CT': (-8.2, -22.2), 'CG/GC': (-10.6, -27.2), 'GC/CG': (-9.8, -24.4), 
   'GG/CC': (-8.0, -19.9)} 
    
   # SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440 
#   DNA_NN4 = { 
#       'init': (0.2, -5.7), 'init_A/T': (2.2, 6.9), 'init_G/C': (0, 0), 
#       'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0), 
#       'sym': (0, -1.4), 
#       'AA/TT': (-7.6, -21.3), 'AT/TA': (-7.2, -20.4), 'TA/AT': (-7.2, -20.4), 
#       'CA/GT': (-8.5, -22.7), 'GT/CA': (-8.4, -22.4), 'CT/GA': (-7.8, -21.0), 
#       'GA/CT': (-8.2, -22.2), 'CG/GC': (-10.6, -27.2), 'GC/CG': (-9.8, -24.4), 
#       'GG/CC': (-8.0, -19.0)} 


    
###
def complement_calculator(seq):
    seq = seq.upper()
    base_dict = {
        'A':'T', 'T':'A', 'G':'C', 'C':'G'
    }
    forward_complement = ''
    for i in range(0, len(seq)):
        forward_complement += base_dict[seq[i]]
    #print seq
    #print forward_complement
    return forward_complement
###

###
def primer_test(position, mutation):
    print my_seq[position-12:position+12]
    fp = my_seq[position-12:position-1] + gencode_preferred[mutation] + my_seq[position+2:position+12]
    rp = ''#complement_calculator(my_seq[position-28:position-13])#[::-1]
    #fp = my_seq[position-10:position]
    #rp = complement_calculator(my_seq[position-28:position-10])[::-1]
    print fp, rp
    return fp, rp
###

#primer_test(833,'D')

###
def GC_calculator(seq):
    seq = seq.upper()
    GC_count = 0
    for i in range(0, len(seq)-1):
        if seq[i] == 'G' or seq[i] == 'C':
            GC_count += 1
    return float(GC_count)/float(len(seq))
###
###
def library_mutation_input():
    mutation_df = pd.DataFrame(columns = ['Mutation', 'Forward', 'Reverse', 'Tm'])
    #mutation_df.loc[0] = ['Q221G', 'ACTGTGCCAGTG', 'ACGTGTGCAGGTG', 69]
    print mutation_df
    lib_mutations = 'GALEQKFYRWTM'
    lib_input = pd.read_excel('batch_input_test.xlsx', header=0)
    for mutation in lib_input.Residue:
        if mutation == 'EMPTY':
            pass
        if len(mutation) == 4:
            wt_AA = mutation[0]
            AA_num = mutation[1:4]
        elif len(mutation) == 3:
            wt_AA = mutation[0]
            AA_num = mutation[1:3]
        position = (int(AA_num) * 3) - 3
        if wt_AA == codontable[my_seq[position:position+3]]:
            pass
        elif wt_AA != codontable[my_seq[position:position+3]]:
            print 'error'
            break
        for i in range(15):
            for letter in lib_mutations:
                if wt_AA == letter:
                    output_row = ['EMPTY']
                else:
                    output_row = [wt_AA + AA_num + letter, 'primer1', 'primer2', 'Tm']
                    #print output_row
                    mutation_df.loc[i] = output_row
    print mutation_df.head()
    return
###
###
def single_mutation_input():
    single_input = pd.read_excel('batch_input_test.xlsx', header=0)
    for mutation in single_input.Mutation:
        if len(mutation) == 5:
            wt_AA = mutation[0]
            AA_num = mutation[1:4]
            new_AA = mutation[-1]
        elif len(mutation) == 4:
            wt_AA = mutation[0]
            AA_num = mutation[1:3]
            new_AA = mutation[-1]
    return
###

###
def melt_temp_calc(forward):
    # Sets primer and template to global values at start of program
    primer = primer_seq.upper()
    template = temp_seq.upper()
    
    # General initiation value
    delta_h = dna_NN['init'][0]
    delta_s = dna_NN['init'][1]
    d_h = 0
    d_s = 1

    
    # Values for primer:template with no G/C(all A/T) or at least one G/C
    if GC_calculator(primer) == 0:
        delta_h += dna_NN['init_allA/T'][d_h]
        delta_s += dna_NN['init_allA/T'][d_s]
    else:
        delta_h += dna_NN['init_oneG/C'][d_h]
        delta_s += dna_NN['init_oneG/C'][d_s]

    # Values for 5' T (and/or 3' A)
    if primer[0] == 'T':
        delta_h += dna_NN['init_5T/A'][d_h]
        delta_s += dna_NN['init_5T/A'][d_s]
    if primer[-1] == 'A':
        delta_h += dna_NN['init_5T/A'][d_h]
        delta_s += dna_NN['init_5T/A'][d_s]

    # Values for terminal A/T or G/C
    ends = primer[0] + primer[-1]
    AT = ends.count('A') + ends.count('T')
    GC = ends.count('G') + ends.count('C')
    delta_h += dna_NN['init_A/T'][d_h] * AT
    delta_s += dna_NN['init_A/T'][d_s] * AT
    delta_h += dna_NN['init_A/T'][d_h] * GC
    delta_s += dna_NN['init_A/T'][d_s] * GC
    #if (primer[0:1] + "/" + template[0:1]) == 'A/T' or (primer[0:1] + "/" + template[0:1]) == 'T/A':
    #    delta_h += dna_NN['init_A/T'][0]
    #    delta_s += dna_NN['init_A/T'][1]
    #if (primer[0:1] + "/" + template[0:1]) == 'G/C' or (primer[0:1] + "/" + template[0:1]) == 'C/G':
    #    delta_h += dna_NN['init_G/C'][0]
    #    delta_s += dna_NN['init_G/C'][1]

    #for i in range(len(primer)-1):
    #    forward_bp = primer[i:i+2]
    #    reverse_bp = template[i:i+2]
    #    pair = forward_bp + "/" + reverse_bp
    #    if pair in dna_NN:
    #        delta_h += dna_NN[pair][0]
    #        delta_s += dna_NN[pair][1]
    #    if pair in dna_mismatch:
    #        delta_h += dna_mismatch[pair][0]
    #        delta_s += dna_mismatch[pair][1]
    for basenumber in range(len(primer)-1):
        neighbors = primer[basenumber:basenumber+2] + '/' + template[basenumber:basenumber+2]
        if neighbors in dna_mismatch:
            delta_h += dna_mismatch[neighbors][d_h]
            delta_s += dna_mismatch[neighbors][d_s]
        elif neighbors[::-1] in dna_mismatch:
            delta_h += dna_mismatch[neighbors[::-1]][d_h]
            delta_s += dna_mismatch[neighbors[::-1]][d_s]
        elif neighbors in dna_NN:
            delta_h += dna_NN[neighbors][d_h]
            delta_s += dna_NN[neighbors][d_s]
        elif neighbors[::-1] in dna_NN:
            delta_h += dna_NN[neighbors[::-1]][d_h]
            delta_s += dna_NN[neighbors[::-1]][d_s]
        else:
            print 'wtf'

    gas_constant = 1.9872
    oligo_conc = 500e-9
    mg_conc = 0.002
    monovalent_conc = 0.05
    fGC = GC_calculator(primer)
    bp_len = len(primer)
    R_value = math.sqrt(mg_conc)/monovalent_conc


#Lets talk about salt adjustment

#IDT Divalent Salt Adjustment:

#a = 3.92*(0.843 - (0.352*math.sqrt(monovalent_conc)*math.log(monovalent_conc)))
#d = 1.42*(1.279 - 4.03e-3*math.log(monovalent_conc) - 8.03e-3*(math.log(monovalent_conc))**2) 
#g = 8.31*((0.486 - 0.258*math.log(monovalent_conc) + 5.25e-3*(math.log(monovalent_conc)**3))) 

#Tm_salt_adjusted = 1/((1/Tm_1M_Na)+(a-0.911*math.log(mg_conc) + fGC*(6.26+d*math.log(mg_conc))+(1/(2*(bp_len-1))*(-48.2+52.5*math.log(mg_conc)+g*(math.log(mg_conc)**2))))*10e-5) 



#Owczarzy(2008) Divalent Salt Adjustment:

#a = 3.92e-5*(0.843 - (0.352*math.sqrt(monovalent_conc)*math.log(monovalent_conc)))
#b = -9.11e-6
#c = 6.26e-5
#d = 1.42e-5*(1.279 - 4.03e-3*math.log(monovalent_conc) - 8.03e-3*(math.log(monovalent_conc))**2)
#e = -4.82e-4
#f = 5.25e-4
#g = 8.31e-5*((0.486 - 0.258*math.log(monovalent_conc) + 5.25e-3*(math.log(monovalent_conc)**3)))

#Tm_salt_adjusted = 1/((1/Tm_1M_Na) + a + b*math.log(mg_conc) + fGC*(c + d*math.log(mg_conc)) + ((e + f*math.log(mg_conc)) + g*(math.log(mg_conc))**2)/2*(bp_len - 1)) 

    # Equation from SantaLucia 2004: Thermodynamics of DNA structural motifs
    #Tm_1M_Na = (delta_h*1000/(delta_s + gas_constant*math.log(oligo_conc)))-273.15 # For self-complementary duplexes
    Tm_1M_Na = (delta_h*1000/(delta_s + gas_constant*math.log(oligo_conc/4)))-273.15 # For nonself-complementary duplexes (oligo_conc/4)
    print Tm_1M_Na


    # IDT constants
    a = 3.92*(0.843 - (0.352*math.sqrt(monovalent_conc)*math.log(monovalent_conc)))
    d = 1.42*(1.279 - 4.03e-3*math.log(monovalent_conc) - 8.03e-3*(math.log(monovalent_conc))**2)
    g = 8.31*((0.486 - 0.258*math.log(monovalent_conc) + 5.25e-3*(math.log(monovalent_conc)**3)))

    # IDT: https://www.idtdna.com/Calc/Analyzer/Home/Definitions#MeltTemp 
    # From Owczarzy, R. et al., Biochemistry, 47, 5336
    #Tm_salt_adjusted = 1/((1/Tm_1M_Na)+(a-0.911*math.log(mg_conc) + fGC*(6.26+d*math.log(mg_conc))+(1/(2*(bp_len-1))*(-48.2+52.5*math.log(mg_conc)+g*(math.log(mg_conc)**2))))*10e-5) 
    #print Tm_salt_adjusted

    a = 3.92e-5*(0.843 - (0.352*math.sqrt(monovalent_conc)*math.log(monovalent_conc)))
    b = -9.11e-6
    c = 6.26e-5
    d = 1.42e-5*(1.279 - 4.03e-3*math.log(monovalent_conc) - 8.03e-3*(math.log(monovalent_conc))**2)
    e = -4.82e-4
    f = 5.25e-4
    g = 8.31e-5*((0.486 - 0.258*math.log(monovalent_conc) + 5.25e-3*(math.log(monovalent_conc)**3)))
    
    # Owczarzy,2008, Predicting stability.... The above eq from IDT just simplifies constants b, c, e, and f
    Tm_salt_adjusted = 1/((1/Tm_1M_Na) + a + b*math.log(mg_conc) + fGC*(c + d*math.log(mg_conc)) + ((e + f*math.log(mg_conc)) + g*(math.log(mg_conc))**2)/(2*(bp_len - 1)))
    print Tm_salt_adjusted

    return
###

#library_mutation_input()
melt_temp_calc('ACTG')

