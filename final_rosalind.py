#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from typing import List


# ## Task 1: Evolution as a Sequence of Mistakes
# 
# A mutation is simply an error that occurs during the creation or replication of nucleic acids, particularly DNA. Since nucleic acids are crucial for cellular functions, mutations often have a ripple effect throughout the cell. Although mutations are errors, a rare mutation can provide a beneficial trait to the cell. In fact, the long-term effects of evolution are the cumulative result of advantageous microscopic mutations over generations.
# 
# The most basic and common type of nucleic acid mutation is a point mutation, which involves the replacement of one base with another at a single nucleotide. In DNA, this also requires a change in the complementary base.
# 
# DNA strands from different organisms or species are considered homologous if they share a common ancestor, and counting the base differences in homologous strands gives the minimum number of point mutations that could have occurred over their evolutionary history.
# 
# We aim to minimize the number of point mutations separating two species because, based on the principle of parsimony, evolutionary histories should be explained as simply as possible.

# ### Problem:
# Given two strings $s$ and $t$ of equal length, the Hamming distance between $s$ and $t$, denoted $d_{\mathrm{H}}(s, t)$, is the number of corresponding symbols that differ in $s$ and $t$
# 
# <span style="color: green;">Given</span>: Two DNA strings $s$ and $t$ of equal length (not exceeding 1 kbp).
# 
# <span style="color: green;">Return</span>: The Hamming distance $d_{\mathrm{H}}(s, t)$.

# <span style="color: blue;">Sample dataset</span>:
# 
# GAGCCTACTAACGGGAT \
# CATCGTAATGACGGCCT
# 
# <span style="color: blue;">Sample output</span>: \
# 7

# In[635]:


def HammingDistance(s: str, t: str) -> int:
    #check that s and t are strings, know that they are equal lengths
    s = str(s)
    t = str(t)

    #initialize hammond_dist at 0
    hammond_dist = 0

    #iterate through each bp in 's', if it does not match t increase hammond_dist by 1
    for bp in range(len(s)):
        if s[bp] != t[bp]:
            hammond_dist += 1    
    
    return hammond_dist


# In[637]:


HammingDistance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")


# ## Task 2: Translating RNA into Protein
# 
# Just as nucleic acids are polymers of nucleotides, proteins are chains of smaller molecules called amino acids; 20 amino acids commonly appear in every species. How are proteins created? The genetic code, discovered throughout the course of a number of ingenious experiments in the late 1950s, details the translation of an RNA molecule called messenger RNA (mRNA) into amino acids for protein creation. The apparent difficulty in translation is that somehow 4 RNA bases must be translated into a language of 20 amino acids; in order for every possible amino acid to be created, we must translate 3-nucleobase strings (called codons) into amino acids. Note that there are 4<sup>3</sup>=64
#  possible codons, so that multiple codons may encode the same amino acid. Two special types of codons are the start codon (AUG), which codes for the amino acid methionine always indicates the start of translation, and the three stop codons (UAA, UAG, UGA), which do not code for an amino acid and cause translation to end.

# ### Problem:
# The 20 commonly occurring amino acids are abbreviated by using 20 letters from the English alphabet (all letters except for B, J, O, U, X, and Z). Protein strings are constructed from these 20 symbols. Henceforth, the term genetic string will incorporate protein strings along with DNA strings and RNA strings.
# 
# The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.
# 
# <span style="color: green;">Given</span>:  An RNA string $s$ corresponding to a strand of mRNA (of length at most 10 kbp).
# 
# <span style="color: green;">Return</span>: The protein string encoded by $s$.

# <span style="color: blue;">Sample dataset</span>:
# 
# AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA
# 
# <span style="color: blue;">Sample output</span>: \
# MAMAPRTEINSTRING

# In[109]:


CODON_TABLE = {
    'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S',
    'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C',
    'UGC': 'C', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H',
    'CAA': 'Q', 'CAG': 'Q', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T',
    'ACG': 'T', 'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGU': 'S',
    'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
    'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAU': 'D',
    'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}


def Translation(s: str) -> str:

    #ensure s is a string
    s = str(s).upper()

    #initialize protein_string to hold AA sequence
    protein_seq = ''

    #iterate through 's' in steps of 3 (nucleotides / codon)
    #trim off last 3 nucleotides to account for STOP codon not being translated
    for i in range(0, len(s) - 3, 3):

        #each codon = 3 nucleotides from s
        codon = s[i:i+3]

        #assign aa_id using CODON_TABLE, assign "?" if no matches
        aa_id = CODON_TABLE.get(codon, "?")

        #add aa_id to protein_seq
        protein_seq += aa_id
    
    return protein_seq


# In[111]:


Translation("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")


# ## Task 3: Finding a Motif in DNA
# 
# Discovering the same DNA segment in the genomes of two different organisms strongly suggests that it serves a similar function in both. Such shared DNA segments are referred to as motifs. However, the presence of multiple repeated DNA sequences, known as repeats, complicates the situation. These repeats occur far more frequently than random chance would predict, highlighting the structured nature of genomes

# ### Problem:
# Given two strings $s$ and $t$, $t$ is a substring of $s$. If $t$ is contained as a contiguous collection of symbols in $s$ (as a result, $t$ must be no longer than $s$).
# 
# The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position $i$ of $s$ is denoted by $s[i]$.
# 
# A substring of $s$ can be represented as $s[j:k]$, where $j$ and $k$ represent the starting and ending positions of the substring in $s$; for example, if $s$ = "AUGCUUCAGAAAGGUCUUACG", then $s[2:5]$ = "UGCU".
# 
# The location of a substring $s[j:k]$ is its beginning position $j$; note that $t$ will have multiple locations in $s$ if it occurs more than once as a substring of $s$.
# 
# <span style="color: green;">Given</span>: Two DNA strings $s$ and $t$ (each of length at most 1 kbp).
# 
# <span style="color: green;">Return</span>: All locations of $t$ as a substring of $s$.
# 

# <span style="color: blue;">Sample dataset</span>:
# 
# GATATATGCATATACTT \
# ATAT
# 
# <span style="color: blue;">Sample output</span>: \
# 2 4 10

# In[631]:


def FindingMotif(s: str, t: str) -> List[int]:
    
    #ensure that the inputs are uppercase strings
    s = str(s).upper()
    t = str(t).upper()

    #initialize an empty list to store the locations
    t_location = []

    #iterate through s searching for occurences of 't' and adding to t_location
    for i in range(0, len(s)):
        #use startswith instead of find to get all occurences of t
        #add 1 to account for 0-indexing
        t_location = [i+1 for i in range(len(s)) if s.startswith(t, i)]

    #print t_location cleanly
    return print(*t_location)


# In[633]:


FindingMotif("GATATATGCATATACTT","ATAT")


# ## Task 4: RNA Splicing
# 
# In the nucleus, an enzyme called RNA polymerase (RNAP) starts transcription by breaking the bonds between DNA's complementary bases. It then creates pre-mRNA using one DNA strand as a template, adding complementary RNA bases, with uracil replacing thymine. The other DNA strand, known as the coding strand, is nearly identical to the RNA strand, except for thymine being replaced by uracil. As RNAP progresses, the separated DNA strands quickly rejoin. Pre-mRNA is then processed by removing non-coding segments (introns) and joining coding segments (exons) through a process called splicing, carried out by the spliceosome. Exons together form the gene's coding region, responsible for protein production.

# ### Problem:
# After identifying the exons and introns of an RNA string, we only need to delete the introns and concatenate the exons to form a new string ready for translation.
# 
# The RNA codon table dictates the details regarding the encoding of specific codons into the amino acid alphabet.
# 
# <span style="color: green;">Given</span>:  A DNA string $s$ (of length at most 1 kbp) and a collection of substrings of $s$
#  acting as introns. All strings are given in FASTA format.
# 
# <span style="color: green;">Return</span>: A protein string resulting from transcribing and translating the exons of $s$
# . (Note: Only one solution will exist for the dataset provided)

# <span style="color: blue;">Sample dataset</span>:
# 
# &gt;Pseudo_DNA \
# ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG \
# &gt;Pseudo_intron1 \
# ATCGGTCGAA \
# &gt;Pseudo_intron2\
# ATCGGTCGAGCGTGT
# 
# <span style="color: blue;">Sample output</span>: \
# MVYIADKQHVASREAYGHMFKVCA

# In[627]:


def RNASplicing(s: str, introns: List[str]) -> str:

    #check inputs
    if type(s) != str:
        return print("This is not a sequence. Please input a seqence as a string and retry.")

    if type(introns) != list:
        return print("This is not a list of introns. Please input the introns as a list and retry.")

    for pseud_intron in introns:
        if len(pseud_intron) > len(s):
            return print("One of these introns is longer than the sequence.")
            
    #convert any lowercase inputs
    s = s.upper()
    introns = [intron.upper() for intron in introns]
    
    #initialiaze empty string 'dna' and 'introns' to store seqs without FASTA
    dna = ''
    intron = ''
    
    #get rid of FASTA headers or lines beginning with '>'
    #for DNA
    if s.startswith(">"):
        dna = s.split("\n", 1)[1]
    #for introns
    intron = [intron.split("\n", 1)[1] if intron.startswith(">") else intron for intron in introns]

    #remove intron list from dna
    #initialize empty string to hold exon dna (pre-mRNA processing)
    dna_exon = dna
    #remove each pseud_intron from dna, store in dna_exon
    for pseud_intron in intron:
        dna_exon = dna_exon.replace(pseud_intron, "")

    #transcribe, replace T with U to get mRNA
    mrna = ''
    mrna = dna_exon.replace("T", "U")

    #translate into protein sequence using Translation function
    prot_seq = Translation(mrna)  
    
    return prot_seq


# In[629]:


#stored sample dataset for clarity with function
s = ">Pseudo_DNA\nATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG"
introns = [">Pseudo_intron1\nATCGGTCGAA", ">Pseudo_intron2\nATCGGTCGAGCGTGT"]

RNASplicing(s, introns)


# ## Task 5: Finding a Shared Motif

# ### Problem:
# A common substring of a collection of strings is a substring of every member of the collection. We say that a common substring is a longest common substring if there does not exist a longer common substring. For example, "CG" is a common substring of "ACGTACGT" and "AACCGTATA", but it is not as long as possible; in this case, "CGTA" is a longest common substring of "ACGTACGT" and "AACCGTATA".
# 
# Note that the longest common substring is not necessarily unique; for a simple example, "AA" and "CC" are both longest common substrings of "AACC" and "CCAA".
# 
# <span style="color: green;">Given</span>: A collection of $k$ (k≤100) DNA strings of length at most 1 kbp each in FASTA format.
# 
# <span style="color: green;">Return</span>: A longest common substring of the collection. (If multiple solutions exist, you may return any single solution)

# <span style="color: blue;">Sample dataset</span>:
# 
# &gt;seq_1 \
# GATTACA \
# &gt;seq_2 \
# TAGACCA \
# &gt;seq_3 \
# ATACA
# 
# <span style="color: blue;">Sample output</span>: \
# AC

# In[623]:


def LongestCommonSubstring(k: List[str]) -> str:

    #check inputs
    if type(k) != list:
        return print("Please input the sequences in list format and retry.")

    #convert any lowercase inputs
    fasta_seqs = [seq.upper() for seq in k]

    #remove FASTA headers
    clean_seqs = ''
    clean_seqs = [seq.split("\n", 1)[1] if seq.startswith(">") else seq for seq in fasta_seqs]

    #find the shortest sequence, separate into own variable for iterative clarity later
    short_seq = ''
    for seq in clean_seqs:
        short_seq = min(clean_seqs, key = len)

    #initialize lcs_seq
    lcs_seq = ''
    
    #iterate through each potential substring length (max = short_seq, min = 1)
    for length in range(len(short_seq), 0, -1):
        
        #check each potential substring 
        for start in range(len(short_seq) - length + 1):
            substring = short_seq[start:start + length]
            
            #check if the substring is in all sequences
            if all(substring in seq for seq in clean_seqs):
                #return the first longest common substring found
                return substring
    
    return lcs_seq


# In[625]:


#stored sample dataset for clarity with function
k = [">seq_1\nGATTACA", ">seq_2\nTAGACCA", ">seq_3\nATACA"]

LongestCommonSubstring(k)


# ## Task 6: Finding a Spliced Motif
# 
# In “Finding a Motif in DNA”, we searched for occurrences of a motif as a substring of a larger database genetic string. However, a DNA strand coding for a protein is often interspersed with introns (see “RNA Splicing”), thus we need a way to recognize a motif that has been chopped up into pieces along a chromosome.

# ### Problem:
# A subsequence of a string is a collection of symbols contained in order (though not necessarily contiguously) in the string (e.g., ACG is a subsequence of TATGCTAAGATC). The indices of a subsequence are the positions in the string at which the symbols of the subsequence appear; thus, the indices of ACG in TATGCTAAGATC can be represented by (2, 5, 9).
# 
# As a substring can have multiple locations, a subsequence can have multiple collections of indices, and the same index can be reused in more than one appearance of the subsequence; for example, ACG is a subsequence of AACCGGTT in 8 different ways.
# 
# <span style="color: green;">Given</span>: Two DNA strings $s$ and $t$ (each of length at most 1 kbp) in FASTA format.
# 
# <span style="color: green;">Return</span>: One collection of indices of $s$ in which the symbols of $t$ appear as a subsequence of $s$. If multiple solutions exist, you may return any one.

# <span style="color: blue;">Sample dataset</span>:
# 
# &gt;seq_1 \
# ACGTACGTGACG \
# &gt;seq_2 \
# GTA
# 
# <span style="color: blue;">Sample output</span>: \
# 3 8 10

# In[617]:


def FindingSubsequence(s: str, t: str) -> List[int]:

    #check inputs
    if len(s) < len(t):
        return print("Your query sequence is longer than your sequence. Please input a new query and retry.")

    #remove FASTA headers, store in seq and query_seq respectively (else is just in case not in FASTA)
    if s.startswith(">"):
        seq = s.split("\n", 1)[1]
    else:
        seq = s
    
    if t.startswith(">"):
        query_seq = t.split("\n", 1)[1]
    else:
        query_seq = t

    #initialize empty list to store locations
    subseq = []
    
    #initialize current_position to track main sequence
    current_position = 0
    
    #find each base pair in query_seq in seq
    for bp in query_seq:
        current_position = seq.find(bp, current_position) + 1
        #address if there is no shared base pairs
        if current_position == 0:
            raise ValueError("Subsequence not found.")
        #append current_position to track where query_seq[bp] = seq[bp]
        subseq.append(current_position)
    
    return print(*subseq)


# In[621]:


#stored sample dataset for clarity with function
s = '>seq_1\nACGTACGTGACG'
t = '>seq_2\nGTA'

FindingSubsequence(s, t)


# ## Task 7: Finding a Shared Spliced Motif
# 
# In "Finding a Shared Motif," we explored how to search a database of genetic strings to find the longest common substring, which represented a motif shared by both strings. However, as discussed in "RNA Splicing," coding regions in DNA are interrupted by introns that don't code for proteins.
# 
# Thus, we need to identify shared motifs that are spread across exons, meaning the motifs don't have to be continuous. To represent this, we must use subsequences.

# ### Problem:
# A string $u$ is a common subsequence of strings $s$ and $t$. If the symbols of $u$ appear in order as a subsequence of both $s$ and $t$. For example, "ACTG" is a common subsequence of "AACCTTGG" and "ACACTGTGA".
# 
# Analogously to the definition of longest common substring, $u$ is a longest common subsequence of $s$ and $t$ if there does not exist a longer common subsequence of the two strings. Continuing our above example, "ACCTTG" is a longest common subsequence of "AACCTTGG" and "ACACTGTGA", as is "AACTGG".
# 
# <span style="color: green;">Given</span>: Two DNA strings $s$ and $t$ (each of length at most 1 kbp) in FASTA format.
# 
# <span style="color: green;">Return</span>:  A longest common subsequence of $s$ and $t$. If multiple solutions exist, you may return any one.

# <span style="color: blue;">Sample dataset</span>:
# 
# &gt;seq_1 \
# AACCTTGG \
# &gt;seq_2 \
# ACACTGTGA
# 
# <span style="color: blue;">Sample output</span>: \
# AACTGG

# In[479]:


def LongestCommonSubsequence(s: str, t: str) -> str:
    
    #remove FASTA headers, store in seq and query_seq respectively
    if s.startswith(">"):
        seq1 = s.split("\n", 1)[1]
    else:
        seq1 = s.upper()
    if t.startswith(">"):
        seq2 = t.split("\n", 1)[1]
    else:
        seq2 = t.upper()

    #initialize dynamic programming talbe to store all lengths of common subsequences
    len_seq1, len_seq2 = len(seq1), len(seq2)
    #initialiaze subseq_table with 0s using lengths of sequences
    subseq_table = [[0] * (len_seq2 + 1) for _ in range(len_seq1 + 1)]

    #fill in table, accounting for 0 indexing and extra col/row throughout with +/- 1
    for i in range(1, len_seq1 + 1):
        for j in range(1, len_seq2 + 1):

            #if seq1 and seq2 base pairs match:
            if seq1[i - 1] == seq2[j - 1]: 
                #"mark" corresponding spot in subseq_table
                subseq_table[i][j] = subseq_table[i - 1][j - 1] + 1

            #if seq1 and seq2 base pairs do not match
            else:
                #apply max on values above/left mismatch (this excludes bp from seq1 or seq2 based on score)
                subseq_table[i][j] = max(subseq_table[i - 1][j], subseq_table[i][j - 1])

    #use subseq_table calculations to find the longest common subsequence
    #initialize empty list to store lc_subsequence
    lc_subsequence = []

    #create variables that correspond to base pairs relative to sequence (set to length)
    bp_seq1, bp_seq2 = len_seq1, len_seq2

    #bound the iterations based on size of sequences
    while bp_seq1 > 0 and bp_seq2 > 0:

        #match: if seq1 matches sequence 2 add character to lc_subsequence
        if seq1[bp_seq1 - 1] == seq2[bp_seq2 - 1]:
            lc_subsequence.append(seq1[bp_seq1 - 1])
            #move diagnoally/left to account for indexing in subseq_table
            bp_seq1 -= 1
            bp_seq2 -= 1
            
        #mismatch: if mismatch check value of previous character to see if it is greater than current
        elif subseq_table[bp_seq1 - 1][bp_seq2] > subseq_table[bp_seq1][bp_seq2 - 1]:
            #if greater then go back a to check previous bp in seq1
            bp_seq1 -= 1
        else:
            #if less than then go back a to check previous bp in seq2
            bp_seq2 -= 1

    #have to reverse to get correct order
    return ''.join(reversed(lc_subsequence))


# In[481]:


#stored sample dataset for clarity with function
s = '>seq_1\nAACCTTGG'
t = '>seq_2\nACACTGTGA'

LongestCommonSubsequence(s,t)


# ## Task 8: Two Motifs, One Gene
# 
# In the previous task, we found the longest motif that could have been shared by two genetic strings, allowing for the motif to be split onto multiple exons in the process. As a result, we needed to find a longest common subsequence of the two strings (which extended the problem of finding a longest common substring from “Finding a Shared Motif”).
# 
# In this problem, we consider an inverse problem of sorts in which we are given two different motifs and we wish to interleave them together in such a way as to produce a shortest possible genetic string in which both motifs occur as subsequences.
# 
# Thus, we need to identify shared motifs that are spread across exons, meaning the motifs don't have to be continuous.

# ### Problem:
# A string $s$ is a supersequence of another string $t$ if $s$ contains $t$ as a subsequence. A common supersequence of strings $s$ and $t$ is a string that serves as a supersequence of both $s$ and $t$. For example, "GACCTAGGAACTC" serves as a common supersequence of "ACGTC" and "ATAT". A shortest common supersequence of $s$ and $t$ is a supersequence for which there does not exist a shorter common supersequence. Continuing our example, "ACGTACT" is a shortest common supersequence of "ACGTC" and "ATAT".
# 
# <span style="color: green;">Given</span>: Two DNA strings $s$ and $t$.
# 
# <span style="color: green;">Return</span>: A shortest common supersequence of $s$ and $t$. If multiple solutions exist, you may return any one.

# <span style="color: blue;">Sample dataset</span>:
# 
# ATCTGAT \
# TGCATA
# 
# <span style="color: blue;">Sample output</span>: \
# ATGCATGAT

# In[611]:


def ShortestCommonSupersequence(s: str, t: str) -> str:
    #get longest common subsequence
    lc_subseq = LongestCommonSubsequence(s, t)

    #create indices i,j for seq1, seq2 respectively
    i, j = 0, 0
    #initialize short_cs to hold shortest common supersequence
    short_cs = ''

    #change names for iterative clarity
    seq1 = s.upper()
    seq2 = t.upper()

    #iterate through each base pair in lc_subseq
    for bp in lc_subseq:
        
        #use counter 'i' to stay within range of seq1
        if i < len(seq1):

            #if there is not a match append seq1 base pair, increase counter i
            while seq1[i] != bp:
                short_cs += seq1[i]
                i += 1
            i += 1
            
        #use counter 'j' to stay within range of seq2
        if j < len(seq2):

            #if there is no match append seq2 base pair, increase counter j
            while seq2[j] != bp:
                short_cs += seq2[j]
                j += 1
            j += 1

        #case of a match append bp directly
        short_cs += bp
    
    #account for mismatching lengths, append these directly on
    if i < len(seq1):
        short_cs += seq1[i:]
    if j < len(seq2):
        short_cs += seq2[j:]

    return short_cs


# In[613]:


#stored sample dataset for clarity with function
s = 'ATCTGAT'
t = 'TGCATA'

ShortestCommonSupersequence(s,t)

