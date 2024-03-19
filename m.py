#! /usr/bin/python3

sample = "dna.example.fasta"
test = "dna2.fasta"

def fasta_to_dict(filename=test, verbous=False):
    """
    This function converts fasta file to a python dictionary.
    """
    reads_dict = {}
    with open(filename) as f:
        current_read = ""
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                read_id = line.split(" ")[0]
                current_read = read_id
                reads_dict[current_read] = ""
                # print(f"current_read => {current_read}")
            else: 
                # print(current_read)
                reads_dict[current_read] += line

    if verbous == True:
        for k, v in reads_dict.items():
            print(k, "  =>  ", v)

    print(f"Total number of sequence reads: {len(reads_dict)}")

    return reads_dict

def read_lengths(read_dictionary):
    """This Function tells about length of reads"""
    length_dict = {}

    for k in read_dictionary.keys():
        # print(k)
        length_dict[k] = len(read_dictionary[k])
    
    return length_dict

def analyze_longest_reads(fasta_dict):
    """Analyze Reads"""
    length = 0
    for k in fasta_dict.keys():
        read_len = len(fasta_dict[k])
        if read_len >= length:
            length = read_len
    print(f"Highest sequence read length is {str(length)}")

    highests = []
    for k in fasta_dict.keys():
        read_len = len(fasta_dict[k])
        if read_len == length:
            highests.append(k)

    return length, highests

def analyze_shortest_reads(fasta_dict):
    """Analyze Reads"""
    length = float("inf")
    for k in fasta_dict.keys():
        read_len = len(fasta_dict[k])
        if read_len <= length:
            length = read_len
    print(f"Shortest sequence read length is {str(length)}")

    shortests = []
    for k in fasta_dict.keys():
        read_len = len(fasta_dict[k])
        if read_len == length:
            shortests.append(k)

    return length, shortests

def kmer_frequencies(fasta_dict, k):
    kmer_dict = {}

    for key in fasta_dict.keys():
        read = fasta_dict[key]
        # print(read)
        for i in range(len(read)-int(k)+1):
            # print(i)
            kmer = read[i:i+k]
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict

def most_frequent_kmer(freq):
    length = 0
    for f in freq.keys():
        kmer_f = freq[f]
        if kmer_f >= length:
            length = kmer_f

    most_freqs = []
    for f in freq.keys():
        kmer_f = freq[f]
        if kmer_f == length:
            most_freqs.append({f: freq[f]})

    return length, len(most_freqs), most_freqs
    
# What is the starting position of the longest ORF in reading frame 3 in any of the sequences? The position should indicate the character number where the ORF begins. For instance, the following ORF:
# def analyze_orf(fasta_dict, pos):
#     """
#     This function analyzes ORFs in a fasta file.
#     pos input is actual not translated to python. 
#     Means possible values are like 1, 2, 3 Not 0, 1, 2
#     """
#     if pos > 3 or pos < -3:
#         print("Not a very common frame, edit code to make it work!")
#         return

#     all_reads_orf = {}
#     for seq in fasta_dict.keys():
#         read = fasta_dict[seq]
#         reads_orf = []
#         frames = read[int(pos)-1:]

#         for i in range(0, len(frames), 3):
#             codon = frames[i:i+3]
#             if codon == "ATG":
#                 for j in range(i, len(frames), 3):
#                     ncodon = frames[j:j+3]
#                     if ncodon == "TAA" or ncodon == "TAG" or ncodon == "TGA":
#                         # print(f"ORF found => start {i+1} end {j+3} length {j+3-i}") 
#                         orf = frames[i:j+3]          
#                         # print(orf)  
#                         length_of_orf = j+3-i
#                         item_to_add = {"start": i, "length": length_of_orf}
#                         reads_orf.append(item_to_add)
#         all_reads_orf[seq] = reads_orf

#     return all_reads_orf

def analyze_orf(fasta_dict, pos):
    """
    This function analyzes ORFs in a fasta file.
    pos input is actual not translated to python. 
    Means possible values are like 1, 2, 3 Not 0, 1, 2
    """
    if pos not in [1, 2, 3]:
        print("Invalid frame position. Please provide a value of 1, 2, or 3.")
        return {}

    all_reads_orf = {}
    for seq_id, seq in fasta_dict.items():
        reads_orf = []
        frames = seq[pos - 1:]

        i = 0
        while i < len(frames):
            if frames[i:i+3] == "ATG":
                j = i + 3
                while j < len(frames):
                    codon = frames[j:j+3]
                    if codon in ["TAA", "TAG", "TGA"]:
                        orf_length = j - i + 3
                        reads_orf.append({"start": i + 1, "length": orf_length})
                        break
                    j += 3
                i = j
            else:
                i += 3

        all_reads_orf[seq_id] = reads_orf

    return all_reads_orf

reads_in_py = fasta_to_dict()
# reads_len = read_lengths(reads_in_py)
# longests = analyze_longest_reads(reads_in_py)
# shortests = analyze_shortest_reads(reads_in_py)

# freq_12mer = kmer_frequencies(reads_in_py, 12)
# highest_12mers = most_frequent_kmer(freq_12mer)

# freq_6mer = kmer_frequencies(reads_in_py, 6)
# highest_6mers = most_frequent_kmer(freq_6mer)

# freq_7mer = kmer_frequencies(reads_in_py, 7)
# highest_7mers = most_frequent_kmer(freq_7mer)

# orf3 = analyze_orf(reads_in_py, 3)
orf2 = analyze_orf(reads_in_py, 3)

# print(longests)
# print(shortests)
# print(freq_12mer)
# print(highest_12mers)
# print(highest_6mers)
# print(highest_7mers)
# print(orf3)

def orf2_():
    orf2_length = 0
    for read in orf2.keys():
        print(read)
        print(" =>  ")
        data = orf2[read]
        print()
        for i in data:
            print(i)
            # print(i["length"])
            # if i["length"] >= orf2_length:
            #     orf2_length = i["length"]
            if i["start"] >= orf2_length:
                orf2_length = i["start"]
    return orf2_length
        

print(len(orf2))
print(orf2_())