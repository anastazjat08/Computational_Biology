from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


"""
  The function dynamically generates a set of k-mers.  

  Input:  
    - s – nucleotide sequence  
    - k – length of k-mers  
    - prev_k – previous k value  
      - if 'None', the function creates the initial set, which will be modified in subsequent iterations  
    - all_seq_kmers – dictionary where the key is a k-mer, and the value is a list of sequences in which that k-mer appears  
    
  Output:
    - k_mers – set of new k-mers  
    - all_seq_kmers – modified dictionary
"""
def new_kmers(s, k, prev_k, prev_kmers, all_seq_kmers):

    if k > len(s) or k <= 0:
        return []  # return an empty list for invalid values of k

    k_mers = set()  # set of tuples (k-mer, end position of the k-mer in the sequence)

    # Creating k-mers for the first time
    if prev_k == None:
        for i in range(len(s) - k + 1):
            k_mers.add((s[i:i + k], i+k))


    # Dynamic editing of the previous k-mer set
    else:

        r = k - prev_k

        if r == 0: return prev_kmers
        # extend the previous k-mers
        elif r > 0:
            for prev_mer, end_id in prev_kmers:
                new_mer = prev_mer + s[end_id:end_id+r]
                k_mers.add((new_mer, end_id+r))
        # shorten the previous k-mers
        else:
            last_id = None
            for prev_mer, end_id in prev_kmers:
                new_mer = prev_mer[:len(prev_mer)+r]
                k_mers.add((new_mer, end_id+r))
                last_id = end_id+r

            # Complete the set with k-mers from the end of the sequence
            for i in range(last_id, len(s)-k+1):
                k_mers.add((s[i:i + k], i + k))

    # Complete the all_seq_kmers
    for krotka in k_mers:
        kmer = krotka[0]
        if kmer not in all_seq_kmers.keys():
            all_seq_kmers[kmer] = [s]
        else:
            if s not in all_seq_kmers[kmer]:
                all_seq_kmers[kmer].append(s)


    return (k_mers, all_seq_kmers)

"""
    The function determines the minimum k.
    
    Input:
      fasta file
    
    Output:
      Value of k

"""
def find_k(fasta_file):

  sequences = [] # list of all sequences in the file

# l - length of the shortest sequence in the given set
# determines the maximum length of the probe
  l = float('inf')

# saving the sequence
# determining l
  for record in SeqIO.parse(fasta_file, 'fasta'):
    sequences.append(record.seq)
    l = min(l, len(record.seq))

  left = 1
  right = l
  prev_k = None
  min_k = None

  seq_kmers = {} # dictionary, key is the sequence, value is a list (k-mer, its end in the sequence)

# binary search for k
  while left < right:
      all_seq_kmers = {}
      used_seq = set()
      mid = (left + right )//2

      for seq in sequences:
          if seq not in seq_kmers.keys():
              prev_kmers = None
          else:
              prev_kmers = seq_kmers[seq]

          k_mers, all_seq_kmers = new_kmers(seq, mid, prev_k, prev_kmers, all_seq_kmers)
          seq_kmers[seq] = k_mers

      for kmer in all_seq_kmers:
          if len(all_seq_kmers[kmer]) == 1:
              used_seq.add(list(all_seq_kmers[kmer])[0])


      prev_k = mid

      # checking if a k-mer has been found for all sequences
      if len(used_seq) == len(sequences):

          min_k = mid
          right = mid

      else:
          left = mid + 1

  return(min_k)
