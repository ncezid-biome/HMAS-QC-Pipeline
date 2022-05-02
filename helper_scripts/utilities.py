

def create_fasta_dict(fasta):
    '''
    this method reads a fasta file and convert it into a dictionary, with seq_ID being the key and actual sequence
    being the value

    Note: if there are duplicate seq ID, only the first seq_ID/sequence will be saved. 

    Parameters
    ----------
    fasta: String name of the fasta file

    Returns a dictionary

    '''
    seq_dict = {}
    with open(fasta, 'r') as f:
        for ind, row in enumerate(f.readlines(), start=1):
            if ind%2 == 1:
                last_seqID = row.strip().split()[0][1:] #removes the '>'
            else:
                last_seq = row.strip()
                if last_seqID in seq_dict:
                    print (f"warning ! has a duplicate seq ID {last_seqID}")
                else:
                    seq_dict[last_seqID] = last_seq
    return seq_dict
