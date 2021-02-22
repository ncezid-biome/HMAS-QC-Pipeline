import regex
import argparse
from Bio.Seq import Seq

def check_primer_fasta():
    '''
    python3 check_primer.py -p 191018.oligos -g sal.final.groups -f sal.final.FASTA -o result
    '''

    parser = argparse.ArgumentParser(description = 'check if primer sequences exist in reads (fasta)')
    parser.add_argument('-p', '--primer', metavar='', required=True, help='Specify primer file') # oligos file
    parser.add_argument('-g', '--group', metavar='', required=True, help='Specify group file')
    parser.add_argument('-f', '--fasta', metavar='', required=True, help='Specify fasta reads file')
    parser.add_argument('-o', '--out', metavar='', required=True, help='Specify output file')
    args = parser.parse_args()

    primer_file = args.primer
    fasta_file = args.fasta
    group_file = args.group
    out_put = args.out

    dict_primer = {}
    with open(primer_file) as p:
        for line in p.readlines():
            if 'primer' in line:
                fp = line.split()[1]
                rp = line.split()[2]
                primer_id = line.split()[3]
                dict_primer[primer_id] = (fp,rp)

    dict_seq = {}
    line_counter = 0
    with open(fasta_file) as f:
        for line in f.readlines():
            if line_counter % 2 == 0:
                # seq_id = line.split()[0].replace(':', '_')[1:]  # to get the part of seq_id we need
                seq_id = line.split()[0][1:]
            else:
                dict_seq[seq_id] = line
            line_counter += 1

    dict_group = {}
    with open(group_file) as g:
        counter = 0
        for line in g.readlines():
            seq_id = line.split()[0]
            primer_id = line.split()[1].split('.')[1]
            dict_group[seq_id] = primer_id

            if counter %100 == 0:
                print (f'{seq_id} , {primer_id}')
                counter += 1


    result = []
    counter = 0
    for seq_id in dict_seq:

        seq = dict_seq[seq_id]
        primer_id = dict_group[seq_id]
        fp = dict_primer[primer_id][0]
        r_rp = str(Seq(dict_primer[primer_id][1]).reverse_complement())

        if counter %500000 == 0:
            print (f'fp: {fp} , rp: {dict_primer[primer_id][1]}, r_rp: {r_rp}')
            counter += 1

        f_result = match_primer(fp, seq)
        r_result = match_primer(r_rp, seq)
        if f_result:
            print (f'---Found forward primer {fp} in seq: {seq_id}\n'
                   f'{seq}\n')
            result.append(f'---Found forward primer {fp} in seq: {seq_id}\n'
                          f'{seq}\n')
            # index = seq.index(fp)
            # result.append(' '*index + fp + '\n')
            # print (' '*index + fp + '\n')
            result.append(f'fp: {fp} , matched: {f_result.group()}\n')
            print(f'fp: {fp} , matched: {f_result.group()}\n')
        if r_result:
            print (f'---Found reverse_complement primer {r_rp} in seq: {seq_id}\n'
                   f'{seq}\n')
            result.append(f'---Found reverse_complement primer {r_rp} in seq: {seq_id}\n'
                          f'{seq}\n')
            # index = seq.index(r_rp)
            # result.append(' '*index + r_rp + '\n')
            # print (' '*index + r_rp + '\n')
            result.append(f'r_rp: {r_rp} , matched: {r_result.group()}\n')
            print(f'r_rp: {r_rp} , matched: {r_result.group()}\n')


    with open (out_put, 'w') as f:
        f.write(''.join(result))

def match_primer(primer,seq):
    '''
    this method checks the primer against the seq, trying to find any matches
    use regex fuzzy match {e<=2} allowing 2 mis-match

    param: primer (str)
    param: seq (str)

    return: the first found matches or None if not found
    '''
    primer_list = IUPAC_translation(primer)
    for p in primer_list:
        found = regex.search(f"({p}){{e<=3}}", seq)
        if found:
            return found


def IUPAC_translation(seq):
    '''
    this method translate the passed-in seq into a list of all possible IUPAC coded seq
    param: seq (str)
    return: the result list
    '''
    IUPAC_dict = {'Y':'CT','R':'AG','W':'AT','S':'CG','K':'TG','M':'AC',
                  'D':'AGT','V':'ACG','H':'ACT','B':'CGT'}
    result_list = [seq]
    for key in IUPAC_dict:
        temp_list = []
        to_remove_list = []
        for s in result_list:
            if key in s:
                to_remove_list.append(s)
                for char in IUPAC_dict[key]:
                    temp_list.append(s.replace(key,char))
        result_list.extend(temp_list)
        result_list = list(set(result_list) - set(to_remove_list))

    return result_list


if __name__ == "__main__":
    check_primer_fasta()

