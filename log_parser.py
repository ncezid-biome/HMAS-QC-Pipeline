#!/usr/bin/env python
import re


counter = 0
def replace(matched):
    """
    This helper method substitutes all but the first matched occurrence with a literal string

    Parameters
    ----------
    A regex match object

    Returns
    -------
    either the match object itself or literal string mothur >

    """
    global counter
    if counter == 0:
        counter += 1
        return matched.group(0)
    else:
        counter += 1
        return 'mothur >'

def replace_no_error(matched):
    """
    This helper method substitutes the matched occurrence with a literal string,
    if there is no ERROR in it

    Parameters
    ----------
    A regex match object

    Returns
    -------
    either the match object itself or literal string mothur >
    """
    if ('ERROR' in matched.group(0)):
        return matched.group(0)
    else:
        return 'mothur >'


def parse(file_to_parse):
    """
    This function parses a given file (MOTHUR LOG FILE), makes substitutions based on certain patterns,
    and overwrites the MOTHUR LOG FILE.

    Parameters
    ----------
    the name of the file to be parsed

    Returns
    -------
    None

    """
    with open(file_to_parse, 'r') as f:
        file = f.read()

    #match(non-greedy) any contents between get.current() and mothur > quit()
    pattern_1 = r'mothur > get.current\(\).*?mothur > quit\(\)'

    #match any contents between Linux version (or Windows version) and mothur >'
    pattern_2 = r'(Linux|Windows) version.*?mothur >'

    #match any contents between mothur > set.logfile and 3 or more new lines followed by mothur >'
    pattern_3 = r'mothur > set.logfile.*?[\n]{3,}mothur >'

    file = re.sub(pattern_1, '', file, flags=re.DOTALL|re.I)
    #keep the first occurrence of Linux version etc. information
    file = re.sub(pattern_2,replace,file, flags=re.DOTALL|re.I)
    #replace redundancy if there is no error information in it
    file = re.sub(pattern_3,replace_no_error,file, flags=re.DOTALL|re.I)

    with open(file_to_parse, 'w') as f:
        f.write(file)

if __name__ == "__main__":
    print("This module is called by pipeline.py.  Please run pipeline.py --help for more information")