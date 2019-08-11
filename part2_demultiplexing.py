#!/usr/bin/env python

# Define the problem:
# We have recieved 4 raw fastq output files from an illumina sequencer. It is possible
# for index hopping to take place during the sequencing process. Out goal is to get rid
# of the low coverage reads, the index hopped reads, and any read that has an index that is
# not one of the 24. Storing all of these full records in their respective files. There should be
# 52 files in total. 48 for the forward and reverse index files, 2 for undetermined, 2 for index hopped
# We would only care about the 48 files afterward because those are the reads that pass.

import gzip

COVERAGE_CUTOFF = 31
ALL_INDEXES = [] # I would read all 24 indexes in from the file
R1_RECORD = []
R2_RECORD = []
R3_RECORD = []
R4_RECORD = []


INDEX_DICT_FW = {} # I would store each index as the key and the value as a filepointer to the respective forward file
INDEX_DICT_RV = {} # I would store each index as the key and the value as a filepointer to the respective reverse file


# FUNCTION that converts a single character into a numerical phred score
def convert_phred(letter):                                # method to convert a single character into a phred score
    """Converts a single character into a phred score""" 
    return(ord(letter)-33)                                # ord returns the unicode code for any given unicode character. Then 33 is subracted off because that is the standard in order to avoid the first 33 non-visable characters
# unit test:
#   'A' goes in
#   32 gets returned


# FUNCTION to check if the quality of the given phred score meets my quality cutoff
def coverage_check(phred_score):
    """ calls convert_phred and compares the numerical phred value to my chosen cutoff"""
    # if convert_phred(phred_score) > COVERAGE_CUTOFF
    #   then it passes the cutoff check
    #  **Return True**
    # else
    #   the read does not pass and will be put into the R1 or R2 undetermined file 
    # **Return False**
# unit test:
#   30: fail
#   38: pass


# FUNCTION to calculate the reverse compliment
def reverse_compliment(index):
    """Converts a given string to the reverse compliment string"""
    # convert the given string to the reverse compliment
    # return the converted string
# unit test:
#   TATG goes in 
#   CATA gets returned


# FUNCTION to check if an index is one of the 24 given
def check_index(current_index):
    """checks to see if the current_index is in my list ALL_INDEXES"""
    # I will have a list that contains all 24 of the indexes
    # If current_index is in ALL_INDEXES
    #   current_index passes check
    #  **Return True**
    # else
    #   does not pass and read will be put into the R1 or R2 undetermined file
    #  **Return False**
# unit test:
#   GTAGCGTA    pass
#   GTACCGTA    fail  (not one of the 24 given indexes)

count = 0
with gzip.open('r1test.fastq.gz','rt') as R1, gzip.open('r2test.fastq.gz','rt') as R2, gzip.open('r3test.fastq.gz','rt') as R3, gzip.open('r4test.fastq.gz','rt') as R4:
    while(True):
        count += 1
        # We want to start off by grabbing the first record from each of the files
        #grab current R1 record
        R1_RECORD.append(R1.readline().strip())
        R1_RECORD.append(R1.readline().strip())
        R1_RECORD.append(R1.readline().strip())
        R1_RECORD.append(R1.readline().strip())

        #grab current R2 record
        R2_RECORD.append(R2.readline().strip())
        R2_RECORD.append(R2.readline().strip())
        R2_RECORD.append(R2.readline().strip())
        R2_RECORD.append(R2.readline().strip())

        #grab current R3 record
        R3_RECORD.append(R3.readline().strip())
        R3_RECORD.append(R3.readline().strip())
        R3_RECORD.append(R3.readline().strip())
        R3_RECORD.append(R3.readline().strip())

        #grab current R4 record
        R4_RECORD.append(R4.readline().strip())
        R4_RECORD.append(R4.readline().strip())
        R4_RECORD.append(R4.readline().strip())
        R4_RECORD.append(R4.readline().strip())

        # condition to break out of the loop
        if(R1_RECORD[0] == ""):
            break

        # After we have access to the current record in each file, we then want to test for each of the four possibilities before we can determine if it is dual matched
        # 1. if either of the indexes have Ns, put the whole record in the respective undetermined files
        # 2. if the phred score of any of the index bases is below the coverage cutoff, put the whole record in the respective undetermined files
        # 3. if the forward index or the reverse-compliment of the reverse index do not match the 24 given indexes, put the whole record in the respective undetermined files
        # 4. if both are valid indexes, but are not the reverse compliment of eachother, then the whole record is put in the respective index hopped files
        
        # unreversed_uncomplimented = reverse_compliment(R3_RECORD[1])

        # check 1/3:
        # if(check_index(R2_RECORD[1]) != True or check_index(unreversed_uncomplimented) != True):
        #     One of the indexes had either an N or it just did not match one of the 24 given indexes.
        #     The forward sequence record will be appended to the undetermined_R1.fastq file with the header modified to contain both indices
        #     The reverse sequence record will be appended to the undetermined_R2.fastq file with the header modified to contain both indices
        # else:
        #     Both indexes do not have Ns and match one of the 24 given indexes

        # check 2:
        # for i in range(len(R2_RECORD[1])):
        #     if(coverage_check(R2_RECORD[1][i]) != True):
        #         The forward sequence record will be appended to the undetermined_R1.fastq file with the header modified to contain both indices
        #         The reverse sequence record will be appended to the undetermined_R2.fastq file with the header modified to contain both indices
        #     else:
        #         The current base in this index passes the coverage check
        #     if(coverage_check(R3_RECORD[1][i]) != True):
        #         The forward sequence record will be appended to the undetermined_R1.fastq file with the header modified to contain both indices
        #         The reverse sequence record will be appended to the undetermined_R2.fastq file with the header modified to contain both indices
        #     else:
        #         The current base in this index passes the coverage check

        # check 4:
        # if(R2_RECORD[1] != unreversed_uncomplimented):
        #     this read has been index hopped
        #     The forward sequence record will be appended to the hopped_R1.fastq file with the header modified to contain both indices
        #     The reverse sequence record will be appended to the hopped_R2.fastq file with the header modified to contain both indices
        # else:
        #     There is no index hopping and therefore the two reads are dual matched
        #     The forward sequence record will be appended to the dual matched file that goes with this index (ex. indx1_R1.fastq) with the header modified to contain both indices
        #     The reverse sequence record will be appended to the ual matched file that goes with this index (ex. indx1_R2.fastq) with the header modified to contain both indices




        R1_RECORD = []
        R2_RECORD = []
        R3_RECORD = []
        R4_RECORD = []