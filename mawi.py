#!/usr/bin/env python3
import pandas as pd
import numpy as np
import fnmatch #module for unix style pattern matching
import glob #module is used to retrieve files/pathnames matching a specified pattern
from yattag import Doc, indent
import argparse, hashlib, os, subprocess, sys, time
import gzip
from Bio.Seq import Seq
import re
import regex
parser = argparse.ArgumentParser(prog='mawi.py', formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog="""
        + =================================================================================================================================== +
        |                                                  Match and Write Indexer (MAWI)                     |
        + =================================================================================================================================== +
        """)
parser.add_argument('-s', '--spreadsheet', help='metadata spreadsheet', type=str, required=True)
parser.add_argument('-f', '--folder', help='data folder directory', type=str, required=True)
parser.add_argument('-o', '--output', help='output folder directory', type=str, required=True)
parser.add_argument('-m', '--mismatch', help='character number mismatch allowed in the index', type=int, required=True)
parser.add_argument('-ft', '--fileType', help='file type ( fastq or fasta )', type=int, required=True)
args = parser.parse_args()


def spreadsheet_upload():
    if fnmatch.fnmatch(args.spreadsheet, '*.xls*'):
        metadata_df = pd.read_excel(args.spreadsheet, header=0, sheet_name='Sheet1')
    else:
        print(f'you have used an unsupported spreadsheet: {args.spreadsheet}, please try again')
        exit(1)

    return metadata_df

def find_pattern_with_mismatch(text, pattern, mismatch):
    # Constructing the pattern with allowed mismatches
    pattern_regex = ''.join(['({}){{e<={}}}'.format(pattern,mismatch)])

    # Searching for the pattern in the text with the constructed regex
    for match in regex.findall(pattern_regex, text, overlapped=True):
        return match


def writeToFile(text,sampleName,index_label,index_seq):
    with open(f'{args.output}/{sampleName}.{args.fileType.lower()}', 'a') as file:
        file.write(f'{text}\n')

    with open(f'{args.output}/{sampleName}.{index_label}.txt', 'a') as file:
        file.write(f'{index_seq}\n')

def merge_NanoporeFastq(metadata_df):
    chunk_list = []
    for folder in os.listdir(args.folder):
        if fnmatch.fnmatch(folder, f'*.{args.fileType.lower()}*'):
            if args.fileType.lower() == 'fasta':
                header = '>'
            elif args.fileType == 'fastq':
                header = '@'
            else:
                print('Unsupported file type, please choose allowed value (fasta or fastq)')
            with open(f'{args.folder}/{folder}', 'r') as f:
                 file_content = f.read()
                 results = file_content.split(f'\n{header}')
                 for chunk in results:
                     chunk_list.append(f"{header}{chunk.lstrip(header)}")
    index1_list = []
    index2_list = []
    index1Rev_list= []
    index2Rev_list= []
    chunkList2 = []
    print(len(chunk_list))
    for item in chunk_list:
        for index, row in metadata_df.iterrows():
            index2= Seq(row['index2'])
            index2_reversed = index2.reverse_complement()
            index1 = Seq(row['index'])
            index1_reversed = index1.reverse_complement()
            matches_index1 = find_pattern_with_mismatch(item, str(index1), args.mismatch)
            matches_index2 = find_pattern_with_mismatch(item, str(index2), args.mismatch)
            matches_index1_rev = find_pattern_with_mismatch(item, str(index1_reversed), args.mismatch)
            matches_index2_rev = find_pattern_with_mismatch(item, str(index2_reversed), args.mismatch)
            if f'{item}\n' not in chunkList2:
                if matches_index1:
                    writeToFile(item,row["Sample Name"], 'index1Forward',matches_index1)
                    chunkList2.append(f'{item}\n')
                    index1_list.append(matches_index1)
                if matches_index2 and not matches_index1:
                    writeToFile(item, row["Sample Name"],'index2Forward',matches_index2)
                    chunkList2.append(f'{item}\n')
                    index2_list.append(matches_index2)
                if matches_index1_rev and not (matches_index1 or matches_index2):
                    writeToFile(item, row["Sample Name"],'index1Reverse',matches_index1_rev)
                    chunkList2.append(f'{item}\n')
                    index1Rev_list.append(matches_index1_rev)
                if matches_index2_rev and not (matches_index1 or matches_index2 or matches_index1_rev):
                    writeToFile(item, row["Sample Name"],'index2Reverse',matches_index2_rev)
                    chunkList2.append(f'{item}\n')
                    index2Rev_list.append(matches_index2_rev)



    print(len(chunkList2))


#########################
#      MAIN              #
#########################


if __name__ == '__main__':

    metadata = spreadsheet_upload()
    merge_NanoporeFastq(metadata)
