# Match and Write Indexer - MAWI
MAWI script provide the ability to match indices to sequences ( FASTA or FASTQ) and then writing (or sorting) these matches into different files according to a spreadsheet.

## Getting Started

#### Spreadsheet: an excel spreadsheet with at least three columns titled `Sample Name`, `Index` and `Index2`.
-- `Sample Name: will contain the resulted file name`
-- `Index: will contain the index1 sequence`
-- `Index2: will contain the index2 sequence`

#### Command to run:
` python3 maiw.py -s <path/to/spreadsheet.xlsx> -f <path/to/sequencefilesDirectory> -o <path/to/outputDirectory> -m <number of mismatch allowed> -ft <file type>`

##### Parameters:
`'-s', '--spreadsheet', 'metadata spreadsheet'`
`'-f', '--folder',  'data folder directory'`
`'-o', '--output', help='output folder directory'`
`'-m', '--mismatch', help='character number mismatch allowed in the index'`
`'-ft', '--fileType', help='file type ( fastq or fasta )'`
