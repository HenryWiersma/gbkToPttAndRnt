# GenBank file to Protein Table file (.ptt) and RNA table file (.rnt) parser
This script creates a Protein Table file (.ptt) and a RNA table file (.rnt) from the given GenBank file
Multiple files can be given (or using *) for batch processing

This project is licensed under the terms of the MIT license.

## Install
Download or clone this repository and run the Python script

## Requerements
1. Python 3.2 or higher 
2. BioPython package version 1.43 or higher (with the GenBank parser)

## Usage
The script can be used for a singe file or for a batch of files

### Single file
```
python3 gbkToPttAndRnt.py pathToGenBankFile.gb
```

The script creates a .ptt and a .rnt file with the same name of the input file

### Batch of files
```
python3 gbkToPttAndRnt.py ./*.gb
```

The script creates a .ptt and a .rnt file for every single GenBank file in the batch. 
The name of the files or the same as the orginal GenBank file names

## Multi GenBank file
When multiple GenBank records are present in a single GenBank file (multi-GenBank file), 
a .ptt and a .rnt for every single record will be created.
