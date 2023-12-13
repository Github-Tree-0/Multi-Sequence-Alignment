#main.py

from align import multi_sequence_alignment

#Notes
#1. This program was tested in pycharm ide
#2. program include : main.py, _init_.py, align.py, alignment.py, utils.py
#3. You need to install numpy, pandas, matplotlib, scipy
#4. Output Alignment and score is displayed
#5. a pdf (msa.pdf) is also generated

if __name__ == '__main__':
    sequences = ['GARFIELD-THE-LAST-FAT-CAT',
                 'GARFIELD-THE-FAST-CAT',
                 'GARFIELD-THE-VERY-FAST-CAT',
                 'THE-FAT-CAT']
    print(sequences)
    alignment = multi_sequence_alignment(sequences)
    print(alignment)
    print(alignment.score())
    alignment.plot("msa.pdf")

