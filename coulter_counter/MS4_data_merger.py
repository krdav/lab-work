import numpy as np
from pathlib import Path
import pandas as pd
from functools import reduce
import os
import warnings

### These are the rows to select from the .tab files to merge.
### All other rows will be omitted.
rows = ['File name', 'File ID', 'Sample ID', 'Control mode', 'Acquired', 'Counts above threshold',
        'Coincidence corrected', 'From', 'To', 'Number', 'Mean', 'Median', 'S.D.']
rows_numeric = ['Counts above threshold', 'Coincidence corrected', 'From', 'To', 'Number', 'Mean', 'Median', 'S.D.']

def get_row_mask(df, rows, filename):
    '''
    Return a mask for selecting rows in a dataframe.
    From and To are special cases for this dat.
    Only the last From, To set is desired.
    '''
    mask = list()
    for row in rows:
        idxs = df[df[0] == row].index.values

        if len(idxs) == 1:
            mask.append(idxs[0])
        elif row == 'From' or row == 'To':
            mask.append(max(idxs))
        elif len(idxs) > 1:
            raise Exception('Cannot parse filename: {}\nWith multiple rows with the same name: {}'.format(filename, row))
        elif len(idxs) == 0:
            raise Exception('Cannot parse filename: {}\nUnable to find row with the name: {}'.format(filename, row))

    return(mask)


print('This program will gather data files and merge sample information into a single excel file.')
print('For each sample MS4 software should generate a \'.tab\' and a \'.#m4\'. Please organize all these files into one subfolder under the MS4_data_dump folder.\n')

### Take input:
answer = False
while answer is False:
    folder_name = input('What is the name of the folder?\n')
    data_fname = input('What should be the name of the resulting datafile (excluding file extension)?\n')
    answer = input('You typed folder: \'{}\' and data filename: \'{}\' is that correct? Y/N\n'.format(folder_name, data_fname))
    if answer.upper() == 'Y':
        answer = True
    elif answer.upper() == 'N':
        print('Ok, try again.')
        answer = False
    else:
        print('Did not understand answer: {} try again.'.format(answer))
        answer = False
    if os.path.isdir(folder_name) is False:
    	print('No folder named: {}\nBe sure that the specified folder is a subfolder of the MS4_data_dump folder.'.format(folder_name))
    	answer = False
# This filename needs to have extension for Pandas to know which format to write:
data_fname = data_fname + '.xlsx'
# Generate the path for the output file:
folder_path = Path(folder_name)
data_fname_path = Path(data_fname)
data_fname_target = os.path.join(folder_path, data_fname)


print('\nStripping \':\' and space from the beginning and end of row names.')
print('Selecting the following rows:\n{}\nTo merge from all the \'.tab\' files.'.format(', '.join(rows)))

# Locate the files:
tab_files = list(folder_path.glob('*.tab'))
m4_files = list(folder_path.glob('*.#m4'))
tab_file_names = [f.name for f in tab_files]
m4_file_names = [f.name for f in m4_files]
if len(tab_file_names) == 0:
	raise Exception('No .tab files found in folder: {}'.format(folder_name))
elif len(tab_file_names) != len(m4_file_names):
	print('\n****************\nWarning: Found {} .tab file and {} .#m4 files. Usually, these files come in pairs.\n****************\n'.format(len(tab_file_names), len(m4_file_names)))



print('\nThese are the datafiles to merge:\n{}'.format('\n'.join(tab_file_names)))


# Read all the samples into a list of dataframes:
df_list = []
for fnam in tab_files:
    with fnam.open() as fh:
        df = pd.read_csv(fnam.open(), sep='\t', header=None)
    # Strip colon and space from the row names:
    df[0] = [s.strip(': ') for s in df[0].values]
    # Extract the desired rows:
    mask = get_row_mask(df, rows, fnam.name)
    df = df.loc[mask]
    df_list.append(df)


# Now merge all dataframes based on row name:
df_merged = reduce(lambda  left, right: pd.merge(left, right, on=[0], how='outer'), df_list)

# Swap rows so the Sample ID appears as the first entry:
df_idx = np.array(df_merged.index.values)
Sid_idx = df_merged[df_merged[0] == 'Sample ID'].index.values[0]
df_idx[Sid_idx] = 0
df_idx[0] = Sid_idx
# The redo index:
df_merged = df_merged.loc[df_idx].reset_index(drop=True)

# Transpose the dataframe so each sample is one row and the row names are column names:
df_merged = df_merged.T.reset_index(drop=True)
# Move the row names up into column names:
df_merged.columns = df_merged.loc[0]
df_merged = df_merged.drop(0)

# Convert numeric rows to the correct datatype:
df_merged[rows_numeric] = df_merged[rows_numeric].apply(pd.to_numeric)
# Also convert the date to correct datatype:
df_merged['Acquired'] = pd.to_datetime(df_merged['Acquired'])

# Dump to Excel file:
df_merged.to_excel(data_fname_target, header=True, index=False)

