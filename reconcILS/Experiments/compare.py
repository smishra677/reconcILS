import pandas as pd

# Load the CSV files
df1 = pd.read_csv('./bio_result_ZeroCol_ASTRAL_ML_SCO_MIN27_pruned_fourth.csv')
df2 = pd.read_csv('./bio_result_ZeroCol_ASTRAL_ML_SCO_MIN27_intro_1_pruned.csv')

# Select the columns you want to compare
columns_to_compare = ['Duplication', 'NNI', 'Loss']
replicate_column = 'Replicate'  # assuming the name of the replicate column

# Creating a new DataFrame to hold comparisons and differences
comparison_df = pd.DataFrame()

for column in columns_to_compare:
    comparison_df[f'{column}_df1'] = df1[column]
    comparison_df[f'{column}_df2'] = df2[column]
    comparison_df[f'{column}_differences'] = df1[column] != df2[column]

# Add replicate column from both DataFrames for reference
comparison_df['replicate_df1'] = df1[replicate_column]
comparison_df['replicate_df2'] = df2[replicate_column]

# Filter to find rows where there are differences
differences = comparison_df[comparison_df[[f'{col}_differences' for col in columns_to_compare]].any(axis=1)]

print("Differences:")
print(differences[['replicate_df1', 'replicate_df2'] + [f'{col}_df1' for col in columns_to_compare] + [f'{col}_df2' for col in columns_to_compare]])
