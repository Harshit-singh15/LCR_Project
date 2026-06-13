import pandas as pd

# 1. Define the path to your TSV file
# Replace 'your_file.tsv' with the actual path to your file
file_path = r"Fig8\fig8A_disprot_complexity.tsv"

try:
    # 2. Open and read the TSV file
    # sep='\t' tells pandas that the file is tab-separated
    df = pd.read_csv(file_path, sep="\t")

    # 3. Get unique, sorted values from the 'ec_name' column
    # We check if the column exists first to avoid a KeyError
    if "ec_name" in df.columns:
        sorted_unique_names = sorted(df["ec_name"].unique())

        # 4. Print the result
        print(sorted_unique_names)
    else:
        print("Error: The column 'ec_name' was not found in the TSV file.")
        print(f"Available columns are: {list(df.columns)}")

except FileNotFoundError:
    print(f"Error: The file at '{file_path}' could not be found. Please check the path.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")