import pandas as pd
import sys

# Ensure a filename is provided when running the script
if len(sys.argv) < 2:
    print("Usage: python animal_analysis.py <filename.csv>")
    sys.exit(1)

# Get the filename from command-line arguments
filename = sys.argv[1]

# Load the CSV file
try:
    # df = pd.read_csv(filename, header=None, names=["animalid", "sireid", "damid", "sex", "birthdate"])
    df = pd.read_csv(filename, header=None, names=["animalid", "sireid", "damid", "birthdate"])
except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    sys.exit(1)

# Convert animalid, sireid, and damid to integers for proper comparison
df["animalid"] = df["animalid"].astype(int)
df["sireid"] = df["sireid"].astype(int)
df["damid"] = df["damid"].astype(int)

# Normalize gender values to uppercase (convert 'f' -> 'F' and 'm' -> 'M')
# df["sex"] = df["sex"].str.upper()

# Find all founders (animals with both sireid and damid == 0)
# founders = df[(df["sireid"] == 0) & (df["damid"] == 0)]
# Find all founders (animals with both sireid and damid == 0)
founders = df[(df["sireid"] == -1) & (df["damid"] == -1)]

# Find founders that have offspring
offspring_sires = df[df["sireid"].isin(founders["animalid"])]
offspring_dams = df[df["damid"].isin(founders["animalid"])]

# FIX: Use pd.concat() instead of append()
founders_with_offspring = pd.concat([
    offspring_sires["sireid"], offspring_dams["damid"]
]).dropna().unique()

founders_with_offspring_df = founders[founders["animalid"].isin(founders_with_offspring)]

# Find animals with "01010001" as birthdate that have offspring
# unknown_birthdate_animals = df[df["birthdate"] == "01010001"]
# Find animals with "0000-00-00" as birthdate that have offspring
unknown_birthdate_animals = df[df["birthdate"] == "0000-00-00"]
offspring_from_unknown_birthdate = df[
    df["sireid"].isin(unknown_birthdate_animals["animalid"]) | 
    df["damid"].isin(unknown_birthdate_animals["animalid"])
]
unknown_birthdate_with_offspring_df = unknown_birthdate_animals[
    unknown_birthdate_animals["animalid"].isin(offspring_from_unknown_birthdate["sireid"]) |
    unknown_birthdate_animals["animalid"].isin(offspring_from_unknown_birthdate["damid"])
]

# Print results
# print("\n=== Founders ===")
# print(founders)

# print("\n=== Founders with Offspring ===")
# print(founders_with_offspring_df)

# print("\n=== Animals with Unknown Birthdate and Offspring ===")
# print(unknown_birthdate_with_offspring_df)

with open("output_results.txt", "w") as f:
    f.write("\n=== Founders ===\n")
    f.write(founders.to_string())  # Save DataFrame as text
    
    f.write("\n\n=== Founders with Offspring ===\n")
    f.write(founders_with_offspring_df.to_string())
    
    f.write("\n\n=== Animals with Unknown Birthdate and Offspring ===\n")
    f.write(unknown_birthdate_with_offspring_df.to_string())

print("Results saved to 'output_results.txt'")
