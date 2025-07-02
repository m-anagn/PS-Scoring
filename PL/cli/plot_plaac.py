import pandas as pd
import matplotlib.pyplot as plt

# Replace with your actual filename
input_file = "plotdata.txt"

# Load the data
# Assuming tab-delimited with a header row
df = pd.read_csv(input_file, sep='\t', skiprows=9)

# Check the columns available (uncomment the next line to see)
# print(df.columns)

# Plot scores for each residue number
plt.figure(figsize=(12, 5))
plt.plot(df['AANUM'], df['PLAAC'], label='PLAAC', color='blue')
plt.plot(df['AANUM'], df['HMM.PrD-like'], label='HMM PrD-like', color='purple')
plt.plot(df['AANUM'], df['PAPA'], label='PAPA', color='green')

plt.xlabel('Residue Number')
plt.ylabel('Score')
plt.title('PLAAC Prion-Like Domain Prediction Scores')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('plaac_scores.png', dpi=300)
print("Plot saved as plaac_scores.png")

