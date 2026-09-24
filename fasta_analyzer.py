#FASTA analyzer
import pandas as pd
import matplotlib.pyplot as plt

def read_fasta(fasta):
    sequences = {}
    names = []

    with open(fasta, "r") as file:
        for line in file:
            line = "".join(line.split())

            if line.startswith(">"):
                name = line[1:]
                names.append(name)
                sequences[name] = ""
            else:
                sequences[name] += line

    return names, sequences

names, sequences = read_fasta("sequences.fasta")

#Transform the sequences into a DataFrame
def calculate_length(sequence):
    return len(sequence)
def calculate_gc(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    return round((gc_count / len(sequence)) * 100, 2)

results=[]

for name, sequence in sequences.items():
    length = calculate_length(sequence)
    gc = calculate_gc(sequence)
    results.append({"Gene": name, "Length": length, "GC Content": gc})

df = pd.DataFrame(results)
print(df)

#Operations on the DataFrame
print("Average Length:", round(df["Length"].mean(), 2))
print("Average GC Content:", round(df["GC Content"].mean(), 2))

sorted_df = df.sort_values(by="GC Content", ascending=False)
print(sorted_df)

high_gc = df[df["GC Content"] > 60]
print("Genes with high GC Content:")
print(high_gc)

#Save the DataFrame to a CSV file
df.to_csv("fasta_analysis_results.csv", index=False)

#Plotting the results of GC content
plt.bar(df["Gene"], df["GC Content"])
plt.xlabel("Gene")
plt.ylabel("GC Content (%)")
plt.title("GC Content of Genes")

plt.savefig("gc_content_plot.png")
plt.show()