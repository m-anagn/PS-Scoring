import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def extract_data(input_txt, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    proteins = {}

    with open(input_txt, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) >= 4 and parts[0] == 'ResiduePScore:':
                residue = int(parts[1])
                score = float(parts[3])
                protein_name = parts[-1].replace('>', '')

                if protein_name not in proteins:
                    proteins[protein_name] = []
                proteins[protein_name].append([residue, score])

    csv_files = []
    for protein_name, data in proteins.items():
        df = pd.DataFrame(data, columns=['residue', 'score'])
        output_file = os.path.join(output_dir, f"{protein_name}.csv")
        df.to_csv(output_file, index=False)
        csv_files.append((output_file, protein_name))
        print(f"Saved {output_file}")

    return csv_files

def create_bar_plot(input_file, title, output_file, y_ticks=1, x_ticks=20, y_lim=None):
    data = pd.read_csv(input_file)

    plt.figure(figsize=(10, 6), dpi=300)
    plt.bar(data['residue'], data['score'], color='black', linewidth=0)

    plt.title(title, fontsize=14)
    plt.xlabel('Residue', fontsize=12)
    plt.ylabel('Score', fontsize=12)
    plt.axhline(0, color='black', linestyle='dotted', linewidth=1)
    plt.xlim([min(data['residue']), max(data['residue'])])

    if y_lim:
        plt.ylim(y_lim)
        y_min, y_max = y_lim
    else:
        y_min, y_max = min(data['score']), max(data['score'])

    plt.xticks(ticks=range(0, len(data['residue']), x_ticks))
    plt.yticks(ticks=np.arange(y_min, y_max + y_ticks, y_ticks))

    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Plot saved to {output_file}")

def process_data(input_txt, output_dir, y_ticks=1, x_ticks=20, y_lim=None):
    csv_files = extract_data(input_txt, output_dir)
    plot_dir = os.path.join(output_dir, "plots")
    os.makedirs(plot_dir, exist_ok=True)

    for csv_file, protein_name in csv_files:
        plot_file = os.path.join(plot_dir, f"{protein_name}.png")
        create_bar_plot(csv_file, protein_name, plot_file, y_ticks, x_ticks, y_lim)

def main():
    parser = argparse.ArgumentParser(description="Process protein scoring data from TXT and generate plots.")
    parser.add_argument('input_txt', type=str, help="Path to the input TXT file.")
    parser.add_argument('output_dir', type=str, help="Path to save the output CSV and plots.")
    parser.add_argument('--y_ticks', type=float, default=1.0, help="Interval for y-axis ticks.")
    parser.add_argument('--x_ticks', type=int, default=20, help="Interval for x-axis ticks.")
    parser.add_argument('-y_lim', nargs=2, type=float, metavar=('MIN', 'MAX'),
                        help="Override y-axis limits with minimum and maximum values.")
    args = parser.parse_args()

    y_lim = tuple(args.y_lim) if args.y_lim else None
    process_data(args.input_txt, args.output_dir, args.y_ticks, args.x_ticks, y_lim)

if __name__ == "__main__":
    main()

