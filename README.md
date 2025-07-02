# Phase separation predictive analysis

Phase separation (PS) is a biophysical mechanism that drives the formation of membraneless, organelle-like condensates, that arise when specific macromolecules demix creating concentrated phases. This process relies on multivalent, weak, and reversible interactions such as cation-π, hydrophobic, and electrostatic interactions mediated by intrinsically disordered regions (IDRs). An important model describing the PS behavior of biomolecular condensates in vitro and in vivo is the known as the stickers-and-spacers. Stickers are specific motifs that mediate attractive, reversible interactions enabling the condensates' assembly, while spacers are the residues between stickers which play a crucial regulatory and structural role.

Using predictive computational tools, you can identify motifs that aid PS on your protein of interst, annotate and visualize them. First, make sure to download the FASTA file your protein of interest from [UniProt](https://www.uniprot.org/). Alternatively, you can use FASTA files for FUS, a protein with high PS propensity, and GAI, a more structired protein. Open your terminal and create a new directory using the following command:

```
mkdir yourname
cd /path/to/yourname
```

Clone the github to your directory with this command:

```
git clone 

```
## PLAAC

Let's start with PLAAC! PLAAC (Prion-Like Amino Acid Composition) analyzes protein sequences to    identify candidate prion-like domains (PLDs) using a hidden Markov model (HMM) based on compositional     similarity to known prions. First, we need to move in the directory.

```
cd /data1/projects/pi-vriesendorpb/shared/PS/plaac/cli
```

Make sure to upload your FASTA file in this directory. You can check whether it is uploaded by listing the contents.

```
ls
```

Next, load the necessary prerequisites:

```
module load  Java/11.0.20
module load RStudio-Desktop/2024.04.2+764-gfbf-2023a-Java-11-R-4.4.0
```

Now, you can use the following command to calculate the per-residue scores, which are suitable for plotting. Remember to change the name of the input.fasta with the one you are using.

```
java -jar plaac.jar -i input.fasta -p all > plotdata.txt
```

You can visualise the results using the following commands.

```
pip install pandas matplotlib  
python plot_plaac.py
```
 

## LCD-Composer

LCD-Composer is a customizable computational tool for identifying low-complexity domains (LCDs). By calculating the amino acid composition and linear dispersion of amino acids in a protein sequence, it is designed to output the results in a tab-separated values (.tsv) file. If Microsoft Excel is not set as the default program to open .tsv files, the file must be opened from within Excel.

First, you need to move to the correct directory and make sure to copy your protein FASTA file in this directory. The FUS and GAI files are already present.

```
cd /data1/projects/pi-vriesendorpb/shared/LCD-Composer
```
```
cp /pathto/yourprotein.fasta .
```


Afterwards, you can run this command in your terminal. Make sure to specify the residue or group of residues of interest after ```a``` and the threshold composition after ```-c``` (0-100). It is recommended to scan for G,P,D,E residues but you have the option to check for any specified residues that might promote PS in your protein.

``` 
python LCD-Composer.py input.fasta output -a -c
```

Example:

```
python LCD-Composer.py FUS.fasta output -a DE -c 60
```

You can download your output file and open the document in Excel or use the following command to view your results in your directory.

```
cat output.tsv
```


## Pi-Pi contacts

Utilizing this PS predictive algorithm, individual amino acids are scored based on π interaction frequency. The scoring plots provide a residue-level visualization of their propensity to undergo PS, based on π-interaction frequency. High-scoring residues, which are typically aromatic or positively charged, are enriched in intrinsically disordered regions that can facilitate PS.

Once again, you need to change to the appropriate directory and make sure that your desired FASTA file is present.
 
```
cd /data1/projects/pi-vriesendorpb/shared/PS/SourceCodeS2/
```

Next, run the following command.

```
python elife_phase_separation_predictor.py fasta/input.fasta -residue_scores -output output.txt -overwrite
```

```
python plot.py output.txt "PScore prediction yourprotein" -y_lim -10 15--y_ticks 5 --x_ticks 20
```

Your results will be stored in the "PScore prediction yourprotein" and you can view and download them.

 
## AlphaFold Modelling

By integrating all these results, you can pinpoint the specific motifs responsible for phase separation behavior. The most effective way to visualize these findings is to use the [AlphaFold Server](https://alphafoldserver.com/). Simply paste your protein’s amino acid sequence into the designated field, and once the model is generated, you can examine the relevant domains of interest.




##### If you can this far, you are awesome! We hope this was useful and feel free to let us know if you have any questions.





# References

1.  Alberti, S.; Halfmann, R.; King, O.; Kapila, A.; Lindquist, S. A Systematic Survey Identifies Prions and Illuminates Sequence Features of Prionogenic Proteins. Cell 2009, 137, 146–158, doi:10.1016/j.cell.2009.02.044.
2.  Ginell, G.M.; Holehouse, A.S. An Introduction to the Stickers-and-Spacers Framework as Applied to Biomolecular Condensates. In Phase-Separated Biomolecular Condensates: Methods and Protocols; Zhou, H.-X., Spille, J.-H., Banerjee, P.R., Eds.; Springer US: New York, NY, 2023; pp. 95–116 ISBN 978-1-0716-2663-4.
3.  Cascarina, S.M.; King, D.C.; Osborne Nishimura, E.; Ross, E.D. LCD-Composer: An Intuitive, Composition-Centric Method Enabling the Identification and Detailed Functional Mapping of Low-Complexity Domains. NAR Genomics and Bioinformatics 2021, 3, lqab048, doi:10.1093/nargab/lqab048.
4.  Vernon, R.M.; Chong, P.A.; Tsang, B.; Kim, T.H.; Bah, A.; Farber, P.; Lin, H.; Forman-Kay, J.D. Pi-Pi Contacts Are an Overlooked Protein Feature Relevant to Phase Separation. eLife 7, e31486, doi:10.7554/eLife.31486.
