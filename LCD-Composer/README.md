## New!
The LCD-Composer webserver has been officially released! You can access the server at (**http://lcd-composer.bmb.colostate.edu/**) and the corresponding paper ([***Cascarina and Ross (2022), Bioinformatics***](https://pubmed.ncbi.nlm.nih.gov/36282522/)). We have also added new command-line scripts in the WebserverScripts directory for those interested in running high-throughput analyses for any of the options available on the webserver.

# LCD-Composer
LCD-Composer is a ***l***ow-***c***omplexity ***d***omain ***compo***sition ***s***cann***er*** designed to identify simple or multifaceted low-complexity domains (LCDs) in protein sequences, described in [***Cascarina et al. (2021), NAR Genomics and Bioinformatics***](https://academic.oup.com/nargab/article/3/2/lqab048/6285187).

LCD-Composer identifies LCDs by calculating the amino acid composition and linear dispersion of amino acids at each position in a protein sequence using a scanning window of defined size. For a full description of how the algorithm works (complete with graphical representations of algorithm workflow and extensive parameter testing), see the publication cited above.

LCD-Composer runs on Python 3 and does not require any external dependencies. The latest version of Python can be downloaded from https://www.python.org/downloads/.<br/>

## Basic Usage
    python LCD-Composer.py Sequences_File Results_File [-a AMINO_ACIDS] [-c COMPOSITION] [-w WINDOW_SIZE] [-d DISPERSION] [-i IGNORE_DISPERSION_THRESHOLD] [-v]

***positional arguments:***<br/>
| Positional Argument | Description |
| --- | --- |
| Sequences_File | The name of the file containing the protein sequences you wish to search for LCDs (in FASTA format). |
| Results_File | The name of the results file that you want to create and store the resulting LCD data. |
<br />

***required arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -a AMINO_ACIDS &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | --amino_acids AMINO_ACIDS &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Amino acid(s) of interest (single letter abbreviation, with a single underscore between each amino acid group when desired). |
<br />

***optional arguments:***<br/>
| Shortcut usage | Argument usage | Description |
| --- | --- | --- |
| -c COMPOSITION | --composition COMPOSITION | Percent composition threshold for amino acid(s) of interest (0-100). |
| -w WINDOW_SIZE | --window_size WINDOW_SIZE | Sliding window size. |
| -d DISPERSION | --dispersion DISPERSION | Linear dispersion threshold for amino acid(s) of interest (0.0-1.0). |
| -i IGN_DISP_THRESHOLD | --ignore_dispersion_threshold IGN_DISP_THRESHOLD | Threshold at which to ignore the linear dispersion parameter (0.0-1.0). |
| -v | --verbose | Output verbose data for each protein (per-position composition and dispersion values). |
<br />

***The following tutorial video demonstrates how to run LCD-Composer with default or customizable parameters:***<br/>
https://www.youtube.com/watch?v=0fybpRAMi6k&feature=youtu.be

***If you are completely new to using the command line, the following tutorial video demonstrates a few basic commands that will help you navigate through your file system and run Python scripts:***<br/>
https://www.youtube.com/watch?v=I00rFsnfgeg&feature=youtu.be

### Run LCD-Composer with Default Parameters
LCD-Composer is a Python3 script designed to be run as a stand-alone command-line application. To run LCD-Composer on your sequences of interest, download the LCD-Composer.py script and save to a location containing a FASTA file with your protein sequences of interest. Navigate to that location via command line, and run LCD-Composer with the following command (will use default parameters):

    python LCD-Composer.py Sequences_File Results_File -a <amino_acids>

The "-a" flag is required, and allows you to specify the amino acid(s) that you are interested in for your domain search (see **Detailed Usage and Customizable Parameters** section below for complete description).

NOTE: Make sure to include the file extension in the command above for your file containing FASTA-formatted sequences. FASTA files will often have the file extension ".fa", ".fsa", or ".FASTA", but are sometimes also provided as plain-text files (.txt), which should still work with LCD-Composer. LCD-Composer is designed to output your results in a **t**ab-**s**eparated **v**alues (.tsv) file. This file type was chosen for two main reasons: 1) .tsv files can be easily parsed in downstream computational processing and avoids using comma-delimiters which are often present in FASTA headers, and 2) .tsv files can be opened by Microsoft Excel for the typical user. However, if Microsoft Excel is not set as the default program to open .tsv files, the file must be opened from *within* Excel (i.e. first open Excel, then open the results file from within Excel). Alternatively, you can first change your system settings to open .tsv files with Excel by default. Please note that Excel notoriously re-formats some gene/protein names to date formats automatically: however, this only occurs if the gene/protein name constitutes the entire FASTA record header in the sequence file, and only for a small number of genes/proteins (e.g. SEPT7, MARCH1, etc.).

## Multi-proteome Version of LCD-Composer
We've also created a "MultiProteome" version of LCD-Composer to make it easier to run the same LCD search on multiple proteomes with a single command! This version uses the exact same LCD-identification algorithm as the original LCD-Composer, but allows you to scale your analyses easily to hundreds or thousands of proteomes at once and compare LCD content across proteomes. It will also output basic statistics for each proteome including the total number of LCDs, the number of proteins containing at least one LCD, and the percentage of proteins with at least one LCD.

To run this version, simply download the "LCD-Composer_MultiProteome.py" script, place it in a folder/directory containing your FASTA files of interest, and run via command line like you would for LCD-Composer.py (examples below). This new version contains two additional optional parameters: [-k FILE_KEYWORD] and [-s]. The -k flag allows you to limit your analyses to only FASTA files containing a keyword of your choosing. The -s flag allows you to choose whether you would like to analyze FASTA files in the subdirectories of the current folder/directory. Below are example commands using these flags:

    python LCD-Composer_MultiProteome.py Sequences_File Results_File -a A -k YeastStrain    *(analyzes all FASTA files containing the phrase "YeastStrain" in the current directory)*
    
    python LCD-Composer_MultiProteome.py Sequences_File Results_File -a A -k YeastStrain -s     *(analyzes all FASTA files containing the phrase "YeastStrain" in the current directory and all subdirectories)*

## Detailed Usage and Customizable Parameters
Mulitple LCD-Composer parameters are customizable at runtime for more targeted CompoSer searches. These include:
1. Amino acid(s) of interest
2. Composition threshold(s)
3. Window size
4. Linear dispersion threshold
5. "Verbose" output

In the following sections, we illustrate the usage of each parameter.
<br/>
### Amino Acid(s) of Interest (-a)
LCD-Composer requires users to specify an amino acid or group(s) of amino acids of interest for each run. Amino acid(s) are specified using the "-a" flag followed by a space, followed by the single-letter abbreviation for the amino acid of interest. For example, a search for Q-rich regions (with other LCD-Composer parameters set to default), the command would be:

    python LCD-Composer.py Sequences_File Results_File -a Q

LCD-Composer also permits searches for domains enriched in a specific set of residues. For example, users can search for domains enriched in negatively charged residues (D and/or E) with the following command:

    python LCD-Composer.py Sequences_File Results_File -a DE

### Composition thresholds (-c)
The default composition threshold for LCD-Composer is 40% composition. This means that the specified amino acid(s) of interest must be at least 40% of a given sequence window for further consideration as a domain of interest. Alternative composition thresholds can be specified using the "-c" flag, followed by a space, followed by any value from 0-100 (inclusive). For example, to search for N-rich regions that are at least 60% N, the command would be:

    python LCD-Composer.py Sequences_File Results_File -a N -c 60

LCD-Composer also allows for "multifaceted" composition criteria, where distinct amino acids or groups of amino acids are assigned differnt thresholds. For example, users may be interested in domains that are at least 40% S ***and*** at least 20% A. This can be done by separating the amino acids with an underscore, while also separating distinct composition thresholds by an underscore. The command would be:

    python LCD-Composer.py Sequences_File Results_File -a S_A -c 40_20

This can also be performed with groups of amino acids, by separating groups of amino acids with the underscore delimiter. For example, a search for domains that are at least 50% Q and/or N ***and*** at least 15% Y, the command would be:

    python LCD-Composer.py Sequences_File Results_File -a QN_Y -c 50_15

Notice that any amino acids that are not separated by the underscore will be considered a group, and their combined composition for each window will be considered in the calculation.

### Window size (-w)
By default, LCD-Composer uses a sliding window size of 20 amino acids. To use an alternative window size, use the "-w" flag, followed by a space, followed by any positive integer value. For example, a search for long S-rich domains at least 60 residues in length would be:

    python LCD-Composer.py Sequences_File Results_File -a S -w 60

### Linear dispersion threshold (-d)
The linear dispersion parameter is a normalized measure of the dispersion of the amino acid(s) of interest within each window. Linear dispersion will always be a decimal value from 0-1, with higher values indicating greater linear dispersion (i.e. approaching perfect spacing of the amino acid(s) of interest within a given window sequence). By default, LCD-Composer uses a default linear dispersion threshold of 0.5, with the added caveat that the linear dispersion parameter is ignored if the composition of a window sequences exceeds the midpoint between the composition threshold and 100%. For example, at a composition threshold of 40%, the linear dispersion parameter will be ignored for windows with at least 70% composition corresponding to the amino acid(s) of interest.

To define a new linear dispersion threshold, e.g. 0.7, in a search for N-rich domains, the command would be:

    python LCD-Composer.py Sequences_File Results_File -a N -d 0.7

### "Verbose" output (-v)
Some users may be interested in the composition and linear dispersion values assigned to each position of the protein. This can be achieved by using the "-v" flag (with no other trailing arguments or characters). For example, the per-position composition and dispersion of Q residues can be assessed using the following command:

    python LCD-Composer.py Sequences_File Results_File -a Q -v

NOTE: verbose mode forces CompoSer to perform the complete set of calculations for each position in a protein. Consequently, verbose output runs will be slower than the default LCD-Composer, which only performs full calculations for identified LCDs.

### Combining customizable commands
All LCD-Composer commands can be used in combination for highly selective searches. Below are examples of combined commands, and a brief description of what each search is designed to accomplish:

__Search for P-rich domains that are at least 65% P, at least 50aa long, and with moderate minimum dispersion of 0.6:__

    python LCD-Composer.py Sequences_File Results_File -a P -c 65 -w 50 -d 0.6

__Search for domains that are at least 40% D or E **and** at least 40% K or R (mixed charge domains), that are at least 30aa long, and with moderate minimum dispersion of 0.4:__
    
    python LCD-Composer.py Sequences_File Results_File -a DE_KR -c 40_40 -w 30 -d 0.4

__Search for domains that are at least 30% G **and** at least 15% R or Y (often associated with RGG domains), that are at least 60aa long, and with relatively high minimum dispersion of 0.7:__
    
    python LCD-Composer.py Sequences_File Results_File -a G_RY -c 30_15 -w 60 -d 0.7

__Search for aromatic-rich domains that are at least 25% F, W, or Y; at least 35aa long; with high minimum dispersion of 0.8; and output the per-position values for each protein in "verbose" mode:__
    
    python LCD-Composer.py Sequences_File Results_File -a FWY -c 25 -w 35 -d 0.8 -v

## License info
LCD-Composer is subject to the terms of the GPLv3 license.
