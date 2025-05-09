---
Final: "Giant Panda Associated Gemycircularvirus Sequence Analysis Project"
Author: "Josie Martinez"
Date: "May 1, 2025"
Accession: "MF327571"

---

# Introduction
![image](https://github.com/user-attachments/assets/88f3bc75-1d8a-458e-86d9-55ff5185cb48)

- **Viral classification:**  
  - *ICTV classification*:
    - Realm: Monodnaviria, Kingdom: Shotokuvirae, Phylum: Cressdnaviricota, Class: Repensiviricetes, Order: Geplafuvirales, Family: Genomoviridae, Genus: Gemykibivirus, Species: Gemykibivirus giapa1. ([3](#ncbi))
      
  - *Baltimore classification*:
    - Group II single-stranded DNA(+/-) genome. ([4](#viralzone))
      
- **Physical size:**  
  - The physical size of the gemycircular virus is approximately 20-22 nm, which is smaller than a typical human cell (~10,000 nm) and smaller than SARS-CoV-2 (~120 nm). ([4](#viralzone))
    
- **Shape and envelope:**  
  - The virus exhibits an icosahedral morphology and does not possess an envelope. ([4](#viralzone))

- **Discovery and outbreaks:**  
  - According to a CDC article, gemycircularviruses was found initially in fungi in 2010. There are no known putbreaks since gemycircularviruses are very stable in the environment. ([2](#uch2015))
    
- **Host range:**  
  - This virus infects humans, mammals, birds, fungi and is host-specific. It is suggested that the detection of gemycircularvirus genomes in giant pandas may be consumed from fungi or insects. ([4](#viralzone)) ([1](#zhang2017))

- **Cell entry:**  
  - Entry mechanism is unknown but it could most likely be through endocytosis due to the virus lacking an envelope. ([4](#viralzone))
    
- **Replication strategy:**  
  - Virus attaches to host cell, and genomic DNA penetrates into the cytoplasm.
  - Presumably the ssDNA is converted to dsDNA.
  - Transcription of viral genes, possibly by a host RNA polymerase.
  - Genome replication, possibly by rolling circle
  - Assembly and release of virions.
    ([4](#viralzone))

- **Release mechanism:**  
  - Although there is no confirmed release mechanism, it is assumed that the virus penetrates the host cell via cell lysis. Cell lysis is typically very common in small non-enveloped ssDNA viruses. ([4](#viralzone))

- **Latency:**  
  - There is no direct evidence that gemycircularvirues establish latency. ([2](#uch2015)) But we can hypothesize that it could be possible due to the structure of the genome. 

- **Equilibrium and antigenic shift:**  
  - Because very little is known about the pathlogy of gemycircularviruses, we do not currently know if they reach equilibrium with any host. 

- **Vaccines:**  
  - There are no known vaccines for this virus.

- **Antiviral drugs:**  
  - There are no antiviral drugs indentified.

# Methods

1. **First, I downloaded the viral sequence by accession number.** 

```python
from Bio import Entrez, SeqIO

Entrez.email = "jmart411@uncc.edu" 
accession = "MF327571"

with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
    record = SeqIO.read(handle, "fasta")

SeqIO.write(record, "MF327571.fasta", "fasta")

```
2. **Next, I used the code below which uses Biopython to identify open reading frames (ORF's) in the virus sequence and translate them to proteins.**

```python
from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(sequence, min_length=300):
    orfs = []
    seq_len = len(sequence)
    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            trans = nuc[frame:].translate(to_stop=False)
            trans_len = len(trans)
            aa_start = 0
            while aa_start < trans_len:
                aa_start = trans.find("M", aa_start)
                if aa_start == -1:
                    break
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    break
                if (aa_end - aa_start) * 3 >= min_length:
                    start = frame + aa_start * 3
                    end = frame + aa_end * 3 + 3
                    if strand == -1:
                        start, end = seq_len - end, seq_len - start
                    orf_seq = nuc[start:end] if strand == 1 else nuc.reverse_complement()[start:end]
                    orfs.append(orf_seq.translate())
                aa_start = aa_end + 1
    return orfs

# Load sequence
record = SeqIO.read("MF327571.fasta", "fasta")
orf_proteins = find_orfs(record.seq)

# Save protein sequences
with open("MF327571_proteins.faa", "w") as output:
    for i, protein in enumerate(orf_proteins, 1):
        output.write(f">orf{i}\n{str(protein)}\n")


```
I now have the MF327571_proteins.faa file with the translated proteins. 

3. **I gathered accession codes provided in the spreadsheet along with the sequence lengths.** This is needed in order to make the IQ Tree. 
```python
Entrez.email = "jmart411@uncc.edu"

accession_codes = {
    "Porcine feces-associated gemycircularvirus": "KY214433",
    "Badger associated gemykibivirus 1": "KP263543",
    "Blackbird associated gemykibivirus 1": "KF371633",
    "Black robin associated gemykibivirus 1": "KF371634",
    "Bovine associated gemykibivirus 1": "LK931483",
    "Canine feces-associated gemycircularvirus": "KY214441",
    "Cattle blood-associated gemycircularvirus": "MF669480",
    "Chicken genomovirus mg7_73": "MN379612",
    "Chicken genomovirus mg8_401": "MN379615",
    "Chicken genomovirus mg4_1218": "MN379608",
    "Genomoviridae sp.": "MK032708",
    "Finch associated genomovirus 3": "MK249305",
    "Gopherus associated genomovirus 1": "MF373638",
    "Finch associated genomovirus 2": "MK249239",
    "Finch associated genomovirus 4": "MK249240",
    "Finch associated genomovirus 8": "MK249245",
    "Capybara genomovirus 10": "MK483082",
    "Rabbit associated gemykroznavirus 1": "KF371631",
    "Ostrich associated gemytondvirus 1": "KF371630",
    "Human associated gemyvongvirus 1": "KP974693",
    "Circoviridae sp. ctcc28": "MK012475"
}
def fetch_gemycircularvirus_fasta(accession_codes):
    sequences = {}
    for name, accession in accession_codes.items():
        try:
            with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
                fasta_data = handle.read().strip()
            sequences[name] = fasta_data
            print(f"Retrieved {name} ({accession})")
            time.sleep(0.35)  # to avoid overloading NCBI servers
        except Exception as e:
            print(f"Error retrieving {name} ({accession}): {str(e)}")
            sequences[name] = None
    return sequences

fasta_sequences = fetch_gemycircularvirus_fasta(accession_codes)

with open("gemycircularviruses.fasta", "w") as f:
    for name, fasta in fasta_sequences.items():
        if fasta:
            f.write(fasta + "\n\n")
```
My output shows the sequence lengths of all similar sequences. 
```python
Output: 
Retrieved Porcine feces-associated gemycircularvirus (KY214433)
Retrieved Badger associated gemykibivirus 1 (KP263543)
Retrieved Blackbird associated gemykibivirus 1 (KF371633)
Retrieved Black robin associated gemykibivirus 1 (KF371634)
Retrieved Bovine associated gemykibivirus 1 (LK931483)
Retrieved Canine feces-associated gemycircularvirus (KY214441)
Retrieved Cattle blood-associated gemycircularvirus (MF669480)
Retrieved Chicken genomovirus mg7_73 (MN379612)
Retrieved Chicken genomovirus mg8_401 (MN379615)
Retrieved Chicken genomovirus mg4_1218 (MN379608)
Retrieved Genomoviridae sp. (MK032708)
Retrieved Finch associated genomovirus 3 (MK249305)
Retrieved Gopherus associated genomovirus 1 (MF373638)
Retrieved Finch associated genomovirus 2 (MK249239)
Retrieved Finch associated genomovirus 4 (MK249240)
Retrieved Finch associated genomovirus 8 (MK249245)
Retrieved Capybara genomovirus 10 (MK483082)
Retrieved Rabbit associated gemykroznavirus 1 (KF371631)
Retrieved Ostrich associated gemytondvirus 1 (KF371630)
Retrieved Human associated gemyvongvirus 1 (KP974693)
Retrieved Circoviridae sp. ctcc28 (MK012475)
>>> 
>>> with open("gemycircularviruses.fasta", "w") as f:
...     for name, fasta in fasta_sequences.items():
...         if fasta:
...             f.write(fasta + "\n\n")
... 
2328
2227
2229
2274
2230
2364
2311
2278
2273
2264
2251
2351
2329
2343
2259
2227
2155
2269
2222
2380
2428
>>>
```

4. **Next, I aligned the sequences using the MAFFT slurm script**
```bash
#!/bin/bash
#SBATCH --job-name=mafft_align
#SBATCH --output=mafft_align.out
#SBATCH --error=mafft_align.err
#SBATCH --time=01:00:00
#SBATCH --partition=Centaurus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

module load mafft

mafft --auto all_gemycircularviruses.fasta > aligned_gemycircularviruses.fasta
```
I submitted the slurm job with sbatch mafft.slurm. Once the sequences have been aligned, I submitted another slurm job to make my phylogeny tree. 

5. **IQ TREE**
```bash
#!/bin/bash
#SBATCH --job-name=iqtree_run
#SBATCH --output=iqtree_run.out
#SBATCH --error=iqtree_run.err
#SBATCH --time=02:00:00
#SBATCH --partition=Centaurus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

source ~/miniconda3/bin/activate iqtree_env

iqtree -s aligned_gemycircularviruses.fasta -m MFP -bb 1000 -nt AUTO
```
     

# Results and Discussion


## FigTree
![FinalFigTree](https://github.com/user-attachments/assets/392314bd-a03a-4845-97d3-2e5367afd7a6)

## Hydrophobicity Comparison
![image](https://github.com/user-attachments/assets/c4ad81b6-2cd5-46bf-b687-e66d2d5178ed)

# References Cited
### References

1. <a name="zhang2017"></a>Zhang, Wen et al. “Virome comparisons in wild-diseased and healthy captive giant pandas.” Microbiome vol. 5, 1 90. 7 Aug. 2017, doi: [10.1186/s40168-017-0308-0](https://doi.org/10.1186/s40168-017-0308-0).

2. <a name="uch2015"></a>Uch, R., Fournier, P., Robert, C., Blanc-Tailleur, C., Galicher, V., Barre, R., Biagini, P. (2015). Divergent Gemycircularvirus in HIV-Positive Blood, France. Emerging Infectious Diseases, 21(11), 2096-2098. [https://doi.org/10.3201/eid2111.150486](https://doi.org/10.3201/eid2111.150486).

3. <a name="ncbi"></a>U.S. National Library of Medicine. (n.d.). Giant Panda Associated Gemycircularvirus strain GPGE014, complete Genomic Sequence - NCBI. National Center for Biotechnology Information. [https://www.ncbi.nlm.nih.gov/nuccore/MF327571.1/](https://www.ncbi.nlm.nih.gov/nuccore/MF327571.1/).

4. <a name="viralzone"></a>ViralZone is operated by the Swiss-Prot group of the SIB Swiss Institute of Bioinformatics. (n.d.). Gemycircularvirus (taxid:1542744). ViralZone. [https://viralzone.expasy.org/6717](https://viralzone.expasy.org/6717).


