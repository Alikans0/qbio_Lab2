
# Monday 10/10/2022

## Get familiar with the model
There is 7 published models in the pdb with **GHS-R**:
- [7F83][1]: Crystal Structure of a receptor in Complex with inverse agonist.
- [6KO5][2]: Complex structure of ghrelin receptor with Fab.
- [7F9Z][3]: GHRP-6-bound Ghrelin receptor in complex with Gq.
- [7NA8][4]: Structures of human ghrelin receptor-Gi complexes with ghrelin and a synthetic agonist.
- [7F9Y][5]: ghrelin-bound ghrelin receptor in complex with Gq.
- [7NA7][6]: Structures of human ghrelin receptor-Gi complexes with ghrelin and a synthetic agonist.
- [7W2Z][7]: Cryo-EM structure of the ghrelin-bound human ghrelin receptor-Go complex.

I first aligned all the models in Pymol, and compared them.

## Issues
We are interested in [7F9Z][8] (see [Article][9]).  

This model shows a lot of issues:
- The peptide has **3 Cis peptide bonds**!
- The model used [6KO5][10] as a template, which also has a problem since helix 5 shows a shift of half helix turn.   

### Using VMD, load density CPP4
So when we try to fit [7F9Z][11] in the **density map ** , the peptide does not fit well, and we notice that the **M213** in helix 5 is not in its place.   
The resolution is reported to be **3.20Å**, this is the average value. We can tell that the resolution in the whole model is so much better than in the ligand pocket, suggesting that the resolution is lower for the peptide. So the density map for the peptide is probably for the backbone, explaining the abnormal 3 cis peptide bonds when they tried to fit it with the side chains.  

_Side note:_ The configuration L and D of the lysine in the peptide is so important (switch from agonist to and antagonist respectively)

![][image-1]

## Discussion 
Got introduced to multiple concepts: 
- Calculations with Floating 
- NVT & NPT
- PBD
- Harmonic potential
- Morse potential 
- Cutoff







# Tuesday 10/11/2022

[1]:	https://www.rcsb.org/structure/7F83
[2]:	https://www.rcsb.org/structure/6KO5
[3]:	https://www.rcsb.org/structure/7F9Z
[4]:	https://www.rcsb.org/structure/7NA8
[5]:	https://www.rcsb.org/structure/7F9Y
[6]:	https://www.rcsb.org/structure/7NA7
[7]:	https://www.rcsb.org/structure/7W2Z
[8]:	https://www.rcsb.org/structure/7F9Z
[9]:	https://www.nature.com/articles/s41467-021-25364-2
[10]:	https://www.rcsb.org/structure/6KO5
[11]:	https://www.rcsb.org/structure/7F9Z

[image-1]:	file:///.file/id=6571367.43553643