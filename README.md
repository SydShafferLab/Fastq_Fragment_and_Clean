# Fastq_Fragment_and_Clean
Based on the Qscore how much of my sequence should I keep? The error on next-gen sequencing increases linearly as you read the bases,
therefore it can be dificult to draw a cut off-point. Here we fragment the sequence (user defined length) and check if the fragments mean Qscore
are passing a Qscore threshold (user defined). The sequence is then build up one fragment at a time, once a fragment doesn't pass
the Qscore threshold the remaining fragments are discarded. 

Example:

Seq_1

'____________________________________'


## First fragment the Seq_1


fragment_1 - fragment_2 - fragment_3

'____________ ____________ ___________'

## Second determine if fragments pass threshold
If fragments 1 & 3 pass the Qscore threshold and fragment_2 fails then the modified sequence will be:


Mod Seq_1

'____________'

fragment_1 - fragment_2 - fragment_3


Only fragment_1 is kept. A summary.txt file is outputed with the mean, median and std of modified sequence length and Qscore.
