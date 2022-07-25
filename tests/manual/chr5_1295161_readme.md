## Method

* Output from:
  
    mumbai sam postgresql://mumbai@hpcnode058/570_MELA_mk2_normal \
        chr5 1295161 1295161

  was converted to BAM and indexed using `samtools`. The position was viewed in
  IGV with all read filters turned off in `View --> Preferences --> Alignments`.
  `chr5_1295161_IGV.png` is a screenshot of the chr5:1295161 position counts.

* Output from:

    mumbai pileup postgresql://mumbai@hpcnode058/570_MELA_mk2_normal \
        chr5 1295161 1295161

  was directed to chr5_1295161_pileup.txt.

## Result

  The numbers match exactly, when taking account of the two reads with N at
  the position that are in the SAM format output but aren't counted by the
  current implementation of pileup in `mumbai`.
