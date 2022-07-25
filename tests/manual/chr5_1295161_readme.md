## Method

* Output from:

    ```
    mumbai pileup postgresql://mumbai@hpcnode058/570_MELA_mk2_normal \
        chr5 1295161 1295161
    ```

  was directed to `chr5_1295161_pileup.txt`.

* Output from:
    ```
    mumbai sam postgresql://mumbai@hpcnode058/570_MELA_mk2_normal \
        chr5 1295161 1295161
    ```

  was converted to BAM and indexed using `samtools`. The position was viewed in
  IGV with all read filters turned off in `View --> Preferences --> Alignments`.

  `chr5_1295161_IGV.png` is a screenshot of the chr5:1295161 position counts.


* Output from:

    ```
    samtools mpileup \
        -A \
        -Q 0 \
        -d 30000 \
        -r chr5:1295161-1295161 \
        --ff 0 \
        MELA_mk2_normal_TERT_5_1295161.bam|\
        cut -f 5|\
        sed 's/./&\n/g'|\
        tr '[:lower:]' '[:upper:]'|\
        sort|uniq -c|grep [ACGTN]|sort -n
    ```

  was directed to `chr5_1295161_mpileup.txt`. Note the options required to
  force mpileup to consider all reads:

    ```
    -A        # do not discard anomalous read pairs
    -Q 0      # do not skip bases with low BAQ
    -d 30000  # increase max read depth beyond default of 8000
    -ff 0     # don't filter on any alignment flags (principally DUP)
    ```

## Result

  The numbers match exactly, when taking account of the two reads with N at
  the position that are in the SAM format output but aren't counted by the
  current implementation of pileup in `mumbai`.
