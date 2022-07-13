# on CRL_dev branch using BGZFile but manually constructing AlignedSegments:

```bash
[conradL@hpcnode062 mumbai]$ time LOG_LEVEL=DEBUG mumbai pileup postgresql://mumbai@hpcnode058/570_MELA_mk2_tumour chr7 140453136 140453136
chr7:[140453136-140453136]
b"SELECT bam.filepath, chunk.chunk_beg FROM bam JOIN intv ON intv.bam_id = bam.id JOIN chunk ON intv.bam_id = chunk.bam_id AND intv.ref_id = chunk.ref_id JOIN ref ON chunk.ref_id = ref.id WHERE ref.sn = 'chr7' AND intv.window = 8572 AND chunk.bin_id IN (0, 3, 25, 206, 1656, 13253) AND chunk.chunk_end >= ioffset AND chunk.chunk_beg <= ioffset"
570 bams queried
chr7	140453136	A	  A:48297   T:8533    C:144     G:69     -:23      |:0

real	0m13.149s
user	3m11.943s
sys	0m0.991s
[conradL@hpcnode062 mumbai]$ time LOG_LEVEL=DEBUG mumbai pileup postgresql://mumbai@hpcnode058/570_MELA_mk2_normal chr12 20001191 20001195
chr12:[20001191-20001195]
b"SELECT bam.filepath, chunk.chunk_beg FROM bam JOIN intv ON intv.bam_id = bam.id JOIN chunk ON intv.bam_id = chunk.bam_id AND intv.ref_id = chunk.ref_id JOIN ref ON chunk.ref_id = ref.id WHERE ref.sn = 'chr12' AND intv.window = 1220 AND chunk.bin_id IN (0, 1, 11, 92, 737, 5901) AND chunk.chunk_end >= ioffset AND chunk.chunk_beg <= ioffset"
570 bams queried
chr12	20001191	C	  C:24605    A:136     T:97     |:75     -:39     G:29
chr12	20001192	T	  T:24420    |:156     A:75     G:69     C:36     -:34
chr12	20001193	G	  |:17963   G:6236    A:467     -:36     T:30     C:12
chr12	20001194	A	  A:25517     -:55     |:35      G:8      C:5      T:1
chr12	20001195	A	  A:25534     |:37      C:7      G:7      -:2      T:1

real	0m6.184s
user	1m51.995s
sys	0m0.531s
```

# On main using AlignmentFile and using pysam native AlignedSegment reading:

```bash
[conradL@hpcnode062 mumbai]$ time LOG_LEVEL=DEBUG mumbai pileup postgresql://mumbai@hpcnode058/570_MELA_mk2_tumour chr7 140453136 140453136
chr7:[140453136-140453136]
b"SELECT bam.filepath, chunk.chunk_beg FROM bam JOIN intv ON intv.bam_id = bam.id JOIN chunk ON intv.bam_id = chunk.bam_id AND intv.ref_id = chunk.ref_id JOIN ref ON chunk.ref_id = ref.id WHERE ref.sn = 'chr7' AND intv.window = 8572 AND chunk.bin_id IN (0, 3, 25, 206, 1656, 13253) AND chunk.chunk_end >= ioffset AND chunk.chunk_beg <= ioffset"
570 bams queried
chr7	140453136	A	  A:48297   T:8533    C:144     G:69     -:23      |:0

real	0m2.389s
user	0m28.670s
sys	0m0.468s
[conradL@hpcnode062 mumbai]$ time LOG_LEVEL=DEBUG mumbai pileup postgresql://mumbai@hpcnode058/570_MELA_mk2_normal chr12 20001191 20001195
chr12:[20001191-20001195]
b"SELECT bam.filepath, chunk.chunk_beg FROM bam JOIN intv ON intv.bam_id = bam.id JOIN chunk ON intv.bam_id = chunk.bam_id AND intv.ref_id = chunk.ref_id JOIN ref ON chunk.ref_id = ref.id WHERE ref.sn = 'chr12' AND intv.window = 1220 AND chunk.bin_id IN (0, 1, 11, 92, 737, 5901) AND chunk.chunk_end >= ioffset AND chunk.chunk_beg <= ioffset"
570 bams queried
chr12	20001191	C	  C:24605    A:136     T:97     |:75     -:39     G:29
chr12	20001192	T	  T:24420    |:156     A:75     G:69     C:36     -:34
chr12	20001193	G	  |:17963   G:6236    A:467     -:36     T:30     C:12
chr12	20001194	A	  A:25517     -:55     |:35      G:8      C:5      T:1
chr12	20001195	A	  A:25534     |:37      C:7      G:7      -:2      T:1

real	0m1.392s
user	0m18.528s
sys	0m0.377s
```

# conclusion

It's 4-5 times faster using pysam's reader (which I believe is using htslib directly)
to create AlignedSegments even though we're paying a disk read cost to read the header
to create the AlignmentFile.
