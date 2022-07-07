#!/usr/bin/env python3
"""
Inspect the contents of a .bai file.

See http://samtools.github.io/hts-specs/SAMv1.pdf section 5.2
"""
import struct
import sys

if len(sys.argv) != 2 or (
       len(sys.argv) == 2 and (sys.argv[1] == '-h' or
                               sys.argv[1] == '--help')):
    print(f"USAGE: {sys.argv[0]} BAIFILE")
    sys.exit(1)

bai = sys.argv[1]

SHIFT_AMOUNT = 16

def split_vfp(vfp):
    """
    Splits a long virtual file pointer into (file pointer, block offset)
    """
    start = vfp >> SHIFT_AMOUNT
    return (start, vfp ^ (start << SHIFT_AMOUNT))

with open(bai, 'rb') as fh:
    fh.seek(0)
    magic = fh.read(4).decode('ascii')
    print(f"magic: {magic}")
    n_ref = struct.unpack('i', fh.read(4))[0]
    print(f"n_ref: {n_ref}")
    for r in range(n_ref):
        n_bin = struct.unpack('i', fh.read(4))[0]
        print(f"n_bin: {n_bin}")
        for i in range(n_bin):
            bin_ = struct.unpack('i', fh.read(4))[0]
            print(f"bin: {bin_}")
            n_chunk = struct.unpack('i', fh.read(4))[0]
            print(f"n_chunk: {n_chunk}")
            if bin_ != 37450:
                for c in range(n_chunk):
                    chunk_beg = struct.unpack('l', fh.read(8))[0]
                    chunk_end = struct.unpack('l', fh.read(8))[0]
                    print(f"chunk_beg: {split_vfp(chunk_beg)}")
                    print(f"chunk_end: {split_vfp(chunk_end)}")
            else:
                # metadata pseudo-bin
                ref_beg = struct.unpack('l', fh.read(8))[0]
                ref_end = struct.unpack('l', fh.read(8))[0]
                n_mapped = struct.unpack('l', fh.read(8))[0]
                n_unmapped = struct.unpack('l', fh.read(8))[0]
                print(f"ref_beg: {split_vfp(ref_beg)}")
                print(f"ref_end: {split_vfp(ref_end)}")
                print(f"n_mapped: {n_mapped}")
                print(f"n_unmapped: {n_unmapped}")
            #input("press any key to continue...")
        n_intv = struct.unpack('i', fh.read(4))[0]
        print(f"n_intv: {n_intv}")
        for v in range(n_intv):
            ioffset = struct.unpack('l', fh.read(8))[0]
            print(f"ioffset: {split_vfp(ioffset)}")
    n_no_coor = struct.unpack('l', fh.read(8))[0]
    print(f"n_no_coor: {n_no_coor}")
    assert fh.read() == b'', 'not at EOF'
