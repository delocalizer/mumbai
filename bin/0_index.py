#!/usr/bin/env python3
"""
Create a correct but empty bai file
"""
import struct

with open('0_index', 'wb') as fh:
    # magic bytes
    fh.write('BAI\1'.encode('ascii'))
    # n_ref = 0
    fh.write(struct.pack('i', 0))
