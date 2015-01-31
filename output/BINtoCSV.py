#! /usr/bin/env python

import array
import os
import csv

num_files = 75

def read_in_chunks(file_object, chunk_size):
    read = 0
    num_doubles = os.path.getsize('sim-0.bin')/8

    while True:
        data = array.array('d')
        if num_doubles - read < chunk_size:
            chunk_size = num_doubles - read
        data.read(file_object, chunk_size)
        if not data:
            break
        read += len(data)
        yield data

for file_num in range(0,num_files):

    # read in binary file
    file_name = "sim-"+ str(file_num) + ".bin"
    bin = open(file_name,'rb')
    # array to hold double values
    data = array.array('d')

    # get number of doubles in file
    num_doubles = os.path.getsize('sim-0.bin')/8
    num_tuples = num_doubles/3

    print num_tuples

    # output file
    file_out = "sim-"+ str(file_num) + ".csv"
    output = open(file_out,'wb')
    writer = csv.writer(output)

    for piece in read_in_chunks(bin,chunk_size=3*8*4194304):
        # partion list into tuples
        data = [piece[i:i+3] for i in range(0,len(piece),3)]
        # write csv
        writer.writerows(data)

