#! /usr/bin/env python

import array
import os
import csv
import multiprocessing
import math

num_files = 75
num_threads=8

def processfile(file_name, read_start, num_doubles_to_read):
    # read in binary file
    bin = open(file_name,'rb')

    # array to hold double values
    data = array.array('d')

    # set file cursor position to start
    bin.seek(read_start)

    # read doubles into data
    data.read(bin,num_doubles_to_read);

    # partion list into tuples
    data = [data[i:i+3] for i in range(0,len(data),3)]

    return data

for file_num in range(0,num_files):

    # read in binary file
    file_name = "sim-"+ str(file_num) + ".bin"

    # open file for writing
    file_out = "sim-"+ str(file_num) + ".csv"
    output = open(file_out,'wb')
    writer = csv.writer(output)

    # get number of doubles in file
    num_bytes = os.path.getsize(file_name)
    num_doubles = num_bytes/8
    num_tuples = num_doubles/3

    print num_tuples

    # create thread pool
    pool = multiprocessing.Pool(processes=num_threads)

    max_chunk_size = 3*8*65536000
    num_chunks = int(math.ceil(num_bytes/float(max_chunk_size)))
    if num_chunks==0:
        num_chunks=1

    # num_threads-1 added to ensure we have enough for the last partial chunk
    for i in range(0,num_chunks+num_threads-1, num_threads):
        # store resulting object
        results = []

        # set work for num_threads
        for j in range(0,num_threads):
            start = (i*num_threads+j)*max_chunk_size
            num_to_read = max_chunk_size/8 if start+max_chunk_size/8 >= num_doubles else (num_bytes-start)/8
            if num_to_read > 0:
                result = pool.apply_async(processfile, args=[file_name,start,num_to_read])
                results.append(result)

            # wait for all threads to finish before writing to file
#            pool.wait()

            print 'processed chunks'
            # write csv
            for result in results:
                processfile_result = result.get(timeout=None)
                # need to seek to current position in file first?
                writer.writerows(processfile_result)
            print 'wroterows'
    # close and wait for pool to finish
    pool.close()
    pool.join()
