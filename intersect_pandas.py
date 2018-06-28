#importing pandas
import csv 
import pandas as pd 
import multiprocessing as mp
from functools import partial

#setting the parameters for the two csv's 
names_vntr = ['header', 'tag', 'trid', 'reedid', 'reedstart', 'reedstop','trid2','chr', 'firstindex', 'lastindex', 'subheader']
names_marzie = ['tag', 'chr', 'start', 'end', 'fragments', 'num1', 'num2', 'calulation']
dtypes_marzie  = {'tag': 'str', 'chr':'str', 'start':'int', 'end':'int', 'fragments':'int', 'num1':'int', 'num2':'int', 'calulation':'float'}
dtypes_vntr = {'header':'str', 'tag':'str', 'trid':'int','reedid':'int', 'reedstart':'int', 'reedstop':'int','trid2':'int','chr':'str', 'firstindex':'int', 'lastindex':'int', 'subheader':'str'}

#loading the two csv's
marzie_10x = pd.read_csv('barcodes_mapping_witho_chr.csv',names = names_marzie, dtype = dtypes_marzie,nrows=4000)
vntr = pd.read_csv('sorted_without_chr_vntr.csv',names = names_vntr, dtype = dtypes_vntr,nrows=4000)

# marzie_10x = pd.read_csv('test1intersect.csv', names = ['chr', 'start', 'end', 'tag', 'fragments','num1','num2','calculation'])
# vntr = pd.read_csv('test2intersect.csv', names = ['chr', 'reedstart','reedstop','header','tag','trid','reedid','trid2','firstindex','lastindex','subheader'])
#can read a line at a time
#test1 = pd.read_csv('untitled.csv', names = ['name','num'], dtype = {'name': 'str', 'num':'int'})

#sort by tags and chr 
sorted_vntr_whole = vntr.sort_values(by= ['tag','chr'])
sorted_marzie1_whole = marzie_10x.sort_values(by = ['tag', 'chr'])

#set index 
sorted_marzie1_whole.set_index('tag', inplace=True)

#set dictionary of indices
list_of_index = sorted_marzie1_whole.index.drop_duplicates().tolist()
dictionary_of_tags = {el:0 for el in list_of_index}
del list_of_index[:]


#grouping by chromosome so we can index it 
grouped = sorted_vntr_whole.groupby('chr')

#all the chromosomes as seperate dataframes
chr1 = grouped.get_group('1').reset_index() 
chr2 = grouped.get_group('2').reset_index()
chr3=grouped.get_group('3').reset_index()
chr4=grouped.get_group('4').reset_index()
chr5=grouped.get_group('5').reset_index()
chr6=grouped.get_group('6').reset_index()
chr7=grouped.get_group('7').reset_index()
chr8=grouped.get_group('8').reset_index()
chr9=grouped.get_group('9').reset_index()
chr10=grouped.get_group('10').reset_index()
chr11= grouped.get_group('11').reset_index()
chr12=grouped.get_group('12').reset_index()
chr13=grouped.get_group('13').reset_index()
chr14=grouped.get_group('14').reset_index()
chr15=grouped.get_group('15').reset_index()
chr16=grouped.get_group('16').reset_index()
chr17=grouped.get_group('17').reset_index()
chr18=grouped.get_group('18').reset_index()
chr19=grouped.get_group('19').reset_index()
chr20=grouped.get_group('20').reset_index()
chr21=grouped.get_group('21').reset_index()
chr22=grouped.get_group('22').reset_index()
chrX=grouped.get_group('X').reset_index()
chrY=grouped.get_group('Y').reset_index()


previous_line= sorted_marzie1_whole.iloc[0]
print(previous_line)
print(sorted_marzie1_whole.loc[str(previous_line.name)])
#all the chromosomes
arguments = [chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY]

def intersecting(sorted_vntr, sorted_marzie):
    chromosome=sorted_vntr.chr[0] # get the chromosome name
    #opening the necessary files 
    intersected_dynamic = open('intersected_dynamic_parallel_test%s.csv' % chromosome, 'a') 
    intersected_failed = open('intersected_failed_dynamic_parallel_test%s.csv' % chromosome, 'a')
#     no_tag_found = open('no_tag_found_intersect_dynamic_test_%s.csv' % chromosome, 'a')

    #writing the headers
    intersected_dynamic.write('tag,chr,marzie_start,marzie_stop,reedstart,reedstop,trid,reed_header,vntr_chr' + "\n")
    intersected_failed.write('header,tag,trid,reedid,reedstart,reedstop,vntr_chr,firstindex,lastindex' + "\n")
#     no_tag_found.write('header,tag,trid,reedid,reedstart,reedstop,vntr_chr,firstindex,lastindex' + "\n")

    line1_marzie = None # setting marzie's line  
    for index,row_vntr in sorted_vntr.iterrows(): # for the number of reeds in vntr_file
        tag1 = row_vntr['tag'] # get their tags line by line 
        if str(tag1) in dictionary_of_tags: # if the tag is in the dictionary_of_tags runs in O(1) time
            if index>0: # so it doesn't repeatedly search for tags, make sure index is greater than 0 
                previous_line = sorted_vntr.iloc[index-1] # get the previous reed's tag
                if previous_line.tag == tag1: # if that reeds the same then stick with the current marzie's reeds
                    pass
                else: 
                    line1_marzie=sorted_marzie.loc[[str(tag1)]] # get all of marzie's reeds that have the same tag, enlarging the result to get a
            else:
                line1_marzie=sorted_marzie.loc[[str(tag1)]] # get all of marzie's reeds that have the same tag, enlarging the result to get a
            #dataframe every time
            reedstart_vntr= row_vntr['reedstart'] 
            reedstop_vntr= row_vntr['reedstop']
            chr_vntr= row_vntr['chr']
            mismatched_marzie_reeds=0 
            num_tags = len(line1_marzie)
            for index2,row_marzie in line1_marzie.iterrows(): # for every row that is matching with vntr
                reed_start_marzie = row_marzie['start']
                reed_stop_marzie=row_marzie['end']
                marzie_chr=row_marzie['chr']
                if str(marzie_chr) == str(chr_vntr) and int(reedstart_vntr) in range(int(reed_start_marzie), int(reed_stop_marzie)):
                    #if the chr matches and its in the range of the reed then add it to the file 
                    intersected_dynamic.write(str(row_marzie.name)+','+str(marzie_chr) + ',' + str(reed_start_marzie)+ ',' +str(reed_stop_marzie) + ',' +str(reedstart_vntr) + ',' +str(reedstop_vntr) + ',' +str(row_vntr['trid'])+ ',' +str(row_vntr['header'])+ ',' +str(chr_vntr)+"\n")
                else:
                    mismatched_marzie_reeds = mismatched_marzie_reeds + 1 # will count how many mismatched_tags have been processed so far
                    if mismatched_marzie_reeds == num_tags: # if you cycle through all the tags, then write the reed as mismatched             
                        intersected_failed.write(str(row_vntr['header'])+','+str(tag1) + ',' + str(row_vntr['trid'])+ ',' +str(row_vntr['reedid']) + ',' +str(reedstart_vntr) + ',' +str(reedstop_vntr) + ',' +str(chr_vntr)+ ',' +str(row_vntr['firstindex'])+ ',' +str(row_vntr['lastindex'])+"\n")
        else:
            pass
#         if str(tag1) in dictionary_of_tags: # now constant run time for sure, but a lot of memory to look up if tag is in marzie
#             line1_marzie = sorted_marzie.loc[str(tag1)] # get the lines with the tags from marzies file 
#             matching_chromosomes = line1_marzie
#             reedstart_vntr = line_vntr.reedstart
#             reedstop_vntr = line_vntr.reedstop
#             num_of_tags = len(matching_chromosomes) #number of tags (if there are 5 tags this will be 5)
#             mismatched_marzie_reeds= 0
#             for y in range(num_of_tags): # need to iterate over all the matching tags that show up "[0,1,2,3,4]"
#                 matching_line_marzie= line1_marzie.iloc[y]
#                 reed_start_marzie = matching_line_marzie.start
#                 reed_stop_marzie = matching_line_marzie.end
#                 if str(matching_line_marzie.chr) == str(line_vntr.chr) and int(reedstart_vntr) in range(int(reed_start_marzie),int(reed_stop_marzie)): # check to see if the vntr range is in the marzie range and have the same chr
#                     intersected_dynamic.write(str(matching_line_marzie.name)+','+str(matching_line_marzie.chr) + ',' + str(reed_start_marzie)+ ',' +str(reed_stop_marzie) + ',' +str(reedstart_vntr) + ',' +str(reedstop_vntr) + ',' +str(line_vntr.trid)+ ',' +str(line_vntr.header)+ ',' +str(line_vntr.chr)+"\n")
#                     # if it is in the range and has the same chromosome then add to the intersected data frame
#                 else:
#                     mismatched_marzie_reeds = mismatched_marzie_reeds + 1 # will count how many mismatched_tags have been processed so far
#                     if mismatched_marzie_reeds == num_of_tags: # if you cycle through all the tags              
#                         intersected_failed.write(str(line_vntr.header)+','+str(line_vntr.tag) + ',' + str(line_vntr.trid)+ ',' +str(line_vntr.reedid) + ',' +str(reedstart_vntr) + ',' +str(reedstop_vntr) + ',' +str(line_vntr.chr)+ ',' +str(line_vntr.firstindex)+ ',' +str(line_vntr.lastindex)+"\n")
#         else:
#             pass
        
    intersected_dynamic.close()
    intersected_failed.close()
    #no_tag_found.close()
    
#intersecting(chr1, sorted_marzie1_whole)

#multiprocessing_part 
if __name__ == '__main__': # safetly import the main function
    pool = mp.Pool(processes = 16) # make 16 processes 
    intersect_1=partial(intersecting, sorted_marzie=sorted_marzie1_whole) # can't pass in 2 arguments, so using a partial
    pool.map_async(intersect_1, arguments) # passing in the different chromosomes to worker processes, random
    pool.close() # close
    pool.join() # wait to continue





