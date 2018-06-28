
# coding: utf-8

# ##### The First graph

# In[1]:


import matplotlib.pyplot as plt
import csv


# In[6]:


num_of_tags = []
with open('tags_span_tr_vntr.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    plots.next()
    for row in plots:
        num_of_tags.append(int(row[1]))
count= {}
# so the values gives me the number of trs that have that many tags the keys so 10705:1 means 
# that 1 tr has 10705 tags spanning it 
for x in range(len(num_of_tags)):
    if num_of_tags[x] in count:
        count[num_of_tags[x]] += 1
    else:
        count[num_of_tags[x]] = 1


#plt.scatter(count.keys(), count.values())
plt.bar(count.keys(), count.values(), edgecolor='black')
#plt.hist(count.keys(), bins=200, color='purple')
plt.xlim(0,50)
plt.show()


# In[ ]:


# #counting the number of regions 
# import pandas as pd 
# import csv


# names_10_marzie = ['tag', 'chr','reedstart', 'reedstop','fragments','num1','num2','calculation']
# dtypes_marzie= {'chr':'str','tag':'str','reedstart':'int', 'reedstop':'int', 'fragments':'int', 'num1':'int', 'num2':'str', 'calculation':'float'}

# marzie = pd.read_csv('barcodes_mapping_witho_chr.csv', names = names_10_marzie, dtype = dtypes_marzie, nrows=4000)

# sorted_marzie = marzie.sort_values(by= ['tag'])
# sorted_marzie.set_index('tag', inplace=True)

# # f = open('tag_regions_marzie_2_test .csv', 'w')
# # f.write('tag,num_regions_covered'+'\n')
# dataframe = pd.DataFrame(columns=['tag','num_regions_covered'])

# line_number= 0
# seen_tags = {}
# for index, x in sorted_marzie.iterrows():
#     complete_regions = sorted_marzie.loc[[index]]
#     count = len(complete_regions)
#     if x.name in seen_tags:
#         pass
#     else:
#         seen_tags[x.name] = count
# #         f.write(str(x.name)+','+str(count)+'\n')
#     if index.count == 0:
#         print("Started writing to the file")
#     if index == (sorted_marzie_total_line//4):
#         print("About 1/4 the way there")
#     if index == (sorted_marzie_total_line//2):
#         print("About 1/2 the way there")
#     if index == 3 * (sorted_marzie_total_line//4):
#         print("About 3/4 the way there")
# # f.close()
    


# In[19]:


region_count= {}
#put it all in a dictionary with index, and regions lengths 
with open('barcodes__10x_mapping_with_region_length.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        region_count[plots.line_num] = int(row[8])   

#first plot
plt.figure(figsize=(10,10))
plt.hist(int_values,bins=range(min(int_values), max(int_values), 10000), color='purple')
plt.ylabel("Individual Regions(Non-Unique)")
plt.xlabel("Region Length(10,000)")
plt.savefig('region_length_rows.png')


#second plot
plt.figure(figsize=(10,10))
plt.hist(int_values,  bins=range(min(int_values), max(int_values), 10000), color='purple', edgecolor = 'black')
plt.xlim(0, 300000)
plt.ylabel("Individual Regions(Non-Unique)")
plt.xlabel("Region Length(10,000)")
plt.savefig('region_length_rows_2.png')
#
#plt.hist(int)


# In[ ]:




