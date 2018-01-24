CDF_TABLE=[7,67,247,382]
my_list=[]
for i in range(1,1024):
    def xfunc(x):
        y=(x*i)%1024
        return y
    mod_table=sorted(list(map(xfunc, CDF_TABLE)))
    my_list.append((i,mod_table))
#    print("{0}: {1}".format(i,mod_table))
my_new=sorted(my_list,key=lambda x:x[1][3],reverse=True)

for i in my_new:
    print("{0}: {1}".format(i[0],i[1]))
