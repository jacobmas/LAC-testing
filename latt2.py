K=IntegerModRing(1024)
n=8
f=x**n+1
c=0
g=x**(n/4)+1
#h=x**(n/4)+160*x**(n/4)+250
R=K['x'].quotient(x**(n)+1)
T=K['x'].quotient(x**(n/4)+1)

fact_dict={}
for ww in range(0,2):
    fact_dict={}
    for i in range(0,3**n):
        j=i
        coeff_list=[]
        for l in range(0,n):
            k=j%3
            coeff_list.append(k-1)
            j=(j-k)/3
        vv=R(coeff_list)*(ww)
#        print("vv={0}".format(vv))
        t=T(vv)
        if t==0:
            continue
        y=t.list()
        if len(y)==0:
            print('wtf,t={0}'.format(t))
        if y[0] not in fact_dict:
            fact_dict[y[0]]=1
        else:
            fact_dict[y[0]]=fact_dict[y[0]]+1
    f_list=fact_dict.keys()
    f_list2=sorted(f_list)
    print("ww={0:4d}, fact_dict={1}".format(ww,f_list2))
