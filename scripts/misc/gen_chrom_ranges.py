def chunks(l, n):
    if n<1:
        n=1
    return [l[i:i+n] for i in range(0, len(l), n)]

#Chr III
maximum = 13783700
step = 500000
print '\n'.join(["chrIII:%s-%s" % (x+1,(x+maximum/(maximum/step))) for x in range(0,maximum, maximum/(maximum/step)) if x+step < maximum])

# Chr V

maximum = 20924149
step = 500000
print '\n'.join(["chrV:%s-%s" % (x+1,(x+maximum/(maximum/step))) for x in range(0,maximum, maximum/(maximum/step)) if x+step < maximum])
