from wrapper import Multigene

targets = ['A','B']
multigene = Multigene()
pair = multigene.design(targets) 

#print(pair)
#print(pair.contents)
#print(pair.contents.next is None)
#print("%d %d %0.2f" % (pair.contents.targetId, pair.contents.productSize, pair.contents.freeEnergy))
p = pair
while(p): 
    forward = p.contents.forward
    reverse = p.contents.reverse
    print("%d %d %0.2f" % (p.contents.targetId, p.contents.productSize, p.contents.freeEnergy))
    print("%d %d %d %d %0.2f %0.2f %0.2f %0.2f %0.2f" % (forward.contents.length, 
        forward.contents.position, 
	forward.contents.forward,
	forward.contents.reverse,
	forward.contents.GCContent,
	forward.contents.temperature,
	forward.contents.forwardElongationEfficiency[0],
	forward.contents.reverseElongationEfficiency[0],
	forward.contents.freeEnergy))
    print("%d %d %d %d %0.2f %0.2f %0.2f %0.2f %0.2f" % (reverse.contents.length,
        reverse.contents.position, 
	reverse.contents.forward,
	reverse.contents.reverse,
	reverse.contents.GCContent,
	reverse.contents.temperature,
	reverse.contents.forwardElongationEfficiency[0],
	reverse.contents.reverseElongationEfficiency[0],
	reverse.contents.freeEnergy))
    p = p.contents.next


