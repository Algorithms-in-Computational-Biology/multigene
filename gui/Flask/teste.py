from wrapper import Multigene

#targets = ["ATGGCGGTGGCTTCGACCTCGCCGCTATCCGCCACGGCCCCCTCGCCGCCCGCTCCGGTGTCCGGGTTCCTCGCTCTCCCCGCCCGCCGCGGCTGCGCAACGCGCCTCGGCTCCGCCGCCGCGTGGAGGAGGCTTCGCGTGGAGGCGATCTGGAAGCAGCAGGAGAAGCAGCGGGCAGAGGTGTCCGTCGAGGAACCCGCCCCCGTCAGGGAGGCCGCCGCGCCCCTGGACGGAGTCGGAGCTGACGACCCCATGGTTCCTTCCTCGGACGAGAGCTGGGTGGTCAGGCTCGAGCAGTCGGTCAACATTTTCCTCACGGAATCGGTGATTATACTACTCAATACCGTGTACCGTGATCGGAACTACGCCAGGTTTTTTGTGCTGGAGACGATTGCCAGGGTGCCGTATTTCGCGTTCATATCGGTGCTTCACATGTATGAAACCTTTGGCTGGTGGAGACGAGCTGATTATCTAAAAGTTCACTTTGCGCAGAGCTTGAACGAGTTTCATCATCTCTTGATCATGGAAGAATTGGGTGGCAACGCTATATGGATTGATTGTTTCCTTGCTCGATTTATGGCGTTTTTTTACTACTTCATGACTGTTGCGATGTACATGTTGAGCCCACGAATGGCATATCACTTCTCTGAATGTGTGGAGAGACATGCGTACTCCACCTATGATAAGTTCCTCAAGCTCCATGAAGAGGAATTGAAAACACTACCAGCTCCAGAGGCAGCATTGAACTATTACCTGAATGAGGACCTTTACTTATTTGATGAGTTTCAGACAACAAGAATTCCATGTTCTAGGAGGCCTAAAATAGATAACTTGTATGATGTATTCGTCAATATACGAGATGACGAGGCAGAGCACTGCAAGACAATGAAGGCATGTCAAACACATGGAACTCTTCGTTCTCCTCACTCAATGCCGAACTGCTTAGAAGCTGCTACAGAATGTGTAATACCTGAAAACGATTGTGAAGGTATTGTGGACTGTGTCAAAAAGTCCCTTACAAAGTAA", 
# 			"ATGGCGGTGGCCTCGACCTCGCCGCTGTCCGCCAAGCCCGCCACGGCCCCTTCGCCGCCCGCTCCCGGATCCGGGCTCCTCGCTCTCGGCGTTCGCCGCGCCCCCGCCACTGCCGCGTGGAGGAGGCTCCGCGTGGAGGCGATCAGGACGCAGCGAACGGAGGTGCCCGTCGAGGAGTCCGCCCCCGCCAGGGACGCCGCCGCTGCCGCGCCCCTGGACGGAAACGGAGCCGGAGCGGACGGCTCCGTGGTTCCTTCCTCGGACGACAGCTGGGTTGTCAAGCTCGAGCAGTCGTTCAACATTTTCGCCACGGATTCGGTGATTATGGTACTCAAGGGCGTGTACGGTGATCGGTACTACGCCAGGTTCTTTGCGCTGGAGACGATTGCGAGGGTGCCGTACTTCGCATTCATATCGGTGCTTCACTTGTATGCGACCTTTGGATGGTGGAGACGAGCTGATTACATAAAGGTTCACTTTGCGCAGAGCTGGAACGAGTTCCATCACCTCTTGATCATGGAAGAATTGGGTGGCGACTCTTTGTGGTTTGACTGTTTTCTTGCTCGGTTTATGGCATTCTTTTACTACTTCATGACTGTTGCAATGTACATGCTGAGCCCACGAATGGCATATCACTTTTCCGAATGTGTGGAGAGACATGCATATTCCACCTATGATGAGTTCCTCAAGCTCCATGAAGAGGAATTGAAAAGACTACCAGCTCCAGAGGCAGCATTGAACTATTACATGAATGAGGACCTTTACTTATTCGATGAGTTTCAGGCATCAAGAACTCCAGGTTCTAGGAGGCCTAAAATAGATAACTTATACGATGTATTCGTTAATATACGAGAAGATGAGGCAGAGCACTGCAAGACAATGAAGACCTGTCAAACACATGGAAATCTTCGTTCTCCTCATTCAACGCCGAACTGCTTAGAAGATGATACGGAATGTGTAATACCTGAAAACGACTGTGAAGGTATTGTGGACTGTGTCAAAAAGTCCCTTACAAAGTAA"]
#multigene = Multigene()

#_targets = multigene.design(targets)
#for i in range(len(_targets)):
#	print(_targets[i].id, _targets[i].sequence)

#target = Target("1".encode('utf-8'), "ACTG".encode('utf-8'))
#print(target.id, target.sequence)

#targets = ['AC', 'BD']
#multigene = Multigene()
#pair = multigene.design(targets) 

#print(pair)
#p = pair
#while(p): 
#    forward = p.contents.forward
#    reverse = p.contents.reverse
#    print("%s %d" % (p.contents.target_id, p.contents.product_size)) #  %0.2f, p.contents.free_energy
#    print("%s %d %d %d %d %0.2f %0.2f %0.2f %0.2f" % (forward.contents.sequence,
#       forward.contents.length, 
#        forward.contents.position, 
#	forward.contents.forward,
#	forward.contents.reverse,
#	forward.contents.gc_content,
#	forward.contents.temperature,
#	forward.contents.forward_elongation_efficiency[0],
#	forward.contents.reverse_elongation_efficiency[0])) #  %0.2f, forward.contents.freeEnergy
#    print("%s %d %d %d %d %0.2f %0.2f %0.2f %0.2f" % (reverse.contents.sequence,
#        reverse.contents.length,
#        reverse.contents.position, 
#	reverse.contents.forward,
#	reverse.contents.reverse,
#	reverse.contents.gc_content,
#	reverse.contents.temperature,
#	reverse.contents.forward_elongation_efficiency[0],
#	reverse.contents.reverse_elongation_efficiency[0])) #  %0.2f, reverse.contents.freeEnergy
#    p = p.contents.next

s = "ATCGACTAGCAGATCA"
min_len = 2
max_len = 4
min_prod = 5 
max_prod = 7

for i in range(min_len, max_len + 1):
	#print(i)
	for j in range(len(s)):
		lenght = len(s) - j
		size = max_len + max_prod + max_len
		if (size <= lenght):  
			print((j, i), s[j: j + i], s[j: j + i])
			#print(size, lenght)

