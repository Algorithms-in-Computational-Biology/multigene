import ctypes

#lib = ctypes.cdll.LoadLibrary("./libteste.so")

#lib.getK.argtypes = [ctypes.c_float, ctypes.c_float]
#lib.getK.restype = ctypes.c_double
#k = lib.getK(28.5, 154.6)
#print(k)

class Primer(ctypes.Structure):
	pass

Primer._fields_ = [("targetId", ctypes.c_int), 
		("sequence", ctypes.POINTER(ctypes.c_char)),
		("length", ctypes.c_int),
		("position", ctypes.c_int),
		("forward", ctypes.c_bool),
		("reverse", ctypes.c_bool),
		("GCContent", ctypes.c_float),
		("temperature", ctypes.c_float),
		("freeEnergy", ctypes.c_float),
		("forwardElongationEfficiency", ctypes.POINTER(ctypes.c_float)),
		("reverseElongationEfficiency", ctypes.POINTER(ctypes.c_float)),
		("next", ctypes.POINTER(Primer))]


class Pair(ctypes.Structure):
	pass

Pair._fields_ = [("targetId", ctypes.c_int), 
		("productSize", ctypes.c_int),
		("freeEnergy", ctypes.c_float),
		("forward", ctypes.POINTER(Primer)),
		("reverse", ctypes.POINTER(Primer)),
		("next", ctypes.POINTER(Pair))]

class Multigene:
	def __init__(self):
		self.lib = ctypes.CDLL("./libteste.so")
		self.lib.design.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_int]
		self.lib.design.restype = ctypes.POINTER(Pair)

	def design(self, targets):
		size = len(targets)
		_targets = (ctypes.c_char_p * size)()
		_targets[:] = targets

		return self.lib.design(_targets, size - 1)


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



