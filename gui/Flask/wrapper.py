import ctypes

class Primer(ctypes.Structure):
    pass

Primer._fields_ = [("target_id", ctypes.c_char_p), 
		("sequence", ctypes.c_char_p),
		("length", ctypes.c_int),
		("position", ctypes.c_int),
		("forward", ctypes.c_bool),
		("reverse", ctypes.c_bool),
		("gc_content", ctypes.c_float),
		("temperature", ctypes.c_float),
		#("free_energy", ctypes.c_float),
		("forward_elongation_efficiency", ctypes.POINTER(ctypes.c_float)),
		("reverse_elongation_efficiency", ctypes.POINTER(ctypes.c_float)),
		("next", ctypes.POINTER(Primer))]

class Pair(ctypes.Structure):
    pass

Pair._fields_ = [("target_id", ctypes.c_char_p), 
		("product_size", ctypes.c_int),
		#("free_energy", ctypes.c_float),
		("forward", ctypes.POINTER(Primer)),
		("reverse", ctypes.POINTER(Primer)),
		("next", ctypes.POINTER(Pair))]

class Target(ctypes.Structure):
	_fields_ = [("id", ctypes.c_char_p),
		("sequence", ctypes.c_char_p)]

class Multigene:
	def __init__(self):
		#self.lib = ctypes.CDLL("./libteste.so")
		self.lib = ctypes.CDLL("./multigene.so")
		#self.lib.design.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_int]
		self.lib.design.argtypes = [ctypes.POINTER(Target), ctypes.c_int]
		self.lib.design.restype = ctypes.POINTER(Pair)

	#def design(self, targets):
	#	length = len(targets)
	#	_targets = (ctypes.c_char_p * length)()
	#	for i in range(length):
	#		_targets[i] = targets[i].encode('utf-8')
	#	
	#	return self.lib.design(_targets, length)

	def design(self, targets):
		length = len(targets)
		_targets = (Target * length)()
		for i in range(length):
			_id = str(i).encode('utf-8')
			_sequence = targets[i].encode('utf-8')
			_targets[i] = Target(_id, _sequence)
			#print(_targets[i].id, _targets[i].sequence)

		return self.lib.design(_targets, length)