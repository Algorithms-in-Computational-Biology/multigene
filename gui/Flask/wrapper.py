import ctypes

#lib = ctypes.cdll.LoadLibrary("./libteste.so")

#lib.getK.argtypes = [ctypes.c_float, ctypes.c_float]
#lib.getK.restype = ctypes.c_double
#k = lib.getK(28.5, 154.6)
#print(k)

class Primer(ctypes.Structure):
    pass

Primer._fields_ = [("targetId", ctypes.c_int), 
		("sequence", ctypes.c_char_p),
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
        #self.lib = ctypes.CDLL("./libteste.so")
        self.lib = ctypes.CDLL("./multigene.so")
        self.lib.design.argtypes = [ctypes.POINTER(ctypes.c_char_p), ctypes.c_int]
        self.lib.design.restype = ctypes.POINTER(Pair)

    def design(self, targets):
        length = len(targets)
        _targets = (ctypes.c_char_p * length)()
        for i in range(length):
            _targets[i] = targets[i].encode('utf-8')

        return self.lib.design(_targets, length)
