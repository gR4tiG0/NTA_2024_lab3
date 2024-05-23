from ctypes import CDLL, c_uint64
import time

VERBOSE = True

class DLPSolver:
    def __init__(self):
        self.lib_name = "./libdlpT.so"
        self.lib = self.initLib()

    def initLib(self) -> CDLL:
        lib = CDLL(self.lib_name)

        lib.power.argtypes = [c_uint64, c_uint64, c_uint64]
        lib.power.restype = c_uint64

        lib.inv.argtypes = [c_uint64, c_uint64]
        lib.inv.restype = c_uint64

        lib.factor.argtypes = [c_uint64]
        lib.factor.restype = c_uint64

        lib.isPrime.argtypes = [c_uint64]
        lib.isPrime.restype = c_uint64

        lib.factorizeAndPrint.argtypes = [c_uint64]
        lib.factorizeAndPrint.restype = None

        lib.indexCalculus.argtypes = [c_uint64, c_uint64, c_uint64]
        lib.indexCalculus.restype = c_uint64

        return lib
    
    def indexCalculus(self, a: int, b: int, n: int) -> int:
        return self.lib.indexCalculus(a, b, n)