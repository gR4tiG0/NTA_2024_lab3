#!/usr/bin/env python3
from dlpModule import DLPSolver




def main() -> None:
    a,b,p = 185482,194035,376039 #6
    # a,b,p = 7099845, 28363162, 39370057 #8
    # a,b,p = 272914025769, 139911379996, 454327134073 #12
    #a,b,p = 2276694821408773, 6624933408925346, 8532352737210737 #16
    #a,b,p = 9591297096748371381, 4648727071064622208, 13322254769939007437 #64bit
    s = DLPSolver()
    print(a,b,p)
    res = s.indexCalculus(a,b,p)
    print(res)
    print("pow =",pow(a,res,p))

if __name__ == "__main__":
    main()