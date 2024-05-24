#!/usr/bin/env python3
from dlpModule import DLPSolver
import time

def main() -> None:
    s = DLPSolver()

    try:
        a = int(input('> a = '))
        b = int(input('> b = '))
        p = int(input('> p = '))

    except Exception as e:
       print(e)
    start = time.time()
    res = s.indexCalculus(a, b, p)
    print(f"Time taken: {time.time() - start}s")
    assert pow(a, res, p) == b, "Result not OK"
    print(f"Result: {res}")

if __name__ == "__main__":
    main()
