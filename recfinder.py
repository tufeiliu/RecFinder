#!/usr/bin/env python3
import sys
import os

recfinder_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "src"))
sys.path.insert(0, recfinder_dir)

if __name__ == "__main__":
    from recfinder.run.__main__ import main
    args = sys.argv[1:]
    main(args)