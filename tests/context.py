import sys
import os
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(1, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '/data')))
from src import dssp
from src import common
