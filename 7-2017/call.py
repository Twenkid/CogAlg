from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

#@with PyCallGraph
with PyCallGraph(output=GraphvizOutput()):
    f = open(r"D:\Py\bk\1.py", "rt").read()
    print(f)
    PyCallGraph(f)