from aiida import load_dbenv
load_dbenv()

from aiida.orm.nodes.base import Int
from aiida.engine.workfunction import workfunction as wf


# Define the workfunction
@wf
def sum(a, b):
  return a + b

# Run it with some input
r = sum(Int(4), Int(5))
print(r)

                                    
