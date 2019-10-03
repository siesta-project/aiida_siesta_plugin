from __future__ import absolute_import
from __future__ import print_function
from aiida.orm import Int
from aiida.engine import workfunction as wf
from aiida.cmdline.utils import decorators


# Define the workfunction
@wf
@decorators.with_dbenv()
def sum(a, b):
  return a + b

# Run it with some input
r = sum(Int(4), Int(5))
print(r)
