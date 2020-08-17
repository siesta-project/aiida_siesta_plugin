from .converge_absclass import BasicConverger, SequentialConverger
from .future_iterate import SiestaIterator


class SiestaConverger(BasicConverger, SiestaIterator):
    pass


class SiestaSequentialConverger(SequentialConverger):
    _process_class = SiestaConverger
