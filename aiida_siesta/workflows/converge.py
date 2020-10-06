from ..utils.converge_absclass import BasicConverger, SequentialConverger
from .iterate import SiestaIterator


class SiestaConverger(BasicConverger, SiestaIterator):
    pass


class SiestaSequentialConverger(SequentialConverger):
    _process_class = SiestaConverger
