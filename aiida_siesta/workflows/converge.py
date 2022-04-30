# -*- coding: utf-8 -*-
"""
Cllect the convergers.
"""
from ..utils.converge_absclass import BasicConverger, SequentialConverger
from .iterate import SiestaIterator


class SiestaConverger(BasicConverger, SiestaIterator):
    """
    Only connects the two parent classes.
    """


class SiestaSequentialConverger(SequentialConverger):
    """
    An iterator of convergers.
    """

    _process_class = SiestaConverger
