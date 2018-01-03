#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

def test_sample(sample):
    with open(sample('test.txt'), 'r') as f:
        assert f.read() == 'Test text\n'
