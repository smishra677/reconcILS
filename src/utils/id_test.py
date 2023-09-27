import unittest
import uuid
import itertools

class idmaker1():
    def get_id(self):
        return uuid.uuid1()

class idmaker2():
    newid = itertools.count().__next__
    def __init__(self):
        self.id = idmaker2.newid()

print(idmaker2().id)

class idtests(unittest.TestCase):
    def test_idmaker1_returns_unique(self):
        ids = set(idmaker1().get_id() for i in range(1000))
        self.assertEqual(1000, len(ids))

    def test_idmaker2_returns_unique(self):
        ids = set(idmaker2().id for i in range(1000))
        self.assertEqual(1000, len(ids))