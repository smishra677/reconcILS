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
