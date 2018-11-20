import pytest
from jdna import DoubleLinkedList

@pytest.fixture(scope="function")
def test_str():
    return '^This is a test string for the linked data set.'


@pytest.fixture(scope="function")
def linked_list(test_str):
    return DoubleLinkedList(data=test_str)