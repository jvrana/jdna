import pytest

from jdna.sequence import Sequence


@pytest.fixture(scope="function")
def test_str():
    return "^This is just some test string for testing."


@pytest.fixture(scope="function")
def test_seq(test_str):
    return Sequence(test_str)
