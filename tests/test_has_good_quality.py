import pytest

from src.filtering import has_good_quality


def test_has_good_quality_basic():
    row = {
        'read': 'ACGT',
        'reference': 'ACGT',
        'modules': '0123',
        'quality': 'IIII',  # Quality values of 40
    }
    assert has_good_quality(row, 30, 0, 3) is True


def test_has_good_quality_below_minimum():
    row = {
        'read': 'ACGT',
        'reference': 'ACGT',
        'modules': '0123',
        'quality': '!!!!',  # Quality values of 0
    }
    assert has_good_quality(row, 30, 0, 3) is False


def test_has_good_quality_no_annotated_part():
    row = {
        'read': 'ACGT',
        'reference': 'ACGT',
        'modules': '----',
        'quality': 'IIII',
    }
    assert has_good_quality(row, 30, 1, 2) is True


def test_has_good_quality_with_insertions():
    row = {
        'read': 'A_CGT',
        'reference': 'AC_GT',
        'modules': '01234',
        'quality': 'IIII',
    }
    assert has_good_quality(row, 30, 1, 3) is True


def test_has_good_quality_edge_case_start():
    row = {
        'read': 'ACGT',
        'reference': 'ACGT',
        'modules': '0123',
        'quality': 'I!!!',
    }
    assert has_good_quality(row, 30, 0, 1) is False


def test_has_good_quality_edge_case_end():
    row = {
        'read': 'ACGT',
        'reference': 'ACGT',
        'modules': '0123',
        'quality': '!!!I',
    }
    assert has_good_quality(row, 30, 2, 3) is False


def test_has_good_quality_standard_case():
    row = {
        'read': 'ACGTACGTACAGT',
        'reference': 'ACGTACGTACAGT',
        'modules': '-0111222333--',
        'quality': '!!IIIIIII!!!!',
    }
    assert has_good_quality(row, 1, 0, 3) is False
    assert has_good_quality(row, 1, 1, 3) is False
    assert has_good_quality(row, 0, 1, 3) is True
    assert has_good_quality(row, 1, 1, 2) is True
    assert has_good_quality(row, 1, 0, 2) is False
