"""Layer A random-number bridge: conventions, accounting, replay and offsets."""
import hashlib
import struct

import numpy as np
import pytest

from rnglab.nuclear import NumpyPCG64Source, WordStream, labelled_source, words_to_uniforms
from rnglab.seeds import namespace_seed
from rnglab14.sources import Xoshiro256PlusPlusSource

TEST_NAMESPACE = "rnglab-v020-nuclear-test"


def test_uniform_convention_matches_random_source_random_float():
    reference = Xoshiro256PlusPlusSource(seed=namespace_seed(TEST_NAMESPACE, "streams", "float"))
    expected = [reference.random_float() for _ in range(64)]
    stream = WordStream(Xoshiro256PlusPlusSource(seed=namespace_seed(TEST_NAMESPACE, "streams", "float")))
    assert stream.uniforms(64).tolist() == expected


def test_uniform_extremes():
    words = np.array([0, 2**11 - 1, 2**11, 2**64 - 1], dtype=np.uint64)
    u = words_to_uniforms(words)
    assert u[0] == 0.0 and u[1] == 0.0 and u[2] == 2.0**-53 and u[3] == 1.0 - 2.0**-53


def test_word_stream_accounts_and_hashes_exactly_what_it_hands_out():
    source = labelled_source(TEST_NAMESPACE, "streams", "accounting")
    stream = WordStream(source)
    a, b = stream.raw(3), stream.raw(5)
    assert stream.words_read == 8 and stream.bytes_read == 64 == source.consumed_bytes
    assert source.stream_digest == hashlib.sha256(a + b).hexdigest()
    words = np.frombuffer(a + b, dtype=">u8")
    assert stream.provenance()["source_stream_sha256"] == source.stream_digest
    assert words.size == 8


def test_numpy_pcg64_source_is_the_raw_pcg64_stream():
    seed = namespace_seed(TEST_NAMESPACE, "streams", "raw")
    expected = np.random.PCG64(np.random.SeedSequence(seed)).random_raw(10)
    stream = WordStream(NumpyPCG64Source(seed))
    assert np.array_equal(stream.words(10), expected)


def test_split_reads_give_the_same_byte_stream():
    seed = namespace_seed(TEST_NAMESPACE, "streams", "split")
    whole = NumpyPCG64Source(seed).read_bytes(80)
    source = NumpyPCG64Source(seed)
    pieces = b"".join(source.read_bytes(n) for n in (3, 13, 1, 8, 55))
    assert pieces == whole and source.consumed_bytes == 80


def test_start_word_jumps_without_reading():
    seed = namespace_seed(TEST_NAMESPACE, "streams", "offset")
    full = WordStream(NumpyPCG64Source(seed)).words(12)
    jumped_source = NumpyPCG64Source(seed, start_word=7)
    jumped = WordStream(jumped_source).words(5)
    assert np.array_equal(jumped, full[7:])
    assert jumped_source.consumed_bytes == 40
    assert jumped_source.manifest.metadata["start_word"] == 7


def test_labels_select_streams():
    a1 = WordStream(labelled_source(TEST_NAMESPACE, "streams", "a")).words(4)
    a2 = WordStream(labelled_source(TEST_NAMESPACE, "streams", "a")).words(4)
    b = WordStream(labelled_source(TEST_NAMESPACE, "streams", "b")).words(4)
    assert np.array_equal(a1, a2) and not np.array_equal(a1, b)
    assert labelled_source(TEST_NAMESPACE, "streams", "a").manifest.metadata["seed_label"] == \
        f"{TEST_NAMESPACE}/streams/a"


@pytest.mark.parametrize("bad", [-1, 1.5, True, "3"])
def test_word_counts_are_validated(bad):
    with pytest.raises(ValueError):
        WordStream(labelled_source(TEST_NAMESPACE, "streams", "bad")).raw(bad)


def test_bad_sources_and_seeds_are_rejected():
    with pytest.raises(TypeError):
        WordStream(object())
    with pytest.raises(ValueError):
        NumpyPCG64Source(-5)
    with pytest.raises(ValueError):
        NumpyPCG64Source(5, start_word=-1)


def test_word_decoding_is_big_endian():
    raw = struct.pack(">Q", 0x0123456789ABCDEF)
    assert int(np.frombuffer(raw, dtype=">u8")[0]) == 0x0123456789ABCDEF
