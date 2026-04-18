"""Tests for the DataTypes mixins."""

import pytest
from iSLAT.Modules.DataTypes._mixins import (
    CacheStatsMixin,
    WavelengthRangeMixin,
    ObservableMixin,
    ClassObservableMixin,
)


# ═══════════════════════════════════════════════════════════════════
#  Concrete test classes that compose the mixins with __slots__
# ═══════════════════════════════════════════════════════════════════

class CacheTester(CacheStatsMixin):
    __slots__ = ("_cache_stats",)

    def __init__(self):
        self._init_cache_stats()


class WaveTester(WavelengthRangeMixin):
    __slots__ = ("_wavelength_range", "_invalidated")

    def __init__(self, wr=None):
        self._wavelength_range = wr
        self._invalidated = False

    def _on_wavelength_range_changed(self, old, new):
        self._invalidated = True


class ObservableTester(ObservableMixin):
    __slots__ = ("_callbacks",)

    def __init__(self):
        self._init_callbacks()


class ClassObservableTester(ClassObservableMixin):
    __slots__ = ()
    _class_callbacks: list = []


# ═══════════════════════════════════════════════════════════════════
#  CacheStatsMixin tests
# ═══════════════════════════════════════════════════════════════════

class TestCacheStatsMixin:
    def test_initial_stats(self):
        c = CacheTester()
        assert c.get_cache_stats() == {"hits": 0, "misses": 0, "invalidations": 0}

    def test_record_hit(self):
        c = CacheTester()
        c._record_cache_hit()
        c._record_cache_hit()
        assert c.get_cache_stats()["hits"] == 2

    def test_record_miss(self):
        c = CacheTester()
        c._record_cache_miss()
        assert c.get_cache_stats()["misses"] == 1

    def test_record_invalidation(self):
        c = CacheTester()
        c._record_cache_invalidation()
        assert c.get_cache_stats()["invalidations"] == 1

    def test_get_cache_stats_returns_copy(self):
        c = CacheTester()
        stats = c.get_cache_stats()
        stats["hits"] = 999
        assert c.get_cache_stats()["hits"] == 0

    def test_reset_stats(self):
        c = CacheTester()
        c._record_cache_hit()
        c._record_cache_miss()
        c._record_cache_invalidation()
        c._reset_cache_stats()
        assert c.get_cache_stats() == {"hits": 0, "misses": 0, "invalidations": 0}


# ═══════════════════════════════════════════════════════════════════
#  WavelengthRangeMixin tests
# ═══════════════════════════════════════════════════════════════════

class TestWavelengthRangeMixin:
    def test_initial_none(self):
        w = WaveTester()
        assert w.wavelength_range is None

    def test_initial_value(self):
        w = WaveTester((5.0, 8.0))
        assert w.wavelength_range == (5.0, 8.0)

    def test_setter_triggers_hook(self):
        w = WaveTester((5.0, 8.0))
        w.wavelength_range = (6.0, 10.0)
        assert w._invalidated is True
        assert w.wavelength_range == (6.0, 10.0)

    def test_setter_noop_on_same_value(self):
        w = WaveTester((5.0, 8.0))
        w.wavelength_range = (5.0, 8.0)
        assert w._invalidated is False

    def test_setter_to_none(self):
        w = WaveTester((5.0, 8.0))
        w.wavelength_range = None
        assert w._invalidated is True
        assert w.wavelength_range is None


# ═══════════════════════════════════════════════════════════════════
#  ObservableMixin tests
# ═══════════════════════════════════════════════════════════════════

class TestObservableMixin:
    def test_add_and_notify(self):
        o = ObservableTester()
        results = []
        o.add_callback(lambda x: results.append(x))
        o.notify_callbacks(42)
        assert results == [42]

    def test_remove_callback(self):
        o = ObservableTester()
        results = []
        cb = lambda x: results.append(x)
        o.add_callback(cb)
        o.remove_callback(cb)
        o.notify_callbacks(42)
        assert results == []

    def test_remove_missing_is_silent(self):
        o = ObservableTester()
        o.remove_callback(lambda: None)  # Should not raise

    def test_duplicate_add_ignored(self):
        o = ObservableTester()
        results = []
        cb = lambda: results.append(1)
        o.add_callback(cb)
        o.add_callback(cb)
        o.notify_callbacks()
        assert results == [1]  # Only called once

    def test_callback_exception_does_not_halt(self):
        o = ObservableTester()
        results = []

        def bad_cb():
            raise RuntimeError("boom")

        o.add_callback(bad_cb)
        o.add_callback(lambda: results.append("ok"))
        o.notify_callbacks()
        assert results == ["ok"]


# ═══════════════════════════════════════════════════════════════════
#  ClassObservableMixin tests
# ═══════════════════════════════════════════════════════════════════

class TestClassObservableMixin:
    def setup_method(self):
        ClassObservableTester._class_callbacks.clear()

    def test_add_and_notify(self):
        results = []
        ClassObservableTester.add_class_callback(lambda x: results.append(x))
        ClassObservableTester._notify_class_callbacks(99)
        assert results == [99]

    def test_remove(self):
        results = []
        cb = lambda x: results.append(x)
        ClassObservableTester.add_class_callback(cb)
        ClassObservableTester.remove_class_callback(cb)
        ClassObservableTester._notify_class_callbacks(99)
        assert results == []

    def test_shared_across_instances(self):
        results = []
        ClassObservableTester.add_class_callback(lambda x: results.append(x))
        a = ClassObservableTester()
        b = ClassObservableTester()
        # Both instances see the same class-level callbacks
        a._notify_class_callbacks("from_a")
        b._notify_class_callbacks("from_b")
        assert results == ["from_a", "from_b"]
