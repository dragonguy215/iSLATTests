"""Reusable mixins for the DataTypes module.

Each mixin addresses a single cross-cutting concern that is duplicated
across multiple DataType classes.  All mixins use ``__slots__ = ()`` so
they can be composed with ``__slots__``-based concrete classes — the
concrete class **must** declare the slots that the mixin reads/writes
(documented in each mixin's class docstring).

Usage example::

    class MoleculeLineList(CacheStatsMixin, WavelengthRangeMixin):
        __slots__ = (
            ...
            # Slots required by CacheStatsMixin
            '_cache_stats',
            # Slots required by WavelengthRangeMixin
            '_wavelength_range',
            ...
        )
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional


# ═══════════════════════════════════════════════════════════════════
#  CacheStatsMixin
# ═══════════════════════════════════════════════════════════════════

class CacheStatsMixin:
    """Standardised cache hit/miss/invalidation bookkeeping.

    Required slots on the concrete class
    -------------------------------------
    ``_cache_stats`` : dict
        Initialise to ``{'hits': 0, 'misses': 0, 'invalidations': 0}``
        in the concrete ``__init__``.

    Public API provided
    -------------------
    * ``get_cache_stats() -> dict``
    * ``_record_cache_hit()``
    * ``_record_cache_miss()``
    * ``_record_cache_invalidation()``
    * ``_reset_cache_stats()``
    """

    __slots__ = ()

    # ------------------------------------------------------------------
    # Concrete implementations
    # ------------------------------------------------------------------

    def _init_cache_stats(self) -> None:
        """Initialise the ``_cache_stats`` dict.

        Call this from the concrete ``__init__`` instead of manually
        writing the dict literal.
        """
        self._cache_stats: Dict[str, int] = {  # type: ignore[attr-defined]
            "hits": 0,
            "misses": 0,
            "invalidations": 0,
        }

    def get_cache_stats(self) -> Dict[str, int]:
        """Return a *copy* of the cache performance counters."""
        return dict(self._cache_stats)  # type: ignore[attr-defined]

    def _record_cache_hit(self) -> None:
        self._cache_stats["hits"] += 1  # type: ignore[attr-defined]

    def _record_cache_miss(self) -> None:
        self._cache_stats["misses"] += 1  # type: ignore[attr-defined]

    def _record_cache_invalidation(self) -> None:
        self._cache_stats["invalidations"] += 1  # type: ignore[attr-defined]

    def _reset_cache_stats(self) -> None:
        self._cache_stats.update(hits=0, misses=0, invalidations=0)  # type: ignore[attr-defined]


# ═══════════════════════════════════════════════════════════════════
#  WavelengthRangeMixin
# ═══════════════════════════════════════════════════════════════════

class WavelengthRangeMixin:
    """Getter/setter for a ``wavelength_range`` property with an
    invalidation hook.

    Required slots on the concrete class
    -------------------------------------
    ``_wavelength_range`` : tuple | None

    Hook to override
    ----------------
    ``_on_wavelength_range_changed(old, new)``
        Called **after** ``_wavelength_range`` has been updated but
        **only** when the new value differs from the old one.  The
        default implementation is a no-op.
    """

    __slots__ = ()

    @property
    def wavelength_range(self) -> Optional[tuple]:
        """Active wavelength range ``(lam_min, lam_max)`` in microns, or *None*."""
        return self._wavelength_range  # type: ignore[attr-defined]

    @wavelength_range.setter
    def wavelength_range(self, value: Optional[tuple]) -> None:
        old = self._wavelength_range  # type: ignore[attr-defined]
        if value != old:
            self._wavelength_range = value  # type: ignore[attr-defined]
            self._on_wavelength_range_changed(old, value)

    # Override point -------------------------------------------------
    def _on_wavelength_range_changed(
        self, old: Optional[tuple], new: Optional[tuple]
    ) -> None:
        """React to a wavelength-range change.  Override in subclasses."""


# ═══════════════════════════════════════════════════════════════════
#  ObservableMixin
# ═══════════════════════════════════════════════════════════════════

class ObservableMixin:
    """Simple observer pattern — instance-level callback list.

    Required slots on the concrete class
    -------------------------------------
    ``_callbacks`` : list
        Initialise to ``[]`` in the concrete ``__init__``.

    Public API
    ----------
    * ``add_callback(cb)``
    * ``remove_callback(cb)``
    * ``notify_callbacks(*args, **kwargs)``
    """

    __slots__ = ()

    def _init_callbacks(self) -> None:
        """Initialise the callback list.  Call from concrete ``__init__``."""
        self._callbacks: List[Callable] = []  # type: ignore[attr-defined]

    def add_callback(self, callback: Callable) -> None:
        """Register a callback.  Duplicates are silently ignored."""
        if callback not in self._callbacks:  # type: ignore[attr-defined]
            self._callbacks.append(callback)  # type: ignore[attr-defined]

    def remove_callback(self, callback: Callable) -> None:
        """Unregister a callback.  Missing callbacks are silently ignored."""
        try:
            self._callbacks.remove(callback)  # type: ignore[attr-defined]
        except ValueError:
            pass

    def notify_callbacks(self, *args: Any, **kwargs: Any) -> None:
        """Invoke every registered callback with the given arguments.

        Exceptions raised by a callback are printed but do **not**
        prevent the remaining callbacks from being called.
        """
        for cb in list(self._callbacks):  # type: ignore[attr-defined]
            try:
                cb(*args, **kwargs)
            except Exception as exc:
                print(f"Error in callback {cb!r}: {exc}")


# ═══════════════════════════════════════════════════════════════════
#  ClassObservableMixin
# ═══════════════════════════════════════════════════════════════════

class ClassObservableMixin:
    """Observer pattern at the **class** level (shared by all instances).

    Unlike :class:`ObservableMixin`, the callback list is stored on the
    **class** object itself.  This matches the existing pattern in
    :class:`Molecule` where ``_molecule_parameter_change_callbacks`` is
    a class-level list.

    Required class-level attributes
    --------------------------------
    ``_class_callbacks`` : list
        Declare as ``_class_callbacks: list = []`` in the class body.

    Public API (all classmethods)
    ----------------------------
    * ``add_class_callback(cb)``
    * ``remove_class_callback(cb)``
    * ``_notify_class_callbacks(*args, **kwargs)``
    """

    __slots__ = ()

    @classmethod
    def add_class_callback(cls, callback: Callable) -> None:
        """Register a class-level callback."""
        if callback not in cls._class_callbacks:  # type: ignore[attr-defined]
            cls._class_callbacks.append(callback)  # type: ignore[attr-defined]

    @classmethod
    def remove_class_callback(cls, callback: Callable) -> None:
        """Unregister a class-level callback."""
        try:
            cls._class_callbacks.remove(callback)  # type: ignore[attr-defined]
        except ValueError:
            pass

    @classmethod
    def _notify_class_callbacks(cls, *args: Any, **kwargs: Any) -> None:
        """Invoke every class-level callback."""
        for cb in list(cls._class_callbacks):  # type: ignore[attr-defined]
            try:
                cb(*args, **kwargs)
            except Exception as exc:
                print(f"Error in class callback {cb!r}: {exc}")
