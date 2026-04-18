"""Shared lazy-import helper for pandas.

Every DataType module that needs pandas used to carry its own copy of the
``_get_pandas`` function.  This module consolidates them into a single
canonical implementation so the lazy-import logic lives in exactly one
place.

Usage::

    from ._pandas_import import get_pandas

    def some_method(self):
        pd = get_pandas()
        return pd.DataFrame(...)
"""

from __future__ import annotations

_pd_module = None
_pd_imported = False


def get_pandas():
    """Return the ``pandas`` module, importing it lazily on first call.

    Raises
    ------
    ImportError
        If pandas is not installed.
    """
    global _pd_module, _pd_imported
    if not _pd_imported:
        try:
            import pandas as _pd
            _pd_module = _pd
        except ImportError:
            raise ImportError(
                "pandas is required for table functionality. "
                "Install it with: pip install pandas"
            )
        finally:
            _pd_imported = True
    if _pd_module is None:
        raise ImportError(
            "pandas is required for table functionality. "
            "Install it with: pip install pandas"
        )
    return _pd_module
