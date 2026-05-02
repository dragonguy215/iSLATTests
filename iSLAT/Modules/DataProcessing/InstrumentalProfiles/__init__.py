# -*- coding: utf-8 -*-
"""
Instrumental LSF profile plugin system for iSLAT.

Quick usage::

    from iSLAT.Modules.DataProcessing.InstrumentalProfiles import (
        PROFILE_REGISTRY,
        PROFILE_DISPLAY_NAMES,
        InstrumentalProfile,
        ConstantProfile,
        MiriMrsProfile,
    )
"""

from .profiles import (
    InstrumentalProfile,
    ConstantProfile,
    MiriMrsProfile,
    PROFILE_REGISTRY,
    PROFILE_DISPLAY_NAMES,
)

__all__ = [
    "InstrumentalProfile",
    "ConstantProfile",
    "MiriMrsProfile",
    "PROFILE_REGISTRY",
    "PROFILE_DISPLAY_NAMES",
]
