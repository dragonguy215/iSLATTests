# -*- coding: utf-8 -*-
"""
Tests for the multi-format line-list reader infrastructure.

Covers:
- ReadResult dataclass
- detect_format()
- CsvLineListReader
- SavedLinesReader
- HitranParReader protocol conformance
- MoleculeLineList.from_file
- MoleculeLineList.register_reader / reader registry
- Extra-fields round-trip through MoleculeLineList
- Binary cache with source_format + extra_columns (v3)
- LineListMaker.from_file / from_saved_lines / to_linelist extras
"""

import csv
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

import iSLAT.Constants as c
from iSLAT.Modules.FileHandling.line_list_readers import (
    CORE_FIELD_NAMES,
    CsvLineListReader,
    DEFAULT_READERS,
    LineListReader,
    ReadResult,
    SavedLinesReader,
    detect_format,
)


# ────────────────────────────────────────────────────────────────────
#  Fixtures
# ────────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp_dir():
    """Provide a temporary directory that is cleaned up after the test."""
    d = tempfile.mkdtemp()
    yield Path(d)
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture
def csv_linelist_path(tmp_dir):
    """Create a minimal CSV line-list file."""
    path = tmp_dir / "test_lines.csv"
    rows = [
        {"species": "H2O", "lev_up": "0|10", "lev_low": "0|9",
         "lam": "12.407", "a_stein": "1.05e-2", "e_up": "4586.4",
         "e_low": "3379.0", "g_up": "21", "g_low": "19",
         "xmin": "12.3", "xmax": "12.5"},
        {"species": "H2O", "lev_up": "0|8", "lev_low": "0|7",
         "lam": "14.950", "a_stein": "2.10e-2", "e_up": "3200.0",
         "e_low": "2500.0", "g_up": "17", "g_low": "15",
         "xmin": "14.8", "xmax": "15.1"},
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["species", "lev_up", "lev_low", "lam",
                         "a_stein", "e_up", "e_low", "g_up", "g_low",
                         "xmin", "xmax"],
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.fixture
def saved_lines_path(tmp_dir):
    """Create a minimal saved-lines CSV file."""
    path = tmp_dir / "test_saved.csv"
    fieldnames = [
        "species", "lev_up", "lev_low", "lam",
        "a_stein", "e_up", "e_low", "g_up", "g_low",
        "xmin", "xmax", "Flux_data", "Err_data", "Fit_SN",
        "Flux_fit", "FWHM_fit", "Centr_fit", "Red-chisq",
    ]
    rows = [
        {"species": "H2O", "lev_up": "0|10", "lev_low": "0|9",
         "lam": "12.407", "a_stein": "1.05e-2", "e_up": "4586.4",
         "e_low": "3379.0", "g_up": "21", "g_low": "19",
         "xmin": "12.3", "xmax": "12.5",
         "Flux_data": "1.23e-14", "Err_data": "5e-16",
         "Fit_SN": "24.6", "Flux_fit": "1.20e-14",
         "FWHM_fit": "0.005", "Centr_fit": "12.4071",
         "Red-chisq": "1.02"},
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


@pytest.fixture
def par_file_path(tmp_dir):
    """Create a minimal .par file (enough for detect_format)."""
    path = tmp_dir / "test_data.par"
    path.write_text("# HITRAN H2O; id:1; iso:1;gid:1\n# Molar Mass: 18.0\n")
    return path


# ────────────────────────────────────────────────────────────────────
#  ReadResult
# ────────────────────────────────────────────────────────────────────

class TestReadResult:
    def test_defaults(self):
        r = ReadResult(lines_data=[])
        assert r.lines_data == []
        assert r.partition is None
        assert r.metadata == {}
        assert r.extra_columns == {}

    def test_with_data(self):
        r = ReadResult(
            lines_data=[{"lam": 12.0}],
            metadata={"source_format": "csv"},
            extra_columns={"xmin": ["12.3"]},
        )
        assert len(r.lines_data) == 1
        assert r.metadata["source_format"] == "csv"
        assert "xmin" in r.extra_columns

    def test_frozen(self):
        r = ReadResult(lines_data=[])
        with pytest.raises(AttributeError):
            r.lines_data = [{"lam": 1.0}]


# ────────────────────────────────────────────────────────────────────
#  detect_format
# ────────────────────────────────────────────────────────────────────

class TestDetectFormat:
    def test_par_extension(self, par_file_path):
        assert detect_format(par_file_path) == "hitran"

    def test_csv_plain(self, csv_linelist_path):
        assert detect_format(csv_linelist_path) == "csv"

    def test_csv_saved(self, saved_lines_path):
        assert detect_format(saved_lines_path) == "saved"

    def test_unknown_extension_raises(self, tmp_dir):
        p = tmp_dir / "data.xyz"
        p.write_text("random content")
        with pytest.raises(ValueError, match="Cannot auto-detect"):
            detect_format(p)


# ────────────────────────────────────────────────────────────────────
#  CsvLineListReader
# ────────────────────────────────────────────────────────────────────

class TestCsvLineListReader:
    def test_protocol_conformance(self):
        assert isinstance(CsvLineListReader(), LineListReader)

    def test_read_basic(self, csv_linelist_path):
        reader = CsvLineListReader()
        result = reader.read(csv_linelist_path)

        assert isinstance(result, ReadResult)
        assert len(result.lines_data) == 2
        assert result.partition is None
        assert result.metadata["source_format"] == "csv"

    def test_core_fields_present(self, csv_linelist_path):
        reader = CsvLineListReader()
        result = reader.read(csv_linelist_path)
        line = result.lines_data[0]
        for field in CORE_FIELD_NAMES:
            assert field in line, f"Missing core field: {field}"

    def test_freq_computed_from_lam(self, csv_linelist_path):
        reader = CsvLineListReader()
        result = reader.read(csv_linelist_path)
        line = result.lines_data[0]
        expected_freq = c.SPEED_OF_LIGHT_MICRONS / 12.407
        assert abs(line["freq"] - expected_freq) < 1e6  # within 1 MHz

    def test_extra_columns_extracted(self, csv_linelist_path):
        reader = CsvLineListReader()
        result = reader.read(csv_linelist_path)
        assert "xmin" in result.extra_columns
        assert "xmax" in result.extra_columns
        assert len(result.extra_columns["xmin"]) == 2

    def test_species_in_metadata(self, csv_linelist_path):
        reader = CsvLineListReader()
        result = reader.read(csv_linelist_path)
        assert "species" in result.metadata
        assert "H2O" in result.metadata["species"]

    def test_empty_file(self, tmp_dir):
        p = tmp_dir / "empty.csv"
        p.write_text("species,lam,a_stein\n")
        reader = CsvLineListReader()
        result = reader.read(p)
        assert result.lines_data == []

    def test_nr_defaults_to_index(self, csv_linelist_path):
        reader = CsvLineListReader()
        result = reader.read(csv_linelist_path)
        assert result.lines_data[0]["nr"] == 0
        assert result.lines_data[1]["nr"] == 1


# ────────────────────────────────────────────────────────────────────
#  SavedLinesReader
# ────────────────────────────────────────────────────────────────────

class TestSavedLinesReader:
    def test_protocol_conformance(self):
        assert isinstance(SavedLinesReader(), LineListReader)

    def test_read_basic(self, saved_lines_path):
        reader = SavedLinesReader()
        result = reader.read(saved_lines_path)

        assert isinstance(result, ReadResult)
        assert len(result.lines_data) == 1
        assert result.metadata["source_format"] == "saved"

    def test_core_fields_present(self, saved_lines_path):
        reader = SavedLinesReader()
        result = reader.read(saved_lines_path)
        line = result.lines_data[0]
        for field in CORE_FIELD_NAMES:
            assert field in line, f"Missing core field: {field}"

    def test_fit_columns_in_extras(self, saved_lines_path):
        reader = SavedLinesReader()
        result = reader.read(saved_lines_path)
        extras = result.extra_columns
        # These should be in extras, not in core fields
        for col in ("Flux_data", "Err_data", "Fit_SN", "FWHM_fit",
                     "Centr_fit", "Red-chisq"):
            assert col in extras, f"Missing extra column: {col}"

    def test_species_in_metadata(self, saved_lines_path):
        reader = SavedLinesReader()
        result = reader.read(saved_lines_path)
        assert "species" in result.metadata
        assert "H2O" in result.metadata["species"]


# ────────────────────────────────────────────────────────────────────
#  HitranParReader (protocol check only — actual parsing tested elsewhere)
# ────────────────────────────────────────────────────────────────────

class TestHitranParReader:
    def test_protocol_conformance(self):
        from iSLAT.Modules.FileHandling.line_list_readers import HitranParReader
        assert isinstance(HitranParReader(), LineListReader)


# ────────────────────────────────────────────────────────────────────
#  DEFAULT_READERS registry
# ────────────────────────────────────────────────────────────────────

class TestDefaultReaders:
    def test_contains_hitran(self):
        assert "hitran" in DEFAULT_READERS

    def test_contains_csv(self):
        assert "csv" in DEFAULT_READERS

    def test_contains_saved(self):
        assert "saved" in DEFAULT_READERS

    def test_all_satisfy_protocol(self):
        for name, reader in DEFAULT_READERS.items():
            assert isinstance(reader, LineListReader), \
                f"Reader {name!r} does not satisfy LineListReader protocol"


# ────────────────────────────────────────────────────────────────────
#  MoleculeLineList.from_file
# ────────────────────────────────────────────────────────────────────

class TestMoleculeLineListFromFile:
    def test_from_csv_file(self, csv_linelist_path):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(csv_linelist_path, molecule_id="H2O")
        assert ll.molecule_id == "H2O"
        assert ll.num_lines == 2
        assert ll.source_format == "csv"

    def test_extra_fields_from_csv(self, csv_linelist_path):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(csv_linelist_path)
        assert "xmin" in ll.extra_fields
        assert "xmax" in ll.extra_fields
        assert len(ll.extra_fields["xmin"]) == 2

    def test_from_saved_file(self, saved_lines_path):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(saved_lines_path, molecule_id="H2O")
        assert ll.source_format == "saved"
        assert ll.num_lines == 1
        assert "Flux_data" in ll.extra_fields

    def test_wavelengths_correct(self, csv_linelist_path):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(csv_linelist_path, molecule_id="H2O")
        wvl = ll.get_wavelengths()
        assert len(wvl) == 2
        np.testing.assert_allclose(wvl[0], 12.407, atol=1e-3)
        np.testing.assert_allclose(wvl[1], 14.950, atol=1e-3)

    def test_frequencies_computed(self, csv_linelist_path):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(csv_linelist_path, molecule_id="H2O")
        freqs = ll.get_frequencies()
        expected_0 = c.SPEED_OF_LIGHT_MICRONS / 12.407
        assert abs(freqs[0] - expected_0) < 1e6

    def test_molecule_id_derived_from_species(self, csv_linelist_path):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(csv_linelist_path)
        assert ll.molecule_id == "H2O"

    def test_unknown_format_raises(self, tmp_dir):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        p = tmp_dir / "data.fits"
        p.write_text("SIMPLE = T")
        with pytest.raises(ValueError):
            MoleculeLineList.from_file(p)

    def test_explicit_format_override(self, csv_linelist_path):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(
            csv_linelist_path, molecule_id="H2O", format="csv"
        )
        assert ll.source_format == "csv"
        assert ll.num_lines == 2


# ────────────────────────────────────────────────────────────────────
#  MoleculeLineList.register_reader
# ────────────────────────────────────────────────────────────────────

class TestRegisterReader:
    def test_register_custom_reader(self, tmp_dir):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        class DummyReader:
            def read(self, filepath):
                return ReadResult(
                    lines_data=[{
                        "nr": 0, "lev_up": "", "lev_low": "",
                        "lam": 99.0, "freq": 3e12,
                        "a_stein": 0.1, "e_up": 100.0, "e_low": 50.0,
                        "g_up": 1, "g_low": 1,
                    }],
                    metadata={"source_format": "dummy"},
                )

        MoleculeLineList.register_reader("dummy", DummyReader())

        try:
            p = tmp_dir / "test.csv"  # extension doesn't matter — explicit format
            p.write_text("dummy content")
            ll = MoleculeLineList.from_file(p, format="dummy", molecule_id="test")
            assert ll.num_lines == 1
            np.testing.assert_allclose(ll.get_wavelengths()[0], 99.0)
        finally:
            # Clean up — remove from class-level registry
            MoleculeLineList._reader_registry.pop("dummy", None)


# ────────────────────────────────────────────────────────────────────
#  Extra-fields round-trip
# ────────────────────────────────────────────────────────────────────

class TestExtraFieldsRoundTrip:
    def test_extra_fields_on_init(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        lines_data = [
            {"nr": 0, "lev_up": "", "lev_low": "", "lam": 10.0,
             "freq": 3e13, "a_stein": 0.01, "e_up": 100, "e_low": 50,
             "g_up": 1, "g_low": 1},
        ]
        extras = {"note": ["important"], "Flux_data": ["1.5e-14"]}
        ll = MoleculeLineList(
            molecule_id="X",
            lines_data=lines_data,
            extra_fields=extras,
        )
        assert ll.extra_fields == extras
        assert ll.num_lines == 1

    def test_extra_fields_in_pandas_table(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        lines_data = [
            {"nr": 0, "lev_up": "", "lev_low": "", "lam": 10.0,
             "freq": 3e13, "a_stein": 0.01, "e_up": 100, "e_low": 50,
             "g_up": 1, "g_low": 1},
            {"nr": 1, "lev_up": "", "lev_low": "", "lam": 11.0,
             "freq": 2.7e13, "a_stein": 0.02, "e_up": 200, "e_low": 60,
             "g_up": 3, "g_low": 1},
        ]
        extras = {"note": ["a", "b"], "Flux_data": ["1e-14", "2e-14"]}
        ll = MoleculeLineList(
            molecule_id="X",
            lines_data=lines_data,
            extra_fields=extras,
        )
        df = ll.get_pandas_table(include_extras=True)
        assert "note" in df.columns
        assert "Flux_data" in df.columns
        assert list(df["note"]) == ["a", "b"]

    def test_pandas_table_default_no_extras(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        lines_data = [
            {"nr": 0, "lev_up": "", "lev_low": "", "lam": 10.0,
             "freq": 3e13, "a_stein": 0.01, "e_up": 100, "e_low": 50,
             "g_up": 1, "g_low": 1},
        ]
        extras = {"note": ["x"]}
        ll = MoleculeLineList(
            molecule_id="X", lines_data=lines_data, extra_fields=extras,
        )
        df = ll.get_pandas_table()  # include_extras=False by default
        assert "note" not in df.columns


# ────────────────────────────────────────────────────────────────────
#  Binary cache with v3 (source_format + extra_columns)
# ────────────────────────────────────────────────────────────────────

class TestCacheV3:
    def test_cache_round_trip(self, csv_linelist_path, tmp_dir):
        """from_file → _save_to_cache → _load_from_cache preserves extras."""
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList.from_file(csv_linelist_path, molecule_id="H2O")

        cache_dir = str(tmp_dir / "cache_test")
        assert ll._save_to_cache(cache_dir)

        # Create a fresh instance and load from cache
        ll2 = MoleculeLineList(molecule_id="H2O")
        assert ll2._load_from_cache(cache_dir)

        assert ll2.num_lines == ll.num_lines
        assert ll2.source_format == ll.source_format
        assert set(ll2.extra_fields.keys()) == set(ll.extra_fields.keys())
        np.testing.assert_allclose(
            ll2.get_wavelengths(), ll.get_wavelengths()
        )


# ────────────────────────────────────────────────────────────────────
#  MoleculeLineList new __init__ kwargs
# ────────────────────────────────────────────────────────────────────

class TestMoleculeLineListNewKwargs:
    def test_format_kwarg(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList(molecule_id="test", format="csv")
        assert ll.source_format == "csv"

    def test_extra_fields_kwarg(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList(
            molecule_id="test",
            lines_data=[{"lam": 10.0}],
            extra_fields={"note": ["x"]},
        )
        assert ll.extra_fields == {"note": ["x"]}

    def test_backward_compat_no_new_kwargs(self):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList

        ll = MoleculeLineList(molecule_id="test", lines_data=[{"lam": 10.0}])
        assert ll.extra_fields == {}
        assert ll.source_format is None


# ────────────────────────────────────────────────────────────────────
#  LineListMaker integration
# ────────────────────────────────────────────────────────────────────

class TestLineListMakerIntegration:
    def test_from_file_csv(self, csv_linelist_path):
        from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker

        maker = LineListMaker.from_file(csv_linelist_path, molecule_id="H2O")
        assert len(maker) == 2
        assert "xmin" in maker._df.columns

    def test_from_saved_lines(self, saved_lines_path):
        from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker

        maker = LineListMaker.from_saved_lines(saved_lines_path, molecule_id="H2O")
        assert len(maker) == 1
        assert "Flux_data" in maker._df.columns

    def test_from_file_filter_chain(self, csv_linelist_path):
        from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker

        maker = (
            LineListMaker.from_file(csv_linelist_path, molecule_id="H2O")
            .filter_wavelength(min_val=14.0)
        )
        assert len(maker) == 1

    def test_to_linelist_preserves_extras(self, csv_linelist_path):
        from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker

        maker = LineListMaker.from_file(csv_linelist_path, molecule_id="H2O")
        ll = maker.to_linelist()
        assert "xmin" in ll.extra_fields
        assert "xmax" in ll.extra_fields
        assert ll.molecule_id == "H2O"
        assert ll.source_format == "csv"

    def test_to_linelist_round_trip_data(self, csv_linelist_path):
        from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker

        maker = LineListMaker.from_file(csv_linelist_path, molecule_id="H2O")
        ll = maker.to_linelist()

        assert ll.num_lines == 2
        wvl = ll.get_wavelengths()
        np.testing.assert_allclose(wvl[0], 12.407, atol=1e-3)
        np.testing.assert_allclose(wvl[1], 14.950, atol=1e-3)

    def test_existing_linelist_with_extras(self):
        """LineListMaker(MoleculeLineList_with_extras) includes extras in _df."""
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        from iSLAT.Modules.DataProcessing.LineListMaker import LineListMaker

        lines_data = [
            {"nr": 0, "lev_up": "", "lev_low": "", "lam": 10.0,
             "freq": 3e13, "a_stein": 0.01, "e_up": 100, "e_low": 50,
             "g_up": 1, "g_low": 1},
        ]
        extras = {"note": ["important"]}
        ll = MoleculeLineList(
            molecule_id="X", lines_data=lines_data, extra_fields=extras,
        )
        maker = LineListMaker(ll)
        assert "note" in maker._df.columns
        assert maker._df["note"].iloc[0] == "important"
