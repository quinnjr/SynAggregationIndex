from __future__ import annotations

import sys
import types
from pathlib import Path

import pytest

PLUMA = Path(__file__).resolve().parent.parent / "PluMA"
sys.path.insert(0, str(PLUMA))

stub = types.ModuleType("PyPluMA")
stub._prefix = ""
stub.prefix = lambda: stub._prefix
sys.modules["PyPluMA"] = stub

from SynAggregationIndexPlugin import SynAggregationIndexPlugin  # noqa: E402

EXAMPLE_DIR = Path(__file__).resolve().parent / "example"


def write_csv(tmp_path: Path, name: str, header: str, rows: list[str]) -> None:
    (tmp_path / name).write_text(
        "\n".join(['"","' + header + '"'] + rows) + "\n"
    )


def build_params(tmp_path: Path, oligomer: str = "oligomer.csv",
                  monomer: str = "monomer.csv", log_inputs: str | None = None,
                  zero_guard: str | None = None) -> Path:
    lines = ["oligomer\t" + oligomer, "monomer\t" + monomer]
    if log_inputs is not None:
        lines.append("log_inputs\t" + log_inputs)
    if zero_guard is not None:
        lines.append("zero_guard\t" + zero_guard)
    params = tmp_path / "params.txt"
    params.write_text("\n".join(lines) + "\n")
    return params


def run(params: Path, tmp_path: Path, outname: str = "out.tsv") -> list[str]:
    stub._prefix = str(tmp_path)
    plugin = SynAggregationIndexPlugin()
    plugin.input(str(params))
    plugin.run()
    out = tmp_path / outname
    plugin.output(str(out))
    return out.read_text().splitlines()


def parse_row(line: str) -> tuple[str, float, float, float]:
    sample, idx, o, m = line.split("\t")
    return sample, float(idx), float(o), float(m)


def test_agg_index_linear_inputs(tmp_path: Path) -> None:
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",220'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",250'])
    params = build_params(tmp_path)
    lines = run(params, tmp_path)
    assert lines[0] == "sample\tagg_index\toligomer\tmonomer"
    sample, idx, o, m = parse_row(lines[1])
    assert sample == "PD_001"
    assert idx == pytest.approx(220 / (220 + 250))
    assert o == 220
    assert m == 250


def test_log_inputs_true_delogs_before_combining(tmp_path: Path) -> None:
    # log2(221) ~ 7.788..., log2(251) ~ 7.971... just use round numbers:
    # store log2(x+1) values so that 2**v - 1 recovers the linear values.
    import math
    lin_o, lin_m = 220.0, 250.0
    log_o = math.log2(lin_o + 1.0)
    log_m = math.log2(lin_m + 1.0)
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",' + repr(log_o)])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",' + repr(log_m)])
    params = build_params(tmp_path, log_inputs="true")
    lines = run(params, tmp_path)
    sample, idx, o, m = parse_row(lines[1])
    assert o == pytest.approx(lin_o)
    assert m == pytest.approx(lin_m)
    assert idx == pytest.approx(lin_o / (lin_o + lin_m))


def test_log_inputs_false_leaves_values_as_is_even_if_log_scale(
        tmp_path: Path) -> None:
    # Demonstrates documented behavior: with the default log_inputs=false,
    # values are used verbatim (no de-logging), even if they happen to be
    # on a log2 scale (e.g. straight from SynProteinFilter's default output).
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",7.788'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",7.971'])
    params = build_params(tmp_path)  # log_inputs defaults to false
    lines = run(params, tmp_path)
    sample, idx, o, m = parse_row(lines[1])
    assert o == pytest.approx(7.788)
    assert m == pytest.approx(7.971)
    assert idx == pytest.approx(7.788 / (7.788 + 7.971))


def test_zero_guard_used_when_denominator_near_zero(tmp_path: Path) -> None:
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",0'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",0'])
    params = build_params(tmp_path, zero_guard="1e-6")
    lines = run(params, tmp_path)
    sample, idx, o, m = parse_row(lines[1])
    # denom clamped to zero_guard = 1e-6; index = 0 / 1e-6 = 0
    assert idx == pytest.approx(0.0)


def test_zero_guard_default_applies_when_not_specified(tmp_path: Path) -> None:
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",0'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",0'])
    params = build_params(tmp_path)  # zero_guard defaults to 1e-6
    lines = run(params, tmp_path)
    sample, idx, o, m = parse_row(lines[1])
    assert idx == pytest.approx(0.0)


def test_intersection_only_samples_present_in_both(tmp_path: Path) -> None:
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",220', '"PD_002",265', '"ONLY_IN_OLIGOMER",100'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",250', '"PD_002",310', '"ONLY_IN_MONOMER",50'])
    params = build_params(tmp_path)
    lines = run(params, tmp_path)
    samples = [line.split("\t")[0] for line in lines[1:]]
    assert samples == ["PD_001", "PD_002"]
    assert "ONLY_IN_OLIGOMER" not in samples
    assert "ONLY_IN_MONOMER" not in samples


def test_no_shared_samples_raises_value_error(tmp_path: Path) -> None:
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",220'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"CTRL_001",170'])
    params = build_params(tmp_path)
    with pytest.raises(ValueError, match="no shared samples"):
        run(params, tmp_path)


def test_quoted_sample_ids(tmp_path: Path) -> None:
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",220'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",250'])
    params = build_params(tmp_path)
    lines = run(params, tmp_path)
    sample = lines[1].split("\t")[0]
    assert sample == "PD_001"
    assert '"' not in sample


def test_output_tsv_format(tmp_path: Path) -> None:
    write_csv(tmp_path, "oligomer.csv", "alpha_synuclein_oligomer",
              ['"PD_001",220', '"CTRL_001",14'])
    write_csv(tmp_path, "monomer.csv", "alpha_synuclein_monomer",
              ['"PD_001",250', '"CTRL_001",170'])
    params = build_params(tmp_path)
    lines = run(params, tmp_path)
    assert lines[0] == "sample\tagg_index\toligomer\tmonomer"
    assert len(lines) == 3
    for line in lines[1:]:
        fields = line.split("\t")
        assert len(fields) == 4


def test_example_demo_data_reproduces_readme_worked_numbers(
        tmp_path: Path) -> None:
    # example/demo_oligomer.csv and example/demo_monomer.csv ship with the
    # repo; the README's worked example table (rounded to 2 decimals) is
    # derived directly from these files.
    params = build_params(
        tmp_path,
        oligomer=str(EXAMPLE_DIR / "demo_oligomer.csv"),
        monomer=str(EXAMPLE_DIR / "demo_monomer.csv"),
    )
    lines = run(params, tmp_path)
    rows = {}
    for line in lines[1:]:
        sample, idx, o, m = parse_row(line)
        rows[sample] = (idx, o, m)

    idx, o, m = rows["PD_001"]
    assert round(idx, 2) == 0.47
    assert o == 220
    assert m == 250

    idx, o, m = rows["CTRL_001"]
    assert round(idx, 2) == 0.08
    assert o == 14
    assert m == 170
