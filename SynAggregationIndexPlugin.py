import math
import os

import PyIO
import PyPluMA


class SynAggregationIndexPlugin:

    def input(self, filename):
        self.parameters = PyIO.readParameters(filename)
        prefix = PyPluMA.prefix()
        self.oligomer_path = os.path.join(prefix, self.parameters["oligomer"])
        self.monomer_path = os.path.join(prefix, self.parameters["monomer"])
        self.zero_guard = float(self.parameters.get("zero_guard", "1e-6"))
        self.log_inputs = self.parameters.get("log_inputs", "false").strip().lower() == "true"

    def run(self):
        olig = _read_single_feature_csv(self.oligomer_path)
        mono = _read_single_feature_csv(self.monomer_path)
        shared = [s for s in olig if s in mono]
        if not shared:
            raise ValueError(
                "SynAggregationIndex: no shared samples between "
                + self.oligomer_path + " and " + self.monomer_path
            )
        self.rows = []
        for sample in shared:
            o, m = olig[sample], mono[sample]
            if self.log_inputs:
                o = max(2.0 ** o - 1.0, 0.0)
                m = max(2.0 ** m - 1.0, 0.0)
            denom = (o + m) if (o + m) > self.zero_guard else self.zero_guard
            idx = o / denom
            self.rows.append((sample, idx, o, m))

    def output(self, filename):
        with open(filename, "w") as out:
            out.write("sample\tagg_index\toligomer\tmonomer\n")
            for sample, idx, o, m in self.rows:
                out.write(
                    sample + "\t" + _fmt(idx) + "\t" + _fmt(o) + "\t" + _fmt(m) + "\n"
                )


def _read_single_feature_csv(path):
    values = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split(",")]
            if len(parts) < 2:
                continue
            sample = parts[0].strip('"')
            try:
                values[sample] = float(parts[1])
            except ValueError:
                continue
    return values


def _fmt(v):
    if math.isnan(v) or math.isinf(v):
        return "NA"
    return repr(v)
