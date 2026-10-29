"""Physical alignment ownership and conservative graph-only source contexts."""

from collections import Counter, defaultdict
from dataclasses import dataclass
import re
import numpy as np
from .graph_io import reverse_complement
from .path_context import context_runs, occurrence_seeds, extend_terminal_context


@dataclass(frozen=True)
class Alignment:
    identity: str
    path: str
    strand: str
    query_start: int
    query_end: int
    ref_start: int
    ref_end: int
    runs: tuple

    @property
    def index(self):
        runs = list(self.runs)
        ends = []
        end = -1
        for q, stop, ref in runs:
            end = max(end, stop)
            ends.append(end)
        return runs, [r[0] for r in runs], ends


def read_alignments(path, minimum_mapq=20):
    """Keep primary record identity; never connect different PAF records implicitly."""
    result = defaultdict(list)
    with open(path) as stream:
        for number, line in enumerate(stream, 1):
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                raise ValueError("Malformed PAF row")
            tags = {item[:2]: item[5:] for item in fields[12:]}
            if int(fields[11]) < minimum_mapq or tags.get("tp") != "P":
                continue
            if "cg" not in tags or fields[4] not in ("+", "-"):
                raise ValueError("PAF needs primary, oriented --eqx CIGAR alignments")
            q = int(fields[2]) if fields[4] == "+" else int(fields[1]) - int(fields[3])
            qstart = q
            ref = int(fields[7])
            runs = []
            pieces = re.findall(r"(\d+)([MIDNSHP=X])", tags["cg"])
            if "".join(n + op for n, op in pieces) != tags["cg"] or any(
                op == "M" for _, op in pieces
            ):
                raise ValueError("Exact =/X CIGAR operations are required")
            for number_text, op in pieces:
                n = int(number_text)
                if op == "=":
                    runs.append((q, q + n, ref))
                if op in "MIS=X":
                    q += n
                if op in "MDN=X":
                    ref += n
            expected_q = int(fields[3]) if fields[4] == "+" else int(fields[1]) - int(fields[2])
            if q != expected_q or ref != int(fields[8]):
                raise ValueError("PAF CIGAR does not consume its declared intervals")
            result[fields[0]].append(
                Alignment(
                    f"M{number:08d}",
                    fields[0],
                    fields[4],
                    qstart,
                    q,
                    int(fields[7]),
                    ref,
                    tuple(runs),
                )
            )
    return result


def unique_oriented(word, reference, query):
    reverse = reverse_complement(word)
    if set(word) - set("ACGT") or word == reverse:
        return False
    return all(
        seq.find(word) >= 0 and seq.find(word, seq.find(word) + 1) < 0 and reverse not in seq
        for seq in (reference, query)
    )


def mapped_chains(
    records, reference, nodes, ref_occurrences, walk, offsets, query, minimum_anchor, mapping_k
):
    chains = []
    for record in records:
        if record.strand != "+":
            continue
        seeds = occurrence_seeds(walk, offsets, ref_occurrences, nodes, record.index)
        seeds = [
            (r, q, n, node)
            for r, q, n, node in seeds
            if record.query_start <= q
            and q + n <= record.query_end
            and record.ref_start <= r
            and r + n <= record.ref_end
        ]
        runs = context_runs(reference, query, seeds, minimum_anchor)
        runs = extend_terminal_context(reference, query, runs, record.index, mapping_k)
        if runs:
            chains.append(dict(identity=record.identity, kind="mapped", runs=runs, record=record))
    return chains


def graph_chains(
    nodes,
    refwalk,
    ref_offsets,
    walk,
    query,
    reference,
    minimum_anchor=31,
    maximum_gap=1000,
    places=None,
):
    """All reference-sharing steps constrain graph-only contexts, including repeats.

    Unlike mapped components, unaligned graph contexts require bounded gaps and
    unique sequence anchors in both orientations. Repeated/opposite steps break
    the chain; they cannot be skipped to manufacture a convenient long deletion.
    """
    if places is None:
        places = defaultdict(list)
        for i, node in enumerate(refwalk):
            places[abs(node)].append(i)
    counts = Counter(abs(node) for node in walk)
    chunks = []
    run = []
    q = 0
    for node in walk:
        length = len(nodes[abs(node)])
        targets = places.get(abs(node), [])
        if not targets:
            q += length
            continue
        valid = (
            len(targets) == 1
            and counts[abs(node)] == 1
            and refwalk[targets[0]] == node
            and not (set(nodes[abs(node)]) - set("ACGT"))
        )
        r = ref_offsets[targets[0]] if valid else -1
        if run and (
            not valid
            or r < run[-1][0] + run[-1][2]
            or max(r - run[-1][0] - run[-1][2], q - run[-1][1] - run[-1][2]) > maximum_gap
            or set(query[run[-1][1] + run[-1][2] : q]) - set("ACGT")
        ):
            chunks.append(run)
            run = []
        if valid:
            if run and (r, q) == (run[-1][0] + run[-1][2], run[-1][1] + run[-1][2]):
                old = run[-1]
                run[-1] = (old[0], old[1], old[2] + length, old[3] + [(r, q, length, node)])
            else:
                run.append((r, q, length, [(r, q, length, node)]))
        q += length
    if run:
        chunks.append(run)
    out = []
    for chunk in chunks:
        boundary = [
            i
            for i, (r, q, n, payload) in enumerate(chunk)
            if n >= minimum_anchor
            and unique_oriented(query[q : q + minimum_anchor], reference, query)
            and unique_oriented(query[q + n - minimum_anchor : q + n], reference, query)
        ]
        if boundary:
            out.append(chunk[boundary[0] : boundary[-1] + 1])
    return out


def intervals(mask):
    edges = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(np.int8)))
    return list(zip(edges[::2].tolist(), edges[1::2].tolist()))


def clip_runs(runs, left, right, reference, query, anchor):
    """Cut only inside unique exact matches, never through an edit."""
    clipped = []
    for r, q, n, payload in runs:
        lo = max(left, q)
        hi = min(right, q + n)
        if lo < hi:
            clipped.append((r + lo - q, lo, hi - lo, [(r + lo - q, lo, hi - lo, None)]))
    valid = [
        i
        for i, (r, q, n, p) in enumerate(clipped)
        if n >= anchor
        and unique_oriented(query[q : q + anchor], reference, query)
        and unique_oriented(query[q + n - anchor : q + n], reference, query)
    ]
    return clipped[valid[0] : valid[-1] + 1] if valid else []


def plan_path(graph, name, records, minimum_anchor=31, mapping_k=19, maximum_gap=1000):
    """Build nonoverlapping physical query sources before any variant is aligned."""
    from .graph_io import offsets

    sequence = graph.spell(name)
    length = len(sequence)
    ref_occ = graph.ref_occurrences
    mapped = []
    oriented = {}
    for strand in ("+", "-"):
        walk = graph.walks[name] if strand == "+" else [-n for n in reversed(graph.walks[name])]
        query = sequence if strand == "+" else reverse_complement(sequence)
        oriented[strand] = (walk, offsets(graph.nodes, walk), query)
        # Coordinates of each record are already expressed on its oriented query.
        local = [
            Alignment(
                r.identity, r.path, "+", r.query_start, r.query_end, r.ref_start, r.ref_end, r.runs
            )
            for r in records
            if r.strand == strand
        ]
        for block in mapped_chains(
            local,
            graph.reference,
            graph.nodes,
            ref_occ,
            walk,
            oriented[strand][1],
            query,
            minimum_anchor,
            mapping_k,
        ):
            block["strand"] = strand
            mapped.append(block)
    owners = np.zeros(length, dtype=np.int16)
    for block in mapped:
        runs = block["runs"]
        lo = runs[0][1]
        hi = runs[-1][1] + runs[-1][2]
        if block["strand"] == "-":
            lo, hi = length - hi, length - lo
        owners[lo:hi] += 1
    accepted = []
    occupied = np.zeros(length, dtype=bool)
    for block in mapped:
        strand = block["strand"]
        walk, off, query = oriented[strand]
        runs = block["runs"]
        lo = runs[0][1]
        hi = runs[-1][1] + runs[-1][2]
        allowed = (owners if strand == "+" else owners[::-1]) == 1
        for left, right in intervals(allowed[lo:hi]):
            left += lo
            right += lo
            selected = (
                runs
                if (left, right) == (lo, hi)
                else clip_runs(runs, left, right, graph.reference, query, mapping_k)
            )
            if not selected:
                continue
            item = dict(block, runs=selected)
            accepted.append(item)
            x = selected[0][1]
            y = selected[-1][1] + selected[-1][2]
            if strand == "-":
                x, y = length - y, length - x
            occupied[x:y] = True
    # Graph-only contexts occupy genuinely unclaimed territory. Areas of mapper
    # ownership disagreement stay unresolved, even if a graph-only seed exists.
    available = ~occupied & (owners < 2)
    supplements = []
    for strand, (walk, off, query) in oriented.items():
        for j, runs in enumerate(
            graph_chains(
                graph.nodes,
                graph.refwalk,
                graph.ref_offsets,
                walk,
                query,
                graph.reference,
                minimum_anchor,
                maximum_gap,
                graph.ref_places,
            )
        ):
            lo = runs[0][1]
            hi = runs[-1][1] + runs[-1][2]
            mask = available if strand == "+" else available[::-1]
            for left, right in intervals(mask[lo:hi]):
                selected = clip_runs(
                    runs, left + lo, right + lo, graph.reference, query, minimum_anchor
                )
                if selected:
                    supplements.append(
                        dict(
                            identity=f"G{j:06d}{strand}",
                            kind="graph",
                            runs=selected,
                            strand=strand,
                            record=None,
                        )
                    )
    supplement_counts = np.zeros(length, dtype=np.int16)
    for block in supplements:
        rs = block["runs"]
        lo = rs[0][1]
        hi = rs[-1][1] + rs[-1][2]
        if block["strand"] == "-":
            lo, hi = length - hi, length - lo
        supplement_counts[lo:hi] += 1
    for block in supplements:
        rs = block["runs"]
        lo = rs[0][1]
        hi = rs[-1][1] + rs[-1][2]
        if block["strand"] == "-":
            lo, hi = length - hi, length - lo
        if np.all(supplement_counts[lo:hi] == 1):
            accepted.append(block)
            occupied[lo:hi] = True
    return accepted, intervals(~occupied), oriented
