#!/usr/bin/env python3
"""
paf_liftover.py
Lift a base-pair coordinate from one HG002 haplotype to its homologous
chromosome using a pre-generated minimap2 PAF (with cg:Z CIGAR).

Chains were built as:  query = chrNN_p (PATERNAL),  target = chrNN_m (MATERNAL)
So:  p -> m  is query  -> target
     m -> p  is target -> query

Usage (lib):
    lo = PafLiftover('/.../chains/chr15.paf')
    res = lo.lift_point(80262400, src='p')   # -> {'pos':..., 'strand':+1/-1} or None
    span = lo.lift_span(start_bp, end_bp, src='p')  # -> dict with m_start,m_end,strand
"""
import re
import os

_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def _parse_cigar(cg):
    """Split a cg:Z CIGAR string into (length, operation) pairs."""
    return [(int(n), op) for n, op in _CIGAR_RE.findall(cg)]


class PafBlock:
    """One alignment block from the PAF: its query and target spans, strand and CIGAR.

    __slots__ because a whole-chromosome PAF holds many thousands of these.
    """
    __slots__ = ('qstart', 'qend', 'tstart', 'tend', 'strand', 'cigar')

    def __init__(self, qstart, qend, tstart, tend, strand, cigar):
        self.qstart = qstart
        self.qend = qend
        self.tstart = tstart
        self.tend = tend
        self.strand = strand  # '+' or '-'
        self.cigar = cigar    # list of (len, op)


class PafLiftover:
    """Coordinate lift between the two haplotypes of one chromosome, backed by a minimap2 PAF.

    The chains were built with query = chrNN_p (paternal) and target = chrNN_m (maternal), so
    src='p' lifts query to target and src='m' goes the other way. An absent or unreadable PAF
    leaves the object empty rather than raising, and every lift then returns None.
    """

    def __init__(self, paf_path):
        self.path = paf_path
        self.blocks = []
        self.qlen = 0          # query chromosome length (chrNN_p, PATERNAL)
        self.tlen = 0          # target chromosome length (chrNN_m, MATERNAL)
        if paf_path and os.path.exists(paf_path):
            self._load(paf_path)

    def _load(self, paf_path):
        """Read the PAF into PafBlock objects and record both chromosome lengths."""
        with open(paf_path, encoding='utf-8') as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 12:
                    continue
                if not self.qlen:
                    self.qlen = int(f[1])
                    self.tlen = int(f[6])
                qstart, qend = int(f[2]), int(f[3])
                strand = f[4]
                tstart, tend = int(f[7]), int(f[8])
                cg = None
                for tag in f[12:]:
                    if tag.startswith('cg:Z:'):
                        cg = tag[5:]
                        break
                if cg is None:
                    continue
                self.blocks.append(
                    PafBlock(qstart, qend, tstart, tend, strand, _parse_cigar(cg)))

    # ---- core: lift a query-bp position to target-bp ----
    def _q2t(self, blk, qp):
        """Map query position qp -> target position within block blk."""
        if blk.strand == '+':
            q, t = blk.qstart, blk.tstart
            for n, op in blk.cigar:
                if op in 'M=X':
                    if q <= qp < q + n:
                        return t + (qp - q), 1
                    q += n
                    t += n
                elif op == 'I':            # insertion in query (no target)
                    if q <= qp < q + n:
                        return t, 1        # snap to current target base
                    q += n
                elif op in 'DN':           # deletion / skip in query (target only)
                    t += n
            return None
        else:  # '-' : query consumed from qend downwards as target goes up
            q, t = blk.qend, blk.tstart
            for n, op in blk.cigar:
                if op in 'M=X':
                    if q - n <= qp < q:
                        return t + (q - 1 - qp), -1
                    q -= n
                    t += n
                elif op == 'I':
                    if q - n <= qp < q:
                        return t, -1
                    q -= n
                elif op in 'DN':
                    t += n
            return None

    def _t2q(self, blk, tp):
        """Map target position tp -> query position within block blk."""
        if blk.strand == '+':
            q, t = blk.qstart, blk.tstart
            for n, op in blk.cigar:
                if op in 'M=X':
                    if t <= tp < t + n:
                        return q + (tp - t), 1
                    q += n
                    t += n
                elif op == 'I':
                    q += n
                elif op in 'DN':
                    if t <= tp < t + n:
                        return q, 1
                    t += n
            return None
        else:
            q, t = blk.qend, blk.tstart
            for n, op in blk.cigar:
                if op in 'M=X':
                    if t <= tp < t + n:
                        return q - 1 - (tp - t), -1
                    q -= n
                    t += n
                elif op == 'I':
                    q -= n
                elif op in 'DN':
                    if t <= tp < t + n:
                        return q, -1
                    t += n
            return None

    def lift_point(self, pos, src='p'):
        """src='p' lifts p->m (query->target); src='m' lifts m->p (target->query)."""
        best = None
        for blk in self.blocks:
            if src == 'p':
                if blk.qstart <= pos < blk.qend:
                    r = self._q2t(blk, pos)
                    if r is not None:
                        return {'pos': r[0], 'strand': r[1]}
            else:
                if blk.tstart <= pos < blk.tend:
                    r = self._t2q(blk, pos)
                    if r is not None:
                        return {'pos': r[0], 'strand': r[1]}
        return best

    def lift_anchors(self, start_bp, end_bp, src='p', n=9):
        """Lift n evenly-spaced anchors across [start,end].
        Returns list of dicts {src, pos, strand} for anchors that mapped
        (in source order). Empty list if none map."""
        if end_bp < start_bp:
            start_bp, end_bp = end_bp, start_bp
        out = []
        for i in range(n):
            a = int(start_bp + (end_bp - start_bp) * i / (n - 1))
            r = self.lift_point(a, src=src)
            if r is not None:
                out.append({'src': a, 'pos': r['pos'], 'strand': r['strand']})
        return out

    @staticmethod
    def _median(xs):
        """Median of a list, None if empty. Used to combine the anchor lifts of a span."""
        s = sorted(xs)
        m = len(s)
        if m == 0:
            return None
        return s[m // 2] if m % 2 else 0.5 * (s[m // 2 - 1] + s[m // 2])

    def lift_span(self, start_bp, end_bp, src='p', n=9, tol_factor=3.0, tol_min=200000):
        """Lift n anchors across [start,end] and return the mapped span, with
        outlier rejection.

        Anchors in repeat / segmental-duplication regions can be picked up by a
        small paralogous alignment block far away on the homolog; a single such
        outlier would blow up a naive min/max span.  We therefore keep only the
        anchors within a window (median +/- max(tol_factor*region_w, tol_min))
        of the median mapped position, and compute the span from those.

        Returns dict:
          start, end       mapped span (bp) after outlier rejection
          strand           +1 / -1 (majority of kept anchors)
          n_anchor         anchors tried
          n_ok             anchors that lifted at all
          n_kept           anchors kept after outlier rejection
          n_outlier        anchors dropped as outliers
          wide             True if span still > tol_factor * region width
                           (i.e. likely unreliable even after filtering)
        or None if nothing maps."""
        if end_bp < start_bp:
            start_bp, end_bp = end_bp, start_bp
        region_w = max(end_bp - start_bp, 1)
        anchors = [start_bp + (end_bp - start_bp) * i / (n - 1) for i in range(n)]
        mapped = []   # (pos, strand)
        for a in anchors:
            r = self.lift_point(int(a), src=src)
            if r is not None:
                mapped.append((r['pos'], r['strand']))
        if not mapped:
            return None

        med = self._median([p for p, _ in mapped])
        tol = max(region_w * tol_factor, tol_min)
        kept = [(p, s) for (p, s) in mapped if abs(p - med) <= tol]
        if not kept:                      # everything far from its own median: keep all
            kept = mapped

        pos_k = [p for p, _ in kept]
        strand = 1 if sum(s for _, s in kept) >= 0 else -1
        span_w = max(pos_k) - min(pos_k)
        return {
            'start': min(pos_k),
            'end': max(pos_k),
            'strand': strand,
            'n_anchor': n,
            'n_ok': len(mapped),
            'n_kept': len(kept),
            'n_outlier': len(mapped) - len(kept),
            'wide': span_w > region_w * tol_factor,
        }


if __name__ == '__main__':
    import sys
    paf, pos = sys.argv[1], int(sys.argv[2])
    src = sys.argv[3] if len(sys.argv) > 3 else 'p'
    lo = PafLiftover(paf)
    print('blocks:', len(lo.blocks))
    print('point:', lo.lift_point(pos, src=src))
