"""
GUI for diagram algebras: Gram matrix + dual basis visualization.
"""

import sys
import math
import os
from typing import List, Tuple
import re

import sympy
try:
    from PyQt6 import QtCore, QtGui, QtWidgets
except Exception:
    from PyQt5 import QtCore, QtGui, QtWidgets

    class _QtEnumShim:
        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)

    QtCore.Qt.Orientation = _QtEnumShim(
        Horizontal=QtCore.Qt.Horizontal,
        Vertical=QtCore.Qt.Vertical,
    )
    QtCore.Qt.AlignmentFlag = _QtEnumShim(
        AlignCenter=QtCore.Qt.AlignCenter,
    )
    QtCore.Qt.AspectRatioMode = _QtEnumShim(
        KeepAspectRatio=QtCore.Qt.KeepAspectRatio,
    )
    QtCore.Qt.TransformationMode = _QtEnumShim(
        SmoothTransformation=QtCore.Qt.SmoothTransformation,
    )
    QtWidgets.QHeaderView.ResizeMode = _QtEnumShim(
        Fixed=QtWidgets.QHeaderView.Fixed,
    )
    QtWidgets.QAbstractItemView.EditTrigger = _QtEnumShim(
        NoEditTriggers=QtWidgets.QAbstractItemView.NoEditTriggers,
    )
    QtWidgets.QAbstractItemView.SelectionBehavior = _QtEnumShim(
        SelectRows=QtWidgets.QAbstractItemView.SelectRows,
    )
    QtWidgets.QAbstractItemView.SelectionMode = _QtEnumShim(
        ExtendedSelection=QtWidgets.QAbstractItemView.ExtendedSelection,
    )


from algebras.partition import PartitionAlgebra
from algebras.brauer import BrauerAlgebra
from algebras.walled_brauer import WalledBrauerAlgebra
from algebras.symmetric import SymmetricGroupAlgebra
from algebras.half_partition import HalfPartitionAlgebra


class DiagramRenderer:
    def __init__(self, k: int):
        self.k = k
        self.width = 39 + 12 * k
        self.height = 26 + 8 * k
        self.label_height = int(self.height * 0.30)
        body_height = self.height - self.label_height
        self.top_y = int(self.label_height + body_height * 0.30)
        self.bottom_y = int(self.label_height + body_height * 0.75)
        self.margin_x = int(self.width * 0.12)
        self.node_radius = max(6, int(min(self.width, self.height) * 0.045))

    def _positions(self):
        positions = {}
        if self.k <= 1:
            xs = [self.width // 2]
        else:
            span = self.width - 2 * self.margin_x
            xs = [self.margin_x + i * span / (self.k - 1) for i in range(self.k)]
        for i in range(1, self.k + 1):
            positions[i] = QtCore.QPointF(xs[i - 1], self.top_y)
            positions[-i] = QtCore.QPointF(xs[i - 1], self.bottom_y)
        return positions

    def _blocks_for_element(self, elem):
        if hasattr(elem, "blocks"):
            return [set(b) for b in elem.blocks]
        if isinstance(elem, tuple):
            blocks = []
            for i in range(1, self.k + 1):
                blocks.append({i, -elem[i - 1]})
            return blocks
        return []

    def render_pixmap(self, elem, label: str) -> QtGui.QPixmap:
        pixmap = QtGui.QPixmap(self.width, self.height)
        pixmap.fill(QtGui.QColor("white"))
        painter = QtGui.QPainter(pixmap)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)

        font = QtGui.QFont()
        font.setPointSize(max(10, int(self.label_height * 0.45)))
        font.setBold(True)
        painter.setFont(font)
        painter.setPen(QtGui.QColor("black"))
        painter.drawText(
            QtCore.QRectF(0, 0, self.width, self.label_height),
            QtCore.Qt.AlignmentFlag.AlignCenter,
            label,
        )

        positions = self._positions()
        blocks = self._blocks_for_element(elem)
        colors = []
        for i in range(max(1, len(blocks))):
            hue = int(360 * i / max(1, len(blocks)))
            colors.append(QtGui.QColor.fromHsv(hue, 160, 200))

        for idx, block in enumerate(blocks):
            color = colors[idx % len(colors)]
            brush = QtGui.QBrush(color)
            painter.setBrush(brush)
            painter.setPen(QtGui.QPen(QtGui.QColor("black"), 1))
            for node in block:
                pt = positions[node]
                painter.drawEllipse(
                    QtCore.QPointF(pt.x(), pt.y()),
                    self.node_radius,
                    self.node_radius,
                )

        painter.end()
        return pixmap


def leading_degree(expr, d_symbol):
    expr = sympy.together(expr)
    if expr == 0:
        return None
    num, den = sympy.fraction(expr)

    def degree(poly):
        if poly == 0:
            return -sympy.oo
        if not poly.has(d_symbol):
            return 0
        try:
            return sympy.Poly(poly, d_symbol).degree()
        except Exception:
            return None

    deg_num = degree(num)
    deg_den = degree(den)
    if deg_num in (None, -sympy.oo) or deg_den in (None, -sympy.oo):
        return None
    return deg_num - deg_den


def _safe_float(val):
    try:
        f = float(val)
    except Exception:
        return None
    if math.isnan(f) or math.isinf(f):
        return None
    return f


def degree_color(deg, min_deg, max_deg):
    if deg is None or min_deg is None or max_deg is None:
        return QtGui.QColor(255, 255, 255)
    for x in (deg, min_deg, max_deg):
        if x is sympy.nan or (hasattr(x, "is_nan") and x.is_nan):
            return QtGui.QColor(255, 255, 255)
    if min_deg == max_deg:
        t = 0.5
    else:
        t = (deg - min_deg) / (max_deg - min_deg)
    t = _safe_float(sympy.N(t))
    if t is None:
        return QtGui.QColor(255, 255, 255)
    t = max(0.0, min(1.0, t))
    light = QtGui.QColor(220, 235, 255)
    dark = QtGui.QColor(70, 120, 210)
    r = int(light.red() + (dark.red() - light.red()) * t)
    g = int(light.green() + (dark.green() - light.green()) * t)
    b = int(light.blue() + (dark.blue() - light.blue()) * t)
    return QtGui.QColor(r, g, b)


def numeric_color(val, min_val, max_val):
    if val is None or min_val is None or max_val is None:
        return QtGui.QColor(255, 255, 255)
    if math.isnan(val) or math.isinf(val):
        return QtGui.QColor(255, 255, 255)
    max_abs = max(abs(min_val), abs(max_val))
    if max_abs == 0:
        t = 0.0
    else:
        t = val / max_abs
    if math.isnan(t) or math.isinf(t):
        return QtGui.QColor(255, 255, 255)
    t = max(-1.0, min(1.0, t))
    blue = QtGui.QColor(70, 120, 210)
    white = QtGui.QColor(245, 245, 245)
    red = QtGui.QColor(210, 80, 80)
    if t < 0:
        u = abs(t)
        r = int(white.red() + (blue.red() - white.red()) * u)
        g = int(white.green() + (blue.green() - white.green()) * u)
        b = int(white.blue() + (blue.blue() - white.blue()) * u)
    else:
        u = t
        r = int(white.red() + (red.red() - white.red()) * u)
        g = int(white.green() + (red.green() - white.green()) * u)
        b = int(white.blue() + (red.blue() - white.blue()) * u)
    return QtGui.QColor(r, g, b)


def _is_negative_coeff(coeff):
    if hasattr(coeff, "is_negative") and coeff.is_negative is True:
        return True
    if hasattr(coeff, "is_positive") and coeff.is_positive is True:
        return False
    try:
        return bool(coeff < 0)
    except Exception:
        return str(coeff).startswith("-")


class _DSU:
    def __init__(self, n: int):
        self.p = list(range(n))
        self.r = [0] * n

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1


def _partition_join_block_count(blocks_a, blocks_b, element_to_idx):
    """
    Join of two set partitions as the connected components of the graph where
    we connect any two elements that lie in the same block of either partition.
    Returns the number of blocks (connected components) in the join.
    """
    n = len(element_to_idx)
    dsu = _DSU(n)

    def union_block(block):
        if not block:
            return
        it = iter(block)
        first = next(it, None)
        if first is None:
            return
        f = element_to_idx[first]
        for x in it:
            dsu.union(f, element_to_idx[x])

    for blk in blocks_a:
        union_block(blk)
    for blk in blocks_b:
        union_block(blk)

    roots = {dsu.find(i) for i in range(n)}
    return len(roots)


def format_term(coeff, label):
    if coeff == 1:
        return False, label
    if coeff == -1:
        return True, label
    neg = _is_negative_coeff(coeff)
    coeff_abs = -coeff if neg else coeff
    return neg, f"({coeff_abs}) * {label}"


def format_expr(expr):
    s = str(expr)
    s = s.replace(" ", "")

    def sup(num_str):
        sup_map = {
            "0": "⁰",
            "1": "¹",
            "2": "²",
            "3": "³",
            "4": "⁴",
            "5": "⁵",
            "6": "⁶",
            "7": "⁷",
            "8": "⁸",
            "9": "⁹",
            "-": "⁻",
        }
        return "".join(sup_map.get(ch, ch) for ch in num_str)

    s = re.sub(r"d\*\*(-?\d+)", lambda m: f"d{sup(m.group(1))}", s)
    s = re.sub(r"(\d)\*d", r"\1d", s)
    s = s.replace("*d", "d")
    return s


def compute_order(matrix, d_symbol, ascending=True):
    n = len(matrix)
    scores = []
    for i in range(n):
        degs = []
        for j in range(n):
            deg = leading_degree(matrix[i][j], d_symbol)
            if deg is not None:
                degs.append(deg)
        if degs:
            score = sum(degs) / len(degs)
        else:
            score = 0
        scores.append((score, i))
    scores.sort(key=lambda x: x[0], reverse=not ascending)
    return [idx for _, idx in scores]


def compute_order_columns_then_rows(matrix, d_symbol, ascending=True):
    if not matrix:
        return [], []
    # column order from column degree averages
    cols = list(zip(*matrix))
    col_order = compute_order(cols, d_symbol, ascending=ascending)
    # row order from rows after column permutation
    reordered = [[row[j] for j in col_order] for row in matrix]
    row_order = compute_order(reordered, d_symbol, ascending=ascending)
    return row_order, col_order


def _numeric_value(expr):
    return _safe_float(sympy.N(expr))


def compute_order_numeric(matrix, ascending=True):
    n = len(matrix)
    scores = []
    for i in range(n):
        vals = []
        for j in range(len(matrix[i])):
            val = _numeric_value(matrix[i][j])
            if val is not None:
                vals.append(val)
        if vals:
            score = sum(vals) / len(vals)
        else:
            score = 0
        scores.append((score, i))
    scores.sort(key=lambda x: x[0], reverse=not ascending)
    return [idx for _, idx in scores]


def compute_order_columns_then_rows_numeric(matrix, ascending=True):
    if not matrix:
        return [], []
    cols = list(zip(*matrix))
    col_order = compute_order_numeric(cols, ascending=ascending)
    reordered = [[row[j] for j in col_order] for row in matrix]
    row_order = compute_order_numeric(reordered, ascending=ascending)
    return row_order, col_order


def _normalize_shape(shape):
    if isinstance(shape, list):
        shape = tuple(shape)
    if (
        isinstance(shape, tuple)
        and len(shape) == 2
        and isinstance(shape[1], int)
    ):
        part = shape[0]
        if isinstance(part, list):
            part = tuple(part)
        if isinstance(part, tuple):
            return part
    if isinstance(shape, list):
        return tuple(shape)
    return shape


def _partition_label(part):
    part = _normalize_shape(part)
    if part is None:
        return "()"
    if isinstance(part, tuple):
        return "(" + ",".join(str(x) for x in part) + ")"
    return str(part)


def _irrep_label(irrep):
    if (
        isinstance(irrep, tuple)
        and len(irrep) == 2
        and isinstance(irrep[0], tuple)
        and isinstance(irrep[1], int)
    ):
        return f"{_partition_label(irrep[0])}; l={irrep[1]}"
    if isinstance(irrep, tuple) and len(irrep) == 2 and all(
        isinstance(x, tuple) for x in irrep
    ):
        return f"{_partition_label(irrep[0])} | {_partition_label(irrep[1])}"
    return _partition_label(irrep)


def _path_lex_key(path):
    def _node_label(node):
        if isinstance(node, tuple) and len(node) == 2:
            if isinstance(node[1], int):
                node = node[0]
            elif isinstance(node[0], tuple) and isinstance(node[1], tuple):
                return f"{_partition_label(node[0])}|{_partition_label(node[1])}"
        if isinstance(node, list):
            node = tuple(node)
        return _partition_label(node)

    return "->".join(_node_label(step) for step in reversed(path))


def render_young_diagram(part, cell=12, pad=4):
    rows = list(part)
    if not rows:
        rows = [0]
    width = max(rows) if rows else 0
    pix_w = pad * 2 + max(1, width) * cell
    pix_h = pad * 2 + max(1, len(rows)) * cell
    pixmap = QtGui.QPixmap(pix_w, pix_h)
    pixmap.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(pixmap)
    painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
    pen = QtGui.QPen(QtGui.QColor(40, 40, 40), 1)
    painter.setPen(pen)
    painter.setBrush(QtGui.QBrush(QtGui.QColor(230, 230, 230)))
    for r, row_len in enumerate(rows):
        for c in range(row_len):
            x = pad + c * cell
            y = pad + r * cell
            painter.drawRect(x, y, cell, cell)
    painter.end()
    return pixmap


def render_bipartition(left, right, cell=12, pad=4, gap=8):
    pix_left = render_young_diagram(left, cell=cell, pad=pad)
    pix_right = render_young_diagram(right, cell=cell, pad=pad)
    height = max(pix_left.height(), pix_right.height())
    width = pix_left.width() + gap + pix_right.width()
    pixmap = QtGui.QPixmap(width, height)
    pixmap.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(pixmap)
    y_left = (height - pix_left.height()) // 2
    y_right = (height - pix_right.height()) // 2
    painter.drawPixmap(0, y_left, pix_left)
    painter.drawPixmap(pix_left.width() + gap, y_right, pix_right)
    painter.end()
    return pixmap


def _render_shape_pixmap(shape, cell=10, pad=3):
    shape = _normalize_shape(shape)
    if isinstance(shape, tuple) and len(shape) == 2 and all(
        isinstance(x, tuple) for x in shape
    ):
        return render_bipartition(shape[0], shape[1], cell=cell, pad=pad)
    return render_young_diagram(shape, cell=cell, pad=pad)


def _text_pixmap(text, font):
    metrics = QtGui.QFontMetrics(font)
    width = metrics.horizontalAdvance(text) + 4
    height = metrics.height() + 2
    pixmap = QtGui.QPixmap(width, height)
    pixmap.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(pixmap)
    painter.setFont(font)
    painter.setPen(QtGui.QColor("black"))
    painter.drawText(2, metrics.ascent() + 1, text)
    painter.end()
    return pixmap


def _join_pixmaps_h(pixmaps, gap=4):
    if not pixmaps:
        return QtGui.QPixmap(1, 1)
    width = sum(p.width() for p in pixmaps) + gap * (len(pixmaps) - 1)
    height = max(p.height() for p in pixmaps)
    out = QtGui.QPixmap(width, height)
    out.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(out)
    x = 0
    for p in pixmaps:
        y = (height - p.height()) // 2
        painter.drawPixmap(x, y, p)
        x += p.width() + gap
    painter.end()
    return out


def _join_pixmaps_v(pixmaps, gap=4):
    if not pixmaps:
        return QtGui.QPixmap(1, 1)
    width = max(p.width() for p in pixmaps)
    height = sum(p.height() for p in pixmaps) + gap * (len(pixmaps) - 1)
    out = QtGui.QPixmap(width, height)
    out.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(out)
    y = 0
    for p in pixmaps:
        x = (width - p.width()) // 2
        painter.drawPixmap(x, y, p)
        y += p.height() + gap
    painter.end()
    return out


def _render_path_pixmap(path, cell=10, pad=3):
    if not path:
        return QtGui.QPixmap(1, 1)
    font = QtGui.QFont()
    font.setPointSize(10)
    arrow = _text_pixmap("→", font)
    pieces = []
    for idx, shape in enumerate(path):
        pieces.append(_render_shape_pixmap(shape, cell=cell, pad=pad))
        if idx != len(path) - 1:
            pieces.append(arrow)
    return _join_pixmaps_h(pieces, gap=4)


def _render_matrix_unit_row_label(irrep, path_i, path_j):
    font = QtGui.QFont()
    font.setPointSize(10)
    font.setBold(True)
    label_font = QtGui.QFont()
    label_font.setPointSize(10)
    label_font.setBold(False)

    rho_label = _text_pixmap("rho =", label_font)
    i_label = _text_pixmap("i =", label_font)
    j_label = _text_pixmap("j =", label_font)

    rho_pic = _render_shape_pixmap(irrep, cell=12, pad=3)
    i_pic = _render_path_pixmap(path_i, cell=9, pad=2)
    j_pic = _render_path_pixmap(path_j, cell=9, pad=2)

    row1 = _join_pixmaps_h([rho_label, rho_pic], gap=6)
    row2 = _join_pixmaps_h([i_label, i_pic], gap=6)
    row3 = _join_pixmaps_h([j_label, j_pic], gap=6)

    width = max(row1.width(), row2.width(), row3.width())
    height = row1.height() + row2.height() + row3.height() + 6
    out = QtGui.QPixmap(width, height)
    out.fill(QtGui.QColor("white"))
    painter = QtGui.QPainter(out)
    y = 0
    for row in (row1, row2, row3):
        painter.drawPixmap(0, y, row)
        y += row.height() + 3
    painter.end()
    return out


class DiagramHeaderView(QtWidgets.QHeaderView):
    def __init__(self, orientation, parent=None):
        super().__init__(orientation, parent)
        self._pixmaps = None
        self._mode = "label"

    def set_mode(self, mode):
        self._mode = mode
        self.viewport().update()

    def set_pixmaps(self, pixmaps):
        self._pixmaps = pixmaps
        self.viewport().update()

    def paintSection(self, painter, rect, logicalIndex):
        if self._mode == "diagram" and self._pixmaps is not None:
            if 0 <= logicalIndex < len(self._pixmaps):
                pixmap = self._pixmaps[logicalIndex]
                painter.save()
                painter.fillRect(rect, QtGui.QColor(245, 245, 245))
                target = rect.adjusted(4, 4, -4, -4)
                scaled = pixmap.scaled(
                    target.size(),
                    QtCore.Qt.AspectRatioMode.KeepAspectRatio,
                    QtCore.Qt.TransformationMode.SmoothTransformation,
                )
                x = target.x() + (target.width() - scaled.width()) // 2
                y = target.y() + (target.height() - scaled.height()) // 2
                painter.drawPixmap(x, y, scaled)
                painter.restore()
                return
        super().paintSection(painter, rect, logicalIndex)


class AlgebraGui(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        os.environ.setdefault("REPTHEORY_USE_SAGE", "1")
        self.setWindowTitle("Diagram Algebra Viewer")

        self.algebra = None
        self.renderer = None
        self._gram_S_cache_key = None
        self._gram_S_cache_matrix = None

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        control_panel = QtWidgets.QVBoxLayout()
        layout.addLayout(control_panel)

        row1 = QtWidgets.QHBoxLayout()
        control_panel.addLayout(row1)

        self.algebra_selector = QtWidgets.QComboBox()
        self.algebra_selector.addItems(
            ["partition", "half_partition", "brauer", "walled_brauer", "symmetric"]
        )
        row1.addWidget(QtWidgets.QLabel("Algebra:"))
        row1.addWidget(self.algebra_selector)

        self.k_spin = QtWidgets.QSpinBox()
        self.k_spin.setMinimum(1)
        self.k_spin.setMaximum(10)
        self.k_spin.setValue(2)
        row1.addWidget(QtWidgets.QLabel("k:"))
        row1.addWidget(self.k_spin)

        self.r_spin = QtWidgets.QSpinBox()
        self.r_spin.setMinimum(0)
        self.r_spin.setMaximum(10)
        self.r_spin.setValue(1)
        self.s_spin = QtWidgets.QSpinBox()
        self.s_spin.setMinimum(0)
        self.s_spin.setMaximum(10)
        self.s_spin.setValue(1)
        row1.addWidget(QtWidgets.QLabel("r:"))
        row1.addWidget(self.r_spin)
        row1.addWidget(QtWidgets.QLabel("s:"))
        row1.addWidget(self.s_spin)

        self.d_spin = QtWidgets.QSpinBox()
        self.d_spin.setMinimum(-1000)
        self.d_spin.setMaximum(100000)
        self.d_spin.setValue(5)
        row1.addWidget(QtWidgets.QLabel("d:"))
        row1.addWidget(self.d_spin)

        self.symbolic_d_checkbox = QtWidgets.QCheckBox("Symbolic d")
        self.symbolic_d_checkbox.setChecked(False)
        row1.addWidget(self.symbolic_d_checkbox)

        self.load_button = QtWidgets.QPushButton("Load")
        row1.addWidget(self.load_button)
        row1.addStretch()

        row2 = QtWidgets.QHBoxLayout()
        control_panel.addLayout(row2)

        self.display_selector = QtWidgets.QComboBox()
        self.display_selector.addItems(
            ["gram", "gram_S", "gram_S_fourier", "dual", "irreps", "matrix_units"]
        )
        row2.addWidget(QtWidgets.QLabel("Show:"))
        row2.addWidget(self.display_selector)

        self.order_only_checkbox = QtWidgets.QCheckBox("Show order only (Big-O)")
        row2.addWidget(self.order_only_checkbox)

        self.label_selector = QtWidgets.QComboBox()
        self.label_selector.addItems(["diagram", "label"])
        row2.addWidget(QtWidgets.QLabel("Row/Col Labels:"))
        row2.addWidget(self.label_selector)
        row2.addStretch()

        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        layout.addWidget(splitter)

        self.gram_table = QtWidgets.QTableWidget()
        self.gram_table.setHorizontalHeader(DiagramHeaderView(QtCore.Qt.Orientation.Horizontal, self.gram_table))
        self.gram_table.setVerticalHeader(DiagramHeaderView(QtCore.Qt.Orientation.Vertical, self.gram_table))
        self.gram_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.gram_table.horizontalHeader().setSectionsMovable(True)
        self.gram_table.verticalHeader().setSectionsMovable(True)
        splitter.addWidget(self.gram_table)

        self.dual_table = QtWidgets.QTableWidget()
        self.dual_table.setHorizontalHeader(DiagramHeaderView(QtCore.Qt.Orientation.Horizontal, self.dual_table))
        self.dual_table.setVerticalHeader(DiagramHeaderView(QtCore.Qt.Orientation.Vertical, self.dual_table))
        self.dual_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        splitter.addWidget(self.dual_table)

        self.units_table = QtWidgets.QTableWidget()
        self.units_table.setHorizontalHeader(
            DiagramHeaderView(QtCore.Qt.Orientation.Horizontal, self.units_table)
        )
        self.units_table.setVerticalHeader(
            DiagramHeaderView(QtCore.Qt.Orientation.Vertical, self.units_table)
        )
        self.units_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.units_table.setSelectionBehavior(
            QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows
        )
        self.units_table.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
        )
        splitter.addWidget(self.units_table)

        self.units_angle_label = QtWidgets.QLabel(
            "Angle: select two rows (shift-click)."
        )
        self.units_angle_label.setVisible(False)
        control_panel.addWidget(self.units_angle_label)

        self.units_warning_label = QtWidgets.QLabel("")
        self.units_warning_label.setVisible(False)
        control_panel.addWidget(self.units_warning_label)

        self.matrix_group = QtWidgets.QGroupBox("Irrep Matrix Viewer")
        self.matrix_group_layout = QtWidgets.QVBoxLayout(self.matrix_group)
        self.matrix_group_layout.setSizeConstraint(
            QtWidgets.QLayout.SizeConstraint.SetMinAndMaxSize
        )
        controls = QtWidgets.QHBoxLayout()
        self.matrix_group_layout.addLayout(controls)

        self.irrep_selector = QtWidgets.QComboBox()
        controls.addWidget(QtWidgets.QLabel("Irrep:"))
        controls.addWidget(self.irrep_selector)

        self.element_selector = QtWidgets.QComboBox()
        self.element_selector.addItems(
            ["generators", "all basis elements", "dual basis"]
        )
        controls.addWidget(QtWidgets.QLabel("Elements:"))
        controls.addWidget(self.element_selector)

        self.prev_btn = QtWidgets.QPushButton("Prev")
        self.next_btn = QtWidgets.QPushButton("Next")
        controls.addWidget(self.prev_btn)
        controls.addWidget(self.next_btn)
        controls.addStretch()

        self.element_list = QtWidgets.QListWidget()
        self.matrix_group_layout.addWidget(self.element_list)

        self.matrix_table = QtWidgets.QTableWidget()
        self.matrix_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.matrix_group_layout.addWidget(self.matrix_table)

        self.matrix_norm_label = QtWidgets.QLabel("Frobenius norm: ")
        self.matrix_group_layout.addWidget(self.matrix_norm_label)

        self.irrep_scroll = QtWidgets.QScrollArea()
        self.irrep_scroll.setWidgetResizable(True)
        self.irrep_container = QtWidgets.QWidget()
        self.irrep_layout = QtWidgets.QVBoxLayout(self.irrep_container)
        self.irrep_layout.addWidget(self.matrix_group)
        self.irrep_layout.addStretch()
        self.irrep_scroll.setWidget(self.irrep_container)
        splitter.addWidget(self.irrep_scroll)

        self.algebra_selector.currentTextChanged.connect(self._update_param_visibility)
        self.load_button.clicked.connect(self.load_algebra)
        self.display_selector.currentTextChanged.connect(self._refresh_visibility)
        self.order_only_checkbox.stateChanged.connect(self._refresh_tables)
        self.label_selector.currentTextChanged.connect(self._refresh_tables)
        self.symbolic_d_checkbox.toggled.connect(self._refresh_d_controls)
        self.irrep_selector.currentTextChanged.connect(self._reset_matrix_list)
        self.element_selector.currentTextChanged.connect(self._reset_matrix_list)
        self.element_list.currentTextChanged.connect(self._show_matrix_for_label)
        self.prev_btn.clicked.connect(lambda: self._change_page(-1))
        self.next_btn.clicked.connect(lambda: self._change_page(1))
        self.units_table.itemSelectionChanged.connect(self._update_units_angle)
        self._update_param_visibility()
        self._refresh_d_controls()

    def _update_param_visibility(self):
        alg = self.algebra_selector.currentText()
        if alg == "walled_brauer":
            self.k_spin.hide()
            self.r_spin.show()
            self.s_spin.show()
        elif alg == "symmetric":
            self.k_spin.show()
            self.r_spin.hide()
            self.s_spin.hide()
        else:
            self.k_spin.show()
            self.r_spin.hide()
            self.s_spin.hide()

    def _refresh_d_controls(self):
        is_symbolic = self.symbolic_d_checkbox.isChecked()
        self.d_spin.setEnabled(not is_symbolic)
        self.order_only_checkbox.setVisible(is_symbolic)
        if not is_symbolic:
            self.order_only_checkbox.setChecked(False)

    def _d_is_symbolic(self):
        if self.algebra is not None:
            return getattr(self.algebra, "is_symbolic_d", True)
        return self.symbolic_d_checkbox.isChecked()

    def load_algebra(self):
        alg = self.algebra_selector.currentText()
        symbolic_d = self.symbolic_d_checkbox.isChecked()
        d_value = None if symbolic_d else self.d_spin.value()
        if alg == "partition":
            self.algebra = PartitionAlgebra(self.k_spin.value(), d=d_value, symbolic_d=symbolic_d)
        elif alg == "half_partition":
            self.algebra = HalfPartitionAlgebra(self.k_spin.value(), d=d_value, symbolic_d=symbolic_d)
        elif alg == "brauer":
            self.algebra = BrauerAlgebra(
                self.k_spin.value(), d=d_value, symbolic_d=symbolic_d
            )
        elif alg == "walled_brauer":
            self.algebra = WalledBrauerAlgebra(
                self.r_spin.value(),
                self.s_spin.value(),
                d=d_value,
                symbolic_d=symbolic_d,
            )
        else:
            self.algebra = SymmetricGroupAlgebra(
                self.k_spin.value(), d=d_value, symbolic_d=symbolic_d
            )

        self.renderer = DiagramRenderer(self.algebra.k)
        self._populate_gram_matrix()
        self._populate_dual_basis()
        self._populate_irreps()
        self._refresh_d_controls()
        self._refresh_visibility()

    def _refresh_visibility(self):
        mode = self.display_selector.currentText()
        if mode == "gram":
            self.gram_table.setVisible(True)
            self.dual_table.setVisible(False)
            self.units_table.setVisible(False)
            self.irrep_scroll.setVisible(False)
            self.units_angle_label.setVisible(False)
            self.units_warning_label.setVisible(False)
        elif mode == "gram_S":
            self.gram_table.setVisible(True)
            self.dual_table.setVisible(False)
            self.units_table.setVisible(False)
            self.irrep_scroll.setVisible(False)
            self.units_angle_label.setVisible(False)
            self.units_warning_label.setVisible(False)
            self._populate_gram_matrix_S()
        elif mode == "gram_S_fourier":
            self.gram_table.setVisible(True)
            self.dual_table.setVisible(False)
            self.units_table.setVisible(False)
            self.irrep_scroll.setVisible(False)
            self.units_angle_label.setVisible(False)
            self.units_warning_label.setVisible(True)
            self._populate_gram_matrix_S_fourier()
        elif mode == "dual":
            self.gram_table.setVisible(False)
            self.dual_table.setVisible(True)
            self.units_table.setVisible(False)
            self.irrep_scroll.setVisible(False)
            self.units_angle_label.setVisible(False)
            self.units_warning_label.setVisible(False)
        elif mode == "matrix_units":
            self.gram_table.setVisible(False)
            self.dual_table.setVisible(False)
            self.units_table.setVisible(True)
            self.irrep_scroll.setVisible(False)
            self.units_angle_label.setVisible(True)
            self.units_warning_label.setVisible(True)
            self._populate_matrix_units()
        elif mode == "irreps":
            self.gram_table.setVisible(False)
            self.dual_table.setVisible(False)
            self.units_table.setVisible(False)
            self.irrep_scroll.setVisible(True)
            self.units_angle_label.setVisible(False)
            self.units_warning_label.setVisible(False)

    def _refresh_tables(self):
        if self.algebra is None:
            return
        mode = self.display_selector.currentText()
        if mode == "gram":
            self._populate_gram_matrix()
        elif mode == "gram_S":
            self._populate_gram_matrix_S()
        elif mode == "gram_S_fourier":
            self._populate_gram_matrix_S_fourier()
        elif mode == "dual":
            self._populate_dual_basis()
        elif mode == "matrix_units":
            self._populate_matrix_units()

    def _compute_gram_matrix_S(self):
        """
        Compute (and cache) the diagram-basis Gram matrix for the inner product
          <D,E>_S = d^{-cc(D)/2 - cc(E)/2 + cc(D join E)}.
        Returns an n×n Python list-of-lists of SymPy expressions.
        """
        alg = self.algebra
        if alg is None:
            return []
        key = (alg.algebra_id, alg.k, str(alg.d), bool(getattr(alg, "is_symbolic_d", False)))
        if key == self._gram_S_cache_key and self._gram_S_cache_matrix is not None:
            return self._gram_S_cache_matrix

        basis = alg.basis
        n = alg.dim
        if n == 0:
            return []

        # Determine ground set once (any basis element contains all elements).
        ground = sorted({x for blk in basis[0].blocks for x in blk}, key=lambda t: (abs(t), t))
        element_to_idx = {x: i for i, x in enumerate(ground)}

        # Precompute cc(D) and blocks as plain Python lists (speed).
        block_lists = [[list(blk) for blk in elem.blocks] for elem in basis]
        cc = [len(blks) for blks in block_lists]

        M = [[sympy.S.Zero for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                join_cc = _partition_join_block_count(
                    block_lists[i], block_lists[j], element_to_idx
                )
                exp = (
                    sympy.Rational(-cc[i], 2)
                    + sympy.Rational(-cc[j], 2)
                    + sympy.Integer(join_cc)
                )
                M[i][j] = sympy.simplify(alg.d ** exp)

        self._gram_S_cache_key = key
        self._gram_S_cache_matrix = M
        return M

    def _populate_gram_matrix(self):
        alg = self.algebra
        basis = alg.basis
        labels = [alg.label_of(b) for b in basis]
        matrix = alg.gram_matrix
        n = alg.dim

        self.gram_table.clear()
        self.gram_table.setRowCount(n)
        self.gram_table.setColumnCount(n)

        degrees = []
        values = []
        for i in range(n):
            for j in range(n):
                if self._d_is_symbolic():
                    deg = leading_degree(matrix[i, j], alg.d)
                    if deg is not None:
                        degrees.append(deg)
                else:
                    val = _numeric_value(matrix[i, j])
                    if val is not None:
                        values.append(val)
        min_deg = min(degrees) if degrees else None
        max_deg = max(degrees) if degrees else None
        min_val = min(values) if values else None
        max_val = max(values) if values else None

        matrix_list = [[matrix[i, j] for j in range(n)] for i in range(n)]
        if self._d_is_symbolic():
            row_order, col_order = compute_order_columns_then_rows(
                matrix_list, alg.d, ascending=True
            )
        else:
            row_order, col_order = compute_order_columns_then_rows_numeric(
                matrix_list, ascending=True
            )

        for row_idx, i in enumerate(row_order):
            for col_idx, j in enumerate(col_order):
                val = matrix[i, j]
                display_text = self._format_cell(val, alg.d)
                item = QtWidgets.QTableWidgetItem(display_text)
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                if self._d_is_symbolic():
                    color = degree_color(
                        leading_degree(val, alg.d), min_deg, max_deg
                    )
                else:
                    color = numeric_color(
                        _numeric_value(val), min_val, max_val
                    )
                item.setBackground(QtGui.QBrush(color))
                self.gram_table.setItem(row_idx, col_idx, item)

        pixmaps = [self.renderer.render_pixmap(b, labels[i]) for i, b in enumerate(basis)]
        pixmaps_ordered = [pixmaps[idx] for idx in col_order]
        pixmaps_rows = [pixmaps[idx] for idx in row_order]
        show_mode = self.label_selector.currentText()
        hheader = self.gram_table.horizontalHeader()
        vheader = self.gram_table.verticalHeader()
        hheader.set_mode(show_mode)
        vheader.set_mode(show_mode)
        hheader.set_pixmaps(pixmaps_ordered if show_mode == "diagram" else None)
        vheader.set_pixmaps(pixmaps_rows if show_mode == "diagram" else None)
        for pos, idx in enumerate(col_order):
            if show_mode == "label":
                item = QtWidgets.QTableWidgetItem(labels[idx])
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.gram_table.setHorizontalHeaderItem(pos, item)
            else:
                self.gram_table.setHorizontalHeaderItem(pos, QtWidgets.QTableWidgetItem())

        for pos, idx in enumerate(row_order):
            if show_mode == "label":
                vitem = QtWidgets.QTableWidgetItem(labels[idx])
                vitem.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.gram_table.setVerticalHeaderItem(pos, vitem)
            else:
                self.gram_table.setVerticalHeaderItem(pos, QtWidgets.QTableWidgetItem())

        if pixmaps and show_mode == "diagram":
            hheader = self.gram_table.horizontalHeader()
            vheader = self.gram_table.verticalHeader()
            hheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
            vheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
            hheader.setFixedHeight(pixmaps_ordered[0].height())
            vheader.setFixedWidth(pixmaps_rows[0].width())
            hheader.setMinimumSectionSize(pixmaps_ordered[0].width())
            vheader.setMinimumSectionSize(pixmaps_rows[0].height())
            hheader.setDefaultSectionSize(pixmaps_ordered[0].width())
            vheader.setDefaultSectionSize(pixmaps_rows[0].height())
        elif show_mode == "label":
            self.gram_table.horizontalHeader().setDefaultSectionSize(120)
            self.gram_table.verticalHeader().setDefaultSectionSize(28)
        self.gram_table.resizeColumnsToContents()
        self.gram_table.resizeRowsToContents()

    def _populate_gram_matrix_S(self):
        """
        Display Gram matrix on the diagram basis with respect to
          <D,E>_S = d^{-cc(D)/2 - cc(E)/2 + cc(D join E)}.

        Here cc(.) is the number of blocks of the set partition, and "join"
        is the join in the lattice of set partitions (transitive closure of the
        union of equivalence relations).
        """
        alg = self.algebra
        basis = alg.basis
        labels = [alg.label_of(b) for b in basis]
        n = alg.dim

        M = self._compute_gram_matrix_S()

        self.gram_table.clear()
        self.gram_table.setRowCount(n)
        self.gram_table.setColumnCount(n)

        degrees = []
        values = []
        for i in range(n):
            for j in range(n):
                if self._d_is_symbolic():
                    deg = leading_degree(M[i][j], alg.d)
                    if deg is not None:
                        degrees.append(deg)
                else:
                    val = _numeric_value(M[i][j])
                    if val is not None:
                        values.append(val)
        min_deg = min(degrees) if degrees else None
        max_deg = max(degrees) if degrees else None
        min_val = min(values) if values else None
        max_val = max(values) if values else None

        # IMPORTANT: in this view we always use the same order for rows and columns.
        # We compute a symmetric score per basis index from both the row and column.
        scores = []
        if self._d_is_symbolic():
            for i in range(n):
                degs = []
                for j in range(n):
                    d1 = leading_degree(M[i][j], alg.d)
                    if d1 is not None:
                        degs.append(d1)
                    d2 = leading_degree(M[j][i], alg.d)
                    if d2 is not None:
                        degs.append(d2)
                score = (sum(degs) / len(degs)) if degs else 0
                scores.append((score, i))
        else:
            for i in range(n):
                vals = []
                for j in range(n):
                    v1 = _numeric_value(M[i][j])
                    if v1 is not None:
                        vals.append(v1)
                    v2 = _numeric_value(M[j][i])
                    if v2 is not None:
                        vals.append(v2)
                score = (sum(vals) / len(vals)) if vals else 0
                scores.append((score, i))
        scores.sort(key=lambda x: x[0], reverse=False)
        order = [idx for _, idx in scores]
        row_order = order
        col_order = order

        for row_idx, i in enumerate(row_order):
            for col_idx, j in enumerate(col_order):
                val = M[i][j]
                display_text = self._format_cell(val, alg.d)
                item = QtWidgets.QTableWidgetItem(display_text)
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                if self._d_is_symbolic():
                    color = degree_color(
                        leading_degree(val, alg.d), min_deg, max_deg
                    )
                else:
                    color = numeric_color(
                        _numeric_value(val), min_val, max_val
                    )
                item.setBackground(QtGui.QBrush(color))
                self.gram_table.setItem(row_idx, col_idx, item)

        pixmaps = [self.renderer.render_pixmap(b, labels[i]) for i, b in enumerate(basis)]
        pixmaps_ordered = [pixmaps[idx] for idx in col_order]
        pixmaps_rows = [pixmaps[idx] for idx in row_order]
        show_mode = self.label_selector.currentText()
        hheader = self.gram_table.horizontalHeader()
        vheader = self.gram_table.verticalHeader()
        hheader.set_mode(show_mode)
        vheader.set_mode(show_mode)
        hheader.set_pixmaps(pixmaps_ordered if show_mode == "diagram" else None)
        vheader.set_pixmaps(pixmaps_rows if show_mode == "diagram" else None)

        for pos, idx in enumerate(col_order):
            if show_mode == "label":
                item = QtWidgets.QTableWidgetItem(labels[idx])
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.gram_table.setHorizontalHeaderItem(pos, item)
            else:
                self.gram_table.setHorizontalHeaderItem(pos, QtWidgets.QTableWidgetItem())

        for pos, idx in enumerate(row_order):
            if show_mode == "label":
                vitem = QtWidgets.QTableWidgetItem(labels[idx])
                vitem.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.gram_table.setVerticalHeaderItem(pos, vitem)
            else:
                self.gram_table.setVerticalHeaderItem(pos, QtWidgets.QTableWidgetItem())

        if pixmaps and show_mode == "diagram":
            hheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
            vheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
            hheader.setFixedHeight(pixmaps_ordered[0].height())
            vheader.setFixedWidth(pixmaps_rows[0].width())
            hheader.setMinimumSectionSize(pixmaps_ordered[0].width())
            vheader.setMinimumSectionSize(pixmaps_rows[0].height())
            hheader.setDefaultSectionSize(pixmaps_ordered[0].width())
            vheader.setDefaultSectionSize(pixmaps_rows[0].height())
        elif show_mode == "label":
            self.gram_table.horizontalHeader().setDefaultSectionSize(120)
            self.gram_table.verticalHeader().setDefaultSectionSize(28)
        self.gram_table.resizeColumnsToContents()
        self.gram_table.resizeRowsToContents()

    def _populate_gram_matrix_S_fourier(self):
        """
        Gram matrix for the Fourier basis (matrix units E_{ij}^rho) with respect to <.,.>_S.

        We interpret <x,y>_S by bilinear extension from the diagram basis using the
        diagram-basis S-Gram matrix.
        """
        alg = self.algebra
        if alg is None:
            return

        if self._d_is_symbolic():
            self.units_warning_label.setText(
                "gram_S_fourier: numeric d only (symbolic is too slow)."
            )
            self.gram_table.clear()
            self.gram_table.setRowCount(0)
            self.gram_table.setColumnCount(0)
            return

        n = alg.dim
        if n == 0:
            return
        if n > 120:
            self.units_warning_label.setText(
                f"gram_S_fourier: dim={n} is large; this view is disabled to avoid freezing."
            )
            self.gram_table.clear()
            self.gram_table.setRowCount(0)
            self.gram_table.setColumnCount(0)
            return

        self.units_warning_label.setText("")

        S = self._compute_gram_matrix_S()
        units = alg.matrix_units
        unit_labels = alg.matrix_unit_labels
        if len(units) != n or len(unit_labels) != n:
            self.units_warning_label.setText(
                "gram_S_fourier: matrix units unavailable; try clearing cache/*.json and reload."
            )
            self.gram_table.clear()
            self.gram_table.setRowCount(0)
            self.gram_table.setColumnCount(0)
            return

        # Fixed, shared order (irrep_idx, i, j) on both axes.
        unit_label_keys = []
        for irrep_idx, paths in enumerate(alg.bratteli_paths):
            for i in range(len(paths)):
                for j in range(len(paths)):
                    unit_label_keys.append((irrep_idx, i, j))
        order = sorted(range(n), key=lambda idx: unit_label_keys[idx])
        row_order = order
        col_order = order

        # Convert S and units to float matrices/vectors for speed.
        S_float = [[float(sympy.N(S[i][j])) for j in range(n)] for i in range(n)]
        U = [[0.0] * n for _ in range(n)]
        for alpha, u in enumerate(units):
            for a, ca in u.items():
                U[alpha][a] = float(sympy.N(ca))

        # Compute G = U * S * U^T in O(n^3).
        G = [[0.0] * n for _ in range(n)]
        for i in range(n):
            # r = U[i] * S
            r = [0.0] * n
            ui = U[i]
            for a in range(n):
                ca = ui[a]
                if ca == 0.0:
                    continue
                Sa = S_float[a]
                for b in range(n):
                    r[b] += ca * Sa[b]
            for j in range(n):
                uj = U[j]
                s = 0.0
                for b in range(n):
                    if r[b] != 0.0 and uj[b] != 0.0:
                        s += r[b] * uj[b]
                G[i][j] = s

        self.gram_table.clear()
        self.gram_table.setRowCount(n)
        self.gram_table.setColumnCount(n)

        values = []
        for i in range(n):
            for j in range(n):
                values.append(G[i][j])
        min_val = min(values) if values else None
        max_val = max(values) if values else None

        for row_pos, i in enumerate(row_order):
            for col_pos, j in enumerate(col_order):
                val = sympy.Float(G[i][j])
                display_text = self._format_cell(val, alg.d)
                item = QtWidgets.QTableWidgetItem(display_text)
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                color = numeric_color(G[i][j], min_val, max_val)
                item.setBackground(QtGui.QBrush(color))
                self.gram_table.setItem(row_pos, col_pos, item)

        # Use text labels (matrix unit names) on both axes; keep same order.
        hheader = self.gram_table.horizontalHeader()
        vheader = self.gram_table.verticalHeader()
        hheader.set_mode("label")
        vheader.set_mode("label")
        hheader.set_pixmaps(None)
        vheader.set_pixmaps(None)

        for pos, idx in enumerate(col_order):
            item = QtWidgets.QTableWidgetItem(unit_labels[idx])
            item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
            self.gram_table.setHorizontalHeaderItem(pos, item)

        for pos, idx in enumerate(row_order):
            item = QtWidgets.QTableWidgetItem(unit_labels[idx])
            item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
            self.gram_table.setVerticalHeaderItem(pos, item)

        self.gram_table.horizontalHeader().setDefaultSectionSize(160)
        self.gram_table.verticalHeader().setDefaultSectionSize(48)
        self.gram_table.resizeColumnsToContents()
        self.gram_table.resizeRowsToContents()

    def _populate_dual_basis(self):
        alg = self.algebra
        basis = alg.basis
        labels = [alg.label_of(b) for b in basis]
        dual = alg.dual_basis

        degrees = []
        values = []
        for row in dual:
            for coeff in row.values():
                if self._d_is_symbolic():
                    deg = leading_degree(coeff, alg.d)
                    if deg is not None:
                        degrees.append(deg)
                else:
                    val = _numeric_value(coeff)
                    if val is not None:
                        values.append(val)
        min_deg = min(degrees) if degrees else None
        max_deg = max(degrees) if degrees else None
        min_val = min(values) if values else None
        max_val = max(values) if values else None

        n = alg.dim
        matrix_list = [[dual[i].get(j, sympy.S.Zero) for j in range(n)] for i in range(n)]
        if self._d_is_symbolic():
            row_order, col_order = compute_order_columns_then_rows(
                matrix_list, alg.d, ascending=False
            )
        else:
            row_order, col_order = compute_order_columns_then_rows_numeric(
                matrix_list, ascending=False
            )
        # Group rows by irrep (rho) and keep relative order within each group.
        grouped_order = []
        start = 0
        for paths in alg.bratteli_paths:
            block_size = len(paths) ** 2
            block = set(range(start, start + block_size))
            grouped_order.extend([idx for idx in row_order if idx in block])
            start += block_size
        row_order = grouped_order
        self.dual_table.clear()
        self.dual_table.setRowCount(n)
        self.dual_table.setColumnCount(n)

        pixmaps = [self.renderer.render_pixmap(b, labels[i]) for i, b in enumerate(basis)]
        pixmaps_ordered = [pixmaps[idx] for idx in col_order]
        pixmaps_rows = [pixmaps[idx] for idx in row_order]
        show_mode = self.label_selector.currentText()
        hheader = self.dual_table.horizontalHeader()
        vheader = self.dual_table.verticalHeader()
        hheader.set_mode(show_mode)
        vheader.set_mode(show_mode)
        hheader.set_pixmaps(pixmaps_ordered if show_mode == "diagram" else None)
        vheader.set_pixmaps(pixmaps_rows if show_mode == "diagram" else None)
        for pos, idx in enumerate(col_order):
            if show_mode == "label":
                col_item = QtWidgets.QTableWidgetItem(labels[idx])
                col_item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.dual_table.setHorizontalHeaderItem(pos, col_item)
            else:
                self.dual_table.setHorizontalHeaderItem(pos, QtWidgets.QTableWidgetItem())

        for pos, idx in enumerate(row_order):
            if show_mode == "label":
                row_item = QtWidgets.QTableWidgetItem(f"{labels[idx]}*")
                row_item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.dual_table.setVerticalHeaderItem(pos, row_item)
            else:
                self.dual_table.setVerticalHeaderItem(pos, QtWidgets.QTableWidgetItem())

        if pixmaps and show_mode == "diagram":
            hheader = self.dual_table.horizontalHeader()
            vheader = self.dual_table.verticalHeader()
            hheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
            vheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
            hheader.setFixedHeight(pixmaps_ordered[0].height())
            vheader.setFixedWidth(pixmaps_rows[0].width())
            hheader.setMinimumSectionSize(pixmaps_ordered[0].width())
            vheader.setMinimumSectionSize(pixmaps_rows[0].height())
            hheader.setDefaultSectionSize(pixmaps_ordered[0].width())
            vheader.setDefaultSectionSize(pixmaps_rows[0].height())
        elif show_mode == "label":
            self.dual_table.horizontalHeader().setDefaultSectionSize(120)
            self.dual_table.verticalHeader().setDefaultSectionSize(28)

        for row_pos, i in enumerate(row_order):
            row = dual[i]
            for col_pos, j in enumerate(col_order):
                coeff = row.get(j, sympy.S.Zero)
                display_text = self._format_cell(coeff, alg.d)
                item = QtWidgets.QTableWidgetItem(display_text)
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                if self._d_is_symbolic():
                    color = degree_color(
                        leading_degree(coeff, alg.d), min_deg, max_deg
                    )
                else:
                    color = numeric_color(
                        _numeric_value(coeff), min_val, max_val
                    )
                item.setBackground(QtGui.QBrush(color))
                self.dual_table.setItem(row_pos, col_pos, item)

        self.dual_table.resizeColumnsToContents()
        self.dual_table.resizeRowsToContents()

    def _populate_matrix_units(self):
        alg = self.algebra
        basis = alg.basis
        labels = [alg.label_of(b) for b in basis]
        units = alg.matrix_units
        unit_labels = alg.matrix_unit_labels

        degrees = []
        values = []
        for row in units:
            for coeff in row.values():
                if self._d_is_symbolic():
                    deg = leading_degree(coeff, alg.d)
                    if deg is not None:
                        degrees.append(deg)
                else:
                    val = _numeric_value(coeff)
                    if val is not None:
                        values.append(val)
        min_deg = min(degrees) if degrees else None
        max_deg = max(degrees) if degrees else None
        min_val = min(values) if values else None
        max_val = max(values) if values else None

        n = alg.dim
        self.units_table.clear()
        self.units_table.setRowCount(n)
        self.units_table.setColumnCount(n)

        if len(units) != n:
            self.units_warning_label.setText(
                "Warning: matrix units are not available yet. "
                "Try clearing cache/*.json and reload."
            )
            self.units_table.clearContents()
            return

        matrix_list = [[units[i].get(j, sympy.S.Zero) for j in range(n)] for i in range(n)]
        if self._d_is_symbolic():
            row_order, col_order = compute_order_columns_then_rows(
                matrix_list, alg.d, ascending=False
            )
        else:
            row_order, col_order = compute_order_columns_then_rows_numeric(
                matrix_list, ascending=False
            )
        unit_label_keys = []
        for irrep_idx, paths in enumerate(alg.bratteli_paths):
            for i in range(len(paths)):
                for j in range(len(paths)):
                    unit_label_keys.append((irrep_idx, i, j))
        row_order = sorted(range(n), key=lambda idx: unit_label_keys[idx])
        self._units_row_order = row_order
        self._units_col_order = col_order
        self.units_table.clearSelection()
        if self._d_is_symbolic():
            self.units_angle_label.setText(
                "Angle/norm: available only for numeric d."
            )
        else:
            self.units_angle_label.setText(
                "Angle: select two rows (shift-click) or one row for norm."
            )
        zero_rows = [i for i, row in enumerate(units) if not row]
        if zero_rows:
            self.units_warning_label.setText(
                "Warning: some matrix units are zero. This can happen if the "
                "Gram matrix is singular at this d or if cached data is stale. "
                "Try a different d or clear cache/*.json and reload."
            )
        else:
            self.units_warning_label.setText("")

        pixmaps = [self.renderer.render_pixmap(b, labels[i]) for i, b in enumerate(basis)]
        pixmaps_ordered = [pixmaps[idx] for idx in col_order]
        unit_label_pixmaps = []
        for irrep, paths in zip(alg.irreps, alg.bratteli_paths):
            for i in range(len(paths)):
                for j in range(len(paths)):
                    unit_label_pixmaps.append(
                        _render_matrix_unit_row_label(irrep, paths[i], paths[j])
                    )
        unit_label_pixmaps_ordered = [
            unit_label_pixmaps[idx] for idx in row_order
        ]
        show_mode = self.label_selector.currentText()
        hheader = self.units_table.horizontalHeader()
        vheader = self.units_table.verticalHeader()
        hheader.set_mode(show_mode)
        vheader.set_mode("diagram")
        hheader.set_pixmaps(pixmaps_ordered if show_mode == "diagram" else None)
        vheader.set_pixmaps(unit_label_pixmaps_ordered)

        for pos, idx in enumerate(col_order):
            if show_mode == "label":
                col_item = QtWidgets.QTableWidgetItem(labels[idx])
                col_item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.units_table.setHorizontalHeaderItem(pos, col_item)
            else:
                self.units_table.setHorizontalHeaderItem(pos, QtWidgets.QTableWidgetItem())

        for pos, idx in enumerate(row_order):
            row_item = QtWidgets.QTableWidgetItem(unit_labels[idx])
            row_item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
            self.units_table.setVerticalHeaderItem(pos, row_item)

        max_label_height = None
        if unit_label_pixmaps_ordered:
            max_label_height = max(p.height() for p in unit_label_pixmaps_ordered)
        if pixmaps and show_mode == "diagram":
            hheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
            hheader.setFixedHeight(pixmaps_ordered[0].height())
            hheader.setMinimumSectionSize(pixmaps_ordered[0].width())
            hheader.setDefaultSectionSize(pixmaps_ordered[0].width())
            if unit_label_pixmaps_ordered:
                max_height = max(p.height() for p in unit_label_pixmaps_ordered)
                max_width = max(p.width() for p in unit_label_pixmaps_ordered)
                vheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
                vheader.setDefaultSectionSize(max_height)
                vheader.setFixedWidth(max_width)
            else:
                vheader.setDefaultSectionSize(28)
        elif show_mode == "label":
            self.units_table.horizontalHeader().setDefaultSectionSize(120)
            if unit_label_pixmaps_ordered:
                max_height = max(p.height() for p in unit_label_pixmaps_ordered)
                max_width = max(p.width() for p in unit_label_pixmaps_ordered)
                vheader.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Fixed)
                vheader.setDefaultSectionSize(max_height)
                vheader.setFixedWidth(max_width)
            else:
                self.units_table.verticalHeader().setDefaultSectionSize(28)

        for row_pos, i in enumerate(row_order):
            row = units[i]
            for col_pos, j in enumerate(col_order):
                coeff = row.get(j, sympy.S.Zero)
                display_text = self._format_cell(coeff, alg.d)
                item = QtWidgets.QTableWidgetItem(display_text)
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                if self._d_is_symbolic():
                    color = degree_color(
                        leading_degree(coeff, alg.d), min_deg, max_deg
                    )
                else:
                    color = numeric_color(
                        _numeric_value(coeff), min_val, max_val
                    )
                item.setBackground(QtGui.QBrush(color))
                self.units_table.setItem(row_pos, col_pos, item)

        self.units_table.resizeColumnsToContents()
        self.units_table.resizeRowsToContents()
        if max_label_height is not None:
            for row_pos in range(n):
                self.units_table.setRowHeight(row_pos, max_label_height)

    def _update_units_angle(self):
        if self.algebra is None:
            return
        if not getattr(self, "_units_row_order", None):
            return
        if self._d_is_symbolic():
            self.units_angle_label.setText(
                "Angle/norm: available only for numeric d."
            )
            return
        selected = self.units_table.selectionModel().selectedRows()
        if len(selected) == 0:
            self.units_angle_label.setText(
                "Angle: select two rows (shift-click) or one row for norm."
            )
            return
        if len(selected) > 2:
            self.units_angle_label.setText(
                f"Angle: select one or two rows (currently {len(selected)})."
            )
            return

        row_pos_a = selected[0].row()
        row_a = self._units_row_order[row_pos_a]

        units = self.algebra.matrix_units
        unit_labels = self.algebra.matrix_unit_labels
        col_order = self._units_col_order

        vec_a = [units[row_a].get(j, sympy.S.Zero) for j in col_order]
        norm_sq_a = sympy.simplify(sum(a * a for a in vec_a))
        norm_a = sympy.sqrt(norm_sq_a)

        if len(selected) == 1:
            norm_display = norm_sq_a
            approx_text = ""
            if norm_display.is_number:
                # keep exact form if possible, but also show a float approx
                norm_display = sympy.nsimplify(norm_display, rational=True)
                approx_val = _safe_float(sympy.N(norm_sq_a))
                if approx_val is not None:
                    approx_text = f" ({approx_val:.10f})"

            norm_text = f"Norm^2: {format_expr(norm_display)}{approx_text}."
            label_a = unit_labels[row_a]
            self.units_angle_label.setText(f"{label_a} — {norm_text}")
            return

        row_pos_b = selected[1].row()
        row_b = self._units_row_order[row_pos_b]
        vec_b = [units[row_b].get(j, sympy.S.Zero) for j in col_order]
        dot = sum(a * b for a, b in zip(vec_a, vec_b))
        norm_b = sympy.sqrt(sum(b * b for b in vec_b))

        if norm_a == 0 or norm_b == 0:
            angle_text = "Angle: undefined (zero vector)."
        else:
            cos_theta = dot / (norm_a * norm_b)
            angle_expr = sympy.acos(cos_theta)
            angle_val = _safe_float(sympy.N(angle_expr))
            if angle_val is not None:
                deg_val = angle_val * 180 / math.pi
                angle_text = (
                    f"Angle: {angle_val:.6g} rad ({deg_val:.6g}°)."
                )
            else:
                cos_val = _safe_float(sympy.N(cos_theta))
                if cos_val is not None:
                    angle_text = f"Angle: acos({cos_val:.6g})."
                else:
                    angle_text = f"Angle: acos({format_expr(cos_theta)})."

        label_a = unit_labels[row_a]
        label_b = unit_labels[row_b]
        self.units_angle_label.setText(
            f"{label_a} vs {label_b} — {angle_text}"
        )

    def _clear_layout(self, layout):
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)

    def _populate_irreps(self):
        alg = self.algebra
        irreps = alg.irreps
        paths = alg.bratteli_paths

        self._clear_layout(self.irrep_layout)
        self.irrep_layout.addWidget(self.matrix_group)

        self.irrep_selector.clear()
        self._matrix_labels = []
        self._matrix_page = 0
        for irrep in irreps:
            self.irrep_selector.addItem(_irrep_label(irrep))

        for irrep, irrep_paths in zip(irreps, paths):
            group = QtWidgets.QGroupBox()
            group_layout = QtWidgets.QVBoxLayout(group)

            header = QtWidgets.QHBoxLayout()
            if isinstance(irrep, tuple) and len(irrep) == 2 and all(
                isinstance(x, tuple) for x in irrep
            ):
                pixmap = render_bipartition(irrep[0], irrep[1], cell=12)
                label_text = _irrep_label(irrep)
            else:
                pixmap = render_young_diagram(_normalize_shape(irrep), cell=12)
                label_text = _irrep_label(irrep)

            label = QtWidgets.QLabel()
            label.setPixmap(pixmap)
            header.addWidget(label)
            header_text = QtWidgets.QLabel(
                f"{label_text}   dim={len(irrep_paths)}"
            )
            header_text.setStyleSheet("font-weight: bold;")
            header.addWidget(header_text)
            header.addStretch()
            group_layout.addLayout(header)

            for path in irrep_paths:
                row = QtWidgets.QHBoxLayout()
                for idx, step in enumerate(path):
                    pix = _render_shape_pixmap(step, cell=10, pad=2)
                    step_label = QtWidgets.QLabel()
                    step_label.setPixmap(pix)
                    row.addWidget(step_label)
                    if idx != len(path) - 1:
                        row.addWidget(QtWidgets.QLabel("→"))
                row.addStretch()
                group_layout.addLayout(row)

            self.irrep_layout.addWidget(group)

        self.irrep_layout.addStretch()
        self._refresh_matrix_list()

    def _matrix_elements(self):
        if self.algebra is None:
            return []
        mode = self.element_selector.currentText()
        if mode == "dual basis":
            try:
                mats = self.algebra.irrep_matrices
                if not mats:
                    return []
                return [f"{label}*" for label in sorted(mats[0].keys())]
            except Exception:
                return []
        if mode == "generators":
            try:
                mats = self.algebra.irrep_matrices
                if not mats:
                    return []
                return sorted(mats[0].keys())
            except Exception:
                return []
        labels = [self.algebra.label_of(b) for b in self.algebra.basis]
        if mode == "dual basis":
            return [f"{label}*" for label in labels]
        return labels

    def _dual_irrep_matrix_for_label(self, irrep_idx, label):
        alg = self.algebra
        if alg is None:
            return None
        labels = [alg.label_of(b) for b in alg.basis]
        base_label = label[:-1] if label.endswith("*") else label
        try:
            basis_idx = labels.index(base_label)
        except ValueError:
            return None
        coeffs = alg.dual_basis[basis_idx]
        d_rho = len(alg.bratteli_paths[irrep_idx])
        mat = sympy.zeros(d_rho, d_rho)
        for j_idx, coeff in coeffs.items():
            if coeff == 0:
                continue
            mat += coeff * alg.irrep_matrix_for_label(irrep_idx, labels[j_idx])
        return mat

    def _refresh_matrix_list(self):
        self.element_list.clear()
        labels = self._matrix_elements()
        self._matrix_labels = labels
        page_size = 5
        start = self._matrix_page * page_size
        end = start + page_size
        for label in labels[start:end]:
            self.element_list.addItem(label)
        if self.element_list.count() > 0:
            self.element_list.setCurrentRow(0)
        else:
            self.matrix_table.clear()
            self.matrix_table.setRowCount(0)
            self.matrix_table.setColumnCount(0)

    def _reset_matrix_list(self):
        self._matrix_page = 0
        self._refresh_matrix_list()

    def _change_page(self, delta):
        page_size = 5
        if not self._matrix_labels:
            return
        max_page = (len(self._matrix_labels) - 1) // page_size
        self._matrix_page = max(0, min(max_page, self._matrix_page + delta))
        self._refresh_matrix_list()

    def _show_matrix_for_label(self, label):
        if not label or self.algebra is None:
            return
        try:
            irrep_idx = self.irrep_selector.currentIndex()
            mode = self.element_selector.currentText()
            if mode == "dual basis":
                mat = self._dual_irrep_matrix_for_label(irrep_idx, label)
                if mat is None:
                    raise ValueError("Unknown dual basis element.")
            else:
                mat = self.algebra.irrep_matrix_for_label(irrep_idx, label)
        except Exception:
            self.matrix_table.clear()
            self.matrix_table.setRowCount(0)
            self.matrix_table.setColumnCount(0)
            self.matrix_norm_label.setText("Frobenius norm: ")
            return

        n = mat.rows
        paths = self.algebra.bratteli_paths[irrep_idx]
        order = sorted(range(n), key=lambda idx: _path_lex_key(paths[idx]))
        self.matrix_table.clear()
        self.matrix_table.setRowCount(n)
        self.matrix_table.setColumnCount(n)
        for i in range(n):
            for j in range(n):
                item = QtWidgets.QTableWidgetItem(format_expr(mat[order[i], order[j]]))
                item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                item.setForeground(QtGui.QBrush(QtGui.QColor("black")))
                self.matrix_table.setItem(i, j, item)
        norm_sq = sympy.simplify(
            sum(mat[i, j] * mat[i, j] for i in range(n) for j in range(n))
        )
        norm_text = format_expr(sympy.sqrt(norm_sq))
        self.matrix_norm_label.setText(f"Frobenius norm: {norm_text}")
        self.matrix_table.resizeColumnsToContents()
        self.matrix_table.resizeRowsToContents()
        table_width = (
            self.matrix_table.verticalHeader().width()
            + sum(self.matrix_table.columnWidth(i) for i in range(n))
            + self.matrix_table.frameWidth() * 2
        )
        table_height = (
            self.matrix_table.horizontalHeader().height()
            + sum(self.matrix_table.rowHeight(i) for i in range(n))
            + self.matrix_table.frameWidth() * 2
        )
        target_width = table_width
        target_height = table_height
        if n > 1 and self.irrep_scroll.isVisible():
            viewport = self.irrep_scroll.viewport().size()
            max_width = max(200, viewport.width() - 20)
            max_height = max(200, viewport.height() - 20)
            target_width = min(target_width, max_width)
            target_height = min(target_height, max_height)
        self.matrix_table.setMinimumSize(target_width, target_height)
        self.matrix_table.setMaximumSize(target_width, target_height)
        self.matrix_table.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Fixed,
            QtWidgets.QSizePolicy.Policy.Fixed,
        )
        self.matrix_table.updateGeometry()
        self.matrix_group_layout.invalidate()
        self.matrix_group.adjustSize()
        self.irrep_container.adjustSize()
        self.irrep_scroll.widget().adjustSize()

    def _format_cell(self, val, d_symbol):
        if not (self.order_only_checkbox.isChecked() and self._d_is_symbolic()):
            return format_expr(val)
        if val == 0:
            return "0"
        deg = leading_degree(val, d_symbol)
        if deg is None:
            return "0"
        if deg == 0:
            return "O(1)"
        if deg == 1:
            return "O(d)"
        if deg == -1:
            return "O(1/d)"
        if deg < 0:
            return f"O(1/d^{abs(deg)})"
        return f"O(d^{deg})"


def main():
    app = QtWidgets.QApplication(sys.argv)
    win = AlgebraGui()
    win.resize(900, 700)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
