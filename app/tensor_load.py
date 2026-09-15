"""Safe 3x3 stress-history expressions for boundary-condition preparation."""
from __future__ import annotations

import ast
import math
import numpy as np

ALLOWED_NAMES = {"t", "f", "pi", "normal_mean", "normal_amp", "shear_mean", "shear_amp"}
ALLOWED_FUNCS = {"sin": math.sin, "cos": math.cos}


def default_tensor_expressions():
    return "normal_mean+normal_amp*sin(2*pi*f*t),shear_mean+shear_amp*cos(2*pi*f*t),0;shear_mean+shear_amp*cos(2*pi*f*t),0,0;0,0,0"


def _validate(node):
    if isinstance(node, ast.Expression): return _validate(node.body)
    if isinstance(node, ast.Constant):
        try:
            if type(node.value) in (int, float) and math.isfinite(float(node.value)): return
        except OverflowError:
            pass
        raise ValueError("tensor_expression")
    if isinstance(node, ast.Name):
        if node.id not in ALLOWED_NAMES and node.id not in ALLOWED_FUNCS: raise ValueError("tensor_expression")
        return
    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name) or node.func.id not in ALLOWED_FUNCS or node.keywords: raise ValueError("tensor_expression")
        if len(node.args) != 1: raise ValueError("tensor_expression")
        return _validate(node.args[0])
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)): return _validate(node.operand)
    if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
        _validate(node.left); _validate(node.right); return
    raise ValueError("tensor_expression")


def compile_tensor_matrix(text: str):
    if len(text) > 12000: raise ValueError('tensor_expression_too_large')
    rows = [row.strip() for row in str(text).split(";")]
    if len(rows) != 3: raise ValueError("tensor_shape")
    expressions = [[item.strip() for item in row.split(",")] for row in rows]
    if any(len(row) != 3 or any(not item for item in row) for row in expressions): raise ValueError("tensor_shape")
    compiled = []
    for row in expressions:
        out = []
        for item in row:
            tree = ast.parse(item, mode="eval")
            if sum(1 for _ in ast.walk(tree)) > 256: raise ValueError('tensor_expression_too_large')
            _validate(tree)
            # Stress arithmetic is real floating point, not unbounded integer
            # exponentiation supplied by an untrusted setup/AI draft.
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant): node.value = float(node.value)
            out.append(compile(tree, "<tensor-load>", "eval"))
        compiled.append(out)
    return tuple(tuple(row) for row in compiled)


def evaluate_tensor(compiled, *, t, frequency, normal_mean, normal_amplitude, shear_mean, shear_amplitude):
    values = {"t": float(t), "f": float(frequency), "pi": math.pi,
              "normal_mean": float(normal_mean), "normal_amp": float(normal_amplitude),
              "shear_mean": float(shear_mean), "shear_amp": float(shear_amplitude)}
    scope = {"__builtins__": {}, **ALLOWED_FUNCS}
    try:
        result = np.array([[eval(expr, scope, values) for expr in row] for row in compiled], dtype=float)
    except (ArithmeticError, TypeError) as exc:
        raise ValueError('tensor_finite_real') from exc
    if not np.isfinite(result).all(): raise ValueError("tensor_finite")
    return result


def evaluate_tensor_text(text, **kwargs):
    return evaluate_tensor(compile_tensor_matrix(text), **kwargs)
