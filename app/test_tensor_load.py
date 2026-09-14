import numpy as np
import pytest
from app.tensor_load import compile_tensor_matrix, evaluate_tensor, default_tensor_expressions


def test_default_tensor_has_sine_normal_and_cosine_shear():
    result = evaluate_tensor(compile_tensor_matrix(default_tensor_expressions()), t=.25,
        frequency=1, normal_mean=10, normal_amplitude=20, shear_mean=3, shear_amplitude=4)
    np.testing.assert_allclose(result, [[30, 3, 0], [3, 0, 0], [0, 0, 0]])


def test_custom_tensor_is_three_by_three_and_finite():
    result = evaluate_tensor(compile_tensor_matrix('t,0,0;0,cos(pi*f*t),0;0,0,7'),
        t=0, frequency=2, normal_mean=0, normal_amplitude=0, shear_mean=0, shear_amplitude=0)
    np.testing.assert_allclose(result, [[0,0,0],[0,1,0],[0,0,7]])


@pytest.mark.parametrize('expr', ['t,0;0,0', 't,unknown,0;0,0,0;0,0,0',
                                  '__import__("os"),0,0;0,0,0;0,0,0',
                                  'sin(t,0),0,0;0,0,0;0,0,0'])
def test_tensor_expression_rejects_bad_shape_names_or_calls(expr):
    with pytest.raises((ValueError, SyntaxError)):
        compile_tensor_matrix(expr)
