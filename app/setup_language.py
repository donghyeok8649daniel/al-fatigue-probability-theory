"""AFT 1: declarative, mesh-bound surface-load language. No code execution."""
import json
import re
from .surface_setup import encode, decode, mesh_id


def emit(mesh, loads, correction=None):
    lines = ['AFT 1', f'mesh {mesh_id(mesh)}', 'units mm MPa model_time']
    for index, load in enumerate(loads, 1):
        lines += [f'load L{index}', 'faces '+','.join(map(str, load.face_indices)),
                  f'frequency {load.frequency!r}',
                  'parameters '+','.join(repr(x) for x in (load.normal_mean_mpa,
                    load.normal_amplitude_mpa, load.shear_mean_mpa, load.shear_amplitude_mpa)),
                  'stress '+load.tensor_expression, 'end']
    if correction is not None:
        lines += ['balance '+','.join(map(str, correction[0]))]
    return '\n'.join(lines)+'\n'


def compile_setup(text, mesh):
    if len(text) > 4_000_000: raise ValueError('AFT file too large')
    lines = [(i, s.strip()) for i, s in enumerate(text.splitlines(), 1)
             if s.strip() and not s.lstrip().startswith('#')]
    if len(lines) < 3 or [s for _, s in lines[:3]] != [
            'AFT 1', f'mesh {mesh_id(mesh)}', 'units mm MPa model_time']:
        raise ValueError('AFT header, mesh or units mismatch')
    loads, names, row, seen, balance = [], set(), None, set(), None
    for number, line in lines[3:]:
        command, _, argument = line.partition(' ')
        try:
            if command == 'load' and row is None:
                if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,63}', argument) or argument in names:
                    raise ValueError('invalid/duplicate load name')
                names.add(argument); row = {}; seen = set()
            elif command == 'end' and row is not None and not argument:
                if seen != {'faces', 'frequency', 'parameters', 'stress'}:
                    raise ValueError('incomplete load')
                loads.append(row); row = None
            elif command == 'balance' and row is None and balance is None:
                balance = [int(x) for x in argument.split(',')]
            elif row is not None and command not in seen:
                seen.add(command)
                if command == 'faces': row['face_indices'] = [int(x) for x in argument.split(',')]
                elif command == 'frequency': row['frequency'] = float(argument)
                elif command == 'parameters':
                    values = [float(x) for x in argument.split(',')]
                    if len(values) != 4: raise ValueError('four parameters required')
                    row.update(zip(('normal_mean_mpa','normal_amplitude_mpa','shear_mean_mpa','shear_amplitude_mpa'), values))
                elif command == 'stress': row['tensor_expression'] = argument
                else: raise ValueError('unknown command')
            else: raise ValueError('unexpected or duplicate command')
        except (ValueError, TypeError) as exc:
            raise ValueError(f'AFT line {number}: {exc}') from exc
    if row is not None: raise ValueError('missing end')
    for row in loads: row.update(region='custom', area_mm2=0)
    value = json.loads(encode(mesh, [], None))
    value.update(loads=loads, correction_faces=balance)
    return decode(json.dumps(value, allow_nan=False), mesh)


def analyze_setup(mesh, loads, correction_faces):
    """Offline deterministic diagnostics, NOT AI or a solver certification."""
    return dict(schema='aft.analysis/1', load_count=len(loads),
        loaded_face_count=len({i for load in loads for i in load.face_indices}),
        mesh_sha256=mesh_id(mesh), correction_requested=correction_faces is not None,
        spatial_solver_connected=False, physical_time_calibrated=False,
        material_validation='unvalidated', crack_resolution='not_assessed',
        status='validated_surface_setup_only')
