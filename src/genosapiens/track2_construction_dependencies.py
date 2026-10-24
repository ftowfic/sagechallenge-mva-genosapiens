"""Conservative Python dependency closure for unchanged construction semantics."""
import ast
from pathlib import Path

from genosapiens.track2_operator_inputs import require, small_json
from genosapiens.track2_published_genotype import canonical, sha256

# These owners publish and authenticate bytes. They do not choose scientific
# inputs or calculate results. Each artifact owner checks its storage contract.
STORAGE_DEPENDENCIES = frozenset({
    ('genosapiens.track2_scientific_storage', '*'),
    ('genosapiens.track2_scientific_reader', '*'),
    ('genosapiens.track2_scientific_records', '*'),
    ('genosapiens.track2_scientific_completion', '*'),
    ('genosapiens.track2_discovery_outcome', '*'),
    ('genosapiens.track2_run_state', '*'),
    ('genosapiens.track2_construction_dependencies', '*'),
    ('genosapiens.track2_discovery', 'scientific_identity'),
    ('genosapiens.track2_discovery', '_read'),
    ('genosapiens.track2_adjudication', 'reviewed_input_sha256'),
})

# Only these former storage entries may be absent from a current construction
# dependency closure. All surviving dependencies and every input still match.
RETIRED_CONSTRUCTION_STORAGE = frozenset(
    'genosapiens.'+module+':'+symbol for module, symbols in {
        'track2_construction_dependencies': ('construction_dependencies', 'python_dependencies'),
        'track2_run_state': ('RunState', 'SCHEMA', 'STAGE_SCHEMA', 'atomic_record', 'files'),
        'wgs_run_manifest': ('MAX_BYTES', 'RunManifestError', '_require', 'atomic_current'),
        'track2_adjudication': ('reviewed_input_sha256',),
    }.items() for symbol in symbols)


def construction_binding_matches(saved, current):
    if not isinstance(saved, dict) or not isinstance(current, dict):
        return False
    if {k:v for k,v in saved.items() if k != 'dependencies'} != {
            k:v for k,v in current.items() if k != 'dependencies'}:
        return False
    before, after = saved.get('dependencies'), current.get('dependencies')
    if not isinstance(before, dict) or not isinstance(after, dict):
        return False
    return (before.keys()-after.keys() <= RETIRED_CONSTRUCTION_STORAGE
        and all(key in before and before[key] == value for key,value in after.items()))


def python_dependencies(package, roots, *, ignored=(), explicit_empty_ast=False):
    """Bind reachable Python definitions, the resource contract, and input schema.

    An imported module without a selected symbol binds its complete source.
    External libraries remain bound by the unchanged runtime image.
    """
    package = Path(package)
    pending = list(roots)
    visited, result, parsed = set(), {}, {}
    while pending:
        module, symbol = pending.pop()
        # Successful construction evidence supplies the scientific products.
        # Current authorization, storage, and integrity checks run again during
        # continuation. Presentation and disk-capacity guards do not select alleles.
        if (module, symbol) in ignored or (module, '*') in ignored:
            continue
        key = module+':'+symbol
        if key in visited:
            continue
        visited.add(key)
        path = package/'src'/Path(*module.split('.')).with_suffix('.py')
        if not path.exists():
            path = package/'src'/Path(*module.split('.'))/'__init__.py'
        require(path.is_file(), 'Construction dependency source is missing.')
        if module not in parsed:
            parsed[module] = ast.parse(path.read_bytes())
        tree = parsed[module]
        definitions, imports = {}, {}
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                definitions[node.name] = node
            elif isinstance(node, (ast.Assign, ast.AnnAssign)):
                for child in ast.walk(node):
                    if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
                        definitions[child.id] = node
            elif isinstance(node, (ast.Import, ast.ImportFrom)):
                for alias in node.names:
                    imports[alias.asname or alias.name.split('.')[0]] = (node, alias)
        selected = tree if symbol == '*' else definitions.get(symbol)
        if selected is None and symbol in imports:
            selected = imports[symbol][0]
        require(selected is not None, 'Construction dependency symbol is missing.')
        # Python 3.13 omits empty lists by default. Annotation recovery crosses
        # the host/worker boundary and must retain the Python 3.12 spelling.
        import sys
        options = {'show_empty': True} if explicit_empty_ast and sys.version_info >= (3, 13) else {}
        result[key] = sha256(ast.dump(selected, include_attributes=False, **options).encode())
        nodes = list(ast.walk(selected))
        # Bind configuration read by reachable owners, including the clinical
        # session contract and validation schemas. Code hashes alone cannot
        # detect changes to these runtime inputs.
        for node in nodes:
            if (isinstance(node, ast.Constant) and isinstance(node.value, str)
                    and node.value.startswith(('config/', 'schemas/'))
                    and node.value.endswith('.json')):
                relative=Path(node.value)
                require('..' not in relative.parts, 'A dependency path escapes the package.')
                resources=sorted(package.glob(node.value)) if '*' in node.value else [package/relative]
                require(resources and all(p.is_file() for p in resources),
                    'A required stage configuration is missing.')
                for resource in resources:
                    result['file:'+resource.relative_to(package).as_posix()]=sha256(resource.read_bytes())
        # Include imports used inside a selected function, plus referenced globals.
        referenced = {n.id for n in nodes if isinstance(n, ast.Name)}
        import_nodes = [(n, [a for a in n.names if (a.asname or a.name.split('.')[0]) in referenced or a.name == '*'])
                        for n in nodes if isinstance(n, (ast.Import, ast.ImportFrom))]
        import_nodes += [(imports[name][0], [imports[name][1]]) for name in referenced if name in imports]
        if symbol != '*':
            pending.extend((module, name) for name in referenced if name in definitions and name != symbol)
        for node, aliases in import_nodes:
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith('genosapiens.'):
                pending.extend((node.module, alias.name) for alias in aliases)
            elif isinstance(node, ast.ImportFrom) and node.module == 'genosapiens':
                pending.extend(('genosapiens.'+alias.name, '*') for alias in aliases)
            elif isinstance(node, ast.Import):
                pending.extend((alias.name, '*') for alias in aliases if alias.name.startswith('genosapiens.'))
    return dict(sorted(result.items()))


def construction_dependencies(package):
    package = Path(package)
    result = python_dependencies(package, [('genosapiens.track2_operator', 'construct')], ignored=STORAGE_DEPENDENCIES | {
        ('genosapiens.track2_diagnosis_package', 'render_diagnosis_dossier'),
        ('genosapiens.track2_operator', 'storage'),
        ('genosapiens.track2_operator', 'record_source'),
        ('genosapiens.track2_storage', 'ControlledStore'),
        ('genosapiens.track2_storage', 'output_storage'),
        ('genosapiens.bounded_wgs_config', 'DiskLimits'),
    })
    resources = small_json(package/'resources.json', maximum=128*1024*1024)
    def normalized(value):
        if isinstance(value, str):
            return value.replace(str(package)+'/config/', '$PACKAGE/config/')
        if isinstance(value, dict):
            return {k:normalized(v) for k,v in value.items()}
        if isinstance(value, list):
            return [normalized(v) for v in value]
        return value
    result['resources'] = sha256(canonical(normalized(resources)))
    result['public_evidence'] = sha256((package/'public-evidence/evidence.json').read_bytes())
    for path in sorted((package/'schemas/track2').glob('*diagnosis*')):
        result[str(path.relative_to(package))] = sha256(path.read_bytes())
    return dict(sorted(result.items()))
