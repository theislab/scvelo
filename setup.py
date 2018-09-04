# Only use this for installing directly from github
# until this is ready: https://github.com/takluyver/flit/issues/158

from pathlib import Path
from flit import common, inifile
import re, setuptools


fields_same_name = ['name', 'version', 'author', 'author_email', 'description', 'classifiers']
fields_diff_name = {'description': 'summary', 'url': 'home_page'}


def setuptools_requirement_format(req):
    m = re.match(r'(?P<name>[\w\-.]+)\s*(\((?P<vers>.*?)\))?(;(?P<env>.*))?', req)
    name, vers = m.group('name', 'vers')
    if vers:
        return name + ' ' + vers
    else:
        return name


def setup(ini_path=Path('pyproject.toml'), package_data=None):
    ini_info = inifile.read_pkg_ini(ini_path)
    module = common.Module(ini_info['module'], ini_path.parent)
    metadata = common.make_metadata(module, ini_info)

    kwargs = {}
    for field in fields_same_name:
        val = getattr(metadata, field, None)
        if val:
            kwargs[field] = val

    for st_field, metadata_field in fields_diff_name.items():
        val = getattr(metadata, metadata_field, None)
        if val:
            kwargs[st_field] = val

    if module.is_package:
        kwargs['packages'] = setuptools.find_packages(include=[module.name+'*'])
    else:
        kwargs['py_modules'] = [module.name]

    if metadata.requires_dist:
        kwargs['install_requires'] = [setuptools_requirement_format(req) for req in metadata.requires_dist]
        print(kwargs['install_requires'])

    if ini_info['scripts']:
        kwargs['entry_points'] = {'console_scipts':
            ['{} = {}:{}'.format(k, mod, func)
             for k, (mod, func) in ini_info['scripts'].items()]
        }

    if package_data is not None:
        kwargs['package_data'] = package_data

    # You're on your own from here
    setuptools.setup(**kwargs)


setup()
