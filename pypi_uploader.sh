python setup.py sdist bdist_wheel
twine upload dist/*
rm -rf build/
rm -rf dist/
rm -rf splotch.egg-info/
cd docs
make html
